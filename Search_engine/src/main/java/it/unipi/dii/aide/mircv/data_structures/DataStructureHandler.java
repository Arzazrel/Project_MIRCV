package it.unipi.dii.aide.mircv.data_structures;

import it.unipi.dii.aide.mircv.TextProcessor;
import it.unipi.dii.aide.mircv.compression.Unary;
import it.unipi.dii.aide.mircv.QueryProcessor;
import it.unipi.dii.aide.mircv.compression.VariableBytes;

import java.io.*;
import java.nio.CharBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import static it.unipi.dii.aide.mircv.data_structures.DocumentElement.*;
import static it.unipi.dii.aide.mircv.data_structures.PartialIndexBuilder.*;
import static it.unipi.dii.aide.mircv.utils.Constants.*;
import static it.unipi.dii.aide.mircv.utils.Logger.*;

/**
 * This class handles the storage and retrieval of data structures used for document indexing.
 */
public final class DataStructureHandler
{
    // -------- start: functions to store into disk --------

    /**
     * Function to store the whole document table into disk.
     */
    static void storeDocumentTableIntoDisk()
    {
        printLoad("\nStoring block offsets into disk...");

        try (RandomAccessFile raf = new RandomAccessFile(DOCTABLE_FILE, "rw");
             FileChannel channel = raf.getChannel())
        {
            MappedByteBuffer buffer = channel.map(FileChannel.MapMode.READ_WRITE, channel.size(), (long) DOCELEM_SIZE * PartialIndexBuilder.documentTable.size());

            if(buffer == null)      // Buffer not created
                return;
            // scan all document elements of the Document Table
            for(DocumentElement de: PartialIndexBuilder.documentTable.values())
            {
                CharBuffer charBuffer = CharBuffer.allocate(DOCNO_DIM);     //allocate bytes for docno

                for (int i = 0; i < de.getDocno().length(); i++)    //put every char into charbuffer
                    charBuffer.put(i, de.getDocno().charAt(i));

                // write docno, docid and doclength into document file
                buffer.put(StandardCharsets.UTF_8.encode(charBuffer));
                buffer.putInt(de.getDocid());
                buffer.putInt(de.getDoclength());

                if(debug)
                    docTable_logger.logInfo(de.toString());
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Function to store offset of the blocks into disk.
     */
    static void storeBlockOffsetsIntoDisk()
    {
        printLoad("\nStoring block offsets into disk...");

        try (
                RandomAccessFile raf = new RandomAccessFile(BLOCKOFFSETS_FILE, "rw");
                FileChannel channel = raf.getChannel();
        ) {
            MappedByteBuffer buffer = channel.map(FileChannel.MapMode.READ_WRITE, 0, (long) LONG_BYTES * dictionaryBlockOffsets.size()); //offset_size (size of dictionary offset) * number of blocks

            if(buffer == null)      // Buffer not created
                return;

            // scan all block and for each one write offset into disk
            for (int i = 0; i < dictionaryBlockOffsets.size(); i++)
            {
                printDebug("OFFSET BLOCK " + i + ": " + dictionaryBlockOffsets.get(i));
                buffer.putLong(dictionaryBlockOffsets.get(i)); //store into file the dictionary offset of the i-th block
            }

            printDebug(dictionaryBlockOffsets.size() + " blocks stored");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Function to store Dictionary and Inverted Index into disk.
     */
    public static void storeIndexAndDictionaryIntoDisk()
    {
        try (
                RandomAccessFile docidFile = new RandomAccessFile(PARTIAL_DOCID_FILE, "rw");
                RandomAccessFile termfreqFile = new RandomAccessFile(PARTIAL_TERMFREQ_FILE, "rw");
                RandomAccessFile dictFile = new RandomAccessFile(PARTIAL_DICTIONARY_FILE, "rw");
                FileChannel docidChannel = docidFile.getChannel();
                FileChannel termfreqChannel = termfreqFile.getChannel();
                FileChannel dictChannel = dictFile.getChannel()
        ) {
            dictionary.sort();              // Sort the dictionary lexicographically
            dictionaryBlockOffsets.add(PARTIAL_DICTIONARY_OFFSET);// update of the offset of the block for the dictionary file

            // iterate through all the terms of the dictionary ordered
            for (String term : dictionary.getTermToTermStat().keySet())
            {
                ArrayList<Posting> posList = invertedIndex.get(term);   // get posting list of the term

                DictionaryElem dictElem = dictionary.getTermStat(term); // create dictionary element for the term
                dictElem.setOffsetTermFreq(INDEX_OFFSET);               // set termFreq offset
                dictElem.setOffsetDocId(INDEX_OFFSET);                  // set docID offset

                //if(term.equals("0000"))
                //    printDebug("term: 0000 " + dictionary.getTermToTermStat().get("0000") +  " block " + (dictionaryBlockOffsets.size()-1) + " size: " + posList.size());

                // Create buffers for docid and termfreq
                MappedByteBuffer buffer_docid = docidChannel.map(FileChannel.MapMode.READ_WRITE, docidChannel.size(), (long) posList.size() * INT_BYTES); // from 0 to number of postings * int dimension
                MappedByteBuffer buffer_termfreq = termfreqChannel.map(FileChannel.MapMode.READ_WRITE, termfreqChannel.size(), (long) posList.size() * INT_BYTES); //from 0 to number of postings * int dimension

                // iterate through all the postings of the posting list
                for (Posting posting : posList)
                {
                    if (buffer_docid == null || buffer_termfreq == null)    // Buffer not created
                        return;

                    buffer_docid.putInt(posting.getDocId());         // write DocID
                    buffer_termfreq.putInt(posting.getTermFreq());   // write TermFrequency

                    //if(term.equals("0000"))
                    //    printDebug("docId " + posting.getDocId() +  " termfreq " + posting.getTermFreq());
                    INDEX_OFFSET += INT_BYTES;
                }

                dictElem.storeDictionaryElemIntoDisk(dictChannel);  // store dictionary entry to disk
            }
            printDebug(dictionary.getTermToTermStat().size() + " terms stored in block " + (dictionaryBlockOffsets.size()-1));

        } catch (IOException ioException) {
            ioException.printStackTrace();
        }
    }

    /**
     * Function to store one posting list of a term into the disk.
     * @param pl                the posting list to be saved
     * @param termfreqChannel   the channel of term frequency file
     * @param docidChannel      the channel of docID file
     */
    public static void storePostingListIntoDisk(ArrayList<Posting> pl, FileChannel termfreqChannel, FileChannel docidChannel)
    {
        int len = pl.size();        // number of postings in the posting list

        // Create buffers for docid and termfreq
        try {
            MappedByteBuffer bufferdocid = docidChannel.map(FileChannel.MapMode.READ_WRITE, docidChannel.size(), (long) len*Integer.BYTES); // from 0 to number of postings * int dimension
            MappedByteBuffer buffertermfreq = termfreqChannel.map(FileChannel.MapMode.READ_WRITE, termfreqChannel.size(), (long) len*Integer.BYTES); //from 0 to number of postings * int dimension

            // scan all posting in the posting list
            for (Posting posting : pl)
            {
                bufferdocid.putInt(posting.getDocId());         // store the docID
                buffertermfreq.putInt(posting.getTermFreq());   // store the term frequency
                if(debug)
                {
                    docId_logger.logInfo(String.valueOf(posting.getDocId()));
                    termFreq_logger.logInfo(String.valueOf(posting.getTermFreq()));
                }
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    // -------- end: functions to store into disk --------

    // -------- start: functions to read from disk --------

    /**
     * Function to read all document table from disk and put it in memory (HashMap documentTable).
     *
     * @param indexBuilding if true indexing is in progress | if false indexing isn't in progress
     * @throws IOException
     */
    public static void readDocumentTableFromDisk(boolean indexBuilding) throws IOException
    {
        printLoad("Loading document table from disk...");

        try (
             RandomAccessFile docTableRaf = new RandomAccessFile(DOCTABLE_FILE, "r");
             FileChannel channel = docTableRaf.getChannel()
        ) {
            DocumentElement de = new DocumentElement();

            // for to read all DocumentElement stored into disk
            for (int i = 0; i < channel.size(); i += DOCELEM_SIZE)
            {
                de.readDocumentElementFromDisk(i, channel); // get the ith DocElem
                if(indexBuilding)
                    PartialIndexBuilder.documentTable.put(de.getDocid(), new DocumentElement(de.getDocno(), de.getDocid(), de.getDoclength()));
                else
                    QueryProcessor.documentTable.put(de.getDocid(), new DocumentElement(de.getDocno(), de.getDocid(), de.getDoclength()));
            }
        }
    }

    /**
     * Function to read offset of the block from disk.
     */
    public static void readBlockOffsetsFromDisk()
    {
        printLoad("\nLoading block offsets from disk...");

        if(!dictionaryBlockOffsets.isEmpty()) //control check
            dictionaryBlockOffsets.clear();

        try (
                RandomAccessFile raf = new RandomAccessFile(BLOCKOFFSETS_FILE, "rw");
                FileChannel channel = raf.getChannel()
        ) {
            MappedByteBuffer buffer = channel.map(FileChannel.MapMode.READ_WRITE, 0 , channel.size());

            if(buffer == null)      // Buffer not created
                return;

            // iterate through all files for #blocks times
            for(int i = 0; i < channel.size()/ LONG_BYTES; i++)
            {
                dictionaryBlockOffsets.add(buffer.getLong());
                buffer.position((i+1)*LONG_BYTES); //skip to position of the data of the next block to read
                printDebug("OFFSET BLOCK " + i + ": " + dictionaryBlockOffsets.get(i));
            }

            printDebug(dictionaryBlockOffsets.size() + " blocks loaded");

        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    /**
     * function to read and return a posting list from disk
     *
     * @param offsetDocId       offset of the DocID
     * @param offsetTermFreq    offset of the Term Frequency
     * @param posting_size      size of the posting list
     * @param docidChannel      file where read DocID values
     * @param termfreqChannel   file where read Term Frequency values
     * @return the posting lists read from disk
     */
    public static ArrayList<Posting> readPostingListFromDisk(long offsetDocId, long offsetTermFreq, int posting_size, FileChannel docidChannel, FileChannel termfreqChannel) {

        ArrayList<Posting> pl = new ArrayList<>();

        try {
            MappedByteBuffer docidBuffer = docidChannel.map(FileChannel.MapMode.READ_ONLY, offsetDocId, (long) posting_size * Integer.BYTES);
            MappedByteBuffer termfreqBuffer = termfreqChannel.map(FileChannel.MapMode.READ_ONLY, offsetTermFreq, (long) posting_size * Integer.BYTES);

            //while nr of postings read are less than the number of postings to read (all postings of the term)
            for (int i = 0; i < posting_size; i++)
            {
                int docid = docidBuffer.getInt();           // read the DocID
                int termfreq = termfreqBuffer.getInt();     // read the TermFrequency
                pl.add(new Posting(docid, termfreq)); // add the posting to the posting list
            }
            return pl;
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }

    // -------- end: functions to read from disk --------

    /**
     * function to store posting list after compression into disk
     *
     * @param pl                posting list to store
     * @param docidChannel      file where store DocID values
     * @param termfreqChannel   file where store Term Frequency values
     * @return Term Frequency and DocID compressed length
     */
    public static int[] storeCompressedPostingIntoDisk(ArrayList<Posting> pl, FileChannel termfreqChannel, FileChannel docidChannel){

        ArrayList<Integer> tf = new ArrayList<>();
        ArrayList<Integer> docid  = new ArrayList<>();
        int[] length = new int[2];
        //number of postings in the posting list
        for(Posting ps : pl) {
            tf.add(ps.getTermFreq());
            docid.add(ps.getDocId());
        }

        byte[] compressedTf = Unary.integersCompression(tf);
        byte[] compressedDocId = VariableBytes.integersCompression(docid);
        // Create buffers for docid and termfreq
        try {
            MappedByteBuffer buffertermfreq = termfreqChannel.map(FileChannel.MapMode.READ_WRITE, termfreqChannel.size(), compressedTf.length); //number of bytes of compressed tfs
            MappedByteBuffer bufferdocid = docidChannel.map(FileChannel.MapMode.READ_WRITE, docidChannel.size(), compressedDocId.length); //number of bytes of compressed docids

            buffertermfreq.put(compressedTf);
            bufferdocid.put(compressedDocId);

            length[0] = compressedTf.length;
            length[1] = compressedDocId.length;
            return length;

        } catch (IOException e) {
            e.printStackTrace();
        }

        return null;
    }

    /**
     * function to store posting list after compression into disk
     *
     * @param offsetDocId       offset from where to start the read of the DocID values
     * @param offsetTermFreq    offset from where to start the read of the Term Frequency values
     * @param termFreqSize      size of the compressed Term Frequency values
     * @param docIdSize         size of the compressed DocID values
     * @param posting_size      posting list size
     * @param docidChannel      file where read DocID values
     * @param termfreqChannel   file where read Term Frequency values
     * @return termfreq and docid compressed length
     */
    public static ArrayList<Posting> readCompressedPostingListFromDisk(long offsetDocId, long offsetTermFreq, int termFreqSize, int docIdSize, int posting_size, FileChannel docidChannel, FileChannel termfreqChannel) {

        ArrayList<Posting> uncompressed = new ArrayList<>();
        byte[] docids = new byte[docIdSize];
        byte[] tf = new byte[termFreqSize];

        try {
            MappedByteBuffer docidBuffer = docidChannel.map(FileChannel.MapMode.READ_ONLY, offsetDocId, docIdSize);
            MappedByteBuffer termfreqBuffer = termfreqChannel.map(FileChannel.MapMode.READ_ONLY, offsetTermFreq, termFreqSize);

            termfreqBuffer.get(tf, 0, termFreqSize);
            docidBuffer.get(docids, 0, docIdSize );

            ArrayList<Integer> uncompressedTf = Unary.integersDecompression(tf, posting_size);
            ArrayList<Integer> uncompressedDocid = VariableBytes.integersDecompression(docids);
            for(int i = 0; i < posting_size; i++)
            {
                //System.out.println("docid: " + uncompressedDocid.get(i)  + " tf: " + uncompressedTf.get(i));
                uncompressed.add(new Posting(uncompressedDocid.get(i), uncompressedTf.get(i))); // add the posting to the posting list
            }
            return uncompressed;
        } catch (IOException e) {
            e.printStackTrace();
        }

        return null;
    }

}
