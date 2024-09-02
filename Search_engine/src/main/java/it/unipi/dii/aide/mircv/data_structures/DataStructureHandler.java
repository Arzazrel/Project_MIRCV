package it.unipi.dii.aide.mircv.data_structures;

import java.io.*;
import java.nio.CharBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import it.unipi.dii.aide.mircv.compression.Unary;
import it.unipi.dii.aide.mircv.QueryProcessor;
import it.unipi.dii.aide.mircv.compression.VariableBytes;
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
                buffer.putDouble(de.getDenomPartBM25());

                if(debug)
                    docTable_logger.logInfo(de.toString());
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Function to calculate and save (in the docTable) the denominator of the BM25 scoring function, so this
     * calculation is done offline, saving time during query execution.
     */
    public static void calcAndStoreDenPartBM25inDocTable()
    {
        long startTime,endTime;             // variables to calculate the execution time
        double denomPartBM25;
        double k = CollectionStatistics.getK();
        double b = CollectionStatistics.getB();
        double avgDocLen = CollectionStatistics.getAvgDocLen();

        printLoad("\nCalculate and storing the denominator part for BM25 in Document Table into disk...");

        // control check: there is a division for avgDocLen -> it must be != 0. If avgDocLen = 0 -> also all docLen = 0
        //                  -> the whole documents in the collection are empty
        if (avgDocLen == 0)
        {
            printError("Error: the denominator for BM25 cannot be calculated and the same function cannot be used for scoring because the collection consists of documents that are all empty.");
            return;                 // exit to the function
        }

        // if docTable is not in memory load it
        if(QueryProcessor.documentTable.isEmpty())
        {
            try
            {
                readDocumentTableFromDisk(false);
            }
            catch (IOException e) {
                e.printStackTrace();
            }
        }

        // if exist delete file
        File docTable = new File(DOCTABLE_FILE);        // documentTable.txt
        if(docTable.exists())
            docTable.delete();

        // calculate
        startTime = System.currentTimeMillis();         // start time for calculate

        try (RandomAccessFile raf = new RandomAccessFile(DOCTABLE_FILE, "rw");
             FileChannel channel = raf.getChannel())
        {
            MappedByteBuffer buffer = channel.map(FileChannel.MapMode.READ_WRITE, channel.size(), (long) DOCELEM_SIZE * QueryProcessor.documentTable.size());

            if(buffer == null)      // Buffer not created
                return;
            // scan all document elements of the Document Table
            for(DocumentElement de: QueryProcessor.documentTable.values())
            {
                // calculate -> k * ((1 - b) + b * (docLen / avgDocLen))
                denomPartBM25 = k * ((1 - b) + b * (de.getDoclength() / avgDocLen));
                de.setDenomPartBM25(denomPartBM25);
                // store
                CharBuffer charBuffer = CharBuffer.allocate(DOCNO_DIM);     //allocate bytes for docno

                for (int i = 0; i < de.getDocno().length(); i++)    //put every char into charbuffer
                    charBuffer.put(i, de.getDocno().charAt(i));

                // write docno, docid and doclength into document file
                buffer.put(StandardCharsets.UTF_8.encode(charBuffer));
                buffer.putInt(de.getDocid());
                buffer.putInt(de.getDoclength());
                buffer.putDouble(de.getDenomPartBM25());

                if(debug)
                    docTable_logger.logInfo(de.toString());
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        endTime = System.currentTimeMillis();           // end time for calculate
        printTime("\nCalculate and stored the denominator part for BM25 in Document Table in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
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
                    QueryProcessor.documentTable.put(de.getDocid(), new DocumentElement(de.getDocno(), de.getDocid(), de.getDoclength(), de.getDenomPartBM25()));
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
                //printDebug("OFFSET BLOCK " + i + ": " + dictionaryBlockOffsets.get(i));
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
                pl.add(new Posting(docid, termfreq));       // add the posting to the posting list
            }
            return pl;
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }

    // -------- end: functions to read from disk --------

    /**
     * Function to store posting list after compression into disk
     *
     * @param pl                posting list to store
     * @param docIDChannel      file where store DocID values
     * @param termFreqChannel   file where store Term Frequency values
     * @return Term Frequency and DocID compressed length
     */
    public static int[] storeCompressedPostingIntoDisk(ArrayList<Posting> pl, FileChannel termFreqChannel, FileChannel docIDChannel)
    {
        ArrayList<Integer> tf = new ArrayList<>();      // arraylist to contain the term freqs of the PL
        ArrayList<Integer> docid  = new ArrayList<>();  // arraylist to contain the DocID of the PL
        int[] length = new int[2];                      // array to contain the values ot the lens (in Bytes) for the compressed lists
        byte[] compressedTf;                            // array for the compressed term freq
        byte[] compressedDocId;                         // array for the compressed DocID

        // get the term freq and the DocID
        for(Posting ps : pl) {
            tf.add(ps.getTermFreq());       // get term freq
            docid.add(ps.getDocId());       // get DocID
        }

        compressedTf = Unary.integersCompression(tf);                              // compression of term freq (Unary)
        compressedDocId = VariableBytes.integersCompression(docid,true);    // compression of DocID (var bytes)
        // Create buffers for docid and termfreq
        try {
            MappedByteBuffer bufferTermFreq = termFreqChannel.map(FileChannel.MapMode.READ_WRITE, termFreqChannel.size(), compressedTf.length); //number of bytes of compressed tfs
            MappedByteBuffer bufferDocID = docIDChannel.map(FileChannel.MapMode.READ_WRITE, docIDChannel.size(), compressedDocId.length);       //number of bytes of compressed docids

            bufferTermFreq.put(compressedTf);       // write the compressed term freq list
            bufferDocID.put(compressedDocId);       // write the compressed DocID list

            length[0] = compressedTf.length;        // save the bytes of TF list
            length[1] = compressedDocId.length;     // save the bytes of DID list
            //printDebug("Store length -> length[0] (termFreq): " + length[0] + " and length[1] (DID): " + length[1]);
            return length;

        } catch (IOException e) {
            e.printStackTrace();
        }

        return null;
    }

    /**
     * Function to read posting list after compression into disk
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
    public static ArrayList<Posting> readCompressedPostingListFromDisk(long offsetDocId, long offsetTermFreq, int termFreqSize, int docIdSize, int posting_size, FileChannel docidChannel, FileChannel termfreqChannel)
    {
        ArrayList<Posting> uncompressed = new ArrayList<>();    // decompressed posting list
        byte[] docids = new byte[docIdSize];                    // array for the compressed DocID list
        byte[] tf = new byte[termFreqSize];                     // array for the compressed TermFreq list

        try {
            MappedByteBuffer docidBuffer = docidChannel.map(FileChannel.MapMode.READ_ONLY, offsetDocId, docIdSize);
            MappedByteBuffer termfreqBuffer = termfreqChannel.map(FileChannel.MapMode.READ_ONLY, offsetTermFreq, termFreqSize);

            termfreqBuffer.get(tf, 0, termFreqSize);    // read term freq list
            docidBuffer.get(docids, 0, docIdSize );     // read DocID list

            ArrayList<Integer> uncompressedTf = Unary.integersDecompression(tf, posting_size);  // decompress term freq
            ArrayList<Integer> uncompressedDocid = VariableBytes.integersDecompression(docids,true);    // decompress DocID
            for(int i = 0; i < posting_size; i++)
            {
                uncompressed.add(new Posting(uncompressedDocid.get(i), uncompressedTf.get(i))); // add the posting to the posting list
            }
            return uncompressed;
        } catch (IOException e) {
            e.printStackTrace();
        }

        return null;
    }

    /**
     * Function to read and uncompress a whole posting list, compressed and divided in skipping block, stored into disk.
     *
     * @param sl                reference to the SkipList related to the PL
     * @param offsetDocId       offset from where to start the read of the DocID values
     * @param offsetTermFreq    offset from where to start the read of the Term Frequency values
     * @param termFreqSize      size of the compressed Term Frequency values
     * @param docIdSize         size of the compressed DocID values
     * @param skipArrLen        the len of the skip info block array
     * @param posting_size      posting list size
     * @param docidChannel      file where read DocID values
     * @param termfreqChannel   file where read Term Frequency values
     * @return termfreq and docid compressed length
     */
    public static ArrayList<Posting> readAndUncompressCompressedAndSkippedPLFromDisk(SkipList sl, long offsetDocId, long offsetTermFreq, int termFreqSize, int docIdSize,int skipArrLen, int posting_size, FileChannel docidChannel, FileChannel termfreqChannel)
    {
        ArrayList<Posting> uncompressed = new ArrayList<>();    // decompressed posting list
        SkipInfo currSkipInfo;      // the instance of skipInfo related to the term
        byte[] docids;              // array for the compressed DocID list
        byte[] tf;                  // array for the compressed TermFreq list
        int currOffsetDID = 0;      // the current offset (at each iteration) for the compressed DID
        int currOffsetTF = 0;       // the current offset (at each iteration) for the compressed TermFreq
        int currDIDSize = 0;        // the current size (at each iteration) for the compressed DID
        int currTFSize = 0;         // the current size (at each iteration) for the compressed TermFreq

        try {
            MappedByteBuffer docidBuffer = docidChannel.map(FileChannel.MapMode.READ_ONLY, offsetDocId, docIdSize);
            MappedByteBuffer termfreqBuffer = termfreqChannel.map(FileChannel.MapMode.READ_ONLY, offsetTermFreq, termFreqSize);

            if (skipArrLen == 1)    // case of only one skipping block
            {
                docids = new byte[docIdSize];       // initialize the byte array for the compressed DID list
                tf = new byte[termFreqSize];        // initialize the byte array for the compressed TF list

                termfreqBuffer.get(tf, 0, termFreqSize);    // read term freq list
                docidBuffer.get(docids, 0, docIdSize);     // read DocID list

                ArrayList<Integer> uncompressedTf = Unary.integersDecompression(tf, posting_size);  // decompress term freq
                ArrayList<Integer> uncompressedDocid = VariableBytes.integersDecompression(docids,true);    // decompress DocID

                for(int i = 0; i < posting_size; i++)
                {
                    uncompressed.add(new Posting(uncompressedDocid.get(i), uncompressedTf.get(i))); // add the posting to the posting list
                }
            }
            else        // case of more than one skipping block
            {
                // read and uncompress all skipping block
                for (int i = 0; i < skipArrLen; i++)
                {
                    currSkipInfo = sl.getSkipBlockInfo(i);  // get the skip info related to the i-th block
                    // calculate the number of byte for each compressed list
                    currDIDSize = (int) (currSkipInfo.getDocIdOffset() - currOffsetDID - offsetDocId);
                    currTFSize = (int) (currSkipInfo.getFreqOffset() - currOffsetTF - offsetTermFreq);

                    // initialize the bytes array for the compressed list
                    docids = new byte[currDIDSize];     // initialize the byte array for the compressed DID list
                    tf = new byte[currTFSize];          // initialize the byte array for the compressed TF list

                    // read the compressed block of the PL
                    docidBuffer.get(docids, 0, currDIDSize);   // read DocID list
                    termfreqBuffer.get(tf, 0, currTFSize);       // read term freq list

                    // update current offset
                    currOffsetDID += currDIDSize;
                    currOffsetTF += currTFSize;

                    // add the uncompressed PL block to the final uncompressed PL
                    ArrayList<Integer> uncompressedTf;
                    if ((i != (skipArrLen-1)) || (posting_size % SKIP_POINTERS_THRESHOLD == 0))        // last skipping block
                        uncompressedTf = Unary.integersDecompression(tf, SKIP_POINTERS_THRESHOLD);  // decompress term freq
                    else                            // is not the last skipping block
                        uncompressedTf = Unary.integersDecompression(tf, (posting_size % SKIP_POINTERS_THRESHOLD));  // decompress term freq
                    ArrayList<Integer> uncompressedDocid = VariableBytes.integersDecompression(docids,true);    // decompress DocID

                    for(int j = 0; j < uncompressedTf.size(); j++)
                    {
                        uncompressed.add(new Posting(uncompressedDocid.get(j), uncompressedTf.get(j))); // add the posting to the posting list
                    }
                }
            }
            return uncompressed;
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }

    /**
     * Function to read the uncompressed termFrequency of one block of the posting list indicates as parameter.
     *
     * @param sl                SkipList instances related to the term of this posting list
     * @param blockIndex        indicates the index of the PL block to read
     * @param offsetTermFreq    offset to read the first termFrequency
     * @param termFreqSize      the number of term frequency to read from disk (length of skip block or less if is the last block)
     * @param skipArrLen        the length of the SkipInfo array (equal to the number of the skipping block in posting list)
     * @param termfreqChannel   the channel from to read the termFrequency
     * @return      a byte array containing the compressed representation for the term frequency (Unary code)
     */
    public static byte[] readCompTFBlockFromDisk(SkipList sl, int blockIndex, long offsetTermFreq, int termFreqSize, int skipArrLen, FileChannel termfreqChannel)
    {
        SkipInfo currSkipInfo;              // SkipInfo related to current (at each iteration) block
        byte[] tf;                          // array for the compressed TermFreq list
        int startOffsetTF;                  // the current offset (at each iteration) for the compressed TermFreq
        int sizeToReadTF = 0;               // the current size (at each iteration) for the compressed TermFreq

        // control check of boundaries
        if ( (blockIndex < 0) || (blockIndex >= skipArrLen))
            return null;

        // set currOffsetTF
        if (blockIndex == 0)
            startOffsetTF = (int) offsetTermFreq;
        else
            startOffsetTF = (int) sl.getSkipBlockInfo(blockIndex-1).getFreqOffset();

        currSkipInfo = sl.getSkipBlockInfo(blockIndex);  // get the skip info related to the blockIndex-th block

        // case of only one skipping block and want to read the first block
        if (skipArrLen == 1)
            sizeToReadTF = termFreqSize;
        else        // case of more than one skipping block
        {
            // calculate the number of byte for each compressed list
            sizeToReadTF = (int) (currSkipInfo.getFreqOffset() - startOffsetTF);
        }

        try {
            MappedByteBuffer termfreqBuffer = termfreqChannel.map(FileChannel.MapMode.READ_ONLY, startOffsetTF, sizeToReadTF);

            tf = new byte[sizeToReadTF];        // initialize the byte array for the compressed TF list
            termfreqBuffer.get(tf, 0, sizeToReadTF);       // read the TF compressed block of the PL
            return tf;

        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }

    /**
     * Function to read the uncompressed DID of one block of the posting list indicates as parameter.
     *
     * @param sl            SkipList instances related to the term of this posting list
     * @param blockIndex    indicates the index of the PL block to read
     * @param offsetDocId   offset to read the first DID
     * @param docIdSize     the number of DID to read from disk (length of skip block or less if is the last block)
     * @param skipArrLen    the length of the SkipInfo array (equal to the number of the skipping block in posting list)
     * @param docidChannel  the channel from to read the DID
     * @return      a byte array containing the compressed representation for the DID (Variable-Byte code)
     */
    public static byte[] readCompDIDBlockFromDisk(SkipList sl, int blockIndex, long offsetDocId, int docIdSize,int skipArrLen, FileChannel docidChannel)
    {
        SkipInfo currSkipInfo;          // SkipInfo related to current (at each iteration) block
        byte[] docids;                  // array for the compressed TermFreq list
        int startOffsetDID;                  // the current offset (at each iteration) for the compressed TermFreq
        int sizeToReadDID = 0;               // the current size (at each iteration) for the compressed TermFreq

        // control check of boundaries
        if ( (blockIndex < 0) || (blockIndex >= skipArrLen))
            return null;

        // set currOffsetTF
        if (blockIndex == 0)
            startOffsetDID = (int) offsetDocId;
        else
            startOffsetDID = (int) sl.getSkipBlockInfo(blockIndex-1).getDocIdOffset();

        currSkipInfo = sl.getSkipBlockInfo(blockIndex);  // get the skip info related to the blockIndex-th block

        // case of only one skipping block and want to read the first block
        if (skipArrLen == 1)
            sizeToReadDID = docIdSize;       // initialize the byte array for the compressed DID list
        else        // case of more than one skipping block
        {
            // calculate the number of byte for each compressed list
            sizeToReadDID = (int) (currSkipInfo.getDocIdOffset() - startOffsetDID);
        }

        try {
            MappedByteBuffer docidBuffer = docidChannel.map(FileChannel.MapMode.READ_ONLY, startOffsetDID, sizeToReadDID);

            docids = new byte[sizeToReadDID];                 // initialize the byte array for the compressed DID list
            docidBuffer.get(docids, 0, sizeToReadDID);  // read the DID compressed block of the PL
            return docids;

        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }
}
