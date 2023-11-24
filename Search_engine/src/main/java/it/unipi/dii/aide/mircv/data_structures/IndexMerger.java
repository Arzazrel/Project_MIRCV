package it.unipi.dii.aide.mircv.data_structures;

import it.unipi.dii.aide.mircv.utils.FileSystem;
import it.unipi.dii.aide.mircv.utils.Logger;

import java.io.*;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.StandardCharsets;
import java.util.*;

import static it.unipi.dii.aide.mircv.data_structures.DataStructureHandler.*;
import static it.unipi.dii.aide.mircv.data_structures.DictionaryElem.getDictElemSize;
import static it.unipi.dii.aide.mircv.data_structures.PartialIndexBuilder.*;
import static it.unipi.dii.aide.mircv.utils.Constants.*;
import static java.lang.Math.min;

/**
 * class to merge the InverteIndex
 */
public final class IndexMerger {
    // Priority queue which will contain the first term (in lexicographic order) of each block. Used for merge and to
    // take from all the blocks the terms in the right order.
    private static final PriorityQueue<TermBlock> pq = new PriorityQueue<>(dictionaryBlockOffsets.size() == 0 ? 1 : dictionaryBlockOffsets.size(), new CompareTerm());

    private IndexMerger() {
        throw new UnsupportedOperationException();
    }

    static int i = 0;       // counter used only for control prints
    /**
     *  function to merge the block of the inverted index
     */
    public static void mergeBlocks() {
        System.out.println("Merging partial files...");                     // print of the merging start

        int nrBlocks = dictionaryBlockOffsets.size();           // dictionary number
        DataStructureHandler.readBlockOffsetsFromDisk();        // get offsets of dictionary blocks from disk
        MappedByteBuffer buffer;
        // array containing the current read offset for each block
        ArrayList<Long> currentBlockOffset = new ArrayList<>(nrBlocks);
        currentBlockOffset.addAll(dictionaryBlockOffsets);

        // var which indicates the steps of 'i' progression print during merge
        System.out.println("Compression " + Flags.isCompressionEnabled());

        // open file and create channels for reading the partial dictionary and index file and write the complete index and dictionary file
        try (
                // open partial files to read the partial dictionary and index
                RandomAccessFile partialDocidFile = new RandomAccessFile(PARTIAL_DOCID_FILE, "rw");
                RandomAccessFile partialTermfreqFile = new RandomAccessFile(PARTIAL_TERMFREQ_FILE, "rw");
                RandomAccessFile partialDictFile = new RandomAccessFile(PARTIAL_DICTIONARY_FILE, "rw");
                // open complete files to write the merged dictionary and merged index
                RandomAccessFile docidFile = new RandomAccessFile(DOCID_FILE, "rw");
                RandomAccessFile termfreqFile = new RandomAccessFile(TERMFREQ_FILE, "rw");
                RandomAccessFile dictFile = new RandomAccessFile(DICTIONARY_FILE, "rw");
                // open skipping file
                RandomAccessFile skipFile = new RandomAccessFile(SKIP_FILE, "rw");

                // FileChannel in input (partial file)
                FileChannel dictChannel = partialDictFile.getChannel();
                FileChannel docidChannel = partialDocidFile.getChannel();
                FileChannel termfreqChannel = partialTermfreqFile.getChannel();
                // FileChannel in output (complete file)
                FileChannel outDictionaryChannel = dictFile.getChannel();
                FileChannel outDocIdChannel = docidFile.getChannel();
                FileChannel outTermFreqChannel = termfreqFile.getChannel();
                // FileChannel in output (skipping)
                FileChannel outSkipChannel = skipFile.getChannel()
        ) {
            // scroll all blocks and add the first term of each block to priority queue
            for(int i = 0; i <  nrBlocks; i++) {
                buffer = dictChannel.map(FileChannel.MapMode.READ_ONLY, currentBlockOffset.get(i), TERM_DIM); //map current block in memory
                String term = StandardCharsets.UTF_8.decode(buffer).toString().split("\0")[0];  // get first term of the block
                pq.add(new TermBlock(term, i));     // add to the priority queue a TermBlock element (term + its blocks number)
            }

            // build temp structures
            DictionaryElem tempDE = new DictionaryElem();       // empty temporary DictionaryELem, contains the accumulated data for each term
            ArrayList<Posting> tempPL = new ArrayList<>();      // empty temporary PostingList, contains the accumulated data for each term
            DictionaryElem currentDE = new DictionaryElem();       // current DictionaryElem, contains the data taken from the queue in the current iteration
            ArrayList<Posting> currentPL;   // current PostingList, contains the data taken from the queue in the current iteration

            TermBlock currentTermBlock;     // var that contain the TermBlock extract from pq in the current iteration
            String term = "";   // var that contain the Term of the TermBlock extract from pq in the current iteration
            int block_id = -1;  // var that contain the blockID of the TermBlock extract from pq in the current iteration

            // Merging the posting list -> SEE NOTE 1
            while(!pq.isEmpty()) {

                // get first element from priority queue
                currentTermBlock = pq.poll();               // get lowest term from priority queue
                assert currentTermBlock != null;
                term = currentTermBlock.getTerm();          // get the term
                block_id = currentTermBlock.getBlock();     // get the blockID

                // If condition to verify if there are other elements -> SEE NOTE 2
                if (currentBlockOffset.get(block_id) + getDictElemSize()  < (block_id == (currentBlockOffset.size()-1) ? dictChannel.size() : dictionaryBlockOffsets.get(block_id +1))){
                    buffer = dictChannel.map(FileChannel.MapMode.READ_ONLY, currentBlockOffset.get(block_id) + getDictElemSize(), TERM_DIM); // get first element of the block
                    String[] t = StandardCharsets.UTF_8.decode(buffer).toString().split("\0");      // get the term of element
                    if (!(t.length == 0))           // control check if term is not empty
                        pq.add(new TermBlock(t[0], block_id));  //add to the priority queue a term block element (term + its blocks number)
                }

                // get current elem of dictionary
                currentDE.readDictionaryElemFromDisk(currentBlockOffset.get(block_id), dictChannel);
                // get current posting list
                currentPL = readPostingListFromDisk(currentDE.getOffsetDocId(), currentDE.getOffsetTermFreq(), currentDE.getDf(), docidChannel, termfreqChannel);

                if (tempDE.getTerm().equals("")) {        // first iteration

                    //set temp variables values with value of the element taken in the current iteration
                    tempDE = currentDE;
                    tempDE.setOffsetTermFreq(outTermFreqChannel.size());
                    tempDE.setOffsetDocId(outDocIdChannel.size());
                    tempPL = currentPL;

                } else {                                // is not the first iteration

                    // same term found (respect the previous iteration), temporary structures update
                    if (currentDE.getTerm().equals(tempDE.getTerm())) {

                        // update DictionaryElem
                        tempDE.addCf(currentDE.getCf());        // update Cf
                        tempDE.addDf(currentDE.getDf());        // update Df

                        assert tempPL != null;
                        tempPL.addAll(currentPL);               // add all new postings

//                        if(currentDE.getMaxTf() > tempDE.getMaxTf())       //update maxTf
//                            tempDE.setMaxTf(currentDE.getMaxTf());
                    }
                    else{    // different term found, write to disk the complete data of the previous term

                        Flags.setConsiderSkippingBytes(true);

                        //update DocID and Term Frequency offset ( equal to the end of the files)
                        tempDE.setOffsetTermFreq(outTermFreqChannel.size());
                        tempDE.setOffsetDocId(outDocIdChannel.size());

//                        tempDE.computeIdf();
//                        tempDE.computeMaxTFIDF();


                        assert tempPL != null;

                        int lenPL = tempPL.size();
                        int[] tempCompressedLength = new int[2];
                        if(lenPL >= SKIP_POINTERS_THRESHOLD) {
                            int skipInterval = (int) Math.ceil(Math.sqrt(lenPL));        // one skip every rad(docs)
                            int nSkip = 0;

                            for(int i = 0; i < lenPL; i += skipInterval) {
                                List<Posting> subPL = tempPL.subList(i, min(i + skipInterval, lenPL));
                                ArrayList<Posting> tempSubPL = new ArrayList<>(subPL);
                                if (Flags.isCompressionEnabled()) {
                                    int[] compressedLength = DataStructureHandler.storeCompressedPostingIntoDisk(tempSubPL, outTermFreqChannel, outDocIdChannel);//store index with compression - unary compression for termfreq
                                    assert compressedLength != null;
                                    tempCompressedLength[0] += compressedLength[0];
                                    tempCompressedLength[1] += compressedLength[1];
                                    SkipInfo sp = new SkipInfo(subPL.get(subPL.size()-1).getDocId(), outDocIdChannel.size(), outTermFreqChannel.size());
                                    sp.storeSkipInfoToDisk(outSkipChannel);
                                } else {
                                    storePostingListIntoDisk(tempSubPL, outTermFreqChannel, outDocIdChannel);  // write InvertedIndexElem to disk
                                    SkipInfo sp = new SkipInfo(subPL.get(subPL.size()-1).getDocId(), outDocIdChannel.size(),  outTermFreqChannel.size());
                                    sp.storeSkipInfoToDisk(outSkipChannel);
                                }
                                nSkip++;

                            }
                            tempDE.setSkipArrLen(nSkip);
                            tempDE.setSkipOffset(outSkipChannel.size());
                            if(Flags.isCompressionEnabled()) {
                                tempDE.setTermFreqSize(tempCompressedLength[0]);
                                tempDE.setDocIdSize(tempCompressedLength[1]);
                            }
                        }
                        else {
                            if(Flags.isCompressionEnabled()){
                                int[] compressedLength = DataStructureHandler.storeCompressedPostingIntoDisk(tempPL, outTermFreqChannel, outDocIdChannel);//store index with compression - unary compression for termfreq
                                assert compressedLength != null;
                                tempDE.setTermFreqSize(compressedLength[0]);
                                tempDE.setDocIdSize(compressedLength[1]);
                            }
                            else
                                storePostingListIntoDisk(tempPL, outTermFreqChannel, outDocIdChannel);  // write InvertedIndexElem to disk
                        }
                        tempDE.storeDictionaryElemIntoDisk(outDictionaryChannel);

                        Flags.setConsiderSkippingBytes(false);

                        //set temp variables values
                        tempDE = currentDE;
                        tempPL = currentPL;
                    }
                }
                currentDE = new DictionaryElem();
                i++;

                // update the offset of next element to read from the block read in this iteration
                currentBlockOffset.set(block_id, currentBlockOffset.get(block_id) + getDictElemSize());
            }

            printDebug("Merge ended, total number of iterations (i) is: " + i);
//            delete_tempFiles();                                                                       !!!!!!!!!!!!!!!!

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * class to define TermBlock. The priority queue contains instances of TermBlock
     */
    private static class TermBlock {
        String term;    // string of the term related to TermBlock
        int block;      // reference to the id of the block in which are the data

        // constructor with parameters
        public TermBlock(String term, int block) {
            this.term = term;
            this.block = block;
        }

        public String getTerm() {
            return term;
        }

        public int getBlock() {
            return block;
        }

        @Override
        public String toString() {
            return "TermBlock{" +
                    "term='" + term + '\'' +
                    ", block=" + block +
                    '}';
        }
    }

    /**
     * class to compare the block, allows the order of the priority queue
     */
    private static class CompareTerm implements Comparator<TermBlock> {
        @Override
        public int compare(TermBlock tb1, TermBlock tb2) {
            // comparing terms
            int termComparison = tb1.getTerm().compareTo(tb2.getTerm());
            // if the terms are equal, compare by block number
            if (termComparison == 0) {
                return Integer.compare(tb1.getBlock(), tb2.getBlock());
            }

            return termComparison;
        }
    }
}

/*
 * NOTE
 * 1) While explanation:
 *    As long as the priority queue is not empty, extract the first term (in lexicographic order) and merge it.
 *    After each extraction, the new first term of the block (from which it was previously taken) is taken and
 *    put in the priority queue.
 * 2) If condition explanation:
 *    If there are other elements to be processed in the block identified by block_id (block containing the
 *    term taken from the queue in the current iteration) take the next term and add it to the priority queue.
 *    If condition divided whether the block considered is the last one or not
 *    last block -> check if reading one more element is less than file size
 *    not the last block -> check if reading one more element of that block is less than next block start
 *    if condition is satisfied -> read new element
 */
