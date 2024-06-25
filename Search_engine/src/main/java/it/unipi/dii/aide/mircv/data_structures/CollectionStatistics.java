package it.unipi.dii.aide.mircv.data_structures;

import it.unipi.dii.aide.mircv.QueryProcessor;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.*;

import static it.unipi.dii.aide.mircv.utils.Constants.*;
import static it.unipi.dii.aide.mircv.utils.Logger.collStats_logger;

/**
 * class to contain the statistics of the collection
 */
public final class CollectionStatistics
{
    private static HashMap<Integer, Long> termFreqTable = new HashMap<>();   // hash table TermFreqValue to related occurrence
    final static int mostFreqPos = 5;   // indicates how many term Freq gets
    public final static int COLLSTATS_SIZE = INT_BYTES * 7 + DOUBLE_BYTES * 5 + LONG_BYTES * mostFreqPos + DOUBLE_BYTES * mostFreqPos + INT_BYTES * mostFreqPos;   // Size in bytes
    private static int nDocs;           // number of documents in the collection
    private static double totDocLen;    // sum of the all document length in the collection

    // variable to BM25
    private static double avgDocLen;    // average document length (used in BM25 scoring function)
    private static double k = 1.2;      // typical values between 1,2 and 2
    private static double b = 0.75;     // typical values around 0.75
    private static int emptyDocs = 0;   // number of empty docs in the collection
    private static int minLenDoc = 0;   // len of the shortest doc in the collection
    private static int maxLenDoc = 0;   // len of the longest doc in the collection
    private static int maxTermFreq = 0; // max termFreq in the collection
    private static int maxPLLength = 0;     //
    private static int minPLLength = 0;     //
    private static double avgPLLength = 0;  //
    private static long[] mostTFOcc = new long[mostFreqPos];         //
    private static double[] mostTFPerc = new double[mostFreqPos];    //
    private static int[] mostTF = new int[mostFreqPos];              //


    private CollectionStatistics() {
        throw new UnsupportedOperationException();
    }

    // ---- start -- set and get methods ----
    public static int getNDocs() {
        return nDocs;
    }

    public static void setNDocs(int nDocs) {
        CollectionStatistics.nDocs = nDocs;
    }

    public static double getTotDocLen() {
        return totDocLen;
    }

    public static void setTotDocLen(double totDocLen) {
        CollectionStatistics.totDocLen = totDocLen;
    }

    public static double getAvgDocLen() {
        return avgDocLen;
    }

    public static void setAvgDocLen(double avgDocLen) {
        CollectionStatistics.avgDocLen = avgDocLen;
    }

    public static double getK() {
        return k;
    }

    public static void setK(double k) {
        CollectionStatistics.k = k;
    }

    public static double getB() {
        return b;
    }

    public static void setB(double b) {
        CollectionStatistics.b = b;
    }

    public static int getEmptyDocs() {
        return emptyDocs;
    }

    public static void setEmptyDocs(int emptyDocs) {
        CollectionStatistics.emptyDocs = emptyDocs;
    }

    public static int getMinLenDoc() {
        return minLenDoc;
    }

    public static void setMinLenDoc(int minlenDoc) {
        CollectionStatistics.minLenDoc = minlenDoc;
    }

    public static int getMaxLenDoc() {
        return maxLenDoc;
    }

    public static void setMaxLenDoc(int maxLenDoc) {
        CollectionStatistics.maxLenDoc = maxLenDoc;
    }

    public static int getMaxTermFreq() {
        return maxTermFreq;
    }

    public static void setMaxTermFreq(int maxTermFreq) {
        CollectionStatistics.maxTermFreq = maxTermFreq;
    }

    public static int getMaxPLLength() {
        return maxPLLength;
    }

    public static void setMaxPLLength(int maxPLLength) {
        CollectionStatistics.maxPLLength = maxPLLength;
    }

    public static int getMinPLLength() {
        return minPLLength;
    }

    public static void setMinPLLength(int minPLLength) {
        CollectionStatistics.minPLLength = minPLLength;
    }

    public static double getAvgPLLength() {
        return avgPLLength;
    }

    public static void setAvgPLLength(double avgPLLength) {
        CollectionStatistics.avgPLLength = avgPLLength;
    }

    // ---- end -- set and get methods ----

    // -- start -- function to manage the hash table for the termfreq occurrence
    public static HashMap<Integer, Long> getTermFreqTable() { return termFreqTable; }

    /**
     * Function to add termfreq and their occurrence in the collection.
     * If the term freq is not present in the table is the first occurrence and value is set to 1. If the term freq is
     * already present in the table the value of occurrence is incremented by 1.
     *
     * @param termFreq  is the termFreq value
     */
    public static void addTFOccToTermFreqTable(int termFreq)
    {
        long currTFOcc = 0;  // indicate the current Term Freq occurrence related to the TermFreq value passed as parameter

        if (termFreqTable.containsKey(termFreq))    // the term is already present in the hashtable
        {
            currTFOcc = termFreqTable.get(termFreq);    // get current occurrence value
            termFreqTable.put(termFreq, (currTFOcc+1)); // update occurrence value
        }
        else            // term is not present in the hashtable
            termFreqTable.put(termFreq, 1L);            // insert occurrence value
    }

    /**
     * Function to compute the statistics related to the term Freq occurrence in the collection.
     */
    public static void computeTermFreqOccStatistics ()
    {
        PriorityQueue<CollectionStatistics.TermFreqOccBlock> resPQ = new PriorityQueue<>(new CollectionStatistics.CompareTFOTerm());
        ArrayList<Long> ordTFOccList = new ArrayList<>();      // contain scores of all docs
        ArrayList<Integer> ordTFList = new ArrayList<>();      // contain scores of all docs
        TermFreqOccBlock currBlock;
        long totFTOcc = 0;          // var to indicates the total number of term freq occurrences in the collection

        if (termFreqTable.isEmpty())
            return;

        // from hash map to ordered array list
        for (Map.Entry<Integer, Long> entry : termFreqTable.entrySet())
        {
            resPQ.add(new CollectionStatistics.TermFreqOccBlock(entry.getKey(), entry.getValue()));     // add to priority queue
            totFTOcc += entry.getValue();       // update total value
        }

        // get the five most frequent values
        for (int i=0; i < mostFreqPos; i++)
        {
            if (resPQ.isEmpty())    // the priority queue is empty
            {
                mostTFOcc[i] = 0;
                mostTFPerc[i] = 0;
                mostTF[i] = 0;
            }
            else
            {
                currBlock = resPQ.poll();                           // get the head block
                mostTFOcc[i] = currBlock.getTermFreqOcc();          // set occurrence
                mostTF[i] = currBlock.getTermFreqValue();           // set termFreqValue
                mostTFPerc[i] = (double) (mostTFOcc[i] * 100) / totFTOcc ;      // compute the percentage
            }
        }
    }

    // -- end -- function to manage the hash table for the termfreq occurrence
    /**
     * Function that check if there is the 'collectionStatistics.txt' file in "/resources" folder
     *
     * @return  true -> there is
     *          false -> there isn't
     */
    public static boolean isThereStatsFile()
    {
        File docStats = new File(STATS_FILE);        // define file
        return docStats.exists();
    }

    // function to read the collection statistics from disk
    public static void readCollectionStatsFromDisk()
    {
        printLoad("Loading collection statistics from disk...");

        try (
                RandomAccessFile statsRAF = new RandomAccessFile(new File(STATS_FILE), "rw")
        ) {
            ByteBuffer statsBuffer = ByteBuffer.allocate(COLLSTATS_SIZE);   // bytes to read from disk
            statsRAF.getChannel().position(0);

            statsRAF.getChannel().read(statsBuffer);            // Read flag values from file
            statsBuffer.rewind();                               // Move to the beginning of file for reading

            // Get collection statistic values from buffer
            int nDocs = statsBuffer.getInt();               // read number of documents in the collection
            double totDocLen = statsBuffer.getDouble();     // read sum of the all document length in the collection
            double avgDocLen = statsBuffer.getDouble();     // read average document length
            double k = statsBuffer.getDouble();             // read k parameter
            double b = statsBuffer.getDouble();             // read b parameter
            int emptyDocs = statsBuffer.getInt();           // read number of empty docs in the collection
            int minLenDoc = statsBuffer.getInt();           // read len of the shortest doc in the collection
            int maxLenDoc = statsBuffer.getInt();           // read len of the longest doc in the collection
            int maxTermFreq = statsBuffer.getInt();         // read max termFreq in the collection
            int maxPLLength = statsBuffer.getInt();         // read the max posting list len in the collection
            int minPLLength = statsBuffer.getInt();         // read the min posting list len in the collection
            double avgPLLength = statsBuffer.getDouble();   // read the avg posting list len in the collection
            // read the array values
            for (int i=0; i < mostFreqPos; i++)
                mostTFOcc[i] = statsBuffer.getLong();       // read TermFreqOcc (raw value)
            for (int i=0; i < mostFreqPos; i++)
                mostTFPerc[i] = statsBuffer.getDouble();    // read TermFreqOcc (percentage value)
            for (int i=0; i < mostFreqPos; i++)
                mostTF[i] = statsBuffer.getInt();           // read TermFreq value

            // Set collection statistics values with values read
            setNDocs(nDocs);
            setTotDocLen(totDocLen);
            setAvgDocLen(avgDocLen);
            setK(k);
            setB(b);
            setEmptyDocs(emptyDocs);
            setMinLenDoc(minLenDoc);
            setMaxLenDoc(maxLenDoc);
            setMaxTermFreq(maxTermFreq);
            setMaxPLLength(maxPLLength);
            setMinPLLength(minPLLength);
            setAvgPLLength(avgPLLength);

        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

    // function to store the collection statistics into disk
    public static void storeCollectionStatsIntoDisk()
    {
        printLoad("Storing collection statistics into disk...");

        try (
                RandomAccessFile docStats = new RandomAccessFile(STATS_FILE, "rw");
                FileChannel channel = docStats.getChannel()
        ) {
            MappedByteBuffer buffer = channel.map(FileChannel.MapMode.READ_WRITE, 0, COLLSTATS_SIZE); // integer size * number of int to store (1) + double size * number of double to store (1)

            buffer.putInt(nDocs);           // write total number of document in collection
            buffer.putDouble(totDocLen);    // write sum of the all document length in the collection
            buffer.putDouble(avgDocLen);    // write average document length
            buffer.putDouble(k);            // write k parameter
            buffer.putDouble(b);            // write b parameter
            buffer.putInt(emptyDocs);       // write number of empty docs in the collection
            buffer.putInt(minLenDoc);       // write len of the shortest doc in the collection
            buffer.putInt(maxLenDoc);       // write len of the longest doc in the collection
            buffer.putInt(maxTermFreq);     // write max termFreq in the collection
            buffer.putInt(maxPLLength);     // write max posting list len in the collection
            buffer.putInt(minPLLength);     // write min posting list len in the collection
            buffer.putDouble(avgPLLength);  // write avg posting list len in the collection
            // write the array values
            for (int i=0; i < mostFreqPos; i++)
                buffer.putLong(mostTFOcc[i]);       // write TermFreqOcc (raw value)
            for (int i=0; i < mostFreqPos; i++)
                buffer.putDouble(mostTFPerc[i]);    // write TermFreqOcc (percentage value)
            for (int i=0; i < mostFreqPos; i++)
                buffer.putInt(mostTF[i]);           // write TermFreq value

        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

    public static void printCollectionStatistics()
    {
        printDebug("The values of the collection statistics are:");
        printDebug("- number of document in the collection: " + nDocs);
        printDebug("- Documents len parameter:");
        printDebug("-- sum of the length of all document in the collection: " + totDocLen);
        printDebug("-- average length of document in the collection: " + avgDocLen);
        printDebug("-- number of empty document in the collection: " + emptyDocs);
        printDebug("-- len of the shortest doc in the collection: " + minLenDoc);
        printDebug("-- len of the longest doc in the collection: " + maxLenDoc);
        printDebug("- BM25 parameter:");
        printDebug("-- k parameter: " + k);
        printDebug("-- b parameter: " + b);
        printDebug("- Posting list parameter:");
        printDebug("-- max posting list len: " + maxPLLength);
        printDebug("-- min posting list len: " + minPLLength);
        printDebug("-- avg posting list len: " + avgPLLength);
        printDebug("- Term Frequency values parameter:");
        printDebug("-- the max term frequency in the collection: " + maxTermFreq);
        printDebug("-- Top " + mostFreqPos + " Term Frequency occurrence list (from the most common to the rarest):");
        for (int i=0; i < mostFreqPos; i++)
            printDebug("---- pos " + (i+1) + " -> TermFreqValue: " + mostTF[i] + " , occurrence: " + mostTFOcc[i] + " (row value) , occurrence: " + mostTFPerc[i] + "%");
    }

    /**
     * class to define termUpperBoundPostingList. The priority queue contains instances of termUpperBoundPostingList
     * representing ...
     */
    private static class TermFreqOccBlock
    {
        int termFreqValue;      //
        long termFreqOcc;       //

        // constructor with parameters
        public TermFreqOccBlock(int termFreqValue, long termFreqOcc)
        {
            this.termFreqValue = termFreqValue;
            this.termFreqOcc = termFreqOcc;
        }

        // get methods
        public int getTermFreqValue() {
            return termFreqValue;
        }

        public long getTermFreqOcc() {
            return termFreqOcc;
        }

        @Override
        public String toString() {
            return "PB{" +
                    "Term Freq value = '" + termFreqValue + '\'' +
                    ", term Freq occurrence =" + termFreqOcc +
                    '}';
        }
    }

    /**
     * Class to compare the block, allows the order of the priority queue
     * The order is ascending order according to term upper bound. In case of equal term upper bound values will be made
     * order according to query position of the term. Term upper bound and position are sorted in ascending order.
     */
    private static class CompareTFOTerm implements Comparator<CollectionStatistics.TermFreqOccBlock>
    {
        @Override
        public int compare(CollectionStatistics.TermFreqOccBlock tfb1, CollectionStatistics.TermFreqOccBlock tfb2)
        {
            // comparing terms
            int scoreComparison = -Long.compare(tfb1.getTermFreqOcc(), tfb2.getTermFreqOcc());
            // if the term upper bound are equal, compare by position
            if (scoreComparison == 0)
            {
                // return order by both term upper bound and position (the TUB of the two blocks is equal)
                return Integer.compare(tfb1.getTermFreqValue(), tfb2.getTermFreqValue());
            }

            return scoreComparison;     // return order only by score (the score of the two blocks is different)
        }
    }
}
