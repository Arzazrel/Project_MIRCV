package it.unipi.dii.aide.mircv.data_structures;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;

import static it.unipi.dii.aide.mircv.utils.Constants.*;
import static it.unipi.dii.aide.mircv.utils.Logger.collStats_logger;

/**
 * class to contain the statistics of the collection
 */
public final class CollectionStatistics
{
    public final static int COLLSTATS_SIZE = INT_BYTES * 5 + DOUBLE_BYTES * 4;   // Size in bytes
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


    private CollectionStatistics() {
        throw new UnsupportedOperationException();
    }


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

            if(debug)
                collStats_logger.logInfo("nDocs: " + getNDocs() + "\ntotDocLen: " + getTotDocLen());

        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

    public static void printCollectionStatistics()
    {
        printDebug("The values of the collection statistics are:");
        printDebug("- number of document in the collection: " + nDocs);
        printDebug("- sum of the length of all document in the collection: " + totDocLen);
        printDebug("- average length of document in the collection: " + avgDocLen);
        printDebug("- number of empty document in the collection: " + emptyDocs);
        printDebug("- len of the shortest doc in the collection: " + minLenDoc);
        printDebug("- len of the longest doc in the collection: " + maxLenDoc);
        printDebug("- the max term frequency in the collection: " + maxTermFreq);
        printDebug("- BM25 parameter:");
        printDebug("-- k parameter: " + k);
        printDebug("-- b parameter: " + b);
    }
}
