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
public final class CollectionStatistics {

    private static int nDocs;          // number of documents in the collection
    private static double totDocLen;   // sum of the all document length in the collection


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
    public static void readCollectionStatsFromDisk() {

        System.out.println("Loading collection statistics from disk...");

        try (
                RandomAccessFile statsRAF = new RandomAccessFile(new File(STATS_FILE), "rw")
        ) {
            ByteBuffer statsBuffer = ByteBuffer.allocate(INT_BYTES + DOUBLE_BYTES);   // bytes to read from disk
            statsRAF.getChannel().position(0);

            statsRAF.getChannel().read(statsBuffer);            // Read flag values from file
            statsBuffer.rewind();                               // Move to the beginning of file for reading

            // Get collection statistic values from buffer
            int nDocs = statsBuffer.getInt();               // read number of documents in the collection
            double totDocLen = statsBuffer.getDouble();     // read sum of the all document length in the collection

            // Set collection statistics values with values read
            setNDocs(nDocs);
            setTotDocLen(totDocLen);

            printDebug("Collection statistics read -> nDocs: " + nDocs + ", totDocLen: " + totDocLen);

        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

    // function to store the collection statistics into disk
    public static void storeCollectionStatsIntoDisk() {
        System.out.println("Storing collection statistics into disk...");

        try (
                RandomAccessFile docStats = new RandomAccessFile(STATS_FILE, "rw");
                FileChannel channel = docStats.getChannel()
        ) {
            MappedByteBuffer buffer = channel.map(FileChannel.MapMode.READ_WRITE, 0, INT_BYTES + DOUBLE_BYTES); // integer size * number of int to store (1) + double size * number of double to store (1)

            buffer.putInt(nDocs);           // write total number of document in collection
            buffer.putDouble(totDocLen);    // write sum of the all document length in the collection

            if(debug) collStats_logger.logInfo("nDocs: " + getNDocs() + "\ntotDocLen: " + getTotDocLen());

        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

}
