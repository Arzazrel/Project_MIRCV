package it.unipi.dii.aide.mircv.data_structures;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;

import static it.unipi.dii.aide.mircv.utils.Constants.*;

/**
 *  Stores the flag values chosen by the users.
 */
public final class Flags
{
    private static int numberOfFlags = 6;       // indicates the number of flags (to be read or write into or from disk)
    private static boolean sws_flag = false;            // true = stop words removal and stemming enabled, false = stop words removal and stemming disabled
    private static boolean compression_flag = false;    // true = compression enabled, false = compression disabled
    private static boolean scoring_flag = false;        // true = scoring enable, false = scoring disable
    private static boolean skip_flag = false;           // true = skipping enable, false = skipping disable
    private static boolean qdPruning_flag = false;      // true = query executed with dynamic pruning algorithm (Max Score), false = query executed with classic DAAT algorithm
    private static boolean deletePartFile_flag = false; // true = delete the partial file after complete the indexing, false = doesn't delete partial file after complete the indexing
    // -- start -- get method
    public static boolean isSwsEnabled() { return sws_flag; }

    public static boolean isCompressionEnabled() { return compression_flag; }

    public static boolean isScoringEnabled() { return scoring_flag; }

    public static boolean considerSkippingBytes() { return skip_flag; }

    public static boolean isDynamicPruningEnabled() { return qdPruning_flag; }

    public static boolean isDeletePartFileEnabled() { return deletePartFile_flag; }

    // -- start -- set method
    public static void setSws(boolean sws_flag) { Flags.sws_flag = sws_flag; }

    public static void setCompression(boolean compression_flag) {
        Flags.compression_flag = compression_flag;
    }

    public static void setScoring(boolean scoring_flag) {
        Flags.scoring_flag = scoring_flag;
    }

    public static void setConsiderSkippingBytes(boolean skip_flag) { Flags.skip_flag = skip_flag; }

    public static void setDynamicPruning(boolean qdPruning_flag) { Flags.qdPruning_flag = qdPruning_flag;}

    public static void setDeletePartFile(boolean deletePartFile_flag) { Flags.deletePartFile_flag = deletePartFile_flag;}

    // -- start -- functions

    /**
     * Function to store the user's choices for the flags.
     */
    public static void storeFlagsIntoDisk()
    {
        printLoad("Storing flags into disk...");    // control print for the user

        try (
            RandomAccessFile raf = new RandomAccessFile(FLAGS_FILE, "rw");
            FileChannel channel = raf.getChannel()
        ) {
            MappedByteBuffer buffer = channel.map(FileChannel.MapMode.READ_WRITE, 0, (long) Integer.BYTES * numberOfFlags); //offset_size (size of dictionary offset) * number of blocks

            buffer.putInt(isSwsEnabled() ? 1 : 0);             // write stop words removal user's choice
            buffer.putInt(isCompressionEnabled() ? 1 : 0);     // write compression user's choice
            buffer.putInt(isScoringEnabled() ? 1 : 0);         // write scoring user's choice
            buffer.putInt(considerSkippingBytes() ? 1 : 0);    // write skipping user's choice
            buffer.putInt(isDynamicPruningEnabled() ? 1 : 0);  // write dynamic pruning user's choice
            buffer.putInt(isDeletePartFileEnabled() ? 1 : 0);  // write delete partial file user's choice

        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

    /**
     * Function to read the user's choices for the flags
     */
    public static void readFlagsFromDisk()
    {
        printLoad("Loading flags from disk...");   // control print for the user

        try (
            RandomAccessFile flagsRaf = new RandomAccessFile(new File(FLAGS_FILE), "rw"))
        {
            ByteBuffer flagsBuffer = ByteBuffer.allocate(Integer.BYTES * numberOfFlags);
            flagsRaf.getChannel().position(0);
            flagsRaf.getChannel().read(flagsBuffer);            // Read flag values from file
            flagsBuffer.rewind();                               // Move to the beginning of file for reading

            // Get flag values from buffer
            int isSwsEnabled = flagsBuffer.getInt();            // read stop words removal user's choice
            int isCompressionEnabled = flagsBuffer.getInt();    // read compression user's choice
            int isScoringEnabled = flagsBuffer.getInt();        // read scoring user's choice
            int skip_flag = flagsBuffer.getInt();               // read skipping user's choice
            int qdPruning_flag = flagsBuffer.getInt();          // read dynamic pruning user's choice
            int deletePartFile_flag = flagsBuffer.getInt();     // read delete partial file user's choice

            // Set flag values with values read
            setSws(isSwsEnabled == 1);                          // set stop words removal user's choice
            setCompression(isCompressionEnabled == 1);          // set compression user's choice
            setScoring(isScoringEnabled == 1);                  // set scoring user's choice
            setConsiderSkippingBytes(skip_flag == 1);           // set skipping user's choice
            setDynamicPruning(qdPruning_flag == 1);             // set dynamic pruning user's choice
            setDeletePartFile(deletePartFile_flag == 1);        // set delete partial file user's choice
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

    /**
     * Function to show the user's choices for the flags (debug mode).
     * Debug mode -> print only if debug is enabled and show only the value of the flags.
     */
    public static void printFlagsDebug()
    {
        printDebug("- stopwords removal flag : " + isSwsEnabled());         // print choice for stopwords removal
        printDebug("- compression flag : " + isCompressionEnabled());       // print choice for compression
        printDebug("- scoring BM25 flag : " + isScoringEnabled());          // print choice for BM25 scoring function
        printDebug("- skipping flag : " + considerSkippingBytes());         // print choice for skipping
        printDebug("- dynamic pruning flag : " + isDynamicPruningEnabled());// print choice for query with dynamic pruning algorithm
        printDebug("- delete partial file flag : " + isDeletePartFileEnabled());// print choice for delete partial file
    }

    /**
     * Function to show the user's choices for the flags (UI mode).
     * Exhaustively shows to the user the values of the flags and their meaning.
     */
    public static void printFlagsUI()
    {
        printUI("Explanation of the flags:");
        printUI("1) stopwords removal   -> if is 'true': the stopwords will be removed from the queries, the inverted index and dictionary.");
        printUI("2) compression         -> if is 'true': the compression of the inverted index will be enabled.");
        printUI("3) scoring BM25        -> if is 'true': the scoring function used will be BM25.");
        printUI("                       -> if is 'false': the scoring function used will be TFIDF.");
        printUI("4) skipping            -> if is 'true': make the skipping in the inverted index.");
        printUI("5) dynamic pruning     -> if is 'true': queries will be executed with dynamic pruning algorithm (WAND).");
        printUI("                       -> if is 'false': queries will be executed with classic DAAT algorithm.");
        printUI("6) delete partial file -> if is 'true': the partial file will be deleted after complete the inverted index.");
        printUI("                       -> if is 'false': the partial file will not be deleted after complete the inverted index.");

        printUI("\nThe user's choices for the flags are:");                         // control print for the user
        printUI("- is stopwords removal enabled : " + isSwsEnabled());              // print choice for stopwords removal
        printUI("- is compression enabled       : " + isCompressionEnabled());      // print choice for compression
        printUI("- is scoring BM25 enabled      : " + isScoringEnabled());          // print choice for BM25 scoring function
        printUI("- is skipping enabled          : " + considerSkippingBytes());     // print choice for skipping
        printUI("- is dynamic pruning enabled   : " + isDynamicPruningEnabled());   // print choice for dynamic pruning algorithm
        printUI("- delete partial file enabled  : " + isDeletePartFileEnabled());   // print choice for delete partial file
    }

    /**
     * Function that check if there is the 'flags.txt' file in "/resources" folder.
     *
     * @return  true -> there is
     *          false -> there isn't
     */
    public static boolean isThereFlagsFile()
    {
        // define file
        File docFlags = new File(FLAGS_FILE);        // flags.txt

        return docFlags.exists();
    }
}
