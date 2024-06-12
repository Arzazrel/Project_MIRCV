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
public final class Flags {

    private static boolean sws_flag = false;            // true = stop words removal enabled, false = stop words removal disabled
    private static boolean compression_flag = false;    // true = compression enabled, false = compression disabled
    private static boolean scoring_flag = false;        // true = scoring enable, false = scoring disable
    private static boolean skip_flag = false;           // true = skipping enable, false = skipping disable

    /*
    private static boolean isSPIMI = false;
    private static boolean isMerge = false;
    private static boolean isQuery = false;
    */

    // -- start -- get method
    public static boolean isSwsEnabled() { return sws_flag; }

    public static boolean isCompressionEnabled() { return compression_flag; }

    public static boolean isScoringEnabled() { return scoring_flag; }

    public static boolean considerSkippingBytes() { return skip_flag; }

    // -- start -- set method
    public static void setSws(boolean sws_flag) { Flags.sws_flag = sws_flag; }

    public static void setCompression(boolean compression_flag) {
        Flags.compression_flag = compression_flag;
    }

    public static void setScoring(boolean scoring_flag) {
        Flags.scoring_flag = scoring_flag;
    }

    public static void setConsiderSkippingBytes(boolean skip_flag) { Flags.skip_flag = skip_flag; }

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
            MappedByteBuffer buffer = channel.map(FileChannel.MapMode.READ_WRITE, 0, (long) Integer.BYTES * 3); //offset_size (size of dictionary offset) * number of blocks

            buffer.putInt(isSwsEnabled() ? 1 : 0);             // write stop words removal user's choice
            buffer.putInt(isCompressionEnabled() ? 1 : 0);     // write compression user's choice
            buffer.putInt(isScoringEnabled() ? 1 : 0);         // write scoring user's choice

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
            ByteBuffer flagsBuffer = ByteBuffer.allocate(12);
            flagsRaf.getChannel().position(0);

            flagsRaf.getChannel().read(flagsBuffer);            // Read flag values from file
            flagsBuffer.rewind();                               // Move to the beginning of file for reading

            // Get flag values from buffer
            int isSwsEnabled = flagsBuffer.getInt();            // read stop words removal user's choice
            int isCompressionEnabled = flagsBuffer.getInt();    // read compression user's choice
            int isScoringEnabled = flagsBuffer.getInt();        // scoring user's choice

            // Set flag values with values read
            setSws(isSwsEnabled == 1);
            setCompression(isCompressionEnabled == 1);
            setScoring(isScoringEnabled == 1);

            printFlags();   // show to the user the flags values read

        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

    /**
     * Function to show the user's choices for the flags.
     */
    public static void printFlags()
    {
        printUI("The user's choices for the flags are:");                   // control print for the user
        printUI("- is stopwords removal enabled : " + isSwsEnabled());      // print choice for stopwords removal
        printUI("- is compression enabled : " + isCompressionEnabled());    // print choice for compression
        printUI("- is scoring BM25 enabled : " + isScoringEnabled());       // print choice for BM25 scoring function
        printUI("- is skipping enabled : " + considerSkippingBytes());      // print choice for skipping
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

//    public static boolean isSPIMI() {
//        return isSPIMI;
//    }
//
//    public static void setIsSPIMI(boolean isSPIMI) {
//        Flags.isSPIMI = isSPIMI;
//    }
}
