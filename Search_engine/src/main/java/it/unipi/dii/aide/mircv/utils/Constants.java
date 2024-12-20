package it.unipi.dii.aide.mircv.utils;

import java.util.ArrayList;

/**
 * Class to define all the constant (name, path, etc...) used in the application.
 */
public final class Constants {

    private Constants() {
        throw new UnsupportedOperationException();
    }

    // -------------------------------- Constants for folder paths -------------------------------------------
    public static final String RES_FOLDER = "src/main/resources/";          // general path for the resources, is the parent folder of the other
    public static final String PARTIAL_FOLDER = RES_FOLDER + "partial/";    // path for the folder containing the partial file generated by SPIMI algorithm
    public static final String MERGED_FOLDER = RES_FOLDER + "merged/";  // path for the folder containing the complete file generated by SPIMI algorithm
    public static final String DEBUG_FOLDER = RES_FOLDER + "debug/";    // path for debug folder
    public static final String UPPERBOUND_FOLDER = RES_FOLDER + "upperBound/";  //  path for the folder containing Upper Bound (used in query optimization)
    public static final String TEST_FOLDER = RES_FOLDER + "test/";      // path for test folder
    // -------------------------------- Constants for file paths -------------------------------------------
    // -- file in RES_FOLDER --
    public final static String COLLECTION_PATH = RES_FOLDER + "collection.tar.gz";
    // file containing a series of test queries and their QID
    public final static String TEST_QUERY_2020_PATH = TEST_FOLDER + "msmarco-test2020-queries.tsv.gz";  // queries of 2020
    public final static String TEST_QUERY_2019_PATH = TEST_FOLDER + "msmarco-test2019-queries.tsv.gz";  // queries of 2019
    public static final String FLAGS_FILE = RES_FOLDER + "flags"; // file in which flags are stored
    public static final String STATS_FILE = RES_FOLDER + "collectionStatistics"; // file in which collection statistics are stored

    // -- file in PARTIAL_FOLDER --
    public static final String PARTIAL_DICTIONARY_FILE = PARTIAL_FOLDER + "partial_dictionary"; // file in which is stored the vocabulary in blocks
    public static final String PARTIAL_DOCID_FILE = PARTIAL_FOLDER + "partial_docId";  // file containing the docId (element of posting list) for each block
    public static final String PARTIAL_TERMFREQ_FILE = PARTIAL_FOLDER + "partial_termFreq";   // file containing the TermFrequency (element of posting list) for each block
    public static final String BLOCKOFFSETS_FILE = PARTIAL_FOLDER + "blocks"; // file containing the offset of each vocabulary block

    // -- file in MERGED_FOLDER --
    public static final String DOCTABLE_FILE = MERGED_FOLDER + "documentTable"; // file in which is stored the document table
    public static final String DICTIONARY_FILE = MERGED_FOLDER + "dictionary"; // file in which is stored the dictionary
    public static final String DOCID_FILE = MERGED_FOLDER + "docId";   // file containing the docId of the InvertedIndex merged
    public static final String TERMFREQ_FILE = MERGED_FOLDER + "termFreq";   // file containing the termFreq of the InvertedIndex merged
    public static final String SKIP_FILE = MERGED_FOLDER + "skipInfo";

    // -- file in UPPERBOUND_FOLDER --
    public static final String TERMUPPERBOUND_FILE = UPPERBOUND_FOLDER + "termsUpperBound"; // file in which terms upper bound are stored
    public static final String TERMFREQWEIGHT_FILE = UPPERBOUND_FOLDER + "termFreqWeight";  // file in which term freq weight are stored

    // -------------------------------- Constants for variable bytes -------------------------------------------

    public static final int INT_BYTES = Integer.BYTES;      // length in Bytes of the integer type
    public static final int LONG_BYTES = Long.BYTES;        // length in Bytes of the long int type
    public static final int DOUBLE_BYTES = Double.BYTES;    // length in Bytes of the double type

    public static final int TERM_DIM = 20;                  // Length of a term (in bytes)
    public static int N_POSTINGS = 0;                       // Number of partial postings to save in the file
    public static final int SKIP_POINTERS_THRESHOLD = 1024; // minimum size of a posting list to activate skipping (if posting list is greater than threshold will be used skipping otherwise skipping will be not used)

    // -------------------------------------- Constants for file offsets ----------------------------------------------

    public static long PARTIAL_DICTIONARY_OFFSET = 0;   // Offset of the terms in the dictionary
    public static long INDEX_OFFSET = 0;                // Offset of the termfreq and docid in index
    public static double MEMORY_THRESHOLD = 0.8;    // Maximum memory usage threshold for a block in the SPIMI algorithm

    // ---------------------------------------- Utilities for debugging -----------------------------------------------

    // variable that stipulates the behaviour for control printouts. If false there will be no printouts, if true there will be all printouts.
    public static final boolean verbose = true;
    public static final boolean debug = false;

    // ---------------------------------------- Utilities for time printing -------------------------------------------

    /**
     * Function to format the execution time in the specified format.
     *
     * @param start execution start time
     * @param end   execution end time
     * @return      a formatted string that show the execution time.
     */
    public static String formatTime(long start, long end)
    {
        long elapsedTime = end - start;
        long seconds = (elapsedTime / 1000) % 60;
        long minutes = (elapsedTime / 1000 / 60) % 60;
        long hours = (elapsedTime / 1000 / 3600);

        return String.format("%02d:%02d:%02d", hours, minutes, seconds);
    }

    /**
     * Function to format the execution time in the specified format.
     *
     * @return      a formatted string that show the execution time.
     */
    public static String formatTime(long time)
    {
        long seconds = (time / 1000) % 60;
        long minutes = (time / 1000 / 60) % 60;
        long hours = (time / 1000 / 3600);

        return String.format("%02d:%02d:%02d", hours, minutes, seconds);
    }

    /**
     * Function to show the size of a file in the specified format.
     *
     * @param nameFile  file name
     * @param fileSize  size (in Bytes) of the file
     * @return          a formatted string that show the file size.
     */
    public static String formatSize(String nameFile, long fileSize)
    {
        double sizeKiloBytes = (double) fileSize / 1024;
        double sizeMegaBytes = (double) fileSize / (1024 * 1024);
        double sizeGigaBytes = (double) fileSize / (1024 * 1024 * 1024);


        return String.format("The size of the file '%s' is: %d Bytes , %.3f KB , %.3f MB , %.3f GB", nameFile, fileSize, sizeKiloBytes, sizeMegaBytes, sizeGigaBytes);
    }

    /**
     * Function to show the size of a file in the specified format.
     *
     * @param size  size (in Bytes) of the file
     * @return          a formatted string that show the file size.
     */
    public static String formatBytes( long size)
    {
        double sizeKiloBytes = (double) size / 1024;
        double sizeMegaBytes = (double) size / (1024 * 1024);
        double sizeGigaBytes = (double) size / (1024 * 1024 * 1024);


        return String.format(" %d Bytes , %.3f KB , %.3f MB , %.3f GB", size, sizeKiloBytes, sizeMegaBytes, sizeGigaBytes);
    }

    // ------------------------------ Utilities for control printing (and color print) ---------------------------------
    // -- terminal colors
    public static final String ANSI_RESET = "\u001B[0m";        // reset the colour of the print
    public static final String ANSI_CYAN = "\u001B[96m";        // UI print
    public static final String ANSI_YELLOW = "\u001B[93m";      // time print
    public static final String ANSI_RED = "\033[0;31m";         // error print
    public static final String ANSI_GREEN = "\u001B[32m";       // load print
    public static final String ANSI_ORANGE = "\u001b[38;5;208m";// file size print
    public static final String ANSI_MAGENTA = "\u001b[35m";     // debug print

    // -- types print
    public static void printDebug(String s){
        if(verbose)
            System.out.println(ANSI_MAGENTA + s + ANSI_RESET);
    }

    public static void printUIMag(String s){
        System.out.println(ANSI_MAGENTA + s + ANSI_RESET);
    }

    public static void printError(String s){
        System.out.println(ANSI_RED + s + ANSI_RESET);
    }

    public static void printUI(String s){
        System.out.println(ANSI_CYAN + s + ANSI_RESET);
    }

    public static void printTime(String s){
        System.out.println(ANSI_YELLOW + s + ANSI_RESET);
    }

    public static void printLoad(String s){
        System.out.println(ANSI_GREEN + s + ANSI_RESET);
    }

    public static void printSize(String s){
        System.out.println(ANSI_ORANGE + s + ANSI_RESET);
    }

    // ------------------------------ Utilities for array list comparison ---------------------------------

    /**
     * Function to compare two array lists of integers.
     *
     * @param al0   the first integer
     * @param al1   the second integer
     * @return      if 'true' the array lists are the same, if 'false' the array lists aren't the same
     */
    public static boolean isTwoIntArrayListEqual(ArrayList<Integer> al0, ArrayList<Integer> al1)
    {
        if (al0.isEmpty() || al1.isEmpty())     // check if are empty
            return false;

        if (al0.size() != al1.size())           // check if are the same size
            return false;

        // check that all element in both array list are the same and in the same order
        for (int i=0; i < al0.size(); i++)
        {
            if (!al0.get(i).equals(al1.get(i)))
                return false;
        }

        return true;
    }

    /**
     * Function to calculate the difference between two integers either as a raw value or as a percentage.
     *
     * @param num0  the first integer
     * @param num1  the second integer
     * @return  a string containing both the raw and percentage difference.
     */
    public static String differenceBetweenTwoInt(int num0, int num1)
    {
        int difference = 0;
        double percentageDiff = 0;
        int max;
        int min;

        if (num0 >= num1)
        {
            max = num0;
            min = num1;
        }
        else
        {
            max = num1;
            min = num0;
        }

        difference = max - min;
        percentageDiff = (double) (difference * 100) /max;

        return (String.format(" difference is: %d Bytes (in percentage) is: %.3f", difference, percentageDiff) +"%");
    }
}
