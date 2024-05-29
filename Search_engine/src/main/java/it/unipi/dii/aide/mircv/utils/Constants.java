package it.unipi.dii.aide.mircv.utils;

public final class Constants {

    private Constants() {
        throw new UnsupportedOperationException();
    }

    public static final String RES_FOLDER = "src/main/resources/";
    public static final String PARTIAL_FOLDER = RES_FOLDER + "partial/";
    public static final String MERGED_FOLDER = RES_FOLDER + "merged/";
    public static final String DEBUG_FOLDER = RES_FOLDER + "debug/";
    public static final String UPPERBOUND_FOLDER = RES_FOLDER + "upperBound/";



    // -------------------------------- Constants for file paths -------------------------------------------

    public final static String COLLECTION_PATH = RES_FOLDER + "collection.tar.gz";
    //public final static String QUERIES_COLLECTION_PATH = RES_FOLDER + "msmarco-test2020-queries.tsv.gz";  // file containing a series of queries and their best DocID result
    public final static String QUERIES_COLLECTION_PATH = RES_FOLDER + "test_queries.txt";  // file containing a series of queries for testing

    public static final String PARTIAL_DICTIONARY_FILE = PARTIAL_FOLDER + "partial_dictionary"; // file in which is stored the vocabulary in blocks
    public static final String PARTIAL_DOCID_FILE = PARTIAL_FOLDER + "partial_docId";  // file containing the docId (element of posting list) for each block
    public static final String PARTIAL_TERMFREQ_FILE = PARTIAL_FOLDER + "partial_termFreq";   // file containing the TermFrequency (element of posting list) for each block

    public static final String DOCTABLE_FILE = MERGED_FOLDER + "documentTable"; // file in which is stored the document table
    public static final String DICTIONARY_FILE = MERGED_FOLDER + "dictionary"; // file in which is stored the dictionary
    public static final String DOCID_FILE = MERGED_FOLDER + "docId";   // file containing the docId of the InvertedIndex merged
    public static final String TERMFREQ_FILE = MERGED_FOLDER + "termFreq";   // file containing the termFreq of the InvertedIndex merged

    public static final String BLOCKOFFSETS_FILE = PARTIAL_FOLDER + "blocks"; // file containing the offset of each vocabulary block
    public static final String FLAGS_FILE = RES_FOLDER + "flags"; // file in which flags are stored
    public static final String STATS_FILE = RES_FOLDER + "collectionStatistics"; // file in which collection statistics are stored

    public static final String SKIP_FILE = MERGED_FOLDER + "skipInfo";

    public static final String TERMUPPERBOUND_FILE = UPPERBOUND_FOLDER + "termsUpperBound";  // file in which terms upper bound are stored

    // -------------------------------- Constants for variable bytes -------------------------------------------

    public static final int INT_BYTES = Integer.BYTES;
    public static final int LONG_BYTES = Long.BYTES;
    public static final int DOUBLE_BYTES = Double.BYTES;

    public static final int TERM_DIM = 20;                      // Length of a term (in bytes)
    public static int N_POSTINGS = 0;                  // Number of partial postings to save in the file

    public static final int SKIP_POINTERS_THRESHOLD = 1024;

    // -------------------------------------- Constants for file offsets ----------------------------------------------

    public static long PARTIAL_DICTIONARY_OFFSET = 0;          // Offset of the terms in the dictionary
    public static long INDEX_OFFSET = 0;               // Offset of the termfreq and docid in index

    public static double MEMORY_THRESHOLD = 0.8;

    // ---------------------------------------- Utilities for debugging -----------------------------------------------

    // variable that stipulates the behaviour for control printouts. If false there will be no printouts, if true there will be all printouts.
    public static final boolean verbose = true;
    public static final boolean debug = false;

    public static String formatTime(long start, long end) {

        long elapsedTime = end - start;
        long seconds = (elapsedTime / 1000) % 60;
        long minutes = (elapsedTime / 1000 / 60) % 60;
        long hours = (elapsedTime / 1000 / 3600);

        return String.format("%02d:%02d:%02d", hours, minutes, seconds);
    }

    // terminal colors
    public static final String ANSI_RESET = "\u001B[0m";        // reset the colour of the print
    public static final String ANSI_CYAN = "\u001B[96m";        // UI print
    public static final String ANSI_YELLOW = "\u001B[93m";      // time print
    public static final String ANSI_RED = "\033[0;31m";         // error print

    public static final String ANSI_GREEN = "\u001B[32m";       // load print
    public static final String ANSI_MAGENTA = "\u001b[35m";     // debug print

    public static void printDebug(String s){
        if(verbose)
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

}
