package it.unipi.dii.aide.mircv;

import it.unipi.dii.aide.mircv.compression.Unary;
import it.unipi.dii.aide.mircv.compression.VariableBytes;
import it.unipi.dii.aide.mircv.data_structures.*;
import it.unipi.dii.aide.mircv.utils.FileSystem;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

import static it.unipi.dii.aide.mircv.QueryProcessor.queryStartControl;
import static it.unipi.dii.aide.mircv.data_structures.CollectionStatistics.readCollectionStatsFromDisk;
import static it.unipi.dii.aide.mircv.data_structures.DataStructureHandler.*;
import static it.unipi.dii.aide.mircv.utils.FileSystem.*;
import static it.unipi.dii.aide.mircv.utils.Constants.*;
import static it.unipi.dii.aide.mircv.data_structures.Flags.*;
import static java.lang.Math.min;

public class Main
{
    public static void main(String[] args) throws IOException
    {
        Scanner sc = new Scanner(System.in);
        long startTime, endTime;                // variables to calculate the execution time

        // while constituting the user interface
        while (true)
        {
            // print of the user interface
            printUI(
                    "\n********** SEARCH ENGINE **********" +
                    "\n\tSelect an option:" +
                    "\n\t- Functionalities:" +
                    "\n\t  i -> build the index" +
                    "\n\t  q -> query mode" +
                    "\n\t  u -> calculate term upper bound" +
                    "\n\t  f -> see or change flags" +
                    "\n\t- Statistics:" +
                    "\n\t  s -> see files size" +
                    "\n\t  c -> see collection statistics" +
                    "\n\t- Test:" +
                    "\n\t  t -> query test mode" +
                    "\n\t  r -> compression test" +
                    "\n\t  x -> exit" +
                    "\n***********************************\n");
            String mode = sc.nextLine();        // take user's choice

            // switch to run user's choice
            switch (mode)
            {
                case "r":       // execute compression and decompression test

                    compressionTest(sc);

                    continue;
                case "i":       // calculate the indexing

                    makeIndexing(true, sc);             // make the inverted index

                    continue;                           // go next while iteration
                case "u":       // load or calculate the Term Upper Bound

                    calculateTUBs();            // calculate TUBs

                    continue;                           // go next while iteration
                case "f":       // see or change the flags
                    // Variables for control of operations to be redone in case of changed flags
                    boolean rebuildIndexing = false;    // var indicating if it's necessary to rebuild the indexing (recompute of TUBs is included)
                    boolean recomputeTUB = false;       // var indicating if it's necessary to recompute the Term Upper Bounds

                    // -- control for file into disk
                    if (!Flags.isThereFlagsFile())
                        printUI("The flags file is not saved in the disk. The flags will have the default value (all flags equal to false).\n");
                    else
                        readFlagsFromDisk();    // read the flags values from the disk

                    printFlagsUI();         // print the flags values with explanation for user

                    // prints for user
                    printUI("\nChanging the value of the flags may cause to redo some operations ( such as rebuilding the inverted index, recalculating the Term Upper Bounds and etc...)");
                    printUI("\nDo you want change the flags value?");
                    boolean yes = false;        // true = user want change any values
                    boolean no = false;         // true = user don't want change any values
                    do
                    {
                        printUI("Type Y for change or N for don't change.");
                        try
                        {
                            String choice = sc.nextLine().toUpperCase();        // take the user's choice
                            // check the user's input
                            if (choice.equals("Y"))
                                yes = true;             // set yes
                            else if (choice.equals("N"))
                                no = true;              // set no

                        } catch (NumberFormatException nfe) {
                            printError("Insert a valid character.");
                        }
                    } while (!(yes || no));  // continues until isConjunctive or isDisjunctive is set

                    if (no)
                        continue;       // the user doesn't want to change the values

                    // -- user wants to change the values
                    boolean sws = false;            // var for user preferences on the removal of stopwords
                    boolean compression = false;    // var for user preferences on the compression
                    boolean scoring = false;        // var for take user preferences on the scoring
                    boolean skipping = false;       // var for take user preferences on the skipping
                    boolean query_eff = false;      // var for take user preferences on the dynamic pruning algorithm (in query execution)
                    boolean deletePartFile = false; // var for take user preference on the delete partial file

                    // take the new values
                    sws = getUserChoice(sc, "stopwords removal");    // take user preferences on the removal of stopwords
                    compression = getUserChoice(sc, "compression");  // take user preferences on the compression
                    scoring = getUserChoice(sc, "scoring");          // take user preferences on the scoring
                    skipping = getUserChoice(sc, "skipping");        // take user preferences on the skipping
                    query_eff = getUserChoice(sc, "dynamic pruning algorithm"); // take user preferences on the dynamic pruning algorithm
                    deletePartFile = getUserChoice(sc, "delete partial file");  // take user preferences for the delete partial file

                    // check the changed values
                    if (isSwsEnabled() != sws)
                    {
                        printDebug("The value of the stopwords removal flag has been changed.");
                        setSws(sws);                            // change the value
                        rebuildIndexing = true;                 // set the value for the recomputing
                    }
                    if (isCompressionEnabled() != compression)
                    {
                        printDebug("The value of the compression flag has been changed.");
                        setCompression(compression);            // change the value
                        rebuildIndexing = true;                 // set the value for the recomputing
                    }
                    if (isScoringEnabled() != scoring)
                    {
                        printDebug("The value of the scoring flag has been changed.");
                        setScoring(scoring);                    // change the value
                        recomputeTUB = true;                    // set the value for the recomputing
                    }
                    if (considerSkippingBytes() != skipping)
                    {
                        printDebug("The value of the skipping flag has been changed.");
                        setConsiderSkippingBytes(skipping);     // change the value
                        rebuildIndexing = true;                 // set the value for the recomputing
                    }
                    if (isDynamicPruningEnabled() != query_eff)
                    {
                        printDebug("The value of the dynamic pruning algorithm flag has been changed.");
                        setDynamicPruning(query_eff);           // change the value
                        if (isDynamicPruningEnabled())
                            recomputeTUB = true;                // set the value for the recomputing
                    }
                    if (isDeletePartFileEnabled() != deletePartFile)
                    {
                        printDebug("The value of the delete partial file flag has been changed.");
                        setDeletePartFile(deletePartFile);      // change the value
                        if (isDeletePartFileEnabled())
                        {
                            delete_tempFiles();                 // delete the partial files of the indexing
                            printLoad("Partial file deleted.");
                        }
                    }

                    storeFlagsIntoDisk();      // store flags into disk

                    // in according to the changed flags values remake index or other
                    if (rebuildIndexing)
                        makeIndexing(false, sc);
                    else if(recomputeTUB)
                        calculateTUBs();

                    continue;                           // go next while iteration
                case "s":       // see files size

                    printFileSize();            // see and show the files size

                    continue;
                case "c":       // see collection statistic

                    CollectionStatistics.readCollectionStatsFromDisk(); // read statistics
                    CollectionStatistics.printCollectionStatistics();   // show collection statistic

                    continue;
                case "q":       // execute a query

                    //Flags.setConsiderSkippingBytes(true);
                    ArrayList<Integer> rankedResults;       // ArrayList that contain the ranked results of query
                    int numberOfResults = 0;    // take the integer entered by users that indicate the number of results wanted for query

                    // control check that all the files and resources required to execute a query are present
                    if (!queryStartControl()) {
                        return;                           // error exit
                    }
                    // read term upper bound
                    if(TermDocUpperBound.termUpperBoundTableIsEmpty())
                    {
                        if(TermDocUpperBound.termUpperBoundFileExist())     // the file already exist
                            TermDocUpperBound.readTermUpperBoundTableFromDisk();
                        else                                                // the file not exist
                            TermDocUpperBound.calculateTermsUpperBound(false);   // calculate term upper bound for each term of dictionary
                    }

                    printUI("Insert query: \n");
                    String query = sc.nextLine();           // take user's query
                    // control check of the query
                    if (query == null || query.isEmpty()) {
                        printError("Error: the query is empty. Please, retry.");
                        continue;                           // go next while iteration
                    }

                    boolean isConjunctive = false;      // true = Conjunctive query  | false = Disjunctive query (default case)
                    boolean isDisjunctive = false;      // true = user type D
                    // do while for choosing Conjunctive(AND) or Disjunctive(OR) query
                    do {
                        printUI("Type C for choosing Conjunctive query or D for choosing Disjunctive queries.");
                        try {
                            String choice = sc.nextLine().toUpperCase();        // take the user's choice
                            // check the user's input
                            if (choice.equals("C")) {
                                isConjunctive = true;           // set isConjunctive
                            } else if (choice.equals("D")) {
                                isDisjunctive = true;           // set isDisjunctive
                            }
                        } catch (NumberFormatException nfe) {
                            printError("Insert a valid character.");
                        }
                    } while (!(isConjunctive || isDisjunctive));  // continues until isConjunctive or isDisjunctive is set

                    int validN = 0;     // indicates if the entered number is valid or not -> 1 = 10 or 20 - 0 = negative number or not a number or another number
                    // do while for choosing the number of results to return
                    do {
                        printUI("Type the number of results to retrieve (10 or 20)");
                        try {
                            numberOfResults = Integer.parseInt(sc.nextLine());    // take the int inserted by user
                            validN = ((numberOfResults == 10) || (numberOfResults == 20)) ? 1 : 0;  // validity check
                        } catch (NumberFormatException nfe) {
                            printError("Insert a valid positive number");
                        }
                    } while (validN == 0);  // continues until a valid number is entered

                    startTime = System.currentTimeMillis();         // start time of execute query

                    // do query and retry the results
                    rankedResults = QueryProcessor.queryManager(query,isConjunctive,numberOfResults);
                    QueryProcessor.printQueryResults(rankedResults);               // print the results of the query

                    endTime = System.currentTimeMillis();           // end time of execute query

                    // shows query execution time
                    printTime("\nQuery executes in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
                    continue;                       // go next while iteration

                case "t":       // query test mode

                    int validNum = 0;     // 1 = positive number - 0 = negative number or not a number
                    int numberOfQueries = 0;    // take the integer entered by users that indicate the number of queries to test

                    // do while for choosing the number of results to return
                    do {
                        printUI("Type the number of queries to test (must be a positive number)");
                        try {
                            numberOfQueries = Integer.parseInt(sc.nextLine());    // take the int inserted by user
                            validNum = (numberOfQueries > 0) ? 1 : 0;               // validity check of the int
                        } catch (NumberFormatException nfe) {
                            printError("Insert a valid positive number");
                        }
                    } while (validNum == 0);  // continues until a valid number is entered

                    QueryProcessor.readQueryFromCollection(numberOfQueries);
                    continue;

                default:
                    return;     // exit to switch, case not valid
            }
        }
    }

    // ----------------------------------------- start : mode(switch) functions -----------------------------------------

    /**
     * Function to making the inverted index (SPIMI, merge and calculate Term Upper Bound)
     *
     * @param takeFlags if true: have to set the flags | if false: have not to set the flags
     * @param sc
     * @throws IOException
     */
    private static void makeIndexing(boolean takeFlags, Scanner sc) throws IOException
    {
        long startTime, endTime;        // variables to calculate the execution time

        file_cleaner();                 // delete all created files
        // take user choice for the flags
        if (takeFlags)
        {
            printFlagsUI();             // print the flags values with explanation for user

            setSws(getUserChoice(sc, "stopwords removal"));    // take user preferences on the removal of stopwords
            setCompression(getUserChoice(sc, "compression"));  // take user preferences on the compression
            setScoring(getUserChoice(sc, "scoring"));          // take user preferences on the scoring
            setConsiderSkippingBytes(getUserChoice(sc, "skipping"));            // take user preferences on the scoring
            setDynamicPruning(getUserChoice(sc, "dynamic pruning algorithm"));  // take user preferences for the dynamic pruning algorithm
            setDeletePartFile(getUserChoice(sc, "delete partial file"));        // take user preferences for the delete partial file

            storeFlagsIntoDisk();       // store Flags
        }
        // do SPIMI Algorithm
        printLoad("\nIndexing...");
        startTime = System.currentTimeMillis();         // start time to SPIMI Algorithm
        PartialIndexBuilder.SPIMIalgorithm();           // do SPIMI
        endTime = System.currentTimeMillis();           // end time of SPIMI algorithm
        printTime("\nSPIMI Algorithm done in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

        // merge blocks into disk
        startTime = System.currentTimeMillis();         // start time to merge blocks
        IndexMerger.mergeBlocks();                      // merge blocks
        endTime = System.currentTimeMillis();           // end time of merge blocks
        printTime("\nBlocks merged in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

        DataStructureHandler.calcAndStoreDenPartBM25inDocTable();   // calculate the partial denominator BM25 for optimization
        QueryProcessor.calcAndStoreTFWeight();                      // calculate and store the log for the TFWeight
        // calculate term upper bound if dynamic pruning algorithm is enabled
        if (!queryStartControl())
            return;                         // error exit
        TermDocUpperBound.calculateTermsUpperBound(true);   // calculate term upper bound for each term of dictionary

        // check if delete or not the partial files
        if (isDeletePartFileEnabled())
        {
            delete_tempFiles();             // delete the partial files of the indexing
            printLoad("Partial file deleted.");
        }
    }

    /**
     * Function to calculate the Term Upper bound for each term in the vocabulary.
     *
     * @throws IOException
     */
    private static void calculateTUBs() throws IOException
    {
        //Flags.setConsiderSkippingBytes(true);
        if (!queryStartControl())
            return;                           // error exit

        DataStructureHandler.calcAndStoreDenPartBM25inDocTable();

        TermDocUpperBound.calculateTermsUpperBound(false);   // calculate term upper bound for each term of dictionary
        //TermDocUpperBound.readTermUpperBoundTableFromDisk();

        TermDocUpperBound.calculateDocsUpperBound();    // calculate doc upper bound for each doc of docTable
    }

    /**
     * Function to test and show the result of compression and decompression. This function make some different test:
     * 1 - the first test is composed by two different test, one for Unary-code and the other for variable-byte code.
     *     Unary-code: in this test is filled and displays a list of positive integers that will be compressed using
     *                 unary code. The compression bits will also be displayed, then the decompression result will be
     *                 displayed.
     *     Variable-byte code: in this test more integers are added to the list and then compressing and then
     *                 decompressing using both Gap and non-Gap modes, and also mixing the two modes to show the results
     *                 of all cases.
     *     These tests serve to make a logical check of the functioning of the two compression algorithms with a simple
     *     and not too long list.
     * 2 - this second test changes depending on how the index has been constructed, i.e. depending on which flags are
     *     enabled or disabled. in all cases the user can enter a term and its posting list will be used for testing.
     *     The two tests that are performed are:
     *     - compression not enabled: in this case the posting lists are not compressed, therefore saved raw on disk.
     *                 The posting list is read from memory and then the TF list is compressed and decompressed using
     *                 unary-code while the DID list is compressed and decompressed (in the various modes) using
     *                 variable-byte code. At the end, it is also checked that the result of the compressions is equal
     *                 to the posting list read from memory and then the percentage and byte space gains of the
     *                 compressed version compared to the un-compressed one are displayed.
     *     - compression and skipping enabled: in this case the posting list is read from memory and then uncompressed
     *                 using different modes (at once or split into blocks) and then the execution time and how many
     *                 bytes were read are shown.
     *
     * @param sc    scanner to get the choice of the user inserted via keyboard
     */
    private static void compressionTest(Scanner sc)
    {
        printDebug("1)Simple compression test, use a simple and short list: ");
        ArrayList<Integer> myNumbers = new ArrayList<Integer>();
        String chosenTerm = "";     // the term whose posting list will be used for second tests
        int listLen = 10;

        // 1° - test
        // -- Unary test
        for (int i = 0; i < listLen; i++)   // insert in myNumber the number form 1 to listLen
            myNumbers.add(i+1);

        printDebug("-- Unary code test -> the test list is: " + myNumbers);
        // compression
        byte[] compressedResult = Unary.integersCompression(myNumbers);     // compress to Unary
        printDebug("Unary compression -- the list passed len: " + myNumbers.size() + " int with size: " + (myNumbers.size()*4) + " Bytes -> after compression the size is: " + compressedResult.length + " Bytes.");
        Unary.printCompressedList(compressedResult);
        // decompression
        myNumbers = Unary.integersDecompression(compressedResult, listLen);
        printDebug("Unary Decompression -- the compressed list len: " + listLen + " int with size: " + compressedResult.length + " Bytes -> after decompression the size is: " + (myNumbers.size()*4) + " Bytes.");
        Unary.printDecompressedList(myNumbers);    // decompress from Unary

        // -- Variable Bytes test
        myNumbers.clear();
        for (int i = 0; i < (listLen/2); i++)
            myNumbers.add(i+1);
        for (int i = 80000; i < 80000 + listLen; i++)
            myNumbers.add(i);
        myNumbers.add(214577);
        printDebug("\n-- Variable Bytes code test -> the test list is: " + myNumbers);
        // compression DGaps = false
        compressedResult = VariableBytes.integersCompression(myNumbers,false);  // compress to Variable Bytes
        printDebug("VarBytes Compression -- Dgaps: " + false + " the list passed len: " + myNumbers.size() + " int with size: " + (myNumbers.size()*4) + " Bytes -> after compression the size is: " + compressedResult.length + " Bytes.");
        VariableBytes.printCompressedList(compressedResult);
        // decompression DGaps = false
        myNumbers = VariableBytes.integersDecompression(compressedResult,false);
        printDebug("VarBytes Decompression -- Dgaps: " + false + " the compressed list with size: " + compressedResult.length + " Bytes -> after decompression the size is: " + (myNumbers.size()*4) + " Bytes.");
        VariableBytes.printDecompressedList(myNumbers);                                // decompress from Variable Bytes
        // compression DGaps = true
        compressedResult = VariableBytes.integersCompression(myNumbers,true);   // compress to Variable Bytes
        printDebug("VarBytes Compression -- Dgaps: " + false + " the list passed len: " + myNumbers.size() + " int with size: " + (myNumbers.size()*4) + " Bytes -> after compression the size is: " + compressedResult.length + " Bytes.");
        VariableBytes.printCompressedList(compressedResult);
        // decompression DGaps = false
        myNumbers = VariableBytes.integersDecompression(compressedResult,false);
        printDebug("VarBytes Decompression -- Dgaps: " + false + " the compressed list with size: " + compressedResult.length + " Bytes -> after decompression the size is: " + (myNumbers.size()*4) + " Bytes.");
        VariableBytes.printDecompressedList(myNumbers);            // decompress from Variable Bytes
        // decompression DGaps = true
        myNumbers = VariableBytes.integersDecompression(compressedResult,true);
        printDebug("VarBytes Decompression -- Dgaps: " + true + " the compressed list with size: " + compressedResult.length + " Bytes -> after decompression the size is: " + (myNumbers.size()*4) + " Bytes.");
        VariableBytes.printDecompressedList(myNumbers);             // decompress from Variable Bytes
        // 2° - test

        printUI("Please enter the term whose posting list will be used for testing.");
        chosenTerm = sc.nextLine().toUpperCase();                    // take the user's choice
        // preprocess of the entered term
        ArrayList<String> procChosenTerm;                       // array list for containing the query term
        try {
            procChosenTerm = TextProcessor.preprocessText(chosenTerm); // Preprocessing of document text
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        // check if the posting list isn't already compressed
        if (!Flags.isThereFlagsFile())
            printError("Error: missing required files.");
        else
        {
            readFlagsFromDisk();        // read flag from disk
            // the second test depends on the enabled flags
            if (!Flags.isCompressionEnabled())          // compression not enabled
            {
                //compPLTest("giacomo");
                printUIMag("The posting lists are not compressed, testing the compression of a posting list.");
                compPLTest(procChosenTerm.get(0));
            }
            else if (Flags.considerSkippingBytes())     // compression enabled
            {
                printUIMag("The posting lists are compressed, reading test of a compressed posting list.");
                compPLTestRead(procChosenTerm.get(0));
            }
        }
    }

    /**
     * Function to test and show the result of compression and decompression on a posting list.
     *
     * @param term  the term related to the posting list to compress and decompress
     */
    private static void compPLTest(String term)
    {
        ArrayList<Integer> termFreqToCompress = new ArrayList<>();  // array list to contain the TF to compress
        ArrayList<Integer> docIDsToCompress =  new ArrayList<>();   // array list to contain the DID to compress
        ArrayList<Integer> decompressTF;        // array list to contain the decompressed TF
        ArrayList<Integer> decomDIDNoDGaps;     // array list to contain the decompressed DID
        ArrayList<Integer> decomDIDDGaps;       // array list to contain the decompressed DID
        ArrayList<Integer> decomDIDDGapsBlock;  // array list to contain the decompressed DID
        ArrayList<Posting> postingList;         // array list for the posting list of the term
        byte[] compressedTermFreq;              // byte array for the compressed term frequency
        byte[] compDocIDWithoutDGaps;           // byte array for the compressed DID without gap
        byte[] compDocIDDGaps;                  // byte array for the compressed DID with gap
        byte[] compDocIDDGapsBlock;             // byte array for the compressed DID with gap and block
        int NotCompressedLen;                   // len in bytes
        int blockSize = SKIP_POINTERS_THRESHOLD;// block len
        long startTime, endTime;                // variables to calculate the execution time

        // -- control for structures in memory - if not load them from disk
        if (!QueryProcessor.dictionary.dictionaryIsSet())
        {
            startTime = System.currentTimeMillis();
            QueryProcessor.dictionary.readDictionaryFromDisk();
            endTime = System.currentTimeMillis();
            printTime( "Dictionary loaded in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
        }

        printUIMag("\n2)More complicated compression test, usa a posting list of a term: " + term);
        // check if the term is in the dictionarys
        if (QueryProcessor.dictionary.getTermToTermStat().containsKey(term))
        {
            // get the posting list of the term passed as parameter
            try(
                    // open complete files to read the postingList
                    RandomAccessFile docidFile = new RandomAccessFile(DOCID_FILE, "rw");
                    RandomAccessFile termfreqFile = new RandomAccessFile(TERMFREQ_FILE, "rw");

                    // FileChannel
                    FileChannel docIdChannel = docidFile.getChannel();
                    FileChannel termFreqChannel = termfreqFile.getChannel()
            ) {
                DictionaryElem de = QueryProcessor.dictionary.getTermToTermStat().get(term);
                startTime = System.currentTimeMillis();
                postingList = readPostingListFromDisk(de.getOffsetDocId(),de.getOffsetTermFreq(),de.getDf(),docIdChannel,termFreqChannel);
                endTime = System.currentTimeMillis();
                printTime( "Posting list loaded in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

            } catch (IOException e) {
                throw new RuntimeException(e);
            }

            // fill the array list for termFreq and DocID
            assert postingList != null;
            for (Posting p : postingList)
            {
                termFreqToCompress.add(p.getTermFreq());
                docIDsToCompress.add(p.getDocId());
            }
            NotCompressedLen = postingList.size() * 4;

            // test of compression and decompression
            printUIMag("- NOT COMPRESSED -> PL size: " + postingList.size() + " TermFreq: " + (postingList.size()*4) + " Bytes , DID: " + (postingList.size()*4) +" Bytes\n");
            // -- compression
            // ---- TF compression
            printLoad("Term Freq Unary compression...");
            startTime = System.currentTimeMillis();
            compressedTermFreq = Unary.integersCompression(termFreqToCompress);     // compress to Unary
            endTime = System.currentTimeMillis();
            printTime( "Compressed in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

            // ---- DID compression
            printUIMag("\nDID VariableBytes (No DGaps) compression...");
            startTime = System.currentTimeMillis();
            compDocIDWithoutDGaps = VariableBytes.integersCompression(docIDsToCompress,false);    // compress to Variable Bytes
            endTime = System.currentTimeMillis();
            printTime( "Compressed in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

            printUIMag("DID VariableBytes (DGaps) compression...");
            startTime = System.currentTimeMillis();
            compDocIDDGaps = VariableBytes.integersCompression(docIDsToCompress,true);    // compress to Variable Bytes
            endTime = System.currentTimeMillis();
            printTime( "Compressed in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

            printUIMag("DID VariableBytes (DGaps with block) compression...");
            startTime = System.currentTimeMillis();
            compDocIDDGapsBlock = VariableBytes.intCompDGapsBlock(docIDsToCompress,blockSize);    // compress to Variable Bytes
            endTime = System.currentTimeMillis();
            printTime( "Compressed in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

            // -- decompression
            // ---- TF decompression
            printUIMag("\nTerm Freq Unary decompression...");
            startTime = System.currentTimeMillis();
            decompressTF = Unary.integersDecompression(compressedTermFreq, postingList.size());    // decompress from Unary
            endTime = System.currentTimeMillis();
            printTime( "Decompressed in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

            // ---- DID decompression
            printUIMag("\nDID VariableBytes (No DGaps) decompression...");
            startTime = System.currentTimeMillis();
            decomDIDNoDGaps = VariableBytes.integersDecompression(compDocIDWithoutDGaps,false);            // decompress from Variable Bytes
            endTime = System.currentTimeMillis();
            printTime( "Decompressed in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

            printUIMag("DID VariableBytes (DGaps) decompression...");
            startTime = System.currentTimeMillis();
            decomDIDDGaps = VariableBytes.integersDecompression(compDocIDDGaps,true);            // decompress from Variable Bytes
            endTime = System.currentTimeMillis();
            printTime( "Decompressed in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

            printUIMag("DID VariableBytes (DGaps block) decompression...");
            startTime = System.currentTimeMillis();
            decomDIDDGapsBlock = VariableBytes.intDecompDGapsBlock(compDocIDDGapsBlock,blockSize);            // decompress from Variable Bytes
            endTime = System.currentTimeMillis();
            printTime( "Decompressed in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

            // verify the correctness of compression and decompression
            // -- verify term freq (unary)
            if (isTwoIntArrayListEqual(termFreqToCompress, decompressTF))
                printUIMag("Successful compression and decompression of TermFreq list.");
            else
                printUIMag("Not successful compression and decompression of TermFreq list.");

            // -- verify DID (var bytes)
            if (isTwoIntArrayListEqual(docIDsToCompress, decomDIDNoDGaps) && isTwoIntArrayListEqual(docIDsToCompress, decomDIDDGaps) && isTwoIntArrayListEqual(docIDsToCompress, decomDIDDGapsBlock))
                printUIMag("Successful compression and decompression of DocID list.");
            else
                printUIMag("Not successful compression and decompression of DocID list.");

            // print statistics result
            printUIMag("Statistics results:");
            printUIMag("- TermFreq not compressed: " + formatBytes(NotCompressedLen));
            printUIMag("--> Unary compression: " + formatBytes(compressedTermFreq.length) + " -> " + differenceBetweenTwoInt(NotCompressedLen, compressedTermFreq.length));
            printUIMag("- DocID not compressed: " + formatBytes(NotCompressedLen));
            printUIMag("--> Var Bytes NoDGaps compression: " + formatBytes(compDocIDWithoutDGaps.length) + " -> " + differenceBetweenTwoInt(NotCompressedLen, compDocIDWithoutDGaps.length));
            printUIMag("--> Var Bytes DGaps compression: " + formatBytes(compDocIDDGaps.length) + " -> " + differenceBetweenTwoInt(NotCompressedLen, compDocIDDGaps.length));
            printUIMag("--> Var Bytes DGaps Block compression: " + formatBytes(compDocIDDGapsBlock.length) + " Bytes -> " + differenceBetweenTwoInt(NotCompressedLen, compDocIDDGapsBlock.length));
        }
        else
            printError("The term isn't in the list.");
    }

    /**
     * Function to test and show the result of reading a compression posting list a block or in one time.
     *
     * @param term  the term related to the posting list to read
     */
    private static void compPLTestRead(String term)
    {
        ArrayList<Posting> postingList;     // array list for the posting list of the term
        SkipList sList;                     // skipList related to the term
        DictionaryElem de;                  // temp dictionary elem
        int byteTFLoaded = 0;               // number of byte loaded form the disk for the compressed TermFreq
        int byteDIDLoaded = 0;              // number of byte loaded form the disk for the compressed DID
        byte[] tf;                          // array for the compressed TermFreq list
        byte[] docids;                      // array for the compressed TermFreq list
        long startTime, endTime;            // variables to calculate the execution time

        // -- control for structures in memory - if not load them from disk
        if (!QueryProcessor.dictionary.dictionaryIsSet())
        {
            startTime = System.currentTimeMillis();
            QueryProcessor.dictionary.readDictionaryFromDisk();
            endTime = System.currentTimeMillis();
            printTime( "Dictionary loaded in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
        }

        printUIMag("\n2)More complicated reading test, usa a posting list of a term: " + term);
        // check if the term is in the dictionarys
        if (QueryProcessor.dictionary.getTermToTermStat().containsKey(term))
        {
            de = QueryProcessor.dictionary.getTermToTermStat().get(term);
            // set the skipList instance
            sList = new SkipList(de.getSkipOffset(), de.getSkipArrLen(),null, de.getDf());
            // get the posting list of the term passed as parameter
            try(
                    // open complete files to read the postingList
                    RandomAccessFile docidFile = new RandomAccessFile(DOCID_FILE, "rw");
                    RandomAccessFile termfreqFile = new RandomAccessFile(TERMFREQ_FILE, "rw");
                    // FileChannel
                    FileChannel docIdChannel = docidFile.getChannel();
                    FileChannel termFreqChannel = termfreqFile.getChannel()
            ) {
                // 1 - read and not uncompress the whole compressed posting list
                startTime = System.currentTimeMillis();
                // read the termFreq of the whole posting list
                int startOffsetTF;                  // the current offset (at each iteration) for the compressed TermFreq
                int sizeToReadTF = 0;               // the current size (at each iteration) for the compressed TermFreq

                startOffsetTF = (int) de.getOffsetTermFreq();
                sizeToReadTF = (int) (sList.getSkipBlockInfo(sList.getSkipArrLen()-1).getFreqOffset() - startOffsetTF);

                MappedByteBuffer termfreqBuffer = termFreqChannel.map(FileChannel.MapMode.READ_ONLY, startOffsetTF, sizeToReadTF);
                tf = new byte[sizeToReadTF];        // initialize the byte array for the compressed TF list
                termfreqBuffer.get(tf, 0, sizeToReadTF);       // read the TF compressed block of the PL

                // read the DID of the whole posting list
                int startOffsetDID;                  // the current offset (at each iteration) for the compressed TermFreq
                int sizeToReadDID = 0;               // the current size (at each iteration) for the compressed TermFreq

                startOffsetDID = (int) de.getOffsetDocId();
                sizeToReadDID = (int) (sList.getSkipBlockInfo(sList.getSkipArrLen()-1).getDocIdOffset() - startOffsetDID);

                MappedByteBuffer docidBuffer = docIdChannel.map(FileChannel.MapMode.READ_ONLY, startOffsetDID, sizeToReadDID);
                docids = new byte[sizeToReadDID];          // initialize the byte array for the compressed DID list
                docidBuffer.get(docids, 0, sizeToReadDID);   // read the DID compressed block of the PL

                endTime = System.currentTimeMillis();
                printUIMag("Byte loaded for TF: " + sizeToReadTF + " Byte loaded for DID: " + sizeToReadDID);
                printTime( "Whole posting list loaded and not uncompressed in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

                // 1.5 - read and uncompress the whole compressed posting list
                startTime = System.currentTimeMillis();
                // read the termFreq of the whole posting list
                startOffsetTF = (int) de.getOffsetTermFreq();
                sizeToReadTF = (int) (sList.getSkipBlockInfo(sList.getSkipArrLen()-1).getFreqOffset() - startOffsetTF);

                MappedByteBuffer termFreqBuffer = termFreqChannel.map(FileChannel.MapMode.READ_ONLY, startOffsetTF, sizeToReadTF);
                tf = new byte[sizeToReadTF];        // initialize the byte array for the compressed TF list
                termFreqBuffer.get(tf, 0, sizeToReadTF);       // read the TF compressed block of the PL

                // read the DID of the whole posting list
                startOffsetDID = (int) de.getOffsetDocId();
                sizeToReadDID = (int) (sList.getSkipBlockInfo(sList.getSkipArrLen()-1).getDocIdOffset() - startOffsetDID);

                MappedByteBuffer docIdBuffer = docIdChannel.map(FileChannel.MapMode.READ_ONLY, startOffsetDID, sizeToReadDID);
                docids = new byte[sizeToReadDID];          // initialize the byte array for the compressed DID list
                docIdBuffer.get(docids, 0, sizeToReadDID);   // read the DID compressed block of the PL

                // uncompress
                int startTF = 0;
                int startDID = 0;
                int numTFComp = 0;
                for(int i = 0; i < sList.getSkipArrLen(); i++)
                {
                    /*
                    if (i != 0)
                    {
                        startTF = (int)sList.getSkipBlockInfo(i-1).getFreqOffset() - (int) de.getOffsetTermFreq();
                        startDID = (int)sList.getSkipBlockInfo(i-1).getDocIdOffset() - (int) de.getOffsetDocId();
                    }
                     */
                    //printUIMag("Iteration: " + i + " startTF: " + startTF + " startDID: " + startDID);

                    //numTFComp = min(SKIP_POINTERS_THRESHOLD, (de.getDf() - (SKIP_POINTERS_THRESHOLD * 0)));
                    numTFComp = min(SKIP_POINTERS_THRESHOLD, (de.getDf() - (SKIP_POINTERS_THRESHOLD * i)));
                    //Arrays.copyOfRange(tf, startTF, (int)sList.getSkipBlockInfo(i).getFreqOffset());
                    //Arrays.copyOfRange(docids, startDID, (int)sList.getSkipBlockInfo(i).getDocIdOffset());
                    //ArrayList<Integer> uncompressedTf = Unary.integersDecompression(Arrays.copyOfRange(tf, startTF, (int)sList.getSkipBlockInfo(0).getFreqOffset()), numTFComp);  // decompress term freq
                    //ArrayList<Integer> uncompressedDocid = VariableBytes.integersDecompression(Arrays.copyOfRange(docids, startDID, (int)sList.getSkipBlockInfo(0).getDocIdOffset()),true);    // decompress DocID
                    ArrayList<Integer> uncompressedTf = Unary.integersDecompression(Arrays.copyOfRange(tf, startTF, (int)sList.getSkipBlockInfo(i).getFreqOffset()), numTFComp);  // decompress term freq
                    ArrayList<Integer> uncompressedDocid = VariableBytes.integersDecompression(Arrays.copyOfRange(docids, startDID, (int)sList.getSkipBlockInfo(i).getDocIdOffset()),true);    // decompress DocID
                }

                endTime = System.currentTimeMillis();
                printUIMag("Byte loaded for TF: " + sizeToReadTF + " Byte loaded for DID: " + sizeToReadDID);
                printTime( "Whole posting list loaded and uncompressed in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

                // 2 - read and uncompress the blocks of the compressed posting list (one call function)
                startTime = System.currentTimeMillis();
                postingList = readAndUncompressCompressedAndSkippedPLFromDisk(sList, de.getOffsetDocId(), de.getOffsetTermFreq(), de.getTermFreqSize(), de.getDocIdSize(), de.getSkipArrLen(), de.getDf(), docIdChannel, termFreqChannel);
                endTime = System.currentTimeMillis();
                printUIMag("Byte loaded for TF: " + de.getTermFreqSize() + " Byte loaded for DID: " + de.getDocIdSize());
                printTime( "The blocks (one call function) of posting list loaded and uncompressed in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

                // 3 - read and not uncompress the blocks of the compressed posting list (more call function)
                startTime = System.currentTimeMillis();
                for(int i = 0; i < sList.getSkipArrLen(); i++)
                {
                    tf = readCompTFBlockFromDisk(sList, i, de.getOffsetTermFreq(), de.getTermFreqSize(), de.getSkipArrLen(), termFreqChannel);
                    byteTFLoaded += tf.length;
                    docids = readCompDIDBlockFromDisk(sList, i, de.getOffsetDocId(), de.getDocIdSize(), de.getSkipArrLen(), docIdChannel);
                    byteDIDLoaded += docids.length;
                }
                endTime = System.currentTimeMillis();
                printUIMag("Byte loaded for TF: " + byteTFLoaded + " Byte loaded for DID: " + byteDIDLoaded);
                printTime( "The blocks (more call function) of posting list loaded and not uncompressed in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

                // 4 - read and uncompress the blocks of the compressed posting list (more call function)
                startTime = System.currentTimeMillis();
                byteTFLoaded = 0;
                byteDIDLoaded = 0;
                for(int i = 0; i < sList.getSkipArrLen(); i++)
                {
                    tf = readCompTFBlockFromDisk(sList, i, de.getOffsetTermFreq(), de.getTermFreqSize(), de.getSkipArrLen(), termFreqChannel);
                    docids = readCompDIDBlockFromDisk(sList, i, de.getOffsetDocId(), de.getDocIdSize(), de.getSkipArrLen(), docIdChannel);
                    byteTFLoaded += tf.length;
                    byteDIDLoaded += docids.length;
                    numTFComp = min(SKIP_POINTERS_THRESHOLD, (de.getDf() - (SKIP_POINTERS_THRESHOLD * i)));
                    if ( !(tf == null) && !(docids == null) )       // decompress the block
                    {
                        ArrayList<Integer> uncompressedTf = Unary.integersDecompression(tf, numTFComp);  // decompress term freq
                        ArrayList<Integer> uncompressedDocid = VariableBytes.integersDecompression(docids,true);    // decompress DocID
                    }
                }
                endTime = System.currentTimeMillis();
                printUIMag("Byte loaded for TF: " + byteTFLoaded + " Byte loaded for DID: " + byteDIDLoaded);
                printTime( "The blocks (more call function) of posting list loaded and uncompressed in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }
        else
            printError("The term isn't in the list.");
    }

    // ----------------------------------------- end : mode(switch) functions -----------------------------------------

    // ------------------------------------------- start : utility functions -------------------------------------------
    /**
     * fucntion to get the choise of the user for options, the options are pass
     *
     * @param sc     scanner to get the choice of the user inserted via keyboard
     * @param option options passed by parameter
     * @return true if the user chooses yes (enter Y), false if the user chooses no (enter N)
     */
    private static boolean getUserChoice(Scanner sc, String option) {
        while (true) {
            printUI("\nType Y or N for " + option + " options");   // print of the option
            String choice = sc.nextLine().toUpperCase();                    // take the user's choice
            // check the user's input
            if (choice.equals("Y")) {
                return true;
            } else if (choice.equals("N")) {
                return false;
            }
        }
    }
    // ------------------------------------------- end : utility functions -------------------------------------------
}


