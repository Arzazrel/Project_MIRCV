package it.unipi.dii.aide.mircv;

import it.unipi.dii.aide.mircv.compression.Unary;
import it.unipi.dii.aide.mircv.compression.VariableBytes;
import it.unipi.dii.aide.mircv.data_structures.*;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

import static it.unipi.dii.aide.mircv.QueryProcessor.queryStartControl;
import static it.unipi.dii.aide.mircv.utils.FileSystem.*;
import static it.unipi.dii.aide.mircv.utils.Constants.*;
import static it.unipi.dii.aide.mircv.data_structures.Flags.*;

public class Main {

    public static void main(String[] args) throws IOException {

        Scanner sc = new Scanner(System.in);
        long startTime, endTime;                // variables to calculate the execution time

        // while constituting the user interface
        while (true) {
            // print of the user interface
            printUI(
                    "\n********** SEARCH ENGINE **********" +
                    "\n\tSelect an option:" +
                    "\n\t  m -> try merge only" +
                    "\n\t  i -> build the index" +
                    "\n\t  d -> offset debug" +
                    "\n\t  u -> calculate term upper bound and document upper bound" +
                    "\n\t  q -> query mode" +
                    "\n\t  t -> query test mode" +
                    "\n\t  x -> exit" +
                    "\n***********************************\n");
            String mode = sc.nextLine();        // take user's choice

            // switch to run user's choice
            switch (mode) {

                case "m":       // per debugging, prova solo il merge
                    delete_mergedFiles();
                    //setCompression(true);  // take user preferences on the compression

                    DataStructureHandler.readBlockOffsetsFromDisk();

                    setSws(getUserChoice(sc, "stopwords removal"));    // take user preferences on the removal of stopwords
                    setCompression(getUserChoice(sc, "compression"));  // take user preferences on the compression
                    setScoring(getUserChoice(sc, "scoring"));          // take user preferences on the scoring
                    storeFlagsIntoDisk();      // store Flags

                    startTime = System.currentTimeMillis();         // start time to merge blocks from disk
                    IndexMerger.mergeBlocks();                      // merge
                    endTime = System.currentTimeMillis();           // end time to merge blocks from disk
                    printTime( "Merged in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
                    continue;                                   // go next while cycle

                case "i":
                    file_cleaner();                             // delete all created files

                    setSws(getUserChoice(sc, "stopwords removal"));    // take user preferences on the removal of stopwords
                    setCompression(getUserChoice(sc, "compression"));  // take user preferences on the compression
                    setScoring(getUserChoice(sc, "scoring"));          // take user preferences on the scoring

                    storeFlagsIntoDisk();      // store Flags
                    // Do SPIMI Algorithm
                    System.out.println("\nIndexing...");
                    startTime = System.currentTimeMillis();         // start time to SPIMI Algorithm
                    PartialIndexBuilder.SPIMIalgorithm();          // do SPIMI
                    endTime = System.currentTimeMillis();           // end time of SPIMI algorithm
                    printTime("\nSPIMI Algorithm done in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

                    // merge blocks into disk
                    startTime = System.currentTimeMillis();         // start time to merge blocks
                    IndexMerger.mergeBlocks();                      // merge blocks
                    endTime = System.currentTimeMillis();           // end time of merge blocks
                    printTime("\nBlocks merged in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

                    // calculate term upper bound and doc upper bound
                    Flags.setConsiderSkippingBytes(true);
                    if (!queryStartControl()) {
                        return;                           // error exit
                    }
                    TermDocUpperBound.calculateTermsUpperBound();   // calculate term upper bound for each term of dictionary

                    continue;                           // go next while iteration
                case "d":

                    queryStartControl();
                    String term = "0000";
                    printDebug(QueryProcessor.dictionary.getTermStat(term).toString());

                    continue;                           // go next while iteration
                case "u":

                    Flags.setConsiderSkippingBytes(true);
                    if (!queryStartControl()) {
                        return;                           // error exit
                    }
                    TermDocUpperBound.calculateTermsUpperBound();   // calculate term upper bound for each term of dictionary

                    TermDocUpperBound.calculateDocsUpperBound();    // calculate doc upper bound for each doc of docTable

                    continue;                           // go next while iteration
                case "q":       // query

                    Flags.setConsiderSkippingBytes(true);
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
                            TermDocUpperBound.readTermUpperBoundTableIntoDisk();
                        else                                                // the file not exist
                            TermDocUpperBound.calculateTermsUpperBound();   // calculate term upper bound for each term of dictionary
                    }

                    printUI("Insert query: \n");
                    String query = sc.nextLine();           // take user's query
                    // control check of the query
                    if (query == null || query.isEmpty()) {
                        printError("Error: the query is empty. Please, retry.");
                        continue;                           // go next while iteration
                    }

                    boolean isConjunctive = false;          // true = Conjunctive query
                    boolean isDisjunctive = false;          // true = Disjunctive query
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

                    int validN = 0;     // 1 = positive number - 0 = negative number or not a number
                    // do while for choosing the number of results to return
                    do {
                        printUI("Type the number of results to retrieve (10 or 20)");
                        try {
                            numberOfResults = Integer.parseInt(sc.nextLine());    // take the int inserted by user
                            validN = (numberOfResults > 0) ? 1 : 0;               // validity check of the int
                        } catch (NumberFormatException nfe) {
                            printError("Insert a valid positive number");
                        }
                    } while (validN == 0);  // continues until a valid number is entered

                    startTime = System.currentTimeMillis();         // start time of execute query

                    // do query and retry the results
                    rankedResults = QueryProcessor.queryManager(query,isConjunctive,isDisjunctive,numberOfResults);
                    QueryProcessor.printQueryResults(rankedResults);               // print the results of the query

                    endTime = System.currentTimeMillis();           // end time of execute query

                    // shows query execution time
                    printTime("\nQuery executes in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
                    continue;                       // go next while iteration

                case "t":       // query test mode

                    int validNum = 0;     // 1 = positive number - 0 = negative number or not a number
                    int numberOfQueries = 0;    // take the integer entered by users that indicate the number of queries to test

                    Flags.setConsiderSkippingBytes(true);

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


    /**
     * function to shows the user the ranked results (DocID) of the query executed
     *
     * @param rankedResults the results returned by the query
     *//*
    private static void printQueryResults(ArrayList<Integer> rankedResults)
    {
        if (rankedResults.size() != 0)      // there are results
        {
            System.out.println(ANSI_CYAN + "Query results:" + ANSI_RESET);
            for (int i = 0; i < rankedResults.size(); i++)
                printUI((i + 1) + " - " + rankedResults.get(i));
        }
        else                                // there aren't results
            printUI("No results found for this query.");
    }*/
}


