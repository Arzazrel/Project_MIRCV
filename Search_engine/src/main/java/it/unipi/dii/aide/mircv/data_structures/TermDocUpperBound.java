package it.unipi.dii.aide.mircv.data_structures;

import it.unipi.dii.aide.mircv.compression.Unary;
import it.unipi.dii.aide.mircv.QueryProcessor;
import it.unipi.dii.aide.mircv.compression.VariableBytes;
import org.apache.commons.io.FileUtils;

import java.io.*;
import java.nio.CharBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;

import static it.unipi.dii.aide.mircv.data_structures.DocumentElement.*;
import static it.unipi.dii.aide.mircv.data_structures.PartialIndexBuilder.*;
import static it.unipi.dii.aide.mircv.utils.Constants.*;
import static it.unipi.dii.aide.mircv.utils.Logger.*;

/**
 * This class perform the calculating of term upper bound.
 */

public class TermDocUpperBound
{
    // Data structures initialization
    static HashMap<Integer, Double> docUpperBoundTable = new HashMap<>();   // hash table DocID to related DocElement
    static HashMap<String, Double> termUpperBoundTable = new HashMap<>();   // hash table Term to related Posting list


    // ---------------- start: calculating functions ----------------

    /**
     * Calculate the term upper bound for each term in the dictionary
     *
     * @param computeStats indicates whether or not to compute statistics on posting list lengths or occurrences of
     *                     term freq values in the collection. It is usually set to true only after the inverted index
     *                     is computed.
     */
    public static void calculateTermsUpperBound(boolean computeStats)
    {
        ArrayList<String> termsList;    // array list for all the term
        long startTime,endTime;         // variables to calculate the execution time
        long termCount = 0;             // counter for the number of term
        double currentTUB = 0;          // contain the Term Upper Bound for the current term
        boolean scoringFunc;            // user's choice about scoring function
        int minPLLen = 10000000;        // contains the len of the shortest posting list of the collection
        int maxPLLen = 0;               // contains the len of the longest posting list of the collection
        long sumPLLen = 0;              // contains the sum of the len of the all posting lists in the collection
        double avgPLLen = 0;            // contains the average len for the posting lists in the collection
        int currPLLen = 0;              // current posting list len (at each iteration)

        // check if already exist the file
        if (termUpperBoundFileExist())
        {
            deleteTermUpperBoundFile();         // delete the file
            printDebug("The termUpperBoundFile already exist.\nThe termUpperBoundFile erased.");    // control print
        }

        if (!CollectionStatistics.getTermFreqTable().isEmpty() && computeStats)
            CollectionStatistics.getTermFreqTable().clear();    // free the hash table, will be fill

        printDebug("Calculating terms upper bound..."); // control print
        startTime = System.currentTimeMillis();            // start time to calculate all term upper bound

        termsList = new ArrayList<>(QueryProcessor.getDictionary().keySet());   // read all the term of the dictionary
        scoringFunc = Flags.isScoringEnabled();             // take user's choice about using scoring function

        // scan all term in the dictionary
        for (String term : termsList)
        {
            currentTUB = QueryProcessor.maxScoreTerm(term,scoringFunc,computeStats); // calculate the term upper bound for the current term
            termUpperBoundTable.put(term,currentTUB);           // add term upper bound in the hashmap
            termCount++;        // update counter

            if (computeStats)       // compute the statistics related to the PostList len
            {
                currPLLen = QueryProcessor.getDictionary().get(term).getDf();   // get posting list len of the current term
                sumPLLen += currPLLen;
                if (currPLLen > maxPLLen)   // set max
                    maxPLLen = currPLLen;
                if (currPLLen < minPLLen)
                    minPLLen = currPLLen;   // set min
            }
        }
        endTime = System.currentTimeMillis();           // end time to calculate all term upper bound
        // shows term upper bound calculation time
        printTime("Calculated all term upper bound( " + termCount + " term) in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

        if (computeStats)
        {
            avgPLLen = (double) sumPLLen / termCount;           // calculate avgPLLen
            CollectionStatistics.setAvgPLLength(avgPLLen);          // set avgPLLen
            CollectionStatistics.setMinPLLength(minPLLen);          // set minPLLen
            CollectionStatistics.setMaxPLLength(maxPLLen);          // set maxPLLen
            CollectionStatistics.computeTermFreqOccStatistics();    // calculate TF occurrence statistics
            CollectionStatistics.storeCollectionStatsIntoDisk();    // save collection statistics
        }

        startTime = System.currentTimeMillis();         // start time to store all term upper bound
        storeTermUpperBoundTableIntoDisk();     // save the hashmap into disk
        endTime = System.currentTimeMillis();           // end time to store all term upper bound
        // shows term upper bound storing time
        printTime("Stored all term upper bound( " + termCount + " term) in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
    }

    /**
     * calculate the document upper bound for each document in the document table
     */
    public static void calculateDocsUpperBound()
    {

    }

    // ---------------- end: calculating functions ----------------


    // ---------------- start: control functions ----------------

    /**
     * function that return if the termUpperBoundTable in memory is set or not
     *
     * @return  true if the termUpperBoundTable is empty
     *          false if the termUpperBoundTable is not empty
     */
    public static boolean termUpperBoundTableIsEmpty()
    {
        return termUpperBoundTable.isEmpty();
    }

    /**
     * function that return if the docUpperBoundTable in memory is set or not
     *
     * @return  true if the docUpperBoundTable is empty
     *          false if the docUpperBoundTable is not empty
     */
    public static boolean docUpperBoundTableIsEmpty()
    {
        return docUpperBoundTable.isEmpty();
    }

    // method to free memory by deleting the information in termUpperBoundTable and docUpperBoundTable
    private static void freeMemory()
    {
        termUpperBoundTable.clear();
        docUpperBoundTable.clear();
    }

    /**
     * method to check if termUpperBoundFile exist or not
     *
     * @return  true if the file exist
     *          false if the file not exist
     */
    public static boolean termUpperBoundFileExist()
    {
        File termUBFile = new File(TERMUPPERBOUND_FILE);
        return termUBFile.exists();
    }

    /**
     * method to delete termUpperBoundFile (if exist)
     */
    private static void deleteTermUpperBoundFile()
    {
        File termUBFile = new File(TERMUPPERBOUND_FILE);
        if(termUBFile.exists())
        {
            try
            {
                FileUtils.delete(termUBFile);            // delete flags file if exist
            }
            catch (IOException e)
            {
                e.printStackTrace();
            }
        }
    }

    /**
     * method to obtain the termUpperBound from termUpperBoundTable of a term passed as parameter
     *
     * @param term          the term
     * @return the term upper bound for the term passed as parameter
     */
    public static double getTermUpperBound(String term)
    {
        if (!termUpperBoundTableIsEmpty())  // check if the hash map is not empty
        {
            if (termUpperBoundTable.get(term) != null)  // check if the term is in the hash map
                return termUpperBoundTable.get(term);       // return the corresponding term upper bound
            else
                return 0;   // return default value
        }
        else
            return 0;   // return default value
    }

    // ---------------- end: control functions ----------------

    // ---------------- start: read/write into disk functions ----------------

    /**
     * function to store the whole termUpperBoundTable into disk
     */
    static void storeTermUpperBoundTableIntoDisk()
    {

        printLoad("Storing terms upper bound into disk...");

        try (RandomAccessFile raf = new RandomAccessFile(TERMUPPERBOUND_FILE, "rw");
             FileChannel channel = raf.getChannel())
        {
            MappedByteBuffer buffer = channel.map(FileChannel.MapMode.READ_WRITE, channel.size(), (long) DOUBLE_BYTES * termUpperBoundTable.size());
            // Buffer not created
            if(buffer == null)
                return;
            // scan all document elements of the termUpperBoundTable
            for(double termUpperBound: termUpperBoundTable.values())
            {
                buffer.putDouble(termUpperBound);       // write termUpperBoundTable into file
                if(debug)
                    printDebug("write the term upper bound: " + termUpperBound);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        //printDebug("Print the hashtable fot the TUB " + termUpperBoundTable);
    }

    /**
     * function to read the whole termUpperBoundTable into disk
     */
    public static void readTermUpperBoundTableFromDisk()
    {
        ArrayList<String> termsList;    // array list for all the term
        double currentTermUB = 0;       // var to contain the term upper bound read from disk
        int count = 0;                  // counter for the position of the terms list
        long startTime,endTime;         // variables to calculate the execution time

        printLoad("Loading terms upper bound into disk...");

        termsList = new ArrayList<>(QueryProcessor.getDictionary().keySet());   // retrieve all the term of the dictionary
        /*
        printDebug("In term upper bound table (before the read from file).");
        printDebug("Position 0 -> term: " + termsList.get(0) + " and TUB: " + termUpperBoundTable.get(termsList.get(0)));
        printDebug("Position 1 -> term: " + termsList.get(1) + " and TUB: " + termUpperBoundTable.get(termsList.get(1)));
        printDebug("Position " + (termsList.size()-3) + " -> term: " + termsList.get(termsList.size()-3) + " and TUB: " + termUpperBoundTable.get(termsList.get(termsList.size()-3)));
        printDebug("Position " + (termsList.size()-2) + " -> term: " + termsList.get(termsList.size()-2) + " and TUB: " + termUpperBoundTable.get(termsList.get(termsList.size()-2)));
        printDebug("Position " + (termsList.size()-1) + " -> term: " + termsList.get(termsList.size()-1) + " and TUB: " + termUpperBoundTable.get(termsList.get(termsList.size()-1)));
        //*/
        termUpperBoundTable.clear();        // free the hash map table

        startTime = System.currentTimeMillis();         // start time to store all term upper bound
        try (
                RandomAccessFile raf = new RandomAccessFile(TERMUPPERBOUND_FILE, "r");
                FileChannel channel = raf.getChannel()
        ) {
            MappedByteBuffer termUBBuffer = channel.map(FileChannel.MapMode.READ_ONLY, 0, channel.size());
            // size control check
            if (termsList.size() != (channel.size()/DOUBLE_BYTES) )
            {
                printError("The number of term upper bound stored in the file into disk is different from the number of terms in the dictionary.");
                return;     // exit from the method
            }
            // for to read all term upper bound stored into disk and put into termUpperBoundTable
            for (int i = 0; i < channel.size(); i += DOUBLE_BYTES)
            {
                currentTermUB = termUBBuffer.getDouble();       // retrieve the current term upper bound
                termUpperBoundTable.put(termsList.get(count), currentTermUB);   // put into hash map table
                count++;        // update the counter for the current term position in arraylist
            }
        } catch (IOException e)
        {
            e.printStackTrace();
        }
        /*
        printDebug("In term upper bound table (after the read from file).");
        printDebug("Position 0 -> term: " + termsList.get(0) + " and TUB: " + termUpperBoundTable.get(termsList.get(0)));
        printDebug("Position 1 -> term: " + termsList.get(1) + " and TUB: " + termUpperBoundTable.get(termsList.get(1)));
        printDebug("Position " + (termsList.size()-3) + " -> term: " + termsList.get(termsList.size()-3) + " and TUB: " + termUpperBoundTable.get(termsList.get(termsList.size()-3)));
        printDebug("Position " + (termsList.size()-2) + " -> term: " + termsList.get(termsList.size()-2) + " and TUB: " + termUpperBoundTable.get(termsList.get(termsList.size()-2)));
        printDebug("Position " + (termsList.size()-1) + " -> term: " + termsList.get(termsList.size()-1) + " and TUB: " + termUpperBoundTable.get(termsList.get(termsList.size()-1)));
        //*/
        endTime = System.currentTimeMillis();           // end time to store all term upper bound
        // shows term upper bound storing time
        printTime("*** Read all term upper bound( " + termsList.size() + " term) in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
    }
    // ---------------- end: read/write into disk functions ----------------
}
