package it.unipi.dii.aide.mircv.data_structures;

import org.apache.commons.io.FileUtils;
import java.io.*;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import it.unipi.dii.aide.mircv.QueryProcessor;
import static it.unipi.dii.aide.mircv.utils.Constants.*;

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

        // check if already exist the file
        if (termUpperBoundFileExist())
        {
            deleteTermUpperBoundFile();         // delete the file
            printDebug("The termUpperBoundFile already exist.\nThe termUpperBoundFile erased.");    // control print
        }

        if (!CollectionStatistics.getTermFreqTable().isEmpty() && computeStats)
            CollectionStatistics.getTermFreqTable().clear();    // free the hash table, will be fill

        startTime = System.currentTimeMillis();             // start time to calculate all term upper bound
        termsList = new ArrayList<>(QueryProcessor.getDictionary().keySet());   // read all the term of the dictionary
        scoringFunc = Flags.isScoringEnabled();             // take user's choice about using scoring function

        printLoad("Calculating terms upper bound for " + termsList.size() + " terms ...");   // control print
        // scan all term in the dictionary
        for (String term : termsList)
        {
            currentTUB = QueryProcessor.maxScoreTerm(term,scoringFunc,computeStats); // calculate the term upper bound for the current term
            termUpperBoundTable.put(term,currentTUB);       // add term upper bound in the hashmap
            termCount++;        // update counter
        }
        endTime = System.currentTimeMillis();           // end time to calculate all term upper bound
        printTime("Calculated all term upper bound( " + termCount + " term) in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

        if (computeStats)
        {
            CollectionStatistics.computeTermFreqOccStatistics();    // calculate TF occurrence statistics
            CollectionStatistics.storeCollectionStatsIntoDisk();    // save collection statistics
        }

        storeTermUpperBoundTableIntoDisk();     // save the hashmap into disk
    }

    /**
     * Calculate the document upper bound for each document in the document table. (for future implementation)
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

    /**
     * Method to free memory by deleting the information in termUpperBoundTable and docUpperBoundTable
     */
    private static void freeMemory()
    {
        termUpperBoundTable.clear();
        docUpperBoundTable.clear();
    }

    /**
     * Method to check if termUpperBoundFile exist or not
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
     * Method to delete termUpperBoundFile (if exist)
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
     * Method to obtain the termUpperBound from termUpperBoundTable of a term passed as parameter
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
     * Function to store the whole termUpperBoundTable into disk
     */
    public static void storeTermUpperBoundTableIntoDisk()
    {
        ArrayList<String> termsList;    // array list for all the term
        long startTime, endTime;        // variables to calculate the execution time

        printLoad("Storing terms upper bound into disk...");

        startTime = System.currentTimeMillis();  // start time to store all term upper bound
        try (FileOutputStream fos = new FileOutputStream(TERMUPPERBOUND_FILE);
             BufferedOutputStream bos = new BufferedOutputStream(fos);
             DataOutputStream dos = new DataOutputStream(bos))
        {
            // read all the term of the dictionary
            termsList = new ArrayList<>(QueryProcessor.getDictionary().keySet());
            Collections.sort(termsList);                                // order the dictionary term list

            // scan all document elements of the termUpperBoundTable
            for (String term : termsList)
            {
                double termUpperBound = termUpperBoundTable.get(term);  // get TUB related to current term
                dos.writeDouble(termUpperBound);                        // write termUpperBoundTable into file
                if (debug)
                    printDebug("write the term upper bound: " + termUpperBound);
            }
            printLoad("Stored " + termUpperBoundTable.size() + " terms upper bound into disk...");

        } catch (IOException e) {
            e.printStackTrace();
        }

        endTime = System.currentTimeMillis();       // end time to store all term upper bound
        printTime("*** Stored all term upper bound (" + termUpperBoundTable.size() + " terms) in " +
                (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
    }

    /**
     * Function to read the whole termUpperBoundTable into disk
     */
    public static void readTermUpperBoundTableFromDisk()
    {
        ArrayList<String> termsList;    // array list for all the term
        double currentTermUB = 0;       // var to contain the term upper bound read from disk
        int count = 0;                  // counter for the position of the terms list
        long startTime, endTime;        // variables to calculate the execution time

        printLoad("Loading terms upper bound from disk...");

        // retrieve all the term of the dictionary
        termsList = new ArrayList<>(QueryProcessor.getDictionary().keySet());
        Collections.sort(termsList);    // order the list
        termUpperBoundTable.clear();    // free the hash map table

        startTime = System.currentTimeMillis();  // start time to load all term upper bound
        try (FileInputStream fis = new FileInputStream(TERMUPPERBOUND_FILE);
             BufferedInputStream bis = new BufferedInputStream(fis);
             DataInputStream dis = new DataInputStream(bis))
        {
            // Size check to ensure that the file is consistent with the number of terms
            long expectedFileSize = (long) termsList.size() * DOUBLE_BYTES;
            if (expectedFileSize != fis.getChannel().size())
            {
                printError("The number of term upper bound entries in the file on disk is different from the number of terms in the dictionary.");
                return;     // exit from the method
            }

            // for to read all term upper bound stored into disk and put into termUpperBoundTable
            for (String term : termsList)
            {
                currentTermUB = dis.readDouble();               // retrieve the current term upper bound
                termUpperBoundTable.put(term, currentTermUB);   // put into hash map table
                count++;                            // update the counter for the current term position in arraylist
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        endTime = System.currentTimeMillis();  // end time to load all term upper bound
        printTime("*** Loaded all term upper bound (" + termsList.size() + " terms) in " +
                (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
    }
    // ---------------- end: read/write into disk functions ----------------
}
