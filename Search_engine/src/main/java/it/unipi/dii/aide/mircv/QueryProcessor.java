package it.unipi.dii.aide.mircv;

import org.apache.commons.compress.compressors.gzip.GzipCompressorInputStream;
import java.io.*;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.StandardCharsets;
import java.util.*;

import it.unipi.dii.aide.mircv.compression.Unary;
import it.unipi.dii.aide.mircv.compression.VariableBytes;
import it.unipi.dii.aide.mircv.data_structures.*;
import it.unipi.dii.aide.mircv.data_structures.Dictionary;
import it.unipi.dii.aide.mircv.utils.FileSystem;

import static java.lang.Math.log10;
import static java.lang.Math.min;

import static it.unipi.dii.aide.mircv.data_structures.CollectionStatistics.readCollectionStatsFromDisk;
import static it.unipi.dii.aide.mircv.data_structures.DataStructureHandler.readCompressedPostingListFromDisk;
import static it.unipi.dii.aide.mircv.data_structures.DataStructureHandler.readPostingListFromDisk;
import static it.unipi.dii.aide.mircv.data_structures.DataStructureHandler.readAndUncompressCompressedAndSkippedPLFromDisk;
import static it.unipi.dii.aide.mircv.data_structures.DataStructureHandler.readCompTFBlockFromDisk;
import static it.unipi.dii.aide.mircv.data_structures.DataStructureHandler.readCompDIDBlockFromDisk;
import static it.unipi.dii.aide.mircv.data_structures.Flags.readFlagsFromDisk;
import static it.unipi.dii.aide.mircv.utils.Constants.*;
import static it.unipi.dii.aide.mircv.utils.Constants.printDebug;

/**
 * Class to manage and execute query.
 */
public final class QueryProcessor
{
    // indicate whether order all or only first "numberOfResults" results from hash table. TEST VARIABLE
    private static boolean orderAllHashMap = false;
    public static HashMap<Integer, DocumentElement> documentTable = new HashMap<>();    // hash table DocID to related DocElement
    static it.unipi.dii.aide.mircv.data_structures.Dictionary dictionary = new Dictionary();    // dictionary in memory
    private static ArrayList<String> termNotInCollection = new ArrayList<>();   // ArrayList that contain the term that are in the query but not in the collection
    static PriorityQueue<QueryProcessor.ResultBlock> resPQ;     // priority queue for the result of scoring function for the best numberOfResults docs

    // var for when both the compression and the skipping are enabled
    // priority queue for the ordered block of DID and Term freq to execute DAAT when the compression and skipping are enabled
    static PriorityQueue<QueryProcessor.orderedDIDBlock> ordDIDPQ = new PriorityQueue<>(new CompareOrdDIDBlock());
    static ArrayList<Posting>[] skipAndCompPLs;  // contains all the posting lists for each term of the query (case of compression and skipping enabled)
    static PriorityQueue<Integer> pqDID = new PriorityQueue<Integer>();    // priority queue for the DID when use max score when compression and skipping are enabled
    static double[] logTermFreq;    // array to contain the precomputed TF score (1 + log(TermFeq)) for each possible value of TermFreq in this collection (117)
    static byte[][] termFreqPL;     // array for the compressed TermFreq list for each terms of the query
    static byte[][] docIDPL;        // array for the compressed TermFreq list for each terms of the query
    // ---------------------------------------------- start: functions -------------------------------------------------

    // ---------------------------------------- start: set and get functions -------------------------------------------

    /**
     * Function which sets one of the global posting lists used in the case of active compression and skipping with the
     * values given as parameters.
     *
     * @param pl    the posting list passed as parameter
     * @param index the index of the posting list to set
     */
    public static void setPLInSkipAndCompPLs(ArrayList<Posting> pl, int index)
    {
        if ((skipAndCompPLs == null) || (index >= skipAndCompPLs.length))
            return;

        skipAndCompPLs[index] = pl;  // set new uncompressed block of the posting list
    }

    /**
     * Function that returns one of the global posting lists used in the case of active compression and skipping.
     *
     * @param index the index of the posting list to get
     * @return  the requested posting list
     */
    public static ArrayList<Posting> getPLInSkipAndCompPLs(int index)
    {
        if ((skipAndCompPLs == null) || (index >= skipAndCompPLs.length))
            return null;

        return skipAndCompPLs[index];
    }

    /**
     * Function to return the dictionary loaded in memory
     */
    public static HashMap<String, DictionaryElem> getDictionary()
    {
        return dictionary.getTermToTermStat();
    }
    // ----------------------------------------- end: set and get functions --------------------------------------------

    /**
     * Function to manage the query request. Prepare and execute the query and return the results.
     *
     * @param query             is the query of the users (in words)
     * @param isConjunctive     indicates whether the query is conjunctive or disjunctive type (default is disjunctive)
     * @param numberOfResults   the number of results to be returned by the query
     * @return  an ArrayList of integer that representing an ordered list of DocIDs
     */
    public static ArrayList<Integer> queryManager(String query, boolean isConjunctive, int numberOfResults)
    {
        ArrayList<Integer> rankedResults = new ArrayList<>();   // ArrayList that contain the ranked results of query
        ArrayList<String> processedQuery;                       // array list for containing the query term

        try{
            // processed the query to obtain the term
            processedQuery = TextProcessor.preprocessText(query); // Preprocessing of document text

            // check if query is empty
            if (processedQuery.isEmpty() || (processedQuery.size() == 1 && processedQuery.get(0).isEmpty()))
            {
                printError("Error: query is empty, please retry.");     // mex of error
                return rankedResults;
            }

            // manage the query execution algorithm depending on the flags selected by the user
            chooseQueryAlg(processedQuery, isConjunctive, numberOfResults);

            rankedResults = getRankedResults(numberOfResults);          // get ranked results

            // check if array list of the terms in query which are not in the dictionary (collection) is empty
            if (!termNotInCollection.isEmpty())
            {
                printError("No match found for the following query terms: " + termNotInCollection);   // print the terms that are not in collection
                termNotInCollection.clear();        // clear ArrayList
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        return rankedResults;           // return the ranked results of query
    }

    /**
     * Function that checks whether all files and resources required to execute the query are available
     *
     * @return  true -> if all checks are passed,then can proceed with the execution of a query
     *          false -> if at least one check is failed (one file missed), then can't proceed with the execution of a query
     */
    public static boolean queryStartControl() throws IOException
    {
        // -- control for file into disk
        if (!FileSystem.areThereAllMergedFiles() ||
                !Flags.isThereFlagsFile() ||
                !CollectionStatistics.isThereStatsFile()) {
            printError("Error: missing required files.");
            return false;
        }

        readFlagsFromDisk();                                    // read flags from disk
        readCollectionStatsFromDisk();                          // read collection statistics from disk
        double avgDocLen = CollectionStatistics.getAvgDocLen(); // get average doc length
        // control check for the correctness of the BM25 calculation. If avgDocLen = 0 -> also all docLen = 0 -> the whole documents in the collection are empty
        if ( Flags.isScoringEnabled() && (avgDocLen == 0) )
        {
            printError("Error: the denominator for BM25 cannot be calculated and the same function cannot be used for scoring because the collection consists of documents that are all empty.");
            return false;       // exit to the function
        }

        // -- control for structures in memory - if not load them from disk
        if (!dictionary.dictionaryIsSet())
        {
            long startTime = System.currentTimeMillis();
            dictionary.readDictionaryFromDisk();
            long endTime = System.currentTimeMillis();
            printTime( "Dictionary loaded in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
        }
        else
            printDebug("Dictionary is already loaded.");

        if(documentTable.isEmpty())
        {
            long startTime = System.currentTimeMillis();
            DataStructureHandler.readDocumentTableFromDisk(false);
            long endTime = System.currentTimeMillis();
            printTime("Document Table loaded in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
        }
        else
            printDebug("Document Table is already loaded.");

        if ((logTermFreq == null) || (logTermFreq.length == 0))      // check for the precalculated value used to compute the score function
            readTFWeightFromDisk();
        else
            printDebug("termFreq array is already loaded.");

        return true;
    }

    /**
     * Function to choose the query execution algorithm depending on the flags selected by the user.
     *
     * @param processedQuery    array list for containing the query term
     * @param isConjunctive     indicates whether the query is conjunctive or disjunctive type (default is disjunctive)
     * @param numberOfResults   indicated the max number of result to return to user
     */
    private static void chooseQueryAlg(ArrayList<String> processedQuery, boolean isConjunctive, int numberOfResults) throws FileNotFoundException
    {
        // take user's choices that affecting the query execution
        boolean scoringFunc = Flags.isScoringEnabled();         // take user's choice about using scoring function
        boolean dynamicPrun = Flags.isDynamicPruningEnabled();  // take user's choice about using dynamic pruning algorithm
        boolean compression = Flags.isCompressionEnabled();     // take user's choice about using compression
        boolean skipping = Flags.considerSkippingBytes();       // take user's choice about using skipping
        boolean wholePLInMem = Flags.isWholePLInMemEnabled();   // take user's choice for whole PL in memory

        if (dynamicPrun)        // dynamic pruning 'true' -> use Max Score
        {
            if (skipping && compression && wholePLInMem)
                DAATAlgMAXSCORESkipAndCompAndWholePL(processedQuery,scoringFunc,isConjunctive,numberOfResults);   // apply DAAT + MaxScore (comp + skipping + wholePLInMem)  to calculate the score of the Docs
            else if (skipping && compression)
                DAATAlgMAXSCORESkipAndComp(processedQuery,scoringFunc,isConjunctive,numberOfResults);   // apply DAAT + MaxScore (comp+skipping)  to calculate the score of the Docs
            else if (skipping)
                DAATAlgMAXSCORESkipping(processedQuery,scoringFunc,isConjunctive,numberOfResults);  // apply DAAT + MaxScore (skipping) to calculate the score of the Docs
            else
                DAATAlgMAXSCORE(processedQuery,scoringFunc,isConjunctive,numberOfResults);  // apply DAAT + MaxScore (no skipping) to calculate the score of the Docs
        }
        else                    // dynamic pruning 'false' -> use DAAT alg
        {
            if (skipping && compression && wholePLInMem)
                DAATAlgCompSkipWholePL(processedQuery,scoringFunc,isConjunctive,numberOfResults);  // apply DAAT with compression, skipping and whole PL in memory
            else if (compression && skipping)
                DAATAlgCompSkip(processedQuery,scoringFunc,isConjunctive,numberOfResults);  // apply DAAT with compression and skipping
            else
                DAATAlgorithm(processedQuery,scoringFunc,isConjunctive,numberOfResults);    // apply DAAT to calculate the score of the Docs
        }
    }

    // ------------------------------------------- START - Execution Query alg -----------------------------------------

    /**
     * Function for apply the Document at a Time algorithm. It used when dynamic pruning is disabled and both
     * compression and skipping are not enabled at the same time.
     *
     * @param scoringFunc       indicates the preference for scoring. if false use TFIDF, if true use BM25.
     * @param processedQuery    array list for containing the query term
     * @param isConjunctive     indicates whether the query is conjunctive or disjunctive type (default is disjunctive)
     * @param numberOfResults   indicated the max number of result to return to user
     */
    private static void DAATAlgorithm(ArrayList<String> processedQuery, boolean scoringFunc , boolean isConjunctive, int numberOfResults) throws FileNotFoundException
    {
        PriorityQueue<QueryProcessor.orderedDIDBlock> postPQ = new PriorityQueue<>(new CompareOrdDIDBlock());
        resPQ = new PriorityQueue<>(numberOfResults, new CompareTerm());    // length equal to the number of results to be returned to the user
        ArrayList<Posting>[] postingLists;  // contains all the posting lists for each term of the query
        ArrayList<String> newProcQuery;     // new processed query after removing term not in the dictionary
        String[] terms;                     // array containing the terms of the query
        double[] IDFweight;                 // array containing the IDF weight for each posting list
        int[] lengthPostingList;            // array containing the length of the posting lists
        int[] postingListsIndex;            // contain the current position index for the posting list of each term in the query
        orderedDIDBlock currElem;           // the current (at each iteration) element poll from ordDIDPQ
        double partialScore = 0;            // var that contain partial score
        int docScoreCalc = 0;               // indicates the number of documents whose score was calculated (0 to number of results requested by the user)
        int pLNotEmpty = 0;                 // contains the number of posting lists related to the query terms that aren't empty
        int procQLen = 0;                   // the number of term in the original preprocessed query (before removing term not in the dictionary)
        double threshold = 0;               // the minimum score of the DID in the priority queue for the results
        int previousDID = 0;                // DID of the previous doc processed in algorithm
        int currentDID = 0;                 // DID of the current doc processed in algorithm
        int currIndexPL = 0;                // indicates the index of the PL from which the current element taken from the PQ originates
        int countSameDID = 0;               // indicates the number of block of PQ are processed with the same DID
        long startTime,endTime;             // variables to calculate the execution time

        // 0 - verify the term of the query and the number of term
        procQLen = processedQuery.size();               // take the size of the processed Query before term removing
        newProcQuery = removeTermNotInDictFromQuery(processedQuery);    // remove term that are not in the dictionary
        pLNotEmpty = newProcQuery.size();               // take the number of term that have Posting List
        terms = new String[pLNotEmpty];                 // set the size

        if((procQLen != pLNotEmpty) && (isConjunctive)) // there is at least one term not in dictionary, if conjunctive -> no results must be returned
            return;         // exit

        startTime = System.currentTimeMillis();         // start time for retrieve all posting lists of the query
        postingLists = retrieveAllPostListsFromQuery(newProcQuery);   // take all posting lists of query terms
        endTime = System.currentTimeMillis();           // end time for retrieve all posting lists of the query
        // shows query execution time
        printTime("*** Retrieved all posting lists in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

        // check the number of posting lists not empty and perform the best choice
        if (pLNotEmpty == 0)    // all terms in the query aren't in the dictionary or empty query
            return;         // exit
        else if (pLNotEmpty == 1)   // there is only 1 postingList (query with one term or query with more term but only one in dictionary)
        {
            DAATOnePostingList(newProcQuery.get(0), postingLists[0], scoringFunc, numberOfResults);   // execute DAAT algorithm
            return;         // exit
        }

        // 0.5 - more postingLists not empty, set of the utilities for DAAT
        postingListsIndex = getPostingListsIndex(postingLists);     // get the index initialized
        for (int i = 0; i < pLNotEmpty; i++)
            terms[i] = newProcQuery.get(i);
        lengthPostingList = retrieveLengthAllPostingLists(terms);   // take the length of each posting list
        IDFweight = calculateIDFWeight(lengthPostingList);          // calculate the IDF weight

        // 1 - take the first posting from the posting lists of the query terms
        for (int i = 0; i < pLNotEmpty; i++)                // add the block to PQ (postPQ)
        {
            postPQ.add(new QueryProcessor.orderedDIDBlock(postingLists[i].get(0).getDocId(),postingLists[i].get(0).getTermFreq(),i));
            ++postingListsIndex[i];         // update the index for the PL
        }

        // 2 - start DAAT
        startTime = System.currentTimeMillis();           // start time of DAAT (comp + skipping)
        previousDID = postPQ.peek().getDocID();
        if (isConjunctive)
        {   // -- start - if - conj --
            while (!postPQ.isEmpty())
            {   // -- start - while fot scan all block --
                currElem = postPQ.poll();               // get the current element
                currentDID = currElem.getDocID();       // get current DID
                currIndexPL = currElem.getIndexPL();    // get the current PL index

                if (currentDID == previousDID)      // sum the current partial score with the previous
                {
                    // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                    if (scoringFunc)
                        partialScore += ScoringBM25(currentDID, currElem.getTermFreq(), IDFweight[currIndexPL]);    // use BM25
                    else
                        partialScore += ScoringTFIDF(currElem.getTermFreq(), IDFweight[currIndexPL]);               // use TFIDF

                    ++countSameDID;     // update counter
                }
                else            // other DID, the previous partial score is final, save it or not.
                {
                    if ( (countSameDID == pLNotEmpty ) && (partialScore != 0))     // conjunctive and the current DID contain all the term in the query
                    {
                        // insert without control into priority queue (is not full)
                        if (docScoreCalc < numberOfResults)
                        {
                            resPQ.add(new QueryProcessor.ResultBlock(previousDID, partialScore));     // add to priority queue
                            ++docScoreCalc;         // increment result in priority queue counter
                            if (docScoreCalc == numberOfResults)
                                threshold = resPQ.peek().getScore();// set threshold value
                        }
                        else if (threshold < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                        {
                            // substitution of the block
                            resPQ.poll();       // remove the first element
                            resPQ.add(new QueryProcessor.ResultBlock(previousDID, partialScore));     // add to priority queue
                            threshold = resPQ.peek().getScore();// set threshold value
                        }
                    }
                    countSameDID = 1;   // reset counter
                    // calculate SCORE (TFIDF or BM25) for this term and currentDID
                    if (scoringFunc)
                        partialScore = ScoringBM25(currentDID, currElem.getTermFreq(), IDFweight[currIndexPL]);     // use BM25
                    else
                        partialScore = ScoringTFIDF(currElem.getTermFreq(), IDFweight[currIndexPL]);               // use TFIDF

                    previousDID = currentDID;       // update previous DID
                }
                // update the block in the pq -> add the next posting from the posting list from which the current block of the priority queue was taken
                if (postingListsIndex[currIndexPL] < lengthPostingList[currIndexPL])
                {
                    // take the next posting and put it into PQ (postPQ)
                    postPQ.add(new QueryProcessor.orderedDIDBlock(postingLists[currIndexPL].get(postingListsIndex[currIndexPL]).getDocId(),postingLists[currIndexPL].get(postingListsIndex[currIndexPL]).getTermFreq(),currIndexPL));
                    ++postingListsIndex[currIndexPL];         // update the index for the PL
                }
            }   // -- end - while fot scan all block --
            // conjunctive and the current DID contain all the term in the query
            if ((countSameDID == pLNotEmpty) && (partialScore != 0))
            {
                // insert without control into priority queue (is not full) or insert all results (orderAllHashMap = true)
                if (docScoreCalc < numberOfResults)
                {
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                }
                else if (threshold < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                {
                    // substitution of the block
                    resPQ.poll();       // remove the first element
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                }
            }
        }   // -- end - if - conj --
        else        // disjunctive
        {   // -- start - else - disj --
            while (!postPQ.isEmpty())
            {   // -- start - if - conj --
                currElem = postPQ.poll();             // get the current element
                currentDID = currElem.getDocID();       // get current DID
                currIndexPL = currElem.getIndexPL();    // get the current PL index

                if (currentDID == previousDID)      // sum the current partial score with the previous
                {
                    // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                    if (scoringFunc)
                        partialScore += ScoringBM25(currentDID, currElem.getTermFreq(), IDFweight[currIndexPL]);     // use BM25
                    else
                        partialScore += ScoringTFIDF(currElem.getTermFreq(), IDFweight[currIndexPL]);               // use TFIDF
                }
                else            // other DID, the previous partial score is final, save it.
                {
                    // save score
                    if (partialScore != 0)
                    {
                        // insert without control into priority queue (is not full) or insert all results (orderAllHashMap = true)
                        if (docScoreCalc < numberOfResults)
                        {
                            resPQ.add(new QueryProcessor.ResultBlock(previousDID, partialScore));     // add to priority queue
                            docScoreCalc++;         // increment result in priority queue counter
                            if (docScoreCalc == numberOfResults)
                                threshold = resPQ.peek().getScore();// set threshold value
                        }
                        else if (threshold < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                        {
                            // substitution of the block
                            resPQ.poll();       // remove the first element
                            resPQ.add(new QueryProcessor.ResultBlock(previousDID, partialScore));     // add to priority queue
                            threshold = resPQ.peek().getScore();// update threshold value
                        }
                    }

                    // calculate SCORE (TFIDF or BM25) for this term and currentDID
                    if (scoringFunc)
                        partialScore = ScoringBM25(currentDID, currElem.getTermFreq(), IDFweight[currIndexPL]);     // use BM25
                    else
                        partialScore = ScoringTFIDF(currElem.getTermFreq(), IDFweight[currIndexPL]);               // use TFIDF

                    previousDID = currentDID;       // update previous DID
                }
                // update the block in the pq -> add the next posting from the posting list from which the current block of the priority queue was taken
                if (postingListsIndex[currIndexPL] < lengthPostingList[currIndexPL])
                {
                    // take the next posting and put it into PQ (postPQ)
                    postPQ.add(new QueryProcessor.orderedDIDBlock(postingLists[currIndexPL].get(postingListsIndex[currIndexPL]).getDocId(),postingLists[currIndexPL].get(postingListsIndex[currIndexPL]).getTermFreq(),currIndexPL));
                    ++postingListsIndex[currIndexPL];         // update the index for the PL
                }
            }   // -- end - while fot scan all block --

            // save the last partial score
            if ((partialScore != 0))
            {
                // insert without control into priority queue (is not full) or insert all results (orderAllHashMap = true)
                if (docScoreCalc < numberOfResults)
                {
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                }
                else if (threshold < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                {
                    // substitution of the block
                    resPQ.poll();       // remove the first element
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                }
            }
        }   // -- end - else - disj --
        endTime = System.currentTimeMillis();           // end time of DAAT (comp + skipping)
        printTime("*** DAAT (PQ version) execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
    }

    /**
     * Function to apply the Document at a Time algorithm when both skipping and compression are enabled. In this method
     * you take the compressed lists and run the DAAT by decompressing only the necessary block from time to time.
     * At the beginning takes the first block for each PL, decompress it and create the element to be placed in the
     * priority queue (ordDIDPQ). When this is done you will extract the elements of the PQ and calculate the score.
     * When there are no more elements inherent at a PL in the PQ it will decompress the next block of it (if any)
     * and insert the values into the PQ.
     * So on until the PQ is empty.
     *
     * @param scoringFunc       indicates the preference for scoring. if false use TFIDF, if true use BM25.
     * @param processedQuery    array list for containing the query term
     * @param isConjunctive     indicates whether the query is conjunctive or disjunctive type (default is disjunctive)
     * @param numberOfResults   indicated the max number of result to return to user
     */
    private static void DAATAlgCompSkip(ArrayList<String> processedQuery, boolean scoringFunc , boolean isConjunctive, int numberOfResults) throws FileNotFoundException
    {
        PriorityQueue<QueryProcessor.orderedDIDBlock> postPQ = new PriorityQueue<>(new CompareOrdDIDBlock());
        resPQ = new PriorityQueue<>(numberOfResults, new CompareTerm());    // length equal to the number of results to be returned to the user
        ArrayList<Posting>[] postingLists;  // contains all the posting lists for each term of the query (only the current block)
        ArrayList<String> newProcQuery;     // new processed query after removing term not in the dictionary
        SkipList[] skipListArray;           // array of the Skip List reference related to the term of query
        String[] terms;                     // array containing the terms of the query
        int pLNotEmpty = 0;                 // contains the number of posting lists related to the query terms that aren't empty
        double[] IDFweight;                 // array containing the IDF weight for each posting list
        int[] lengthPostingList;            // array containing the length of the posting lists
        int[] blockIndex;                   // contain the index of the next block of skipList for each PL
        int[] elemOfBlockDone;              // contain the counter of the element belong to the same PL block
        int[] totalElemOfBlok;              // contain the total number of element of the same block in the PQ
        int[] numOfBlok;                    // contain the number of the block for each PL
        int previousDID = 0;                // DID of the previous doc processed in algorithm
        int currentDID = 0;                 // DID of the current doc processed in algorithm
        int currIndexPL = 0;                // indicates the index of the PL from which the current element taken from the PQ originates
        orderedDIDBlock currElem;           // the current (at each iteration) element poll from ordDIDPQ
        int docScoreCalc = 0;               // indicates the number of documents whose score was calculated (0 to number of results requested by the user)
        double partialScore = 0;            // var that contain partial score
        int procQLen = 0;                   // the number of term in the original preprocessed query (before removing term not in the dictionary)
        int countSameDID = 0;               // indicates the number of block of PQ are processed with the same DID
        double threshold = 0;               // the minimum score of the DID in the priority queue for the results
        long startTime,endTime;             // variables to calculate the execution time

        procQLen = processedQuery.size();               // take the size of the processed Query before term removing
        newProcQuery = removeTermNotInDictFromQuery(processedQuery);    // remove term that are not in the dictionary
        pLNotEmpty = newProcQuery.size();               // take the number of term that have Posting List
        terms = new String[pLNotEmpty];                 // set the size

        if((procQLen != pLNotEmpty) && (isConjunctive)) // there is at least one term not in dictionary, if conjunctive -> no results must be returned
            return;         // exit

        // create the skip List reference related to the term of query
        skipListArray = skipListInitCompAndSkip(newProcQuery, null, false);
        // check the number of posting lists not empty and perform the best choice
        if (pLNotEmpty == 0)        // all terms in the query aren't in the dictionary or empty query
            return;     // exit
        else if (pLNotEmpty == 1)   // there is only 1 postingList (query with one term or query with more term but only one in dictionary)
        {
            // The PL is only one -> read and decompress the whole PL and use the classic optimization method
            startTime = System.currentTimeMillis();         // start time for retrieve the posting lists of the query
            postingLists = retrieveAllUncompPL(newProcQuery, skipListArray); // get the uncompress PL
            endTime = System.currentTimeMillis();           // end time for retrieve the posting lists of the query
            // shows query execution time
            printTime("*** DAAT (comp+skipping) retrieved PL (case 1 PL) in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

            DAATOnePostingList(newProcQuery.get(0), postingLists[0], scoringFunc, numberOfResults);   // execute DAAT algorithm
            return;     // exit
        }
        // -- more postingLists not empty --
        // 0 - take the first block of the PLs and add to postPQ the first postings
        postingLists = new ArrayList[pLNotEmpty];
        for (int i = 0; i < pLNotEmpty; i++)
        {
            retrieveCompBlockOfPLAndUncompress(newProcQuery.get(i), skipListArray, 0, i, postingLists);
            postPQ.add(new QueryProcessor.orderedDIDBlock(postingLists[i].get(0).getDocId(),postingLists[i].get(0).getTermFreq(),i));
        }

        // 1) setting of all is need for the DAAT algorithm
        for (int i = 0; i < pLNotEmpty; i++)             // take the term of the query
            terms[i] = newProcQuery.get(i);
        lengthPostingList = retrieveLengthAllPostingLists(terms);   // take the length of each posting list
        IDFweight = calculateIDFWeight(lengthPostingList);          // calculate the IDF weight
        // set the index block for each PL
        blockIndex = new int[pLNotEmpty];
        Arrays.fill(blockIndex, 1);                 // the first block (with index 0) has already been taken
        // set the counter and max number of elem of the same block in the PQ -> used to manage the load of next PL block from disk
        elemOfBlockDone = new int[pLNotEmpty];
        Arrays.fill(elemOfBlockDone, 0);             // at the beginning the elem taken by PQ are equal to 0
        totalElemOfBlok = new int[pLNotEmpty];
        numOfBlok = new int[pLNotEmpty];
        for (int i = 0; i < pLNotEmpty; i++)            // take the term of the query
        {
            totalElemOfBlok[i] = min(SKIP_POINTERS_THRESHOLD, lengthPostingList[i]);
            numOfBlok[i] = skipListArray[i].getSkipArrLen();
        }

        // 2) start DAAT
        startTime = System.currentTimeMillis();           // start time of DAAT (comp + skipping)
        previousDID = postPQ.peek().getDocID();
        if (isConjunctive)
        {   // -- start - if - conj --
            while (!postPQ.isEmpty())
            {   // -- start - while fot scan all block --
                currElem = postPQ.poll();             // get the current element
                currentDID = currElem.getDocID();       // get current DID
                currIndexPL = currElem.getIndexPL();    // get the current PL index

                if (currentDID == previousDID)      // sum the current partial score with the previous
                {
                    // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                    if (scoringFunc)
                        partialScore += ScoringBM25(currentDID, currElem.getTermFreq(), IDFweight[currIndexPL]);     // use BM25
                    else
                        partialScore += ScoringTFIDF(currElem.getTermFreq(), IDFweight[currIndexPL]);               // use TFIDF

                    ++countSameDID;     // update counter
                }
                else            // other DID, the previous partial score is final, save it or not.
                {
                    if ( (countSameDID == pLNotEmpty ) && (partialScore != 0))     // conjunctive and the current DID contain all the term in the query
                    {
                        // insert without control into priority queue (is not full)
                        if (docScoreCalc < numberOfResults)
                        {
                            resPQ.add(new QueryProcessor.ResultBlock(previousDID, partialScore));     // add to priority queue
                            ++docScoreCalc;         // increment result in priority queue counter
                            if (docScoreCalc == numberOfResults)
                                threshold = resPQ.peek().getScore();// set threshold value
                        }
                        else if (threshold < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                        {
                            // substitution of the block
                            resPQ.poll();       // remove the first element
                            resPQ.add(new QueryProcessor.ResultBlock(previousDID, partialScore));     // add to priority queue
                            threshold = resPQ.peek().getScore();// set threshold value
                        }
                    }
                    countSameDID = 1;   // reset counter
                    // calculate SCORE (TFIDF or BM25) for this term and currentDID
                    if (scoringFunc)
                        partialScore = ScoringBM25(currentDID, currElem.getTermFreq(), IDFweight[currIndexPL]);     // use BM25
                    else
                        partialScore = ScoringTFIDF(currElem.getTermFreq(), IDFweight[currIndexPL]);               // use TFIDF

                    previousDID = currentDID;       // update previous DID
                }
                // update the counter -> when a counter is equal to maximum value of the block load the next block of the related PL
                if (totalElemOfBlok[currIndexPL] == ++elemOfBlockDone[currIndexPL])  // all elem of a block have been taken
                {
                    // load the new block and put elem into PQ
                    retrieveCompBlockOfPLAndUncompress(newProcQuery.get(currIndexPL), skipListArray, blockIndex[currIndexPL], currIndexPL, postingLists);
                    totalElemOfBlok[currIndexPL] = min(SKIP_POINTERS_THRESHOLD, (lengthPostingList[currIndexPL] - (SKIP_POINTERS_THRESHOLD * blockIndex[currIndexPL])));
                    ++blockIndex[currIndexPL];          // increment the counter of block
                    elemOfBlockDone[currIndexPL] = 0;   // reset the counter of elem
                }
                if (blockIndex[currIndexPL] <= numOfBlok[currIndexPL])
                    postPQ.add(new QueryProcessor.orderedDIDBlock(postingLists[currIndexPL].get(elemOfBlockDone[currIndexPL]).getDocId(),postingLists[currIndexPL].get(elemOfBlockDone[currIndexPL]).getTermFreq(),currIndexPL));
            }   // -- end - while fot scan all block --
            // conjunctive and the current DID doesn't contain all the term in the query
            if (countSameDID != pLNotEmpty)
                partialScore = 0;   // reset counter

            // save the last partial score
            if ((partialScore != 0))
            {
                // insert without control into priority queue (is not full) or insert all results (orderAllHashMap = true)
                if (docScoreCalc < numberOfResults)
                {
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                }
                else if (threshold < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                {
                    // substitution of the block
                    resPQ.poll();       // remove the first element
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                }
            }
        }   // -- end - if - conj --
        else        // disjunctive
        {   // -- start - else - disj --
            while (!postPQ.isEmpty())
            {   // -- start - if - conj --
                currElem = postPQ.poll();             // get the current element
                currentDID = currElem.getDocID();       // get current DID
                currIndexPL = currElem.getIndexPL();    // get the current PL index

                if (currentDID == previousDID)      // sum the current partial score with the previous
                {
                    // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                    if (scoringFunc)
                        partialScore += ScoringBM25(currentDID, currElem.getTermFreq(), IDFweight[currIndexPL]);     // use BM25
                    else
                        partialScore += ScoringTFIDF(currElem.getTermFreq(), IDFweight[currIndexPL]);               // use TFIDF
                }
                else            // other DID, the previous partial score is final, save it.
                {
                    // save score
                    if (partialScore != 0)
                    {
                        // insert without control into priority queue (is not full) or insert all results (orderAllHashMap = true)
                        if (docScoreCalc < numberOfResults)
                        {
                            resPQ.add(new QueryProcessor.ResultBlock(previousDID, partialScore));     // add to priority queue
                            docScoreCalc++;         // increment result in priority queue counter
                            if (docScoreCalc == numberOfResults)
                                threshold = resPQ.peek().getScore();// set threshold value
                        }
                        else if (threshold < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                        {
                            // substitution of the block
                            resPQ.poll();       // remove the first element
                            resPQ.add(new QueryProcessor.ResultBlock(previousDID, partialScore));     // add to priority queue
                            threshold = resPQ.peek().getScore();// update threshold value
                        }
                    }

                    // calculate SCORE (TFIDF or BM25) for this term and currentDID
                    if (scoringFunc)
                        partialScore = ScoringBM25(currentDID, currElem.getTermFreq(), IDFweight[currIndexPL]);     // use BM25
                    else
                        partialScore = ScoringTFIDF(currElem.getTermFreq(), IDFweight[currIndexPL]);               // use TFIDF

                    previousDID = currentDID;       // update previous DID
                }
                // update the counter -> when a counter is equal to maximum value of the block load the next block of the related PL
                if (totalElemOfBlok[currIndexPL] == ++elemOfBlockDone[currIndexPL])  // all elem of a block have been taken
                {
                    // load the new block and put elem into PQ
                    retrieveCompBlockOfPLAndUncompress(newProcQuery.get(currIndexPL), skipListArray, blockIndex[currIndexPL], currIndexPL, postingLists);
                    totalElemOfBlok[currIndexPL] = min(SKIP_POINTERS_THRESHOLD, (lengthPostingList[currIndexPL] - (SKIP_POINTERS_THRESHOLD * blockIndex[currIndexPL])));
                    ++blockIndex[currIndexPL];          // increment the counter of block
                    elemOfBlockDone[currIndexPL] = 0;   // reset the counter of elem
                }
                if (blockIndex[currIndexPL] <= numOfBlok[currIndexPL])
                    postPQ.add(new QueryProcessor.orderedDIDBlock(postingLists[currIndexPL].get(elemOfBlockDone[currIndexPL]).getDocId(),postingLists[currIndexPL].get(elemOfBlockDone[currIndexPL]).getTermFreq(),currIndexPL));
            }   // -- end - while fot scan all block --

            // save the last partial score
            if ((partialScore != 0))
            {
                // insert without control into priority queue (is not full) or insert all results (orderAllHashMap = true)
                if (docScoreCalc < numberOfResults)
                {
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                }
                else if (threshold < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                {
                    // substitution of the block
                    resPQ.poll();       // remove the first element
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                }
            }
        }   // -- end - else - disj --
        endTime = System.currentTimeMillis();           // end time of DAAT (comp + skipping)
        printTime("*** DAAT (comp+skipping) V2 execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
    }

    /**
     * Function to apply the Document at a Time algorithm when compression, skipping and wholePLinMem are enabled.
     * In this method we perform the same algorithm as in the case of compression and skipping both enabled with the
     * difference that the compressed PLs will not be loaded from memory in blocks but whole.
     * Then they will be decompressed and processed in blocks.
     * In this case it has the entire PL compressed in memory and current block decompressed for each PL.
     * This brings a saving in memory used compared to the base case (whole PLs decompressed in memory) but worse than
     * the case with skipping and compression enabled (only the current blocks decompressed in memory).
     * This methodology may be useful in the case of small collections (not large enough to occupy memory entirely with
     * query PLs) because it will occupy less memory than the base case but will have less overhead due to the mechanisms
     * in the compression and skipping enabled case.
     * The latter case may be too algorithm-intensive if memory is more than sufficient to handle the base case, leading
     * to performance degradation rather than improvement. in such cases, this method in the middle may be a good
     * compromise between memory usage, advanced techniques and performance.
     *
     * @param scoringFunc       indicates the preference for scoring. if false use TFIDF, if true use BM25.
     * @param processedQuery    array list for containing the query term
     * @param isConjunctive     indicates whether the query is conjunctive or disjunctive type (default is disjunctive)
     * @param numberOfResults   indicated the max number of result to return to user
     */
    private static void DAATAlgCompSkipWholePL(ArrayList<String> processedQuery, boolean scoringFunc , boolean isConjunctive, int numberOfResults) throws FileNotFoundException
    {
        PriorityQueue<QueryProcessor.orderedDIDBlock> postPQ = new PriorityQueue<>(new CompareOrdDIDBlock());
        resPQ = new PriorityQueue<>(numberOfResults, new CompareTerm());    // length equal to the number of results to be returned to the user
        ArrayList<Posting>[] postingLists;  // contains all the posting lists for each term of the query
        ArrayList<String> newProcQuery;     // new processed query after removing term not in the dictionary
        SkipList[] skipListArray;           // array of the Skip List reference related to the term of query
        String[] terms;                     // array containing the terms of the query
        int pLNotEmpty = 0;                 // contains the number of posting lists related to the query terms that aren't empty
        double[] IDFweight;                 // array containing the IDF weight for each posting list
        int[] lengthPostingList;            // array containing the length of the posting lists
        int[] blockIndex;                   // contain the index of the next block of skipList for each PL
        int[] elemOfBlockDone;              // contain the counter of the element belong to the same PL block
        int[] totalElemOfBlok;              // contain the total number of element of the same block in the PQ
        int[] numOfBlok;                    // contain the number of the block for each PL
        int previousDID = 0;                // DID of the previous doc processed in algorithm
        int currentDID = 0;                 // DID of the current doc processed in algorithm
        int currIndexPL = 0;                // indicates the index of the PL from which the current element taken from the PQ originates
        orderedDIDBlock currElem;           // the current (at each iteration) element poll from ordDIDPQ
        int docScoreCalc = 0;               // indicates the number of documents whose score was calculated (0 to number of results requested by the user)
        double partialScore = 0;            // var that contain partial score
        int procQLen = 0;                   // the number of term in the original preprocessed query (before removing term not in the dictionary)
        int countSameDID = 0;               // indicates the number of block of PQ are processed with the same DID
        double threshold = 0;               // the minimum score of the DID in the priority queue for the results
        long startTime,endTime;             // variables to calculate the execution time

        procQLen = processedQuery.size();               // take the size of the processed Query before term removing
        newProcQuery = removeTermNotInDictFromQuery(processedQuery);    // remove term that are not in the dictionary
        pLNotEmpty = newProcQuery.size();               // take the number of term that have Posting List
        terms = new String[pLNotEmpty];                 // set the size

        if((procQLen != pLNotEmpty) && (isConjunctive)) // there is at least one term not in dictionary, if conjunctive -> no results must be returned
            return;         // exit

        // create the skip List reference related to the term of query
        skipListArray = skipListInitCompAndSkip(newProcQuery, null, false);
        // check the number of posting lists not empty and perform the best choice
        if (pLNotEmpty == 0)        // all terms in the query aren't in the dictionary or empty query
            return;     // exit
        else if (pLNotEmpty == 1)   // there is only 1 postingList (query with one term or query with more term but only one in dictionary)
        {
            // The PL is only one -> read and decompress the whole PL and use the classic optimization method
            startTime = System.currentTimeMillis();         // start time for retrieve the posting lists of the query
            postingLists = retrieveAllUncompPL(newProcQuery, skipListArray); // get the uncompress PL
            endTime = System.currentTimeMillis();           // end time for retrieve the posting lists of the query
            // shows query execution time
            printTime("*** DAAT (comp+skipping+wholePLInMem) case 1 PL executed in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

            DAATOnePostingList(newProcQuery.get(0), postingLists[0], scoringFunc, numberOfResults);   // execute DAAT algorithm
            return;     // exit
        }

        // -- more postingLists not empty --
        // 0 - load all compressed PLs
        startTime = System.currentTimeMillis();         // start time for retrieve first block of PLs of the query
        retrieveAllCompPL(newProcQuery);                // load the whole PLs
        endTime = System.currentTimeMillis();           // end time for retrieve first block of PLs of the query
        printTime("*** Load the whole compressed PLs in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
        // 0.5 - take the first block of the PLs
        postingLists = new ArrayList[pLNotEmpty];
        for (int i = 0; i < pLNotEmpty; i++)
        {
            UncompressOneBlockFromCompressedPL(newProcQuery, skipListArray, 0, i,postingLists);
            postPQ.add(new QueryProcessor.orderedDIDBlock(postingLists[i].get(0).getDocId(),postingLists[i].get(0).getTermFreq(),i));
        }

        // 1 - setting of all is need for the DAAT algorithm
        for (int i = 0; i < pLNotEmpty; i++)             // take the term of the query
            terms[i] = newProcQuery.get(i);
        lengthPostingList = retrieveLengthAllPostingLists(terms);   // take the length of each posting list
        IDFweight = calculateIDFWeight(lengthPostingList);          // calculate the IDF weight
        // set the index block for each PL
        blockIndex = new int[pLNotEmpty];
        Arrays.fill(blockIndex, 1);                 // the first block (with index 0) has already been taken
        // set the counter and max number of elem of the same block in the PQ -> used to manage the load of next PL block from disk
        elemOfBlockDone = new int[pLNotEmpty];
        Arrays.fill(elemOfBlockDone, 0);            // at the beginning the elem taken by PQ are equal to 0
        totalElemOfBlok = new int[pLNotEmpty];
        numOfBlok = new int[pLNotEmpty];
        for (int i = 0; i < pLNotEmpty; i++)            // take the term of the query
        {
            totalElemOfBlok[i] = min(SKIP_POINTERS_THRESHOLD, lengthPostingList[i]);
            numOfBlok[i] = skipListArray[i].getSkipArrLen();
        }

        // 2 - start DAAT
        startTime = System.currentTimeMillis();           // start time of DAAT (comp + skipping)
        previousDID = postPQ.peek().getDocID();
        if (isConjunctive)
        {   // -- start - if - conj --
            while (!postPQ.isEmpty())
            {   // -- start - while fot scan all block --
                currElem = postPQ.poll();             // get the current element
                currentDID = currElem.getDocID();       // get current DID
                currIndexPL = currElem.getIndexPL();    // get the current PL index

                if (currentDID == previousDID)      // sum the current partial score with the previous
                {
                    // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                    if (scoringFunc)
                        partialScore += ScoringBM25(currentDID, currElem.getTermFreq(), IDFweight[currIndexPL]);     // use BM25
                    else
                        partialScore += ScoringTFIDF(currElem.getTermFreq(), IDFweight[currIndexPL]);               // use TFIDF

                    ++countSameDID;     // update counter
                }
                else            // other DID, the previous partial score is final, save it or not.
                {
                    if ( (countSameDID == pLNotEmpty ) && (partialScore != 0))     // conjunctive and the current DID contain all the term in the query
                    {
                        // insert without control into priority queue (is not full)
                        if (docScoreCalc < numberOfResults)
                        {
                            resPQ.add(new QueryProcessor.ResultBlock(previousDID, partialScore));     // add to priority queue
                            ++docScoreCalc;         // increment result in priority queue counter
                            if (docScoreCalc == numberOfResults)
                                threshold = resPQ.peek().getScore();// set threshold value
                        }
                        else if (threshold < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                        {
                            // substitution of the block
                            resPQ.poll();       // remove the first element
                            resPQ.add(new QueryProcessor.ResultBlock(previousDID, partialScore));     // add to priority queue
                            threshold = resPQ.peek().getScore();// set threshold value
                        }
                    }
                    countSameDID = 1;   // reset counter
                    // calculate SCORE (TFIDF or BM25) for this term and currentDID
                    if (scoringFunc)
                        partialScore = ScoringBM25(currentDID, currElem.getTermFreq(), IDFweight[currIndexPL]);     // use BM25
                    else
                        partialScore = ScoringTFIDF(currElem.getTermFreq(), IDFweight[currIndexPL]);               // use TFIDF

                    previousDID = currentDID;       // update previous DID
                }
                // update the counter -> when a counter is equal to maximum value of the block load the next block of the related PL
                if (totalElemOfBlok[currIndexPL] == ++elemOfBlockDone[currIndexPL])  // all elem of a block have been taken
                {
                    // load the new block and put elem into PQ
                    UncompressOneBlockFromCompressedPL(newProcQuery, skipListArray, blockIndex[currIndexPL], currIndexPL,postingLists);
                    totalElemOfBlok[currIndexPL] = min(SKIP_POINTERS_THRESHOLD, (lengthPostingList[currIndexPL] - (SKIP_POINTERS_THRESHOLD * blockIndex[currIndexPL])));
                    ++blockIndex[currIndexPL];          // increment the counter of block
                    elemOfBlockDone[currIndexPL] = 0;   // reset the counter of elem
                }
                if (blockIndex[currIndexPL] <= numOfBlok[currIndexPL])
                    postPQ.add(new QueryProcessor.orderedDIDBlock(postingLists[currIndexPL].get(elemOfBlockDone[currIndexPL]).getDocId(),postingLists[currIndexPL].get(elemOfBlockDone[currIndexPL]).getTermFreq(),currIndexPL));
            }   // -- end - while fot scan all block --
            // conjunctive and the current DID doesn't contain all the term in the query
            if (countSameDID != pLNotEmpty)
                partialScore = 0;   // reset counter

            // save the last partial score
            if ((partialScore != 0))
            {
                // insert without control into priority queue (is not full) or insert all results (orderAllHashMap = true)
                if (docScoreCalc < numberOfResults)
                {
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                }
                else if (threshold < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                {
                    // substitution of the block
                    resPQ.poll();       // remove the first element
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                }
            }
        }   // -- end - if - conj --
        else        // disjunctive
        {   // -- start - else - disj --
            while (!postPQ.isEmpty())
            {   // -- start - if - conj --
                currElem = postPQ.poll();             // get the current element
                currentDID = currElem.getDocID();       // get current DID
                currIndexPL = currElem.getIndexPL();    // get the current PL index

                if (currentDID == previousDID)      // sum the current partial score with the previous
                {
                    // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                    if (scoringFunc)
                        partialScore += ScoringBM25(currentDID, currElem.getTermFreq(), IDFweight[currIndexPL]);     // use BM25
                    else
                        partialScore += ScoringTFIDF(currElem.getTermFreq(), IDFweight[currIndexPL]);               // use TFIDF
                }
                else            // other DID, the previous partial score is final, save it.
                {
                    // save score
                    if (partialScore != 0)
                    {
                        // insert without control into priority queue (is not full) or insert all results (orderAllHashMap = true)
                        if (docScoreCalc < numberOfResults)
                        {
                            resPQ.add(new QueryProcessor.ResultBlock(previousDID, partialScore));     // add to priority queue
                            docScoreCalc++;         // increment result in priority queue counter
                            if (docScoreCalc == numberOfResults)
                                threshold = resPQ.peek().getScore();// set threshold value
                        }
                        else if (threshold < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                        {
                            // substitution of the block
                            resPQ.poll();       // remove the first element
                            resPQ.add(new QueryProcessor.ResultBlock(previousDID, partialScore));     // add to priority queue
                            threshold = resPQ.peek().getScore();// update threshold value
                        }
                    }

                    // calculate SCORE (TFIDF or BM25) for this term and currentDID
                    if (scoringFunc)
                        partialScore = ScoringBM25(currentDID, currElem.getTermFreq(), IDFweight[currIndexPL]);     // use BM25
                    else
                        partialScore = ScoringTFIDF(currElem.getTermFreq(), IDFweight[currIndexPL]);               // use TFIDF

                    previousDID = currentDID;       // update previous DID
                }
                // update the counter -> when a counter is equal to maximum value of the block load the next block of the related PL
                if (totalElemOfBlok[currIndexPL] == ++elemOfBlockDone[currIndexPL])  // all elem of a block have been taken
                {
                    // load the new block and put elem into PQ
                    UncompressOneBlockFromCompressedPL(newProcQuery, skipListArray, blockIndex[currIndexPL], currIndexPL,postingLists);
                    totalElemOfBlok[currIndexPL] = min(SKIP_POINTERS_THRESHOLD, (lengthPostingList[currIndexPL] - (SKIP_POINTERS_THRESHOLD * blockIndex[currIndexPL])));
                    ++blockIndex[currIndexPL];          // increment the counter of block
                    elemOfBlockDone[currIndexPL] = 0;   // reset the counter of elem
                }
                if (blockIndex[currIndexPL] <= numOfBlok[currIndexPL])
                    postPQ.add(new QueryProcessor.orderedDIDBlock(postingLists[currIndexPL].get(elemOfBlockDone[currIndexPL]).getDocId(),postingLists[currIndexPL].get(elemOfBlockDone[currIndexPL]).getTermFreq(),currIndexPL));
            }   // -- end - while fot scan all block --

            // save the last partial score
            if ((partialScore != 0))
            {
                // insert without control into priority queue (is not full) or insert all results (orderAllHashMap = true)
                if (docScoreCalc < numberOfResults)
                {
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                }
                else if (threshold < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                {
                    // substitution of the block
                    resPQ.poll();       // remove the first element
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                }
            }
        }   // -- end - else - disj --
        endTime = System.currentTimeMillis();           // end time of DAAT (comp + skipping)
        printTime("*** DAAT (comp+skipping+wholePLInMem) execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
    }

    /**
     * Function for apply the Document at a Time algorithm in the special case of only one posting list not empty
     *
     * @param postingLists      array list for containing the posting lists of the query's term
     * @param scoringFunc       indicates the preference for scoring. If false use TFIDF, if true use BM25.
     * @param processedQuery    array list for containing the query term
     * @param numberOfResults   indicated the max number of result to return to user
     */
    private static void DAATOnePostingList (String processedQuery, ArrayList<Posting> postingLists, boolean scoringFunc, int numberOfResults)
    {
        String[] terms = new String[1]; // array containing the terms of the query
        int[] lengthPostingList;    // array containing the length of the posting lists
        double[] IDFweight;         // array containing the IDF weight for each posting list
        double partialScore = 0;    // var that contain partial score
        double threshold = 0;       // the minimum score of the DID in the priority queue for the results
        int docScoreCalc = 0;       // indicates the number of documents whose score was calculated (0 to number of results requested by the user)
        long startTime,endTime;     // variables to calculate the execution time

        // initialize
        terms[0] = processedQuery;
        lengthPostingList = retrieveLengthAllPostingLists(terms);   // take the length of each posting list
        IDFweight = calculateIDFWeight(lengthPostingList);          // calculate the IDF weight

        startTime = System.currentTimeMillis();         // start time of DAAT
        // optimization -> there is one posting list -> the DID are already sort
        if (scoringFunc)
        {   // start - if - BM25 -
            for (Posting p : postingLists)
            {
                partialScore = ScoringBM25(p.getDocId(),p.getTermFreq(), IDFweight[0]); // calculate SCORE use BM25
                // save score
                if (partialScore != 0)
                {
                    // insert without control into priority queue (is not full) or insert all results (orderAllHashMap = true)
                    if ((docScoreCalc < numberOfResults) || orderAllHashMap)
                    {
                        resPQ.add(new QueryProcessor.ResultBlock(p.getDocId(), partialScore));     // add to priority queue
                        docScoreCalc++;         // increment result in priority queue counter
                        if (docScoreCalc == numberOfResults)
                            threshold = resPQ.peek().getScore();// set threshold value
                    }
                    else if (threshold < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                    {
                        // substitution of the block
                        resPQ.poll();           // remove the first element (old DID with lower score in the queue)
                        resPQ.add(new QueryProcessor.ResultBlock(p.getDocId(), partialScore));     // add to priority queue
                        threshold = resPQ.peek().getScore();    // set threshold value
                    }
                }
            }
        }   // end - if - BM25 -
        else
        {   // start - else - TFIDF -
            for (Posting p : postingLists)
            {
                partialScore = ScoringTFIDF(p.getTermFreq(), IDFweight[0]); // calculate SCORE  use TFIDF
                // save score
                if (partialScore != 0)
                {
                    // insert without control into priority queue (is not full) or insert all results (orderAllHashMap = true)
                    if ((docScoreCalc < numberOfResults) || orderAllHashMap)
                    {
                        resPQ.add(new QueryProcessor.ResultBlock(p.getDocId(), partialScore));  // add to priority queue
                        docScoreCalc++;         // increment result in priority queue counter
                        if (docScoreCalc == numberOfResults)
                            threshold = resPQ.peek().getScore();// set threshold value
                    }
                    else if (threshold < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                    {
                        // substitution of the block
                        resPQ.poll();           // remove the first element (old DID with lower score in the queue)
                        resPQ.add(new QueryProcessor.ResultBlock(p.getDocId(), partialScore));     // add to priority queue
                        threshold = resPQ.peek().getScore();    // set threshold value
                    }
                }
            }
        }   // end - else - TFIDF -
        endTime = System.currentTimeMillis();           // end time of DAAT
        // shows DAAT execution time
        printTime("*** DAAT V.0.5 (only 1 postingList) execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
    }

    /**
     * Function for apply the Document at a Time algorithm with WAND algorithm as dynamic pruning algorithm.
     * Not used in this current version, but can be implemented instead of the MaxScore except in the case of
     * compression and skipping both enabled.
     *
     * @param processedQuery    array list for containing the query term
     * @param scoringFunc       indicates the preference for scoring. if false use TFIDF, if true use BM25.
     * @param isConjunctive     indicates whether the query is conjunctive or disjunctive type (default is disjunctive)
     * @param numberOfResults   indicated the max number of result to return to user
     * @throws FileNotFoundException
     */
    private static void DAATAlgWAND(ArrayList<String> processedQuery, boolean scoringFunc , boolean isConjunctive, int numberOfResults) throws FileNotFoundException
    {
        resPQ = new PriorityQueue<>(numberOfResults, new CompareTerm());    // length equal to the number of results to be returned to the user
        ArrayList<String> newProcQuery;     // new processed query after removing term not in the dictionary
        ArrayList<Posting>[] postingLists;  // contains all the posting lists for each term of the query
        ArrayList<Integer> ordListDID;      // ordered list of the DocID present in the all posting lists of the term present in the query
        String[] terms;                     // array containing the terms of the query
        double[] IDFweight;                 // array containing the IDF weight for each posting list
        int[] lengthPostingList;            // array containing the length of the posting lists
        double[] termUpperBoundList;        // contains all the posting lists for each term of the query
        int[] postingListsIndex ;           // contain the current position index for the posting list of each term in the query
        int[] postingListsIndexWAND;        // contain the current position index for the posting list of each term in the query used by WAND algorithm
        Posting currentP;                   // support var
        double partialScore = 0;            // var that contain partial score
        double threshold = 0;               // var that contain the current threshold for WAND (is the minimum score value to be in the current best result)
        double tempSumTUB = 0;              // contains the sum of the term upper bound for the current doc
        int docScoreCalc = 0;               // indicates the number of documents whose score was calculated (0 to number of results requested by the user)
        boolean resetScore = false;         // used only in conjunctive case. indicates that the score must be set to 0 (the current Doc there aren't all the term of the query)
        int pLNotEmpty = 0;                 // contains the number of posting lists related to the query terms that aren't empty
        boolean calculateScore = true;      // indicate if of the current DID need to calculate the score or pass to the next
        boolean firstWAND = true;           // indicates if is the first iteration in WAND part of the algorithm
        int procQLen = 0;                   // the number of term in the original preprocessed query (before removing term not in the dictionary)
        long startTime,endTime;             // variables to calculate the execution time

        procQLen = processedQuery.size();               // take the size of the processed Query before term removing
        newProcQuery = removeTermNotInDictFromQuery(processedQuery);    // remove term that are not in the dictionary
        pLNotEmpty = newProcQuery.size();               // take the number of term that have Posting List
        terms = new String[pLNotEmpty];                 // set the size

        if((procQLen != pLNotEmpty) && (isConjunctive)) // there is at least one term not in dictionary, if conjunctive -> no results must be returned
            return;         // exit

        startTime = System.currentTimeMillis();         // start time for retrieve all posting lists of the query
        postingLists = retrieveAllPostListsFromQuery(newProcQuery);   // take all posting lists of query terms
        endTime = System.currentTimeMillis();           // end time for retrieve all posting lists of the query
        printTime("\n*** Retrieved all posting lists in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

        // check the number of posting lists not empty and perform the best choice
        if (pLNotEmpty == 0)    // all terms in the query aren't in the dictionary or empty query
            return;     // exit
        else if (pLNotEmpty == 1)   // there is only 1 postingList (query with one term or query with more term but only one in dictionary)
        {
            DAATOnePostingList(newProcQuery.get(0), postingLists[0], scoringFunc, numberOfResults);   // execute DAAT algorithm
            return;     // exit
        }

        // more postingLists not empty, use WAND algorithm
        ordListDID = DIDOrderedListOfQuery(postingLists, isConjunctive);    // take ordered list of DocID
        postingListsIndex = getPostingListsIndex(postingLists);             // get the index initialized
        postingListsIndexWAND = getPostingListsIndex(postingLists);         // get the WAND index initialized
        termUpperBoundList = getPostingListsTermUpperBound(newProcQuery);   // get the term upper bound for each term(postinglist)
        // initialize array for improvement TFIDF and BM25 scoring
        for (int i = 0; i < pLNotEmpty; i++)                     // get query terms
            terms[i] = newProcQuery.get(i);
        lengthPostingList = retrieveLengthAllPostingLists(terms);           // take the length of each posting list
        IDFweight = calculateIDFWeight(lengthPostingList);                  // calculate the IDF weight

        startTime = System.currentTimeMillis();           // start time of DAAT
        // scan all Doc retrieved and calculate score (TFIDF or BM25)
        for (int currentDID : ordListDID)
        {   // -- start - while to scan all doc retrieved --
            partialScore = 0;           // reset var
            resetScore = false;         // set to false
            tempSumTUB = 0;             // set to 0
            calculateScore = true;      // set to true

            // until the pq is not full has to calculate all score, after it uses the WAND algorithm for optimization
            if(docScoreCalc >= numberOfResults)
            {   // -- start - if.0 - WAND execution --
                // WAND algorithm
                if (firstWAND)
                {
                    // update the WAND index with the values of the DAAT index
                    System.arraycopy(postingListsIndex, 0, postingListsIndexWAND, 0, postingLists.length);
                    firstWAND = false;  // reset
                }
                for (int j = 0; j < postingLists.length; j++)   // take and sum all term upper bound for the currentDID
                {
                    // check if the posting lists of j-th isn't at the end AND if the j-th term of the query is present in the doc identify by currentDID
                    if ( (postingLists[j] != null) && (postingListsIndexWAND[j] < postingLists[j].size()) && (postingLists[j].get(postingListsIndexWAND[j]).getDocId() == currentDID))
                    {
                        tempSumTUB += termUpperBoundList[j];    // update the sum of the term upper bound for the currentDID
                        postingListsIndexWAND[j]++;             // update WAND index of current value
                        //printDebug("------ WAND - find doc: " + currentDID + " in PL: " + j + " in pos: " + (postingListsIndexWAND[j]-1) + " with sumTUB: " + tempSumTUB );
                    }
                    else if (isConjunctive) // must take only the document in which there are all term (DID that compare in all posting lists of the terms)
                    {
                        // if all postings in one posting lists have already been seen the next documents in the posting lists cannot contain all the terms in the query
                        if ((postingLists[j] == null) || (postingListsIndexWAND[j] >= postingLists[j].size()))
                        {
                            //printDebug("Query conjunctive, posting list numero: " + j + " finita. Si è in pos: " + postingListsIndex[j] + " su dimensione: " + postingLists[j].size());
                            endTime = System.currentTimeMillis();           // end time of DAAT
                            // shows query execution time
                            printTime("*** DAAT + WAND V.1.5 execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
                            return;             // exit from function
                        }
                        calculateScore = false;       // reset the partial score
                    }
                }

                //printDebug("-- WAND - find doc: " + currentDID + " with sumTUB: " + tempSumTUB + " and threshold: " + threshold );
                // check tempSumTUB if is greater than threshold
                if ( (tempSumTUB <= threshold) || !calculateScore)
                {
                    // update the DAAT index with the values of the WAND index
                    System.arraycopy(postingListsIndexWAND, 0, postingListsIndex, 0, postingLists.length);
                    continue;       // pass to next DID
                }
            }   // -- end - if.0 - WAND execution --

            //printDebug("-- WAND - start calculating - doc: " + currentDID + " with sumTUB: " + tempSumTUB + " and threshold: " + threshold );
            // calculate the score -> case of tempSumTUB > threshold or docScoreCalc < numberOfResults (default case is query Disjunctive)
            for (int j = 0; j < postingLists.length; j++)   // take all values and calculating the scores in the posting related to currentDID
            {
                // check if the posting lists of j-th isn't at the end AND if the j-th term of the query is present in the doc identify by currentDID
                if ( (postingLists[j] != null) && (postingListsIndex[j] < postingLists[j].size()) && (postingLists[j].get(postingListsIndex[j]).getDocId() == currentDID))
                {
                    currentP = postingLists[j].get(postingListsIndex[j]);   // take posting
                    postingListsIndex[j]++;                                 // update index of current value

                    // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                    if (scoringFunc)
                        partialScore += ScoringBM25(currentDID, currentP.getTermFreq(), IDFweight[j]);     // use BM25
                    else
                        partialScore += ScoringTFIDF(currentP.getTermFreq(), IDFweight[j]);     // use TFIDF
                    //printDebug("------ DAAT: find DID: " + currentDID + " term: '" + processedQuery.get(j) + "' in PL: " + j + " in pos: " + (postingListsIndex[j]-1) + " and partialScore: " + partialScore);
                }
                else if (isConjunctive)
                {
                    // must take only the document in which there are all term (DID that compare in all posting lists of the terms)
                    resetScore = true;       // reset the partial score

                    // if all postings in one posting lists have already been seen the next documents in the posting lists cannot contain all the terms in the query
                    if ((postingLists[j] == null) || (postingListsIndex[j] >= postingLists[j].size())) {
                        //printDebug("Query conjunctive, posting list numero: " + j + " finita. Si è in pos: " + postingListsIndex[j] + " su dimensione: " + postingLists[j].size());
                        endTime = System.currentTimeMillis();           // end time of DAAT
                        // shows query execution time
                        printTime("*** DAAT + WAND V.1.5 execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
                        return;             // exit from function
                    }
                }
            }
            // save score
            if ((partialScore != 0) && !resetScore)
            {
                // insert without control into priority queue (is not full) or insert all results (orderAllHashMap = true)
                if ((docScoreCalc < numberOfResults) || orderAllHashMap)
                {
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                    threshold = resPQ.peek().getScore();    // update threshold
                    docScoreCalc++;                         // increment result in priority queue counter
                }
                else if (threshold < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                {
                    // substitution of the block
                    resPQ.poll();                           // remove the first element
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                    threshold = resPQ.peek().getScore();    // update threshold
                    //printDebug("------ Scoring, new threshold: " + threshold);
                }
            }
        }   // -- end - while to scan all doc retrieved --
        endTime = System.currentTimeMillis();           // end time of DAAT
        // shows DAAT execution time
        printTime("*** DAAT + WAND V.1 execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
    }

    /**
     * Function for apply the Document at a Time algorithm with max Score algorithm as dynamic pruning algorithm.
     * This function apply the Max Score in the case of skipping disabled.
     *
     * @param processedQuery    array list for containing the query term
     * @param scoringFunc       indicates the preference for scoring. if false use TFIDF, if true use BM25.
     * @param isConjunctive     indicates whether the query is conjunctive or disjunctive type (default is disjunctive)
     * @param numberOfResults   indicated the max number of result to return to user
     * @throws FileNotFoundException
     */
    private static void DAATAlgMAXSCORE(ArrayList<String> processedQuery, boolean scoringFunc , boolean isConjunctive, int numberOfResults) throws FileNotFoundException
    {
        resPQ = new PriorityQueue<>(numberOfResults, new CompareTerm());    // length equal to the number of results to be returned to the user
        PriorityQueue<Integer> didPQ = new PriorityQueue<Integer>();    // priority queue for the DID
        ArrayList<Posting>[] postingLists;  // contains all the posting lists for each term of the query
        ArrayList<String> newProcQuery;     // new processed query after removing term not in the dictionary
        List<Posting> tempList;             // sublist of the posting for the search in Non-Essential posting list
        double[] IDFweight;                 // array containing the IDF weight for each posting list
        int[] lengthPostingList;            // array containing the length of the posting lists
        int[] postingListsIndex;            // contain the current position index for the posting list of each term in the query
        Posting currentP;                   // support var
        int firstEssPostListIndex = 0;      // indicates the index of the first (current) essential posting list
        double threshold = 0;               // var that contain the current threshold for MaxScore (is the minimum score value to be in the current best result)
        double currentDocUpperBound = 0;    // current document upper bound (used in max score algorithm for early stopping)
        double partialScore = 0;            // var that contain partial score
        int docScoreCalc = 0;               // indicates the number of documents whose score was calculated (0 to number of results requested by the user)
        boolean resetScore = false;         // used only in conjunctive case. indicates that the score must be set to 0 (the current Doc there aren't all the term of the query)
        int pLNotEmpty = 0;                 // contains the number of posting lists related to the query terms that aren't empty
        int range = 100;                    // value for the hop in the search of currentDID in the Non-Essential Posting lists
        int newIndex = 0;                   // contain the index of the current DID searched (if present) in the posting list (in Non-Essential posting lists)
        int procQLen = 0;                   // the number of term in the original preprocessed query (before removing term not in the dictionary)
        int postListCurrDID = 0;            // the DID in the current position of the posting list (used in no-essential PL)
        int currDID = 0;                    // the current value of DID
        long startTime,endTime;             // variables to calculate the execution time

        // 0 - verify the term of the query and the number of term
        procQLen = processedQuery.size();               // take the size of the processed Query before term removing
        newProcQuery = removeTermNotInDictFromQuery(processedQuery);    // remove term that are not in the dictionary
        pLNotEmpty = newProcQuery.size();               // take the number of term that have Posting List

        if((procQLen != pLNotEmpty) && isConjunctive) // there is at least one term not in dictionary, if conjunctive -> no results must be returned
            return;         // exit

        String[] orderedQueryTerm = new String[pLNotEmpty];     // contains the term query (which is in the dictionary) ordered by their TUB
        double[] termUpperBoundList  = new double[pLNotEmpty];  // contains all the term upper bound for each term of the query
        double[] sumTUBList  = new double[pLNotEmpty];          // array containing the sum of TUB, the value at the i-th position is the sum of TUBs from position 0 to (i-1)

        startTime = System.currentTimeMillis();         // start time for retrieve all posting lists of the query
        postingLists = retrieveAllPostingListsMaxScore(newProcQuery,orderedQueryTerm,termUpperBoundList,sumTUBList);   // take all posting lists of query terms
        endTime = System.currentTimeMillis();           // end time for retrieve all posting lists of the query
        // shows query execution time
        printTime("*** Retrieved all posting lists in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

        // check the number of posting lists not empty and perform the best choice
        if (pLNotEmpty == 0)        // all terms in the query aren't in the dictionary or empty query
            return;         // exit
        else if (pLNotEmpty == 1)   // there is only 1 postingList (query with one term or query with more term but only one in dictionary)
        {
            DAATOnePostingList(newProcQuery.get(0), postingLists[0], scoringFunc, numberOfResults);   // execute DAAT algorithm
            return;         // exit
        }

        // 0.5 - more postingLists not empty, set of the utilities for DAAT
        postingListsIndex = getPostingListsIndex(postingLists);                 // get the index initialized
        lengthPostingList = retrieveLengthAllPostingLists(orderedQueryTerm);    // take the length of each posting list
        IDFweight = calculateIDFWeight(lengthPostingList);                      // calculate the IDF weight

        // 1 - take the first posting from the posting lists of the query terms
        for (int i = 0; i < pLNotEmpty; i++)                // add the DID block to PQ (pqDID)
        {
            currDID = postingLists[i].get(0).getDocId();  // take DID in the first position
            if (!didPQ.contains(currDID))
                didPQ.add(currDID);                 // add to priority queue
        }

        // 2) start DAAT + MaxScore algorithm - scan all Doc retrieved and calculate score (TFIDF or BM25)
        int plsLen = pLNotEmpty;                // take the number of posting list
        int currTF = 0;                         // var for the current Term Freq
        int currPLIndex = 0;                    // var for the index of the current PL at the current iteration
        startTime = System.currentTimeMillis();           // start time of DAAT + MAX SCORE
        if (isConjunctive)
        {   // -- start - if - conj --
            while (!didPQ.isEmpty())
            {   // -- start - for 0: DID --
                currDID = didPQ.poll();     // take current DID
                partialScore = 0;           // reset var
                resetScore = false;         // set to false
                // 0 - scan the essential posting lists, default case is query Disjunctive
                for (int j = firstEssPostListIndex; j < plsLen; j++)
                {   // -- start - for 0.1: EPL --
                    currPLIndex = postingListsIndex[j];     // take the current index for the current posting list
                    // check if the posting lists of j-th isn't at the end AND if the j-th term of the query is present in the doc identify by currentDID
                    if (currPLIndex < lengthPostingList[j])
                    {
                        currentP = postingLists[j].get(currPLIndex);
                        if (currentP.getDocId() == currDID)
                        {
                            currTF = currentP.getTermFreq();   // take posting
                            postingListsIndex[j]++;                         // update index of current value
                            // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                            if (scoringFunc)
                                partialScore += ScoringBM25(currDID, currTF, IDFweight[j]);     // use BM25
                            else
                                partialScore += ScoringTFIDF(currTF, IDFweight[j]);               // use TFIDF
                            // in the current position of this PL there is the currentDID -> insert in PQ new DID (the next DID in this PL)
                            if (postingListsIndex[j] < lengthPostingList[j])
                            {
                                int tempDID = postingLists[j].get(postingListsIndex[j]).getDocId();
                                if (!didPQ.contains(tempDID))
                                    didPQ.add(tempDID);                 // add to priority queue
                            }
                        }
                        else    // must take only the document in which there are all term (DID that compare in all posting lists of the terms)
                            resetScore = true;       // reset the partial score
                    }
                    else
                    {
                        endTime = System.currentTimeMillis();           // end time of DAAT
                        printTime("*** DAAT + MAX SCORE (PQ) execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
                        return;             // exit from function
                    }
                }   // -- end - for 0.1: EPL --

                // 1 - scan non essential posting lists
                currentDocUpperBound = partialScore + sumTUBList[firstEssPostListIndex];    // calculate the current DUB
                // check if the doc has no zero possibility to have a score greater than threshold -- SEE NOTE 0 --
                if ( resetScore || (currentDocUpperBound <= threshold))
                    continue;                           // go to next iteration with next DID

                // scan non essential posting lists
                for (int i = 0; i < firstEssPostListIndex; i++)
                {   // -- start - for - NoEPL --
                    resetScore = true;       // reset the partial score
                    currPLIndex = postingListsIndex[i];     // take the current index for the current posting list
                    // check first position
                    if ( currPLIndex < lengthPostingList[i])
                    {
                        postListCurrDID = postingLists[i].get(currPLIndex).getDocId();
                        // check first position
                        if (postListCurrDID > currDID)   // this PL doesn't have the searched DID
                        {
                            currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                        }
                        else        // postListCurrDID <= currentDID
                        {
                            // jumps (of 'range') until it finds the searched DID or arrives in its vicinity
                            while (((postingListsIndex[i] + range) < lengthPostingList[i]) && (postingLists[i].get(postingListsIndex[i]+range).getDocId() <= currDID))
                            {
                                postingListsIndex[i] += range;         // update index of current value
                            }
                            currPLIndex = postingListsIndex[i];     // take the current index for the current posting list
                            // check the size of the posting list
                            if ((currPLIndex + range) < lengthPostingList[i])
                                tempList = postingLists[i].subList(currPLIndex, currPLIndex + range);
                            else
                                tempList = postingLists[i].subList(currPLIndex, lengthPostingList[i] - 1);

                            newIndex = booleanSearch(tempList, currDID);     // get the new index (boolean search)
                            if (newIndex == -1)     // the searched DID isn't in the posting list
                            {
                                postingListsIndex[i]++;                         // update the index with the new value
                                currentDocUpperBound -= termUpperBoundList[i];  // update currentDocUpperBound
                            }
                            else                    // the searched DID there is in the posting list
                            {
                                postingListsIndex[i] += newIndex;       // update the index with the new value
                                // find the searched DID - update the partialScore and the currentDocUpperBound
                                currTF = postingLists[i].get(postingListsIndex[i]).getTermFreq();   // take posting
                                currentDocUpperBound -= partialScore;               // update currentDocUpperBound
                                // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                                if (scoringFunc)
                                    partialScore += ScoringBM25(currDID, currTF, IDFweight[i]);   // use BM25
                                else
                                    partialScore += ScoringTFIDF(currTF, IDFweight[i]);             // use TFIDF
                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                currentDocUpperBound += partialScore;               // update currentDocUpperBound

                                resetScore = false;         // reset the partial score
                                postingListsIndex[i]++;     // update the index of the postin list
                            }
                        }
                        // check if the doc has no zero possibility to have a score greater than threshold
                        if ((currentDocUpperBound <= threshold) || resetScore)
                            break;
                    } // if all postings in one posting lists have already been seen the next documents in the posting lists cannot contain all the terms in the query
                    else
                    {
                        endTime = System.currentTimeMillis();           // end time of DAAT
                        printTime("*** DAAT + MAX SCORE (PQ) execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
                        return;             // exit from function
                    }
                }   // -- end - for - NoEPL --
                // 2 - save score
                // insert without control into priority queue (is not full)
                if ( (docScoreCalc < numberOfResults) && (partialScore != 0) )
                {
                    resPQ.add(new QueryProcessor.ResultBlock(currDID, partialScore));     // add to priority queue
                    docScoreCalc++;                         // increment result in priority queue counter
                    if (docScoreCalc == numberOfResults)
                        threshold = resPQ.peek().getScore();// update threshold
                }
                else if (threshold < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                {
                    // substitution of the block
                    resPQ.poll();                           // remove the first element
                    resPQ.add(new QueryProcessor.ResultBlock(currDID, partialScore));     // add to priority queue
                    threshold = resPQ.peek().getScore();    // update threshold
                    // calculate new essential posting lists and update firstEssPostListIndex
                    firstEssPostListIndex = updateEssentialPositngLists(sumTUBList, threshold, firstEssPostListIndex);
                    //printDebug("-- **** New threshold: " + threshold + " new first essential posting list: " + firstEssPostListIndex);
                }
            }   // -- end - for: DID --
        }   // -- end - if - conj --
        else    // disjunctive
        {   // -- start - else - disj --
            while (!didPQ.isEmpty())
            {   // -- start - for 0: DID --
                currDID = didPQ.poll();     // take current DID
                partialScore = 0;           // reset var
                // 0 - scan the essential posting lists, default case is query Disjunctive
                for (int j = firstEssPostListIndex; j < plsLen; j++)
                {   // -- start - for 0.1: EPL --
                    currPLIndex = postingListsIndex[j];     // take the current index for the current posting list
                    // check if the posting lists of j-th isn't at the end AND if the j-th term of the query is present in the doc identify by currentDID
                    if (currPLIndex < lengthPostingList[j])
                    {
                        currentP = postingLists[j].get(currPLIndex);
                        if (currentP.getDocId() == currDID)
                        {
                            currTF = currentP.getTermFreq();   // take posting
                            postingListsIndex[j]++;                         // update index of current value
                            // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                            if (scoringFunc)
                                partialScore += ScoringBM25(currDID,currTF, IDFweight[j]);     // use BM25
                            else
                                partialScore += ScoringTFIDF(currTF, IDFweight[j]);               // use TFIDF
                            // in the current position of this PL there is the currentDID -> insert in PQ new DID (the next DID in this PL)
                            if (postingListsIndex[j] < lengthPostingList[j])
                            {
                                int tempDID = postingLists[j].get(postingListsIndex[j]).getDocId();
                                if (!didPQ.contains(tempDID))
                                    didPQ.add(tempDID);                 // add to priority queue
                            }
                        }
                    }
                }   // -- end - for 0.1: EPL --

                // 1 - scan non essential posting lists
                currentDocUpperBound = partialScore + sumTUBList[firstEssPostListIndex];    // calculate the current DUB
                // check if the doc has no zero possibility to have a score greater than threshold -- SEE NOTE 0 --
                if (currentDocUpperBound <= threshold)
                    continue;                           // go to next iteration with next DID

                // scan non essential posting lists
                for (int i = 0; i < firstEssPostListIndex; i++)
                {   // -- start - for - NoEPL --
                    currPLIndex = postingListsIndex[i];     // take the current index for the current posting list
                    // check first position
                    if ( currPLIndex < lengthPostingList[i])
                    {
                        postListCurrDID = postingLists[i].get(currPLIndex).getDocId();
                        // check first position
                        if (postListCurrDID > currDID)   // this PL doesn't have the searched DID
                        {
                            currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                        }
                        else        // postListCurrDID <= currentDID
                        {
                            // jumps (of 'range') until it finds the searched DID or arrives in its vicinity
                            while (((postingListsIndex[i] + range) < lengthPostingList[i]) && (postingLists[i].get(postingListsIndex[i]+range).getDocId() <= currDID))
                            {
                                postingListsIndex[i] += range;         // update index of current value
                            }
                            currPLIndex = postingListsIndex[i];     // take the current index for the current posting list
                            // check the size of the posting list
                            if ((currPLIndex + range) < lengthPostingList[i])
                                tempList = postingLists[i].subList(currPLIndex, currPLIndex + range);
                            else
                                tempList = postingLists[i].subList(currPLIndex, lengthPostingList[i] - 1);

                            newIndex = booleanSearch(tempList, currDID);     // get the new index (boolean search)
                            if (newIndex == -1)     // the searched DID isn't in the posting list
                            {
                                postingListsIndex[i]++;                         // update the index with the new value
                                currentDocUpperBound -= termUpperBoundList[i];  // update currentDocUpperBound
                            }
                            else                    // the searched DID there is in the posting list
                            {
                                postingListsIndex[i] += newIndex;       // update the index with the new value
                                // find the searched DID - update the partialScore and the currentDocUpperBound
                                currTF = postingLists[i].get(postingListsIndex[i]).getTermFreq();   // take posting
                                currentDocUpperBound -= partialScore;               // update currentDocUpperBound
                                // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                                if (scoringFunc)
                                    partialScore += ScoringBM25(currDID, currTF, IDFweight[i]);   // use BM25
                                else
                                    partialScore += ScoringTFIDF(currTF, IDFweight[i]);             // use TFIDF
                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                currentDocUpperBound += partialScore;               // update currentDocUpperBound

                                postingListsIndex[i]++;     // update the index of the postin list
                            }
                        }
                        // check if the doc has no zero possibility to have a score greater than threshold
                        if (currentDocUpperBound <= threshold)
                            break;
                    } // if all postings in one posting lists have already been seen the next documents in the posting lists cannot contain all the terms in the query
                }   // -- end - for - NoEPL --
                // 2 - save score
                // insert without control into priority queue (is not full)
                if ( (docScoreCalc < numberOfResults) && (partialScore != 0) )
                {
                    resPQ.add(new QueryProcessor.ResultBlock(currDID, partialScore));     // add to priority queue
                    docScoreCalc++;                         // increment result in priority queue counter
                    if (docScoreCalc == numberOfResults)
                        threshold = resPQ.peek().getScore();// update threshold
                }
                else if (threshold < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                {
                    // substitution of the block
                    resPQ.poll();                           // remove the first element
                    resPQ.add(new QueryProcessor.ResultBlock(currDID, partialScore));     // add to priority queue
                    threshold = resPQ.peek().getScore();    // update threshold
                    // calculate new essential posting lists and update firstEssPostListIndex
                    firstEssPostListIndex = updateEssentialPositngLists(sumTUBList, threshold, firstEssPostListIndex);
                    //printDebug("-- **** New threshold: " + threshold + " new first essential posting list: " + firstEssPostListIndex);
                }
            }   // -- end - for: DID --
        }   // -- end - else - disj --
        endTime = System.currentTimeMillis();           // end time of DAAT
        // shows DAAT execution time
        printTime("*** DAAT + MAX SCORE (PQ) execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
    }

    /**
     * Function for apply the Document at a Time algorithm with max Score algorithm as dynamic pruning algorithm.
     * This function apply the Max Score in the case of skipping enabled and compression disabled.
     *
     * @param processedQuery    array list for containing the query term
     * @param scoringFunc       indicates the preference for scoring. if false use TFIDF, if true use BM25.
     * @param isConjunctive     indicates whether the query is conjunctive or disjunctive type (default is disjunctive)
     * @param numberOfResults   indicated the max number of result to return to user
     */
    private static void DAATAlgMAXSCORESkipping(ArrayList<String> processedQuery, boolean scoringFunc , boolean isConjunctive, int numberOfResults) throws FileNotFoundException
    {
        resPQ = new PriorityQueue<>(numberOfResults, new CompareTerm());    // length equal to the number of results to be returned to the user
        PriorityQueue<Integer> didPQ = new PriorityQueue<Integer>();    // priority queue for the DID
        String[] orderedQueryTerm;          // contains the term query (which is in the dictionary) ordered by their TUB
        double[] termUpperBoundList;        // contains all the term upper bound (TUB) for each term of the query
        double[] sumTUBList;                // array containing the sum of TUB, the value at the i-th position is the sum of TUBs from position 0 to (i-1)
        ArrayList<String> newProcQuery;     // new processed query after removing term not in the dictionary
        ArrayList<Posting>[] postingLists;  // contains all the posting lists for each term of the query
        double[] IDFweight;                 // array containing the IDF weight for each posting list
        int[] lengthPostingList;            // array containing the length of the posting lists
        int[] postingListsIndex;            // contain the current position index for the posting list of each term in the query
        SkipList[] skipListArray;           // array of the Skip List reference related to the term of query
        int firstEssPostListIndex = 0;      // indicates the index of the first (current) essential posting list
        double threshold = 0;               // var that contain the current threshold for MaxScore (is the minimum score value to be in the current best result)
        double currentDocUpperBound = 0;    // current document upper bound (used in max score algorithm for early stopping)
        double partialScore = 0;            // var that contain partial score
        int docScoreCalc = 0;               // indicates the number of documents whose score was calculated (0 to number of results requested by the user)
        boolean resetScore = false;         // used only in conjunctive case. indicates that the score must be set to 0 (the current Doc there aren't all the term of the query)
        int pLNotEmpty = 0;                 // contains the number of posting lists related to the query terms that aren't empty
        int postListCurrDID = 0;            // the DID in the current position of the posting list (used in no-essential PL)
        int procQLen = 0;                   // the number of term in the original preprocessed query (before removing term not in the dictionary)
        int currDID = 0;                    // the current value of DID
        long startTime,endTime;             // variables to calculate the execution time

        // 0 - verify the term of the query and the number of term
        procQLen = processedQuery.size();               // take the size of the processed Query before term removing
        newProcQuery = removeTermNotInDictFromQuery(processedQuery);    // remove term that are not in the dictionary
        pLNotEmpty = newProcQuery.size();               // take the number of term that have Posting List
        orderedQueryTerm = new String[pLNotEmpty];      // set len
        termUpperBoundList  = new double[pLNotEmpty];   // set len
        sumTUBList  = new double[pLNotEmpty];           // set len

        if((procQLen != pLNotEmpty) && (isConjunctive)) // there is at least one term not in dictionary, if conjunctive -> no results must be returned
            return;         // exit

        startTime = System.currentTimeMillis();         // start time for retrieve all posting lists of the query
        postingLists = retrieveAllPostingListsMaxScore(newProcQuery, orderedQueryTerm, termUpperBoundList, sumTUBList);   // take all posting lists of query terms
        endTime = System.currentTimeMillis();           // end time for retrieve all posting lists of the query
        // shows query execution time
        printTime("*** Retrieved all posting lists in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

        // 0.5 - more postingLists not empty, set of the utilities for DAAT
        if (pLNotEmpty == 0)    // all terms in the query aren't in the dictionary or empty query
            return;     // exit
        else if (pLNotEmpty == 1)   // there is only 1 postingList (query with one term or query with more term but only one in dictionary)
        {
            DAATOnePostingList(newProcQuery.get(0), postingLists[0], scoringFunc, numberOfResults);   // execute DAAT algorithm
            return;     // exit
        }

        // more postingLists not empty, use MAXSCORE algorithm
        postingListsIndex = getPostingListsIndex(postingLists);                 // get the index initialized
        lengthPostingList = retrieveLengthAllPostingLists(orderedQueryTerm);    // take the length of each posting list
        IDFweight = calculateIDFWeight(lengthPostingList);                      // calculate the IDF weight
        skipListArray = SetAllSkipList(orderedQueryTerm, postingLists);         // create the skip List reference related to the term of query

        int maxDID = CollectionStatistics.getNDocs();
        BitSet processedDocs = new BitSet(maxDID + 1);      // "+1" to include the max value

        // 1 - take the first posting from the posting lists of the query terms
        for (int i = 0; i < pLNotEmpty; i++)                // add the DID block to PQ (pqDID)
        {
            currDID = postingLists[i].get(0).getDocId();    // take DID in the first position
            if (!processedDocs.get(currDID))
            {
                didPQ.add(currDID);             // add the first DID in the PLs
                processedDocs.set(currDID);     // set the bit related to currDID in the bitset
            }
        }

        // 2 - start DAAT + MaxScore algorithm, scan all Doc retrieved and calculate score (TFIDF or BM25)
        startTime = System.currentTimeMillis();     // start time of DAAT + MAX SCORE
        double newFEPL = sumTUBList[1];         // if this value is exceeded with the threshold, the first essential posting list must be recalculated
        int plsLen = pLNotEmpty;                // take the number of posting list
        int currTF = 0;                         // var for the current Term Freq
        int currPLIndex = 0;                    // var for the index of the current PL at the current iteration
        if (isConjunctive)
        {   // -- start - if conjunctive --
            while (!didPQ.isEmpty())
            {   // -- start - for 0: DID --
                currDID = didPQ.poll();     // take current DID
                partialScore = 0;           // reset var
                resetScore = false;         // set to false
                // 0 - scan the essential posting lists, default case is query Disjunctive
                for (int j = firstEssPostListIndex; j < plsLen; j++)
                {   // -- start - for - EPL --
                    currPLIndex = postingListsIndex[j];     // take the current index for the current posting list
                    // check if the posting lists of j-th isn't at the end AND if the j-th term of the query is present in the doc identify by currentDID
                    if (currPLIndex < lengthPostingList[j])
                    {
                        if (postingLists[j].get(currPLIndex).getDocId() == currDID)
                        {
                            currTF = postingLists[j].get(postingListsIndex[j]).getTermFreq();   // take posting
                            postingListsIndex[j]++;                         // update index of current value
                            // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                            if (scoringFunc)
                                partialScore += ScoringBM25(currDID, currTF, IDFweight[j]);   // use BM25
                            else
                                partialScore += ScoringTFIDF(currTF, IDFweight[j]);             // use TFIDF

                            // in the current position of this PL there is the currentDID -> insert in PQ new DID (the next DID in this PL)
                            if (postingListsIndex[j] < lengthPostingList[j])
                            {
                                int tempDID = postingLists[j].get(postingListsIndex[j]).getDocId();
                                if (!processedDocs.get(tempDID))
                                {
                                    didPQ.add(tempDID);             // add the first DID in the PLs
                                    processedDocs.set(tempDID);     // set the bit related to currDID in the bitset
                                }
                            }
                        }
                        else
                            resetScore = true;       // reset the partial score
                    }
                    else     // must take only the document in which there are all term (DID that compare in all posting lists of the terms)
                    {
                        // if all postings in one posting lists have already been seen the next documents in the posting lists cannot contain all the terms in the query
                        endTime = System.currentTimeMillis();           // end time of DAAT
                        // shows query execution time
                        printTime("*** Max Score PQ(skipping)  execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
                        return;             // exit from function
                    }
                }   // -- end - for - EPL --

                // 1 - scan non essential posting lists
                currentDocUpperBound = partialScore + sumTUBList[firstEssPostListIndex];    // calculate the current DUB
                // check if the doc has no zero possibility to have a score greater than threshold
                if ( resetScore || (currentDocUpperBound <= threshold))
                    continue;                           // go to next iteration with next DID

                for (int i = 0; i < firstEssPostListIndex; i++)
                {   // -- start - for - NoEPL --
                    resetScore = true;       // reset the partial score
                    currPLIndex = postingListsIndex[i];     // take the current index for the current posting list

                    if ( currPLIndex < lengthPostingList[i])
                    {
                        postListCurrDID = postingLists[i].get(currPLIndex).getDocId();
                        // check first position
                        if (postListCurrDID == currDID)
                        {
                            // find the searched DID - update the partialScore and the currentDocUpperBound
                            currTF = postingLists[i].get(currPLIndex).getTermFreq();

                            currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                            currentDocUpperBound -= partialScore;               // update currentDocUpperBound
                            // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                            if (scoringFunc)
                                partialScore += ScoringBM25(currDID, currTF, IDFweight[i]);  // use BM25
                            else
                                partialScore += ScoringTFIDF(currTF, IDFweight[i]);             // use TFIDF
                            currentDocUpperBound += partialScore;               // update currentDocUpperBound

                            resetScore = false;         // reset the partial score
                            postingListsIndex[i]++;     // update the index of the postin list
                        }
                        else if (postListCurrDID > currDID)      // the searched DID is not in this posting list
                        {
                            currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                        }
                        else    // postListCurrDID < currentDID, search in the posting list
                        {
                            postingListsIndex[i] = skipListArray[i].nextGEQ(currDID, currPLIndex);
                            currPLIndex = postingListsIndex[i];
                            // check the index returned by nextGEQ
                            if (currPLIndex >= lengthPostingList[i])  // check for out of bound in case of reaching the end of the list
                            {
                                // the searched DID isn't in the posting list and the posting list is ended
                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                            }
                            else    // searched DID is in the posting list (founded by the boolean search with skipping)
                            {
                                postListCurrDID = postingLists[i].get(currPLIndex).getDocId(); // take the did

                                if (postListCurrDID != currDID)
                                {
                                    // should be always greater than currentDID
                                    currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                }
                                else                    // the searched DID there is in the posting list
                                {
                                    // find the searched DID - update the partialScore and the currentDocUpperBound
                                    currTF = postingLists[i].get(currPLIndex).getTermFreq();

                                    currentDocUpperBound -= partialScore;               // update currentDocUpperBound
                                    currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                    // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                                    if (scoringFunc)
                                        partialScore += ScoringBM25(currDID, currTF, IDFweight[i]);  // use BM25
                                    else
                                        partialScore += ScoringTFIDF(currTF, IDFweight[i]);             // use TFIDF
                                    currentDocUpperBound += partialScore;               // update currentDocUpperBound

                                    resetScore = false;         // reset the partial score
                                    postingListsIndex[i]++;     // update the index of the postin list
                                }
                            }
                        }
                        // check if the doc has no zero possibility to have a score greater than threshold
                        if ( resetScore || (currentDocUpperBound <= threshold) )
                            break;
                    } // if all postings in one posting lists have already been seen the next documents in the posting lists cannot contain all the terms in the query
                    else
                    {
                        endTime = System.currentTimeMillis();           // end time of DAAT
                        // shows query execution time
                        printTime("*** Max Score PQ(skipping) execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
                        return;             // exit from function
                    }
                }   // -- end - for - NoEPL --

                // 2 - save score
                // insert without control into priority queue (is not full)
                if ( (docScoreCalc < numberOfResults) && (partialScore != 0) )
                {
                    resPQ.add(new QueryProcessor.ResultBlock(currDID, partialScore));     // add to priority queue
                    docScoreCalc++;                         // increment result in priority queue counter
                    if (docScoreCalc == numberOfResults)
                        threshold = resPQ.peek().getScore();    // update threshold
                }
                else if (threshold < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                {
                    // substitution of the block
                    resPQ.poll();                           // remove the first element
                    resPQ.add(new QueryProcessor.ResultBlock(currDID, partialScore));     // add to priority queue
                    threshold = resPQ.peek().getScore();    // update threshold
                    // calculate new essential posting lists and update firstEssPostListIndex
                    //firstEssPostListIndex = updateEssentialPositngLists(sumTUBList, threshold, firstEssPostListIndex);
                    if ((threshold > newFEPL) && (firstEssPostListIndex != (plsLen - 1)))
                    {
                        firstEssPostListIndex = updateEssentialPositngLists(sumTUBList, threshold,firstEssPostListIndex);
                        newFEPL = sumTUBList[min(firstEssPostListIndex + 1, pLNotEmpty - 1)];
                    }
                }
            }   // -- end - for: DID --
        }   // -- end - if conjunctive --
        else        // disjunctive
        {   // -- start - else disjunctive --
            while (!didPQ.isEmpty())
            {   // -- start - for 0: DID --
                currDID = didPQ.poll();     // take current DID
                partialScore = 0;           // reset var
                // 0 - can the essential posting lists, default case is query Disjunctive
                for (int j = firstEssPostListIndex; j < plsLen; j++)
                {   // -- start - for - EPL --
                    currPLIndex = postingListsIndex[j];     // take the current index for the current posting list

                    // check if the posting lists of j-th isn't at the end AND if the j-th term of the query is present in the doc identify by currentDID
                    if ( (currPLIndex < lengthPostingList[j]) && (postingLists[j].get(currPLIndex).getDocId() == currDID))
                    {
                        currTF = postingLists[j].get(currPLIndex).getTermFreq();
                        postingListsIndex[j]++;                         // update index of current value
                        // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                        if (scoringFunc)
                            partialScore += ScoringBM25(currDID, currTF , IDFweight[j]);     // use BM25
                        else
                            partialScore += ScoringTFIDF(currTF, IDFweight[j]);               // use TFIDF

                        // in the current position of this PL there is the currentDID -> insert in PQ new DID (the next DID in this PL)
                        if (postingListsIndex[j] < lengthPostingList[j])
                        {
                            int tempDID = postingLists[j].get(postingListsIndex[j]).getDocId();
                            if (!processedDocs.get(tempDID))
                            {
                                didPQ.add(tempDID);             // add the first DID in the PLs
                                processedDocs.set(tempDID);     // set the bit related to currDID in the bitset
                            }
                        }
                    }
                }   // -- end - for - EPL --

                // 1 - scan non essential posting lists
                currentDocUpperBound = partialScore + sumTUBList[firstEssPostListIndex];    // calculate the current DUB

                // update the score
                for (int i = 0; i < firstEssPostListIndex; i++)
                {   // -- start - for - NoEPL --
                    // check if the doc has no zero possibility to have a score greater than threshold
                    if (currentDocUpperBound <= threshold)
                        break;                              // go to next iteration with next DID

                    currPLIndex = postingListsIndex[i];     // take the current index for the current posting list
                    if (currPLIndex < lengthPostingList[i]) // check if the index is not out of posting list bound
                    {
                        postListCurrDID = postingLists[i].get(currPLIndex).getDocId();
                        // check first position
                        if (postListCurrDID == currDID)
                        {
                            // find the searched DID - update the partialScore and the currentDocUpperBound
                            currTF = postingLists[i].get(currPLIndex).getTermFreq();

                            currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                            currentDocUpperBound -= partialScore;               // update currentDocUpperBound
                            // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                            if (scoringFunc)
                                partialScore += ScoringBM25(currDID, currTF, IDFweight[i]);   // use BM25
                            else
                                partialScore += ScoringTFIDF(currTF, IDFweight[i]);             // use TFIDF
                            currentDocUpperBound += partialScore;               // update currentDocUpperBound

                            postingListsIndex[i]++;     // update the index of the postin list
                        }
                        else if (postListCurrDID > currDID)      // the searched DID is not in this posting list
                        {
                            currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                        }
                        else    // postListCurrDID < currentDID, search in the posting list
                        {
                            postingListsIndex[i] = skipListArray[i].nextGEQ(currDID, postingListsIndex[i]);
                            currPLIndex = postingListsIndex[i];     // take the current index for the current posting list

                            // check the index returned by nextGEQ
                            if (currPLIndex >= lengthPostingList[i])  // check for out of bound in case of reaching the end of the list
                            {
                                // the searched DID isn't in the posting list and the posting list is ended
                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                            }
                            else    // searched DID is in the posting list (founded by the boolean search with skipping)
                            {
                                postListCurrDID = postingLists[i].get(currPLIndex).getDocId(); // take the did

                                if (postListCurrDID != currDID)
                                {
                                    // postListCurrDID should be always greater than currentDID
                                    currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                }
                                else                    // the searched DID there is in the posting list
                                {
                                    // find the searched DID - update the partialScore and the currentDocUpperBound
                                    currTF = postingLists[i].get(currPLIndex).getTermFreq();
                                    currentDocUpperBound -= partialScore;               // update currentDocUpperBound
                                    currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                    // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                                    if (scoringFunc)
                                        partialScore += ScoringBM25(currDID, currTF, IDFweight[i]);   // use BM25
                                    else
                                        partialScore += ScoringTFIDF(currTF, IDFweight[i]);             // use TFIDF
                                    currentDocUpperBound += partialScore;               // update currentDocUpperBound

                                    postingListsIndex[i]++;     // update the index of the postin list
                                }
                            }
                        }
                    }
                }   // -- end - for - NoEPL --
                // 2 - save score
                // insert without control into priority queue (is not full)
                if ( (docScoreCalc < numberOfResults) && (partialScore != 0) )
                {
                    resPQ.add(new QueryProcessor.ResultBlock(currDID, partialScore));     // add to priority queue
                    docScoreCalc++;                         // increment result in priority queue counter
                    if (docScoreCalc == numberOfResults)
                        threshold = resPQ.peek().getScore();    // update threshold
                }
                else if (threshold < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                {
                    // substitution of the block
                    resPQ.poll();                           // remove the first element
                    resPQ.add(new QueryProcessor.ResultBlock(currDID, partialScore));     // add to priority queue
                    threshold = resPQ.peek().getScore();    // update threshold
                    // calculate new essential posting lists and update firstEssPostListIndex
                    //firstEssPostListIndex = updateEssentialPositngLists(sumTUBList, threshold, firstEssPostListIndex);
                    if ((threshold > newFEPL) && (firstEssPostListIndex != (plsLen - 1)))
                    {
                        firstEssPostListIndex = updateEssentialPositngLists(sumTUBList, threshold,firstEssPostListIndex);
                        newFEPL = sumTUBList[min(firstEssPostListIndex + 1, pLNotEmpty - 1)];
                        printDebug("new FEPL: " + firstEssPostListIndex + ", threshold: " + threshold + ", new value for new FEPL: " + newFEPL);
                    }
                }
            }   // -- end - for: DID --
        }   // -- start - if conjunctive --
        endTime = System.currentTimeMillis();           // end time of DAAT
        // shows DAAT execution time
        printTime("*** Max Score PQ(skipping) execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
    }

    /**
     * Function for apply the Document at a Time algorithm with max Score algorithm as dynamic pruning algorithm.
     * This function apply the Max Score in the case of both compression and skipping are enabled.
     *
     * @param processedQuery    array list for containing the query term
     * @param scoringFunc       indicates the preference for scoring. If false use TFIDF, if true use BM25.
     * @param isConjunctive     indicates whether the query is conjunctive or disjunctive type (default is disjunctive)
     * @param numberOfResults   indicated the max number of result to return to user
     */
    private static void DAATAlgMAXSCORESkipAndComp(ArrayList<String> processedQuery, boolean scoringFunc , boolean isConjunctive, int numberOfResults) throws FileNotFoundException
    {
        resPQ = new PriorityQueue<>(numberOfResults, new CompareTerm());    // length equal to the number of results to be returned to the user
        PriorityQueue<Integer> didPQ = new PriorityQueue<Integer>();    // priority queue for the DID
        ArrayList<String> newProcQuery;     // new processed query after removing term not in the dictionary
        String[] orderedQueryTerm;          // contains the query term ordered by TermUpperBound (ascending order)
        double[] termUpperBoundList;        // contains all the term upper bound for each term of the query
        double[] sumTUBList;                // array containing the sum of TUB, the value at the i-th position is the sum of TUBs from position 0 to (i-1)
        double[] IDFweight;                 // array containing the IDF weight for each posting list
        int[] lengthPostingList;            // array containing the length of the posting lists
        int[] postingListsIndex;            // contain the current position index for the posting list (value from 0 to SkipBlockLen) of each term in the query
        int[] blockIndex;                   // contain the index of the next block of skipList for each PL
        SkipList[] skipListArray;           // array of the Skip List reference related to the term of query
        Posting currentP;                   // support var
        boolean resetScore = false;         // used only in conjunctive case. indicates that the score must be set to 0 (the current Doc there aren't all the term of the query)
        double threshold = 0;               // var that contain the current threshold for MaxScore (is the minimum score value to be in the current best result)
        double currentDocUpperBound = 0;    // current document upper bound (used in max score algorithm for early stopping)
        double partialScore = 0;            // var that contain partial score
        int firstEssPostListIndex = 0;      // indicates the index of the first (current) essential posting list
        int docScoreCalc = 0;               // indicates the number of documents whose score was calculated (0 to number of results requested by the user)
        int pLNotEmpty = 0;                 // contains the number of posting lists related to the query terms that aren't empty
        int postListCurrDID = 0;            // the DID in the current position of the posting list (used in no-essential PL)
        int currentDID = 0;                 // DID of the current doc processed in algorithm
        int procQLen = 0;                   // the number of term in the original preprocessed query (before removing term not in the dictionary)
        long startTime,endTime;             // variables to calculate the execution time

        procQLen = processedQuery.size();               // take the size of the processed Query before term removing
        newProcQuery = removeTermNotInDictFromQuery(processedQuery);    // remove term that are not in the dictionary
        pLNotEmpty = newProcQuery.size();               // take the number of term that have Posting List
        orderedQueryTerm = new String[pLNotEmpty];      // set len
        termUpperBoundList  = new double[pLNotEmpty];   // set len
        sumTUBList  = new double[pLNotEmpty];           // set len

        if((procQLen != pLNotEmpty) && (isConjunctive)) // there is at least one term not in dictionary, if conjunctive -> no results must be returned
            return;         // exit

        // check the number of posting lists not empty and perform the best choice
        if (pLNotEmpty == 0)        // all terms in the query aren't in the dictionary or empty query
            return;     // exit
        else if (pLNotEmpty == 1)   // there is only 1 postingList (query with one term or query with more term but only one in dictionary)
        {
            // create the skip List reference related to the term of query
            skipListArray = skipListInitCompAndSkip(newProcQuery, null, false);
            // The PL is only one -> read and decompress the whole PL and use the classic optimization method
            startTime = System.currentTimeMillis();         // start time for retrieve the posting lists of the query
            ArrayList<Posting>[] postingLists = retrieveAllUncompPL(newProcQuery, skipListArray); // get the uncompress PL
            endTime = System.currentTimeMillis();           // end time for retrieve the posting lists of the query
            // shows query execution time
            printTime("*** MAX SCORE (comp+skipping) retrieved PL (case 1 PL) in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

            DAATOnePostingList(newProcQuery.get(0), postingLists[0], scoringFunc, numberOfResults);   // execute DAAT algorithm
            return;     // exit
        }

        // -- more postingLists not empty --
        // 0 - set the array of skip list for the query
        setAllUtilsListMAxScoreCompAndSkipping(newProcQuery, orderedQueryTerm, termUpperBoundList, sumTUBList);   // calculate utilities arrays
        // create the skip List reference related to the term of query
        skipListArray = skipListInitCompAndSkip(newProcQuery, orderedQueryTerm, true);

        int maxDID = CollectionStatistics.getNDocs();
        BitSet processedDocs = new BitSet(maxDID + 1);      // "+1" to include the max value

        // 0.5 - take the first block of the PL
        skipAndCompPLs = new ArrayList[pLNotEmpty];
        for (int i = 0; i < pLNotEmpty; i++)
        {
            if(!retrieveCompBlockOfPLAndUncompressMaxScore(orderedQueryTerm[i], skipListArray, 0, i))
                return;     // error
            int currDID = skipAndCompPLs[i].get(0).getDocId();
            if (!processedDocs.get(currDID))
            {
                didPQ.add(currDID);             // add the first DID in the PLs
                processedDocs.set(currDID);     // set the bit related to currDID in the bitset
            }
        }

        blockIndex = new int[pLNotEmpty];       // set the index block for each PL
        Arrays.fill(blockIndex, 1);         // the first block (with index 0) has already been taken
        postingListsIndex = getPostingListsIndex(skipAndCompPLs);               // get the index initialized
        lengthPostingList = retrieveLengthAllPostingLists(orderedQueryTerm);    // take the length of each posting list
        IDFweight = calculateIDFWeight(lengthPostingList);                      // calculate the IDF weight
        double newFEPL = sumTUBList[1];         // if this value is exceeded with the threshold, the first essential posting list must be recalculated

        startTime = System.currentTimeMillis();           // start time of DAAT + MAX SCORE (comp + skipping)
        // MaxScore algorithm - scan all Doc retrieved and calculate score (TFIDF or BM25)
        if (isConjunctive)
        {   // -- start - if conjunctive --
            while (!didPQ.isEmpty())
            {   // -- start - while 0: DID --
                currentDID = didPQ.poll();  // take the current DID
                partialScore = 0;           // reset var
                // 0 - scan the essential posting lists, default case is query Disjunctive
                for (int j = firstEssPostListIndex; j < pLNotEmpty; j++)
                {   // -- start - for 0.1: EPL --
                    if (skipAndCompPLs[j] != null)
                    {   // -- start - if 0.1.1: check PL null --
                        // check if the j-th term of the query is present in the doc identify by currentDID
                        currentP = skipAndCompPLs[j].get(postingListsIndex[j]);              // take posting
                        if (currentP.getDocId() == currentDID)
                        {
                            postingListsIndex[j]++;                         // update index of current value
                            // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                            if (scoringFunc)
                                partialScore += ScoringBM25(currentDID,currentP.getTermFreq(), IDFweight[j]);     // use BM25
                            else
                                partialScore += ScoringTFIDF(currentP.getTermFreq(), IDFweight[j]);               // use TFIDF

                            // boundary check -> check if it needs to upload another block of the posting list
                            if (postingListsIndex[j] >= skipAndCompPLs[j].size())
                            {
                                // load the new block
                                if (retrieveCompBlockOfPLAndUncompressMaxScore(orderedQueryTerm[j], skipListArray, blockIndex[j], j))
                                {
                                    blockIndex[j]++;            // increment the counter of block
                                    postingListsIndex[j] = 0;   // reset the index of the posting list
                                }
                                else    // was the last block, the PL is over
                                {
                                    skipAndCompPLs[j] = null;   // set to null the posting list
                                    continue;
                                }
                            }
                            // add new DID from this PL
                            int currDID = skipAndCompPLs[j].get(postingListsIndex[j]).getDocId();
                            if (!processedDocs.get(currDID))
                            {
                                didPQ.add(currDID);             // add the first DID in the PLs
                                processedDocs.set(currDID);     // set the bit related to currDID in the bitset
                            }
                        }
                        else
                            resetScore = true;       // reset the partial score
                    }   // -- end - if 0.1.1: check PL null --
                    else     // must take only the document in which there are all term (DID that compare in all posting lists of the terms)
                    {
                        endTime = System.currentTimeMillis();           // end time of DAAT
                        printTime("*** MAX SCORE (comp+skipping) V2 execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
                        return;             // exit from function
                    }
                }   // -- end - for 0.1: EPL --

                // Conditions under which analysis of nonessential posting lists can be skipped -- SEE NOTE 0 --
                if ( (partialScore == 0) || resetScore )
                {
                    resetScore = false;
                    continue;       // go to next iteration, the current doc can't be among the top result
                }

                // 1 - scan non essential posting lists
                if (firstEssPostListIndex != 0)
                {
                    currentDocUpperBound = partialScore + sumTUBList[firstEssPostListIndex];    // calculate the current DUB
                    // check if the doc has no zero possibility to have a score greater than threshold
                    if (currentDocUpperBound <= threshold)
                        continue;                           // go to next iteration with next DID

                    // update the score
                    for (int i = 0; i < firstEssPostListIndex; i++)
                    {   // -- start - for: scan NoEPLs --
                        resetScore = true;       // reset the partial score

                        if (skipAndCompPLs[i] != null)
                        {
                            // check if it needs to upload another block of the posting list
                            if (postingListsIndex[i] >= skipAndCompPLs[i].size())
                            {
                                // load the new block
                                if (retrieveCompBlockOfPLAndUncompressMaxScore(orderedQueryTerm[i], skipListArray, blockIndex[i], i))
                                {
                                    blockIndex[i]++;            // increment the counter of block
                                    postingListsIndex[i] = 0;   // reset the index of the posting list
                                } else    // was the last block, the PL is over
                                {
                                    skipAndCompPLs[i] = null;   // set to null the posting list
                                    continue;
                                }
                            }

                            postListCurrDID = skipAndCompPLs[i].get(postingListsIndex[i]).getDocId();
                            // check first position
                            if (postListCurrDID == currentDID)
                            {
                                // find the searched DID - update the partialScore and the currentDocUpperBound
                                currentP = skipAndCompPLs[i].get(postingListsIndex[i]);              // take posting

                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                double currScore = 0;
                                // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                                if (scoringFunc)
                                    currScore += ScoringBM25(currentDID, currentP.getTermFreq(), IDFweight[i]);   // use BM25
                                else
                                    currScore += ScoringTFIDF(currentP.getTermFreq(), IDFweight[i]);             // use TFIDF
                                partialScore += currScore;
                                currentDocUpperBound += currScore;               // update currentDocUpperBound

                                resetScore = false;         // reset the partial score
                                postingListsIndex[i]++;     // update the index of the postin list
                                continue;
                            } else if (postListCurrDID > currentDID)      // the searched DID is not in this posting list
                            {
                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                continue;       // go to the next step
                            } else    // postListCurrDID < currentDID, search in the posting listw
                            {
                                postingListsIndex[i] = skipListArray[i].nextGEQCompSkip(currentDID, postingListsIndex[i], i, orderedQueryTerm[i]);
                            }
                            // check the index returned by nextGEQ
                            if (postingListsIndex[i] >= skipAndCompPLs[i].size())  // check for out of bound in case of reaching the end of the list
                            {
                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                continue;       // go to the next step
                            }

                            // check if the current target has been found or not
                            if (skipAndCompPLs[i].get(postingListsIndex[i]).getDocId() != currentDID)
                            {
                                // should be always greater than currentDID
                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                            } else                    // the searched DID there is in the posting list
                            {
                                // find the searched DID - update the partialScore and the currentDocUpperBound
                                currentP = skipAndCompPLs[i].get(postingListsIndex[i]);              // take posting

                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                double currScore = 0;
                                // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                                if (scoringFunc)
                                    currScore += ScoringBM25(currentDID, currentP.getTermFreq(), IDFweight[i]);   // use BM25
                                else
                                    currScore += ScoringTFIDF(currentP.getTermFreq(), IDFweight[i]);             // use TFIDF
                                partialScore += currScore;
                                currentDocUpperBound += currScore;               // update currentDocUpperBound

                                resetScore = false;         // reset the partial score
                                postingListsIndex[i]++;     // update the index of the postin list
                            }
                        } // if all postings in one posting lists have already been seen the next documents in the posting lists cannot contain all the terms in the query
                        else
                        {
                            endTime = System.currentTimeMillis();           // end time of DAAT
                            printTime("*** MAX SCORE (comp+skipping) V2 execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
                            return;             // exit from function
                        }
                        // check if the doc has no zero possibility to have a score greater than threshold
                        if ((currentDocUpperBound <= threshold) || resetScore)
                            break;
                    }   // -- start - for: scan NoEPLs --
                }
                // 2 - save score
                // insert without control into priority queue (is not full) or insert all results (orderAllHashMap = true)
                if (docScoreCalc < numberOfResults)
                {
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                    docScoreCalc++;                         // increment result in priority queue counter
                    if (docScoreCalc == numberOfResults)
                        threshold = resPQ.peek().getScore();    // update threshold
                }
                else if (threshold < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                {
                    // substitution of the block
                    resPQ.poll();                           // remove the first element
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                    threshold = resPQ.peek().getScore();    // update threshold
                    // calculate new essential posting lists and update firstEssPostListIndex
                    firstEssPostListIndex = updateEssentialPositngLists(sumTUBList, threshold, firstEssPostListIndex);
                    /*if ((threshold > newFEPL) && (firstEssPostListIndex != (pLNotEmpty - 1)))
                    {
                        firstEssPostListIndex = updateEssentialPositngLists(sumTUBList, threshold,firstEssPostListIndex);
                        newFEPL = sumTUBList[min(firstEssPostListIndex + 1, pLNotEmpty - 1)];
                    }*/
                }
            }   // -- end - while 0: DID --
        }   // -- end - if conjunctive --
        else        // disjunctive case
        {   // -- start - if disjunctive --
            while (!didPQ.isEmpty())
            {   // -- start - while 0: DID --
                currentDID = didPQ.poll();  // take the current DID
                partialScore = 0;           // reset var
                // 0 - scan the essential posting lists, default case is query Disjunctive
                for (int j = firstEssPostListIndex; j < pLNotEmpty; j++)
                {   // -- start - for 0.1: EPL --
                    if (skipAndCompPLs[j] != null)
                    {   // -- start - if 0.1.1: check PL null --
                        // check if the j-th term of the query is present in the doc identify by currentDID
                        currentP = skipAndCompPLs[j].get(postingListsIndex[j]);              // take posting
                        if (currentP.getDocId() == currentDID)
                        {
                            postingListsIndex[j]++;                         // update index of current value
                            // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                            if (scoringFunc)
                                partialScore += ScoringBM25(currentDID,currentP.getTermFreq(), IDFweight[j]);     // use BM25
                            else
                                partialScore += ScoringTFIDF(currentP.getTermFreq(), IDFweight[j]);               // use TFIDF

                            // boundary check -> check if it needs to upload another block of the posting list
                            if (postingListsIndex[j] >= skipAndCompPLs[j].size())
                            {
                                // load the new block
                                if (retrieveCompBlockOfPLAndUncompressMaxScore(orderedQueryTerm[j], skipListArray, blockIndex[j], j))
                                {
                                    blockIndex[j]++;            // increment the counter of block
                                    postingListsIndex[j] = 0;   // reset the index of the posting list
                                }
                                else    // was the last block, the PL is over
                                {
                                    skipAndCompPLs[j] = null;   // set to null the posting list
                                    continue;
                                }
                            }
                            // add new DID from this PL
                            int currDID = skipAndCompPLs[j].get(postingListsIndex[j]).getDocId();
                            if (!processedDocs.get(currDID))
                            {
                                didPQ.add(currDID);             // add the first DID in the PLs
                                processedDocs.set(currDID);     // set the bit related to currDID in the bitset
                            }
                        }
                    }   // -- end - if 0.1.1: check PL null --
                }   // -- end - for 0.1: EPL --

                // Conditions under which analysis of nonessential posting lists can be skipped -- SEE NOTE 0 --
                if (partialScore == 0)
                    continue;       // go to next iteration, the current doc can't be among the top result

                // 1 - scan non essential posting lists
                if (firstEssPostListIndex != 0)
                {
                    currentDocUpperBound = partialScore + sumTUBList[firstEssPostListIndex];    // calculate the current DUB
                    // check if the doc has no zero possibility to have a score greater than threshold
                    if (currentDocUpperBound <= threshold)
                        continue;                           // go to next iteration with next DID

                    // update the score
                    for (int i = 0; i < firstEssPostListIndex; i++)
                    {   // -- start - for: scan NoEPLs --
                        if (skipAndCompPLs[i] != null)
                        {
                            // check if it needs to upload another block of the posting list
                            if (postingListsIndex[i] >= skipAndCompPLs[i].size())
                            {
                                // load the new block
                                if (retrieveCompBlockOfPLAndUncompressMaxScore(orderedQueryTerm[i], skipListArray, blockIndex[i], i))
                                {
                                    blockIndex[i]++;            // increment the counter of block
                                    postingListsIndex[i] = 0;   // reset the index of the posting list
                                } else    // was the last block, the PL is over
                                {
                                    skipAndCompPLs[i] = null;   // set to null the posting list
                                    continue;
                                }
                            }

                            postListCurrDID = skipAndCompPLs[i].get(postingListsIndex[i]).getDocId();
                            // check first position
                            if (postListCurrDID == currentDID)
                            {
                                // find the searched DID - update the partialScore and the currentDocUpperBound
                                currentP = skipAndCompPLs[i].get(postingListsIndex[i]);              // take posting

                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                double currScore = 0;
                                // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                                if (scoringFunc)
                                    currScore = ScoringBM25(currentDID, currentP.getTermFreq(), IDFweight[i]);  // use BM25
                                else
                                    currScore = ScoringTFIDF(currentP.getTermFreq(), IDFweight[i]);             // use TFIDF
                                partialScore += currScore;
                                currentDocUpperBound += currScore;               // update currentDocUpperBound

                                postingListsIndex[i]++;     // update the index of the postin list
                                continue;
                            } else if (postListCurrDID > currentDID)      // the searched DID is not in this posting list
                            {
                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                continue;       // go to the next step
                            } else    // postListCurrDID < currentDID, search in the posting listw
                            {
                                postingListsIndex[i] = skipListArray[i].nextGEQCompSkip(currentDID, postingListsIndex[i], i, orderedQueryTerm[i]);
                            }
                            // check the index returned by nextGEQ
                            if (postingListsIndex[i] >= skipAndCompPLs[i].size())  // check for out of bound in case of reaching the end of the list
                            {
                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                continue;       // go to the next step
                            }

                            // check if the current target has been found or not
                            if (skipAndCompPLs[i].get(postingListsIndex[i]).getDocId() != currentDID)
                            {
                                // should be always greater than currentDID
                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                            } else                    // the searched DID there is in the posting list
                            {
                                // find the searched DID - update the partialScore and the currentDocUpperBound
                                currentP = skipAndCompPLs[i].get(postingListsIndex[i]);              // take posting

                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                double currScore = 0;
                                // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                                if (scoringFunc)
                                    currScore = ScoringBM25(currentDID, currentP.getTermFreq(), IDFweight[i]);  // use BM25
                                else
                                    currScore = ScoringTFIDF(currentP.getTermFreq(), IDFweight[i]);             // use TFIDF
                                partialScore += currScore;
                                currentDocUpperBound += currScore;               // update currentDocUpperBound

                                postingListsIndex[i]++;     // update the index of the postin list
                            }
                        }
                        // check if the doc has no zero possibility to have a score greater than threshold
                        if (currentDocUpperBound <= threshold)
                            break;
                    }   // -- start - for: scan NoEPLs --
                }
                // 2 - save score
                // insert without control into priority queue (is not full) or insert all results (orderAllHashMap = true)
                if (docScoreCalc < numberOfResults)
                {
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                    docScoreCalc++;                         // increment result in priority queue counter
                    if (docScoreCalc == numberOfResults)
                        threshold = resPQ.peek().getScore();    // update threshold
                }
                else if (threshold < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                {
                    // substitution of the block
                    resPQ.poll();                           // remove the first element
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                    threshold = resPQ.peek().getScore();    // update threshold
                    // calculate new essential posting lists and update firstEssPostListIndex
                    firstEssPostListIndex = updateEssentialPositngLists(sumTUBList, threshold, firstEssPostListIndex);
                    /*if ((threshold > newFEPL) && (firstEssPostListIndex != (pLNotEmpty - 1)))
                    {
                        firstEssPostListIndex = updateEssentialPositngLists(sumTUBList, threshold,firstEssPostListIndex);
                        newFEPL = sumTUBList[min(firstEssPostListIndex + 1, pLNotEmpty - 1)];
                    }*/
                }
            }   // -- end - while 0: DID --
        }   // -- end - else disjunctive --
        endTime = System.currentTimeMillis();           // end time of DAAT
        printTime("*** MAX SCORE (comp+skipping) V2 execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
    }

    /**
     * Function for apply the Document at a Time algorithm with max Score algorithm as dynamic pruning algorithm.
     * This function apply the Max Score when compression, skipping and wholePLinMem are enabled.
     * In this method we perform the same algorithm as in the case of compression and skipping both enabled with the
     * difference that the compressed PLs will not be loaded from memory in blocks but whole.
     * Then they will be decompressed and processed in blocks.
     * In this case it has the entire PL compressed in memory and current block decompressed for each PL.
     * This brings a saving in memory used compared to the base case (whole PLs decompressed in memory) but worse than
     * the case with skipping and compression enabled (only the current blocks decompressed in memory).
     * This methodology may be useful in the case of small collections (not large enough to occupy memory entirely with
     * query PLs) because it will occupy less memory than the base case but will have less overhead due to the mechanisms
     * in the compression and skipping enabled case.
     * The latter case may be too algorithm-intensive if memory is more than sufficient to handle the base case, leading
     * to performance degradation rather than improvement. in such cases, this method in the middle may be a good
     * compromise between memory usage, advanced techniques and performance.
     *
     * @param processedQuery    array list for containing the query term
     * @param scoringFunc       indicates the preference for scoring. If false use TFIDF, if true use BM25.
     * @param isConjunctive     indicates whether the query is conjunctive or disjunctive type (default is disjunctive)
     * @param numberOfResults   indicated the max number of result to return to user
     */
    private static void DAATAlgMAXSCORESkipAndCompAndWholePL(ArrayList<String> processedQuery, boolean scoringFunc , boolean isConjunctive, int numberOfResults) throws FileNotFoundException
    {
        resPQ = new PriorityQueue<>(numberOfResults, new CompareTerm());    // length equal to the number of results to be returned to the user
        PriorityQueue<Integer> didPQ = new PriorityQueue<Integer>();    // priority queue for the DID
        ArrayList<String> newProcQuery;     // new processed query after removing term not in the dictionary
        String[] orderedQueryTerm;          // contains the query term ordered by TermUpperBound (ascending order)
        double[] termUpperBoundList;        // contains all the term upper bound for each term of the query
        double[] sumTUBList;                // array containing the sum of TUB, the value at the i-th position is the sum of TUBs from position 0 to (i-1)
        double[] IDFweight;                 // array containing the IDF weight for each posting list
        int[] lengthPostingList;            // array containing the length of the posting lists
        int[] postingListsIndex;            // contain the current position index for the posting list (value from 0 to SkipBlockLen) of each term in the query
        int[] blockIndex;                   // contain the index of the next block of skipList for each PL
        SkipList[] skipListArray;           // array of the Skip List reference related to the term of query
        Posting currentP;                   // support var
        boolean resetScore = false;         // used only in conjunctive case. indicates that the score must be set to 0 (the current Doc there aren't all the term of the query)
        double threshold = 0;               // var that contain the current threshold for MaxScore (is the minimum score value to be in the current best result)
        double currentDocUpperBound = 0;    // current document upper bound (used in max score algorithm for early stopping)
        double partialScore = 0;            // var that contain partial score
        int firstEssPostListIndex = 0;      // indicates the index of the first (current) essential posting list
        int docScoreCalc = 0;               // indicates the number of documents whose score was calculated (0 to number of results requested by the user)
        int pLNotEmpty = 0;                 // contains the number of posting lists related to the query terms that aren't empty
        int postListCurrDID = 0;            // the DID in the current position of the posting list (used in no-essential PL)
        int currentDID = 0;                 // DID of the current doc processed in algorithm
        int procQLen = 0;                   // the number of term in the original preprocessed query (before removing term not in the dictionary)
        long startTime,endTime;             // variables to calculate the execution time

        procQLen = processedQuery.size();               // take the size of the processed Query before term removing
        newProcQuery = removeTermNotInDictFromQuery(processedQuery);    // remove term that are not in the dictionary
        pLNotEmpty = newProcQuery.size();               // take the number of term that have Posting List
        orderedQueryTerm = new String[pLNotEmpty];      // set len
        termUpperBoundList  = new double[pLNotEmpty];   // set len
        sumTUBList  = new double[pLNotEmpty];           // set len

        if((procQLen != pLNotEmpty) && (isConjunctive)) // there is at least one term not in dictionary, if conjunctive -> no results must be returned
            return;         // exit

        // check the number of posting lists not empty and perform the best choice
        if (pLNotEmpty == 0)        // all terms in the query aren't in the dictionary or empty query
            return;     // exit
        else if (pLNotEmpty == 1)   // there is only 1 postingList (query with one term or query with more term but only one in dictionary)
        {
            // create the skip List reference related to the term of query
            skipListArray = skipListInitCompAndSkip(newProcQuery, null, false);
            // The PL is only one -> read and decompress the whole PL and use the classic optimization method
            startTime = System.currentTimeMillis();         // start time for retrieve the posting lists of the query
            ArrayList<Posting>[] postingLists = retrieveAllUncompPL(newProcQuery, skipListArray); // get the uncompress PL
            endTime = System.currentTimeMillis();           // end time for retrieve the posting lists of the query
            // shows query execution time
            printTime("*** MAX SCORE (comp+skipping) retrieved PL (case 1 PL) in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

            DAATOnePostingList(newProcQuery.get(0), postingLists[0], scoringFunc, numberOfResults);   // execute DAAT algorithm
            return;     // exit
        }

        // -- more postingLists not empty --
        // 0) take the first block of the PL
        skipAndCompPLs = new ArrayList[pLNotEmpty];
        setAllUtilsListMAxScoreCompAndSkipping(newProcQuery, orderedQueryTerm, termUpperBoundList, sumTUBList);   // calculate utilities arrays
        // create the skip List reference related to the term of query
        skipListArray = skipListInitCompAndSkip(newProcQuery, orderedQueryTerm, true);

        startTime = System.currentTimeMillis();         // start time for retrieve all compressed PLs of the query
        retrieveAllCompPLMaxScore(orderedQueryTerm);    // load the whole compressed PLs
        endTime = System.currentTimeMillis();           // end time for retrieve all compressed PLs of the query
        printTime("*** MAX SCORE (comp + skipping + WPLInMem) retrieved first block of each PLs in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

        blockIndex = new int[pLNotEmpty];    // set the index block for each PL
        Arrays.fill(blockIndex, 1);                 // the first block (with index 0) has already been taken
        postingListsIndex = getPostingListsIndex(skipAndCompPLs);               // get the index initialized
        lengthPostingList = retrieveLengthAllPostingLists(orderedQueryTerm);    // take the length of each posting list
        IDFweight = calculateIDFWeight(lengthPostingList);                      // calculate the IDF weight

        // 0.5 - take the first block of the PLs
        for (int i = 0; i < pLNotEmpty; i++)
        {
            if(!UncompressOneBlockFromCompressedPLMAxScore(orderedQueryTerm, skipListArray, 0, i))
                return;     // error
            int currDID = skipAndCompPLs[i].get(0).getDocId();
            if (!didPQ.contains(currDID))
                didPQ.add(currDID);             // add the first DID in the PLs
        }
        double newFEPL = sumTUBList[1];         // if this value is exceeded with the threshold, the first essential posting list must be recalculated

        startTime = System.currentTimeMillis();     // start time of DAAT + MAX SCORE (comp + skipping + wholePLInMem)
        // MaxScore algorithm - scan all Doc retrieved and calculate score (TFIDF or BM25)
        if (isConjunctive)
        {   // -- start - if conjunctive --
            while (!didPQ.isEmpty())
            {   // -- start - while 0: DID --
                currentDID = didPQ.poll();  // take the current DID
                partialScore = 0;           // reset var
                // 0 - scan the essential posting lists, default case is query Disjunctive
                for (int j = firstEssPostListIndex; j < pLNotEmpty; j++)
                {   // -- start - for 0.1: EPL --
                    if (skipAndCompPLs[j] != null)
                    {   // -- start - if 0.1.1: check PL null --
                        // check if the j-th term of the query is present in the doc identify by currentDID
                        currentP = skipAndCompPLs[j].get(postingListsIndex[j]);              // take posting
                        if (currentP.getDocId() == currentDID)
                        {
                            postingListsIndex[j]++;                         // update index of current value
                            // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                            if (scoringFunc)
                                partialScore += ScoringBM25(currentDID,currentP.getTermFreq(), IDFweight[j]);     // use BM25
                            else
                                partialScore += ScoringTFIDF(currentP.getTermFreq(), IDFweight[j]);               // use TFIDF

                            // boundary check -> check if it needs to upload another block of the posting list
                            if (postingListsIndex[j] >= skipAndCompPLs[j].size())
                            {
                                // load the new block
                                if (UncompressOneBlockFromCompressedPLMAxScore(orderedQueryTerm, skipListArray, blockIndex[j], j))
                                {
                                    blockIndex[j]++;            // increment the counter of block
                                    postingListsIndex[j] = 0;   // reset the index of the posting list
                                }
                                else    // was the last block, the PL is over
                                {
                                    skipAndCompPLs[j] = null;   // set to null the posting list
                                    continue;
                                }
                            }
                            // add new DID from this PL
                            int currDID = skipAndCompPLs[j].get(postingListsIndex[j]).getDocId();
                            if (!didPQ.contains(currDID))
                                didPQ.add(currDID);             // add the first DID in the PLs
                        }
                        else
                            resetScore = true;       // reset the partial score

                    }   // -- end - if 0.1.1: check PL null --
                    else        // must take only the document in which there are all term (DID that compare in all posting lists of the terms)
                    {
                        endTime = System.currentTimeMillis();           // end time of DAAT
                        printTime("*** MAX SCORE (comp + skipping + WPLInMem) execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
                        return;             // exit from function
                    }
                }   // -- end - for 0.1: EPL --

                // Conditions under which analysis of nonessential posting lists can be skipped -- SEE NOTE 0 --
                if ( (partialScore == 0) || resetScore )
                {
                    resetScore = false;
                    continue;       // go to next iteration, the current doc can't be among the top result
                }

                // 1 - scan non essential posting lists
                if (firstEssPostListIndex != 0)
                {
                    currentDocUpperBound = partialScore + sumTUBList[firstEssPostListIndex];    // calculate the current DUB
                    // check if the doc has no zero possibility to have a score greater than threshold
                    if (currentDocUpperBound <= threshold)
                        continue;                           // go to next iteration with next DID

                    // update the score
                    for (int i = 0; i < firstEssPostListIndex; i++)
                    {   // -- start - for: scan NoEPLs --
                        resetScore = true;       // reset the partial score

                        if (skipAndCompPLs[i] != null)
                        {
                            // check if it needs to upload another block of the posting list
                            if (postingListsIndex[i] >= skipAndCompPLs[i].size())
                            {
                                // load the new block
                                if (UncompressOneBlockFromCompressedPLMAxScore(orderedQueryTerm, skipListArray, blockIndex[i], i))
                                {
                                    blockIndex[i]++;            // increment the counter of block
                                    postingListsIndex[i] = 0;   // reset the index of the posting list
                                } else    // was the last block, the PL is over
                                {
                                    skipAndCompPLs[i] = null;   // set to null the posting list
                                    continue;
                                }
                            }

                            postListCurrDID = skipAndCompPLs[i].get(postingListsIndex[i]).getDocId();
                            // check first position
                            if (postListCurrDID == currentDID)
                            {
                                // find the searched DID - update the partialScore and the currentDocUpperBound
                                currentP = skipAndCompPLs[i].get(postingListsIndex[i]);              // take posting

                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                double currScore = 0;
                                // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                                if (scoringFunc)
                                    currScore += ScoringBM25(currentDID, currentP.getTermFreq(), IDFweight[i]);   // use BM25
                                else
                                    currScore += ScoringTFIDF(currentP.getTermFreq(), IDFweight[i]);             // use TFIDF
                                partialScore += currScore;
                                currentDocUpperBound += currScore;               // update currentDocUpperBound

                                resetScore = false;         // reset the partial score
                                postingListsIndex[i]++;     // update the index of the postin list
                                continue;
                            } else if (postListCurrDID > currentDID)      // the searched DID is not in this posting list
                            {
                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                continue;       // go to the next step
                            } else    // postListCurrDID < currentDID, search in the posting listw
                            {
                                postingListsIndex[i] = skipListArray[i].nextGEQCompSkip(currentDID, postingListsIndex[i], i, orderedQueryTerm[i]);
                            }
                            // check the index returned by nextGEQ
                            if (postingListsIndex[i] >= skipAndCompPLs[i].size())  // check for out of bound in case of reaching the end of the list
                            {
                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                continue;       // go to the next step
                            }

                            // check if the current target has been found or not
                            if (skipAndCompPLs[i].get(postingListsIndex[i]).getDocId() != currentDID)
                            {
                                // should be always greater than currentDID
                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                            } else                    // the searched DID there is in the posting list
                            {
                                // find the searched DID - update the partialScore and the currentDocUpperBound
                                currentP = skipAndCompPLs[i].get(postingListsIndex[i]);              // take posting

                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                double currScore = 0;
                                // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                                if (scoringFunc)
                                    currScore += ScoringBM25(currentDID, currentP.getTermFreq(), IDFweight[i]);   // use BM25
                                else
                                    currScore += ScoringTFIDF(currentP.getTermFreq(), IDFweight[i]);             // use TFIDF
                                partialScore += currScore;
                                currentDocUpperBound += currScore;               // update currentDocUpperBound

                                resetScore = false;         // reset the partial score
                                postingListsIndex[i]++;     // update the index of the postin list
                            }
                        } // if all postings in one posting lists have already been seen the next documents in the posting lists cannot contain all the terms in the query
                        else
                        {
                            endTime = System.currentTimeMillis();           // end time of DAAT
                            printTime("*** MAX SCORE (comp + skipping + WPLInMem) execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
                            return;             // exit from function
                        }
                        // check if the doc has no zero possibility to have a score greater than threshold
                        if ((currentDocUpperBound <= threshold) || resetScore)
                            break;
                    }   // -- start - for: scan NoEPLs --
                }
                // 2 - save score
                // insert without control into priority queue (is not full) or insert all results (orderAllHashMap = true)
                if (docScoreCalc < numberOfResults)
                {
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                    docScoreCalc++;                         // increment result in priority queue counter
                    if (docScoreCalc == numberOfResults)
                        threshold = resPQ.peek().getScore();    // update threshold
                }
                else if (threshold < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                {
                    // substitution of the block
                    resPQ.poll();                           // remove the first element
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                    threshold = resPQ.peek().getScore();    // update threshold
                    // calculate new essential posting lists and update firstEssPostListIndex
                    firstEssPostListIndex = updateEssentialPositngLists(sumTUBList, threshold, firstEssPostListIndex);
                    /*if ((threshold > newFEPL) && (firstEssPostListIndex != (pLNotEmpty - 1)))
                    {
                        firstEssPostListIndex = updateEssentialPositngLists(sumTUBList, threshold,firstEssPostListIndex);
                        newFEPL = sumTUBList[min(firstEssPostListIndex + 1, pLNotEmpty - 1)];
                    }*/
                }
            }   // -- end - while 0: DID --
        }   // -- end - if conjunctive --
        else        // disjunctive case
        {   // -- start - if disjunctive --
            while (!didPQ.isEmpty())
            {   // -- start - while 0: DID --
                currentDID = didPQ.poll();  // take the current DID
                partialScore = 0;           // reset var
                // 0 - scan the essential posting lists, default case is query Disjunctive
                for (int j = firstEssPostListIndex; j < pLNotEmpty; j++)
                {   // -- start - for 0.1: EPL --
                    if (skipAndCompPLs[j] != null)
                    {   // -- start - if 0.1.1: check PL null --
                        // check if the j-th term of the query is present in the doc identify by currentDID
                        currentP = skipAndCompPLs[j].get(postingListsIndex[j]);              // take posting
                        if (currentP.getDocId() == currentDID)
                        {
                            postingListsIndex[j]++;                         // update index of current value
                            // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                            if (scoringFunc)
                                partialScore += ScoringBM25(currentDID,currentP.getTermFreq(), IDFweight[j]);     // use BM25
                            else
                                partialScore += ScoringTFIDF(currentP.getTermFreq(), IDFweight[j]);               // use TFIDF

                            // boundary check -> check if it needs to upload another block of the posting list
                            if (postingListsIndex[j] >= skipAndCompPLs[j].size())
                            {
                                // load the new block
                                if (UncompressOneBlockFromCompressedPLMAxScore(orderedQueryTerm, skipListArray, blockIndex[j], j))
                                {
                                    blockIndex[j]++;            // increment the counter of block
                                    postingListsIndex[j] = 0;   // reset the index of the posting list
                                }
                                else    // was the last block, the PL is over
                                {
                                    skipAndCompPLs[j] = null;   // set to null the posting list
                                    continue;
                                }
                            }
                            // add new DID from this PL
                            int currDID = skipAndCompPLs[j].get(postingListsIndex[j]).getDocId();
                            if (!didPQ.contains(currDID))
                                didPQ.add(currDID);             // add the first DID in the PLs
                        }
                    }   // -- end - if 0.1.1: check PL null --
                }   // -- end - for 0.1: EPL --

                // Conditions under which analysis of nonessential posting lists can be skipped -- SEE NOTE 0 --
                if (partialScore == 0)
                    continue;       // go to next iteration, the current doc can't be among the top result

                // 1 - scan non essential posting lists
                if (firstEssPostListIndex != 0)
                {
                    currentDocUpperBound = partialScore + sumTUBList[firstEssPostListIndex];    // calculate the current DUB
                    // check if the doc has no zero possibility to have a score greater than threshold
                    if (currentDocUpperBound <= threshold)
                        continue;                           // go to next iteration with next DID

                    // update the score
                    for (int i = 0; i < firstEssPostListIndex; i++)
                    {   // -- start - for: scan NoEPLs --

                        if (skipAndCompPLs[i] != null)
                        {
                            // check if it needs to upload another block of the posting list
                            if (postingListsIndex[i] >= skipAndCompPLs[i].size())
                            {
                                // load the new block
                                if (UncompressOneBlockFromCompressedPLMAxScore(orderedQueryTerm, skipListArray, blockIndex[i], i))
                                {
                                    blockIndex[i]++;            // increment the counter of block
                                    postingListsIndex[i] = 0;   // reset the index of the posting list
                                } else    // was the last block, the PL is over
                                {
                                    skipAndCompPLs[i] = null;   // set to null the posting list
                                    continue;
                                }
                            }

                            postListCurrDID = skipAndCompPLs[i].get(postingListsIndex[i]).getDocId();
                            // check first position
                            if (postListCurrDID == currentDID)
                            {
                                // find the searched DID - update the partialScore and the currentDocUpperBound
                                currentP = skipAndCompPLs[i].get(postingListsIndex[i]);              // take posting

                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                double currScore = 0;
                                // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                                if (scoringFunc)
                                    currScore += ScoringBM25(currentDID, currentP.getTermFreq(), IDFweight[i]);   // use BM25
                                else
                                    currScore += ScoringTFIDF(currentP.getTermFreq(), IDFweight[i]);             // use TFIDF
                                partialScore += currScore;
                                currentDocUpperBound += currScore;               // update currentDocUpperBound

                                postingListsIndex[i]++;     // update the index of the postin list
                                continue;
                            } else if (postListCurrDID > currentDID)      // the searched DID is not in this posting list
                            {
                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                continue;       // go to the next step
                            } else    // postListCurrDID < currentDID, search in the posting listw
                            {
                                postingListsIndex[i] = skipListArray[i].nextGEQCompSkip(currentDID, postingListsIndex[i], i, orderedQueryTerm[i]);
                            }
                            // check the index returned by nextGEQ
                            if (postingListsIndex[i] >= skipAndCompPLs[i].size())  // check for out of bound in case of reaching the end of the list
                            {
                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                continue;       // go to the next step
                            }

                            // check if the current target has been found or not
                            if (skipAndCompPLs[i].get(postingListsIndex[i]).getDocId() != currentDID)
                            {
                                // should be always greater than currentDID
                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                            } else                    // the searched DID there is in the posting list
                            {
                                // find the searched DID - update the partialScore and the currentDocUpperBound
                                currentP = skipAndCompPLs[i].get(postingListsIndex[i]);              // take posting

                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                double currScore = 0;
                                // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                                if (scoringFunc)
                                    currScore += ScoringBM25(currentDID, currentP.getTermFreq(), IDFweight[i]);   // use BM25
                                else
                                    currScore += ScoringTFIDF(currentP.getTermFreq(), IDFweight[i]);             // use TFIDF
                                partialScore += currScore;
                                currentDocUpperBound += currScore;               // update currentDocUpperBound

                                postingListsIndex[i]++;     // update the index of the postin list
                            }
                        } // if all postings in one posting lists have already been seen the next documents in the posting lists cannot contain all the terms in the query
                        // check if the doc has no zero possibility to have a score greater than threshold
                        if (currentDocUpperBound <= threshold)
                            break;
                    }   // -- start - for: scan NoEPLs --
                }
                // 2 - save score
                // insert without control into priority queue (is not full) or insert all results (orderAllHashMap = true)
                if (docScoreCalc < numberOfResults)
                {
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                    docScoreCalc++;                         // increment result in priority queue counter
                    if (docScoreCalc == numberOfResults)
                        threshold = resPQ.peek().getScore();    // update threshold
                }
                else if (threshold < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                {
                    // substitution of the block
                    resPQ.poll();                           // remove the first element
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                    threshold = resPQ.peek().getScore();    // update threshold
                    // calculate new essential posting lists and update firstEssPostListIndex
                    firstEssPostListIndex = updateEssentialPositngLists(sumTUBList, threshold, firstEssPostListIndex);
                    /*if ((threshold > newFEPL) && (firstEssPostListIndex != (pLNotEmpty - 1)))
                    {
                        firstEssPostListIndex = updateEssentialPositngLists(sumTUBList, threshold,firstEssPostListIndex);
                        newFEPL = sumTUBList[min(firstEssPostListIndex + 1, pLNotEmpty - 1)];
                    }*/
                }
            }   // -- end - while 0: DID --
        }   // -- end - else disjunctive --
        endTime = System.currentTimeMillis();       // end time of DAAT + MAX SCORE (comp + skipping + wholePLInMem)
        printTime("*** MAX SCORE (comp + skipping + WPLInMem) execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
    }

    // -------------------------------------------- END - Execution Query alg ------------------------------------------

    /**
     * Function to calculate the number of the query terms that are not in the dictionary.
     *
     * @param processedQuery    array list for containing the query term
     * @return      return the new array list containing only the term in the dictionary from processedQuery
     */
    private static ArrayList<String> removeTermNotInDictFromQuery(ArrayList<String> processedQuery)
    {
        ArrayList<String> newQuery = new ArrayList<>();     // new processed query after removing term not in the dictionary

        for (String term : processedQuery)      // scan all query terms
        {
            // there is a posting list for the query term == the term is in the collection
            if (dictionary.getTermToTermStat().containsKey(term))
                newQuery.add(term);             // add the term to newQuery
            else
                termNotInCollection.add(term);  // update array list of the term not in collection
        }

        return newQuery;
    }

    // -------- start: scoring function --------

    /**
     * Function to calculate TFIDF for one term and one document. The complete formula is:
     * TFIDF(t,d) = TFWeight * IDFWeight = (1 + log10(tf)) * log10(NDoc/postListLen(t))
     *
     * @param termFreq      term frequency of the term in the document
     * @param IDFweight     the IDF weight for the current term (precomputed for optimization)
     * @return  the TFIDF score for one term and one document. The total score for a document will be the sum of the
     *          result of this function for each term that is both in the document and in the query
     */
    private static Double ScoringTFIDF(int termFreq, double IDFweight)
    {
        double TFweight, scoreTFIDF;     // variables to calculate the TFIDF score value

        if (termFreq == 0)          // control to avoid acess to position -1
            return (double) 0;

        TFweight = logTermFreq[termFreq-1];   // take TF weight from memory
        scoreTFIDF = TFweight * IDFweight;          // calculate TFIDF weight from Tf and IDF weight values
        return scoreTFIDF;
    }

    /**
     * Function to calculate IDF weight for each term(posting list) of the query. ( IDFw(t,d) = log10(NDoc / postListLen(t)) )
     *
     * @param postingListsLength    array containing the length for the posting lists
     * @return array of IDF weight, the i-th values in the array correspond to the i-th posting list(term of the query)
     */
    private static double[] calculateIDFWeight(int[] postingListsLength)
    {
        // for each term the IDF weight = Math.log10((CollectionStatistics.getNDocs() / postListLength));
        double[] IDFweight = new double[postingListsLength.length];       // array containing the IDF weight for each posting list

        for (int i = 0; i < postingListsLength.length; i++)
        {
            // calculate the IDF weight for the i-th posting list(term of the query)
            IDFweight[i] = Math.log10(((double) CollectionStatistics.getNDocs() / postingListsLength[i]));
        }

        return IDFweight;
    }

    /**
     * Function to calculate and store the precomputed TermFreqWeight.
     * Read the max value of the TermFreq in the collection and compute the weight for each value from 1 to max.
     */
    public static void calcAndStoreTFWeight()
    {
        int maxTF = CollectionStatistics.getMaxTermFreq();  // get the max value of TermFreq in the collection
        double TFweight;
        long startTime, endTime;

        printLoad("Calculated and stored all useful values for termFreqWeight...");

        CollectionStatistics.readCollectionStatsFromDisk();
        logTermFreq = new double[CollectionStatistics.getMaxTermFreq()];    // initialize array for the precomputed value

        startTime = System.currentTimeMillis();
        try (
                RandomAccessFile docStats = new RandomAccessFile(TERMFREQWEIGHT_FILE, "rw");
                FileChannel channel = docStats.getChannel()
        ) {
            // double size * number of term freq weight to store
            MappedByteBuffer buffer = channel.map(FileChannel.MapMode.READ_WRITE, 0, (long) DOUBLE_BYTES * maxTF);

            // scan all possible value of TermFreq weight
            for(int i = 1; i <= maxTF; i++)
            {
                TFweight = (1 + Math.log10(i));         // calculate TF weight
                logTermFreq[i-1] = TFweight;            // save TF weight in memory
                // store
                buffer.putDouble(TFweight);             // write termFreqWeight into disk
            }
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
        //printDebug("The maxTF is: " + maxTF + " and the size of termFreqWeightTable is: " + termFreqWeightTable.size());
        endTime = System.currentTimeMillis();
        printTime("TermFreqWeight calculated and stored in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
    }

    /**
     * Function to read from disk the precomputed TermFreqWeight.
     */
    public static void readTFWeightFromDisk()
    {
        long startTime, endTime;
        int tf = 0;

        printLoad("Loading all useful values for termFreqWeight from disk...");

        startTime = System.currentTimeMillis();
        CollectionStatistics.readCollectionStatsFromDisk();
        logTermFreq = new double[CollectionStatistics.getMaxTermFreq()];

        try (
                RandomAccessFile docStats = new RandomAccessFile(TERMFREQWEIGHT_FILE, "rw");
                FileChannel channel = docStats.getChannel()
        ) {
            MappedByteBuffer buffer = channel.map(FileChannel.MapMode.READ_ONLY, 0, channel.size());

            if(buffer == null)      // Buffer not created
                return;

            // for to read all termFreqWeight stored into disk
            for (int i = 0; i < channel.size(); i += DOUBLE_BYTES)
            {
                logTermFreq[tf] = buffer.getDouble();
                tf++;
            }
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
        endTime = System.currentTimeMillis();
        printTime("TermFreqWeight loaded in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
    }

    /**
     * Function to calculate BM25 for one term and one document. The complete formula is:
     * BM25(t,d) = ( tf / (k * ((1 - b) + b * (docLen / avgDocLen)) + termFreq) ) * log10(NDoc/postListLen(t))
     *
     * @param DocID         DocID of the document processed
     * @param termFreq      term frequency of the term in the document
     * @param IDFweight     the IDF weight for the current term (precomputed for optimization)
     * @return  the BM25 score for one term and one document. The total score for a document will be the sum of the
     *          result of this function for each term that is both in the document and in the query
     */
    private static Double ScoringBM25(int DocID, int termFreq, double IDFweight)
    {
        double denominator, scoreBM25;      // variables to calculate the BM25 score value

        denominator = documentTable.get(DocID).getDenomPartBM25() + termFreq;   // take part of denominator from memory
        scoreBM25 = (termFreq / denominator) * IDFweight;      // calculate TFIDF weight from Tf and IDF weight values
        //printDebug("ScoringBM25 - docLen = " + docLen + " denominator = " + denominator + " IDFweight = " + IDFweight + " scoreBM25 = " + scoreBM25);
        return scoreBM25;
    }

    /**
     * Function that given a term return the corresponding term upper bound
     *
     * @param term          the term to query (for obtain the term upper bound)
     * @param scoringFunc   indicates the preference for scoring. if false use TFIDF, if true use BM25.
     * @param computeStats  indicates whether or not to compute statistics on occurrences of term freq values in the
     *                      collection. It is usually set to true only after the inverted index is computed.
     * @return the term upper bound for the term passed as parameter
     */
    public static Double maxScoreTerm(String term, boolean scoringFunc, boolean computeStats)
    {
        ArrayList<Posting>[] postingLists;  // contains all the posting lists for each term of the query
        ArrayList<String> processedQuery;   // array list for containing the query term
        SkipList[] skipListArray;           // array of the Skip List reference related to the term of query
        double maxScore = 0;                // the var containing the max value for the score of a document for the term
        double partialScore = 0;            // var that contain partial score
        double IDFWeight = 0;               //

        processedQuery = new ArrayList<>(1);
        processedQuery.add(term);           // insert the term

        if (!dictionary.getTermToTermStat().containsKey(processedQuery.get(0)))
            return maxScore;

        if(Flags.considerSkippingBytes() && Flags.isCompressionEnabled())
        {
            skipListArray = skipListInitCompAndSkip(processedQuery,null, false);
            postingLists = retrieveAllUncompPL(processedQuery, skipListArray); // get the uncompress PL
        }
        else
            postingLists = retrieveAllPostListsFromQuery(processedQuery);   // take posting list of the term

        // control check for term that is not in the dictionary
        if (postingLists[0] == null)
            return maxScore;            // return 0

        IDFWeight = Math.log10(((double) CollectionStatistics.getNDocs() / postingLists[0].size()));    // calculate IDFweight

        // execute the term query -> optimization -> there is one posting list -> the DID are already sort
        if (scoringFunc)    // use BM25
        {
            for (Posting p : postingLists[0])
            {
                // calculate SCORE (BM25) for this term and currentDID and sum to partial score
                partialScore = ScoringBM25(p.getDocId(), p.getTermFreq(), IDFWeight);

                // save score if is the new current term upper bound
                if ( partialScore > maxScore)
                {
                    maxScore = partialScore;
                }
                // compute term freq statistics
                if (computeStats)
                    CollectionStatistics.addTFOccToTermFreqTable(p.getTermFreq());
            }
        }
        else                // use TFIDF
        {
            int maxTermFreq = 0;    // indicates the current max term frequency value
            int currTermFreq;       // indicates the current value for the term frequency
            for (Posting p : postingLists[0])
            {
                // in this case with TFIDF all the components for calculating the score are constant except the TFIDF
                // values, so to know whether one doc has a higher score than another just see if it has a higher term
                // frequency. Thus, the term frequencies are compared directly and the score is only calculated if there
                // is a new document with termFreq (score) greater than the current maximum.
                currTermFreq = p.getTermFreq();
                // check if is the new current term upper bound
                if ( currTermFreq > maxTermFreq)
                {
                    // calculate SCORE (TFIDF) for this term and currentDID and sum to partial score
                    partialScore = ScoringTFIDF(currTermFreq, IDFWeight);

                    maxTermFreq = currTermFreq; // update the max term freq
                    maxScore = partialScore;    // save new TUB
                }

                // compute term freq statistics
                if (computeStats)
                    CollectionStatistics.addTFOccToTermFreqTable(p.getTermFreq());
            }
        }

        return maxScore;        // return term upper bound
    }
    // ------------------------ end: scoring function ------------------------

    // ------------------------ start: function to retrieve PL of the term in query ------------------------

    /**
     * Function to retrieve all the posting lists for each term of the query passed as parameter.
     * This function is used when compression and skipping are not both enabled.
     *
     * @param processedQuery    ArrayList of the processed terms of the query
     * @return  an array of posting lists (ArrayList of posting). the array has length equal to the number of terms,
     *          and the i-th position in the array contains the posting list of the i-th term in the processedQuery
     */
    private static ArrayList<Posting>[] retrieveAllPostListsFromQuery(ArrayList<String> processedQuery)
    {
        // array of arrayList (posting list) that contain all the posting lists for each term in the query
        ArrayList<Posting>[] postingLists = new ArrayList[processedQuery.size()];
        int iterator = 0;               // iterator for saving posting lists term in correct position

        try(
                // open complete files to read the postingList
                RandomAccessFile docidFile = new RandomAccessFile(DOCID_FILE, "rw");
                RandomAccessFile termfreqFile = new RandomAccessFile(TERMFREQ_FILE, "rw");

                // FileChannel
                FileChannel docIdChannel = docidFile.getChannel();
                FileChannel termFreqChannel = termfreqFile.getChannel()
        ) {
            // take posting list for each term in query
            for (String term : processedQuery)
            {
                //printDebug("DAAT: retrieve posting list of  " + term);
                DictionaryElem de = dictionary.getTermToTermStat().get(term);
                // take the postingList of term
                if(Flags.isCompressionEnabled() && !Flags.considerSkippingBytes())  // if the compression is enabled and the skipping is not enabled
                    postingLists[iterator] = readCompressedPostingListFromDisk(de.getOffsetDocId(),de.getOffsetTermFreq(), de.getTermFreqSize(), de.getDocIdSize(), de.getDf(), docIdChannel, termFreqChannel); //read compressed posting list
                else    // take the postingList of term
                    postingLists[iterator] = readPostingListFromDisk(de.getOffsetDocId(),de.getOffsetTermFreq(),de.getDf(),docIdChannel,termFreqChannel);
                iterator++;                 // update iterator
            }
            return postingLists;

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Function to retrieve all the posting lists for each term of the query passed as parameter and to sort by
     * term upper bound.
     *
     * @param processedQuery        ArrayList of the processed terms of the query
     * @param orderedQueryTerm      array of string to contain the query terms sorted by term upper bound
     * @param termUpperBoundList    array of double to contain the term upper bound of the query terms
     * @param sumTUBList            array of double to contain the sum of the term upper bound of the previous position
     * @return  an array of posting lists (ArrayList of posting). the array has length equal to the number of terms,
     *          the posting lists are sorted by term upper bound
     */
    private static ArrayList<Posting>[] retrieveAllPostingListsMaxScore(ArrayList<String> processedQuery, String[] orderedQueryTerm, double[] termUpperBoundList, double[] sumTUBList)
    {
        int procQLen = processedQuery.size();   // len of the query passed as parameter
        // array of arrayList (posting list) that contain all the posting lists for each term in the query
        ArrayList<Posting>[] postingLists = new ArrayList[procQLen];
        // priority queue for ordering the posting list according to term upper bound
        PriorityQueue<QueryProcessor.TermUpperBoundBlock> pq = new PriorityQueue<>(procQLen, new CompareTUBTerm());
        TermUpperBoundBlock tempTUBblock;   // the current element of the pq taken at each iteration
        double sumTUB = 0;                  // used to calculate the value for the 'sumTUBList'
        int iterator = 0;                   // iterator for saving posting lists term in correct position

        // control check for the lengths
        if ( (orderedQueryTerm.length != procQLen) || (termUpperBoundList.length != procQLen) || (sumTUBList.length != procQLen))
        {
            printError("Error in retrieveAllPostingListsMaxScore: wrong length in orderedQueryTerm or in termUpperBoundList or in sumTUBList.");
            return postingLists;
        }

        // retrieve the term upper bound for each posting lists and put into PQ
        for (int i = 0; i < procQLen; i++)
            pq.add(new QueryProcessor.TermUpperBoundBlock(i, TermDocUpperBound.getTermUpperBound(processedQuery.get(i))));     // add to priority queue

        // extract the ordered posting lists and related terms and insert them in the array of the term
        for (int i = 0; i < procQLen; i++)
        {
            tempTUBblock = pq.poll();    // get block
            orderedQueryTerm[i] = processedQuery.get(tempTUBblock.getTermPosition());   // get the term
            termUpperBoundList[i] = tempTUBblock.getTermUpperBound();                   // get term upper bound
            sumTUBList[i] = sumTUB;                 // get the term upper bound sum of the previous posting lists
            sumTUB += termUpperBoundList[i];
        }

        // take the posting lists for the related term (ordered by term upper bound instead query order)
        try(
                // open complete files to read the postingList
                RandomAccessFile docidFile = new RandomAccessFile(DOCID_FILE, "rw");
                RandomAccessFile termfreqFile = new RandomAccessFile(TERMFREQ_FILE, "rw");

                // FileChannel
                FileChannel docIdChannel = docidFile.getChannel();
                FileChannel termFreqChannel = termfreqFile.getChannel()
        ) {
            // take posting list for each term in query
            for (String term : orderedQueryTerm)
            {
                //printDebug("DAATMaxScore: retrieve posting list of  " + term);
                DictionaryElem de = dictionary.getTermToTermStat().get(term);

                if(Flags.isCompressionEnabled())
                    postingLists[iterator] = readCompressedPostingListFromDisk(de.getOffsetDocId(),de.getOffsetTermFreq(), de.getTermFreqSize(), de.getDocIdSize(), de.getDf(), docIdChannel, termFreqChannel); //read compressed posting list
                else    // take the postingList of term
                    postingLists[iterator] = readPostingListFromDisk(de.getOffsetDocId(),de.getOffsetTermFreq(),de.getDf(),docIdChannel,termFreqChannel);

                iterator++;                 // update iterator
            }
            return postingLists;

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Function to set all the utility arrays for the Max Score algorithm in the case of both compression and
     * skipping are enabled. Set all the array passed as parameters except processedQuery.
     *
     * @param processedQuery        ArrayList of the processed terms of the query
     * @param orderedQueryTerm      array of string to contain the query terms sorted by term upper bound
     * @param termUpperBoundList    array of double to contain the term upper bound of the query terms
     * @param sumTUBList            array of double to contain the sum of the term upper bound of the previous position
     */
    private static void setAllUtilsListMAxScoreCompAndSkipping(ArrayList<String> processedQuery, String[] orderedQueryTerm, double[] termUpperBoundList, double[] sumTUBList)
    {
        int plsSize = processedQuery.size();
        // priority queue for ordering the posting list according to term upper bound
        PriorityQueue<QueryProcessor.TermUpperBoundBlock> pq = new PriorityQueue<>(plsSize, new CompareTUBTerm());
        TermUpperBoundBlock tempTUBblock;   // the current element of the pq taken at each iteration
        double sumTUB = 0;                  // used to calculate the value for the 'sumTUBList'

        // control check
        if ( (orderedQueryTerm.length != plsSize) || (termUpperBoundList.length != plsSize) || (sumTUBList.length != plsSize))
            printError("Error in retrieveAllPostingListsMaxScore: wrong length in orderedQueryTerm or in termUpperBoundList or in sumTUBList.");

        // retrieve the term upper bound for each posting lists and put into PQ
        for (int i = 0; i < plsSize; i++)
            pq.add(new QueryProcessor.TermUpperBoundBlock(i, TermDocUpperBound.getTermUpperBound(processedQuery.get(i))));     // add to priority queue

        // extract the ordered posting lists and related terms and insert them in the array of the term
        for (int i = 0; i < plsSize; i++)
        {
            tempTUBblock = pq.poll();               // get block
            orderedQueryTerm[i] = processedQuery.get(tempTUBblock.getTermPosition());   // get the term
            termUpperBoundList[i] = tempTUBblock.getTermUpperBound();                   // get term upper bound
            sumTUBList[i] = sumTUB;                 // get the term upper bound sum of the previous posting lists
            sumTUB += termUpperBoundList[i];        // update sumTUB
        }
    }

    /**
     * Function that reads a compressed block of a posting list from disk, decompresses it and inserts it into the posting
     * arraylist array passed as a parameter. This function is used when both compression and skipping are enabled.
     *
     * @param term          term related to the posting list
     * @param slArr         the array of the SkipList of the query terms
     * @param blockIndex    the position of the block to load from disk
     * @param indexPL       the index of the posting list
     * @param pls           where to put the uncompressed posting lists
     */
    private static void retrieveCompBlockOfPLAndUncompress(String term, SkipList[] slArr, int blockIndex, int indexPL, ArrayList<Posting>[] pls)
    {
        byte[] tf;                      // array for the compressed TermFreq list
        byte[] docids;                  // array for the compressed TermFreq list
        // control check
        if ( (indexPL >= slArr.length) || (blockIndex >= slArr[indexPL].getSkipArrLen()))
            return;

        try(
                // open complete files to read the postingList
                RandomAccessFile docidFile = new RandomAccessFile(DOCID_FILE, "rw");
                RandomAccessFile termfreqFile = new RandomAccessFile(TERMFREQ_FILE, "rw");
                // FileChannel
                FileChannel docIdChannel = docidFile.getChannel();
                FileChannel termFreqChannel = termfreqFile.getChannel()
        ) {
            DictionaryElem de = dictionary.getTermToTermStat().get(term);
            // get the compressed block (the control check of the block index is in this function)
            tf = readCompTFBlockFromDisk(slArr[indexPL], blockIndex,de.getOffsetTermFreq(), de.getTermFreqSize(), de.getSkipArrLen(), termFreqChannel);
            docids = readCompDIDBlockFromDisk(slArr[indexPL], blockIndex, de.getOffsetDocId(), de.getDocIdSize(), de.getSkipArrLen(), docIdChannel);
            int numTFComp = min(SKIP_POINTERS_THRESHOLD, (de.getDf() - (SKIP_POINTERS_THRESHOLD * blockIndex)));
            if ( (tf == null) || (docids == null) )
                return;
            // decompress the block
            ArrayList<Integer> uncompressedTf = Unary.integersDecompression(tf, numTFComp);  // decompress term freq
            ArrayList<Integer> uncompressedDocid = VariableBytes.integersDecompression(docids,true);    // decompress DocID
            // add the pl
            pls[indexPL] = new ArrayList<Posting>(numTFComp);
            for (int i = 0; i < numTFComp; i++)
                pls[indexPL].add(new Posting(uncompressedDocid.get(i),uncompressedTf.get(i)));     // add arraylist
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Function that reads a compressed block of a posting list from disk, decompresses it and inserts it into the posting
     * arraylist array passed as a parameter. This function is used when both compression and skipping are enabled in
     * MaxScore algorithm.
     *
     * @param term          term related to the posting list
     * @param slArr         the array of the SkipList of the query terms
     * @param blockIndex    the position of the block to load from disk
     * @param indexPL       the index of the posting list
     * @return              'true' if a block has been read, decompressed, and added, 'false' otherwhise
     */
    private static boolean retrieveCompBlockOfPLAndUncompressMaxScore(String term, SkipList[] slArr, int blockIndex, int indexPL)
    {
        byte[] tf;                      // array for the compressed TermFreq list
        byte[] docids;                  // array for the compressed TermFreq list
        // control check
        if ( (indexPL >= slArr.length) || (blockIndex < 0) || (blockIndex >= slArr[indexPL].getSkipArrLen()))
        {
            return false;
        }
        try(
                // open complete files to read the postingList
                RandomAccessFile docidFile = new RandomAccessFile(DOCID_FILE, "rw");
                RandomAccessFile termfreqFile = new RandomAccessFile(TERMFREQ_FILE, "rw");
                // FileChannel
                FileChannel docIdChannel = docidFile.getChannel();
                FileChannel termFreqChannel = termfreqFile.getChannel()
        ) {
            DictionaryElem de = dictionary.getTermToTermStat().get(term);
            // get the compressed block (the control check of the block index is in this function)
            tf = readCompTFBlockFromDisk(slArr[indexPL], blockIndex,de.getOffsetTermFreq(), de.getTermFreqSize(), de.getSkipArrLen(), termFreqChannel);
            docids = readCompDIDBlockFromDisk(slArr[indexPL], blockIndex, de.getOffsetDocId(), de.getDocIdSize(), de.getSkipArrLen(), docIdChannel);
            int numTFComp = min(SKIP_POINTERS_THRESHOLD, (de.getDf() - (SKIP_POINTERS_THRESHOLD * blockIndex)));
            if ( (tf == null) || (docids == null) )
                return false;
            // decompress the block
            ArrayList<Integer> uncompressedTf = Unary.integersDecompression(tf, numTFComp);  // decompress term freq
            ArrayList<Integer> uncompressedDocid = VariableBytes.integersDecompression(docids,true);    // decompress DocID

            // add the block to the related PL
            skipAndCompPLs[indexPL] = new ArrayList<>(numTFComp);   // decompressed posting list
            for (int i = 0; i < numTFComp; i++)             // add the postings of the block to global PLs
                skipAndCompPLs[indexPL].add(new Posting(uncompressedDocid.get(i), uncompressedTf.get(i)));  // add
            slArr[indexPL].setCurrPostList(skipAndCompPLs[indexPL]);  // set the current uncompressed block of PL in SkipList instance
            return true;
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Function to retrieve and uncompress all the posting lists for each term of the query passed as parameter.
     * This function is used when compression and skipping are both enabled.
     *
     * @param processedQuery    ArrayList of the processed terms of the query
     * @param slArr             the array of the SkipInfo related to the query term
     * @return  an array of posting lists (ArrayList of posting). the array has length equal to the number of terms,
     *          and the i-th position in the array contains the posting list of the i-th term in the processedQuery
     */
    private static ArrayList<Posting>[] retrieveAllUncompPL(ArrayList<String> processedQuery, SkipList[] slArr)
    {
        // array of arrayList (posting list) that contain all the posting lists for each term in the query
        ArrayList<Posting>[] postingLists = new ArrayList[processedQuery.size()];
        int iterator = 0;               // iterator for saving posting lists term in correct position

        try(
                // open complete files to read the postingList
                RandomAccessFile docidFile = new RandomAccessFile(DOCID_FILE, "rw");
                RandomAccessFile termfreqFile = new RandomAccessFile(TERMFREQ_FILE, "rw");
                // FileChannel
                FileChannel docIdChannel = docidFile.getChannel();
                FileChannel termFreqChannel = termfreqFile.getChannel()
        ) {
            // take posting list for each term in query
            for (String term : processedQuery)
            {
                DictionaryElem de = dictionary.getTermToTermStat().get(term);
                // read adn uncompress the whole posting list related to term
                postingLists[iterator] = readAndUncompressCompressedAndSkippedPLFromDisk(slArr[iterator], de.getOffsetDocId(), de.getOffsetTermFreq(), de.getTermFreqSize(), de.getDocIdSize(), de.getSkipArrLen(), de.getDf(), docIdChannel, termFreqChannel);
                iterator++;                 // update iterator
            }
            return postingLists;

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Function to retrieve the compressed all the posting lists for each term of the query passed as parameter.
     * This function is used when compression, skipping and wholePLInMem are enabled.
     * It reads all the posting lists for the query terms in their entirety and places them in memory (in global byte
     * arrays 'termFreqPL' and 'DocIDPL').
     *
     * @param processedQuery    ArrayList of the processed terms of the query
     */
    private static void retrieveAllCompPL(ArrayList<String> processedQuery)
    {
        int queryLen = processedQuery.size();   // number of term in the query
        int iterator = 0;                       // iterator for saving posting lists term in correct position
        long startOffsetTF;                     // the current offset (at each iteration) for the compressed TermFreq
        int sizeToReadTF = 0;                   // the current size (at each iteration) for the compressed TermFreq
        long startOffsetDID;                    // the current offset (at each iteration) for the compressed TermFreq
        int sizeToReadDID = 0;                  // the current size (at each iteration) for the compressed TermFreq

        termFreqPL = new byte[queryLen][];  // set the array
        docIDPL = new byte[queryLen][];     // set the array

        try(
                // open complete files to read the postingList
                RandomAccessFile docidFile = new RandomAccessFile(DOCID_FILE, "rw");
                RandomAccessFile termfreqFile = new RandomAccessFile(TERMFREQ_FILE, "rw");
                // FileChannel
                FileChannel docIdChannel = docidFile.getChannel();
                FileChannel termFreqChannel = termfreqFile.getChannel()
        ) {
            // take posting list for each term in query
            for (String term : processedQuery)
            {
                DictionaryElem de = dictionary.getTermToTermStat().get(term);
                // read the termFreq of the whole posting list
                startOffsetTF = (int) de.getOffsetTermFreq();
                sizeToReadTF = de.getTermFreqSize();

                MappedByteBuffer termFreqBuffer = termFreqChannel.map(FileChannel.MapMode.READ_ONLY, startOffsetTF, sizeToReadTF);
                termFreqPL[iterator] = new byte[sizeToReadTF];        // initialize the byte array for the compressed TF list
                termFreqBuffer.get(termFreqPL[iterator], 0, sizeToReadTF);       // read the TF compressed block of the PL
                // read the DID of the whole posting list
                startOffsetDID = (int) de.getOffsetDocId();
                sizeToReadDID = de.getDocIdSize();

                MappedByteBuffer docIdBuffer = docIdChannel.map(FileChannel.MapMode.READ_ONLY, startOffsetDID, sizeToReadDID);
                docIDPL[iterator] = new byte[sizeToReadDID];          // initialize the byte array for the compressed DID list
                docIdBuffer.get(docIDPL[iterator], 0, sizeToReadDID);   // read the DID compressed block of the PL

                iterator++;     // update the iterator
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Function to retrieve the compressed all the posting lists for each term of the query passed as parameter.
     * This function is used when compression, skipping and wholePLInMem are enabled.
     * It reads all the posting lists for the query terms in their entirety and places them in memory (in global byte
     * arrays 'termFreqPL' and 'DocIDPL').
     *
     * @param orderedQuery  Array of the processed and ordered terms of the query
     */
    private static void retrieveAllCompPLMaxScore(String[] orderedQuery)
    {
        int queryLen = orderedQuery.length;   // number of term in the query
        int iterator = 0;                       // iterator for saving posting lists term in correct position
        long startOffsetTF;                     // the current offset (at each iteration) for the compressed TermFreq
        int sizeToReadTF = 0;                   // the current size (at each iteration) for the compressed TermFreq
        long startOffsetDID;                    // the current offset (at each iteration) for the compressed TermFreq
        int sizeToReadDID = 0;                  // the current size (at each iteration) for the compressed TermFreq

        termFreqPL = new byte[queryLen][];  // set the array
        docIDPL = new byte[queryLen][];     // set the array

        try(
                // open complete files to read the postingList
                RandomAccessFile docidFile = new RandomAccessFile(DOCID_FILE, "rw");
                RandomAccessFile termfreqFile = new RandomAccessFile(TERMFREQ_FILE, "rw");
                // FileChannel
                FileChannel docIdChannel = docidFile.getChannel();
                FileChannel termFreqChannel = termfreqFile.getChannel()
        ) {
            // take posting list for each term in query
            for (String term : orderedQuery)
            {
                DictionaryElem de = dictionary.getTermToTermStat().get(term);
                // read the termFreq of the whole posting list
                startOffsetTF = (int) de.getOffsetTermFreq();
                sizeToReadTF = de.getTermFreqSize();

                MappedByteBuffer termFreqBuffer = termFreqChannel.map(FileChannel.MapMode.READ_ONLY, startOffsetTF, sizeToReadTF);
                termFreqPL[iterator] = new byte[sizeToReadTF];        // initialize the byte array for the compressed TF list
                termFreqBuffer.get(termFreqPL[iterator], 0, sizeToReadTF);       // read the TF compressed block of the PL
                // read the DID of the whole posting list
                startOffsetDID = (int) de.getOffsetDocId();
                sizeToReadDID = de.getDocIdSize();

                MappedByteBuffer docIdBuffer = docIdChannel.map(FileChannel.MapMode.READ_ONLY, startOffsetDID, sizeToReadDID);
                docIDPL[iterator] = new byte[sizeToReadDID];          // initialize the byte array for the compressed DID list
                docIdBuffer.get(docIDPL[iterator], 0, sizeToReadDID);   // read the DID compressed block of the PL

                iterator++;     // update the iterator
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Function used when wholePLInMem is enabled. Fetches a block from the selected compressed posting list.
     * Compressed posting lists are saved in global byte arrays, one for DID values (docIDPL) and one for term
     * frequency values (termFreqPL).
     * The decompressed list block will be added to the posting list arraylist array passed as a parameter.
     *
     * @param processedQuery    ArrayList of the processed terms of the query
     * @param slArr             the array of the SkipInfo related to the query term
     * @param blockIndex        the index of the block wanted
     * @param termIndex         the index of the PL from which to take the compressed block
     * @param pls               where to put the uncompressed posting lists
     */
    private static void UncompressOneBlockFromCompressedPL(ArrayList<String> processedQuery, SkipList[] slArr, int blockIndex, int termIndex, ArrayList<Posting>[] pls)
    {
        int startTF = 0;    // the start offset of the term freq value for the selected block of the PL
        int startDID = 0;   // the start offset of the DID value for the selected block of the PL
        int sizeTF = 0;     // number of bytes to read of the term freq value for the selected block of the PL
        int sizeDID = 0;    // number of bytes to read of the DID value for the selected block of the PL
        int numTFComp;      // the number of termFreq value in the selected block of the PL

        // control check
        if ((blockIndex < 0) || (termIndex < 0) || (termIndex >= processedQuery.size()) || (blockIndex >= slArr[termIndex].getSkipArrLen()))
            return;     // exit

        DictionaryElem de = dictionary.getTermToTermStat().get(processedQuery.get(termIndex));  // get dictionary elem
        if (blockIndex > 0)    // compute start offset
        {
            startTF = (int)slArr[termIndex].getSkipBlockInfo(blockIndex-1).getFreqOffset() - (int) de.getOffsetTermFreq();
            startDID = (int)slArr[termIndex].getSkipBlockInfo(blockIndex-1).getDocIdOffset() - (int) de.getOffsetDocId();
        }
        sizeTF = ((int)slArr[termIndex].getSkipBlockInfo(blockIndex).getFreqOffset() - (int) de.getOffsetTermFreq());
        sizeDID = ((int)slArr[termIndex].getSkipBlockInfo(blockIndex).getDocIdOffset() - (int) de.getOffsetDocId());
        numTFComp = min(SKIP_POINTERS_THRESHOLD, (de.getDf() - (SKIP_POINTERS_THRESHOLD * blockIndex)));    // number to read
        // uncompress
        ArrayList<Integer> uncompressedTf = Unary.integersDecompression(Arrays.copyOfRange(termFreqPL[termIndex], startTF, sizeTF), numTFComp);  // decompress term freq
        ArrayList<Integer> uncompressedDocid = VariableBytes.integersDecompression(Arrays.copyOfRange(docIDPL[termIndex], startDID, sizeDID),true);    // decompress DocID
        // add the pl
        pls[termIndex] = new ArrayList<Posting>(numTFComp);
        for (int i = 0; i < numTFComp; i++)
            pls[termIndex].add(new Posting(uncompressedDocid.get(i),uncompressedTf.get(i)));     // add arraylist
    }

    /**
     * Function used when wholePLInMem is enabled. Fetches a block from the selected compressed posting list.
     * Compressed posting lists are saved in global byte arrays, one for DID values (docIDPL) and one for term
     * frequency values (termFreqPL).
     * The decompressed list block will be added to the global posting lists.
     *
     * @param processedQuery    ArrayList of the processed terms of the query
     * @param slArr             the array of the SkipInfo related to the query term
     * @param blockIndex        the index of the block wanted
     * @param termIndex         the index of the PL from which to take the compressed block
     * @return                  'true' if a block has been read, decompressed, and added, 'false' otherwhise
     */
    private static boolean UncompressOneBlockFromCompressedPLMAxScore(String[] processedQuery, SkipList[] slArr, int blockIndex, int termIndex)
    {
        int startTF = 0;    // the start offset of the term freq value for the selected block of the PL
        int startDID = 0;   // the start offset of the DID value for the selected block of the PL
        int sizeTF = 0;     // number of bytes to read of the term freq value for the selected block of the PL
        int sizeDID = 0;    // number of bytes to read of the DID value for the selected block of the PL
        int numTFComp;      // the number of termFreq value in the selected block of the PL

        // control check
        if ((blockIndex < 0) || (termIndex < 0) || (termIndex >= processedQuery.length) || (blockIndex >= slArr[termIndex].getSkipArrLen()))
            return false;     // exit

        DictionaryElem de = dictionary.getTermToTermStat().get(processedQuery[termIndex]);  // get dictionary elem
        if (blockIndex > 0)    // compute start offset
        {
            startTF = (int)slArr[termIndex].getSkipBlockInfo(blockIndex-1).getFreqOffset() - (int) de.getOffsetTermFreq();
            startDID = (int)slArr[termIndex].getSkipBlockInfo(blockIndex-1).getDocIdOffset() - (int) de.getOffsetDocId();
        }
        sizeTF = ((int)slArr[termIndex].getSkipBlockInfo(blockIndex).getFreqOffset() - (int) de.getOffsetTermFreq());
        sizeDID = ((int)slArr[termIndex].getSkipBlockInfo(blockIndex).getDocIdOffset() - (int) de.getOffsetDocId());
        numTFComp = min(SKIP_POINTERS_THRESHOLD, (de.getDf() - (SKIP_POINTERS_THRESHOLD * blockIndex)));    // number to read
        // uncompress
        ArrayList<Integer> uncompressedTf = Unary.integersDecompression(Arrays.copyOfRange(termFreqPL[termIndex], startTF, sizeTF), numTFComp);  // decompress term freq
        ArrayList<Integer> uncompressedDocid = VariableBytes.integersDecompression(Arrays.copyOfRange(docIDPL[termIndex], startDID, sizeDID),true);    // decompress DocID
        // add the block to the related PL
        skipAndCompPLs[termIndex] = new ArrayList<>(numTFComp);    // decompressed posting list
        for (int i = 0; i < numTFComp; i++)             // add the postings of the block to global PLs
            skipAndCompPLs[termIndex].add(new Posting(uncompressedDocid.get(i), uncompressedTf.get(i)));  // add
        slArr[termIndex].setCurrPostList(skipAndCompPLs[termIndex]);  // set the current uncompressed block of PL in SkipList instance
        return true;
    }
    // ------------------------ end: function to retrieve PL of the term in query ------------------------

    // ------------------------ start: utilities function ------------------------

    /**
     * Function to elaborate all docs and related scores to obtain the ranked list of results
     *
     * @param numResults    number of result(DocID) to return to the user
     * @return  an ordered ArrayList that represent the top numResults results for the query
     */
    private static ArrayList<Integer> getRankedResults(int numResults)
    {
        ArrayList<Integer> rankedResults = new ArrayList<>(numResults); // array list to contain the top "numResults" docs
        String currDocNO;           // indicates the DocNO of the current document in the result (top 'numResults' docs)
        long startTime, endTime;    // variables to calculate the execution time

        if (numResults <= 0)        // control check
            return rankedResults;

        startTime = System.currentTimeMillis();         // start time of hash map ordering

        QueryProcessor.ResultBlock currentResPQ;        // var that contain the resultBlock extract from pq in the current iteration
        while(!resPQ.isEmpty())                         // control if the priority queue for results is empty
        {
            currentResPQ = resPQ.poll();                                    // take the lowest element (score and DID)
            currDocNO = documentTable.get(currentResPQ.getDID()).getDocno();// take the DocNo related to the DID
            //printDebug("DID: " + currentResPQ.getDID() + ", DocNO: " + currDocNO + ", Score: " + currentResPQ.getScore());
            try{
                rankedResults.add(Integer.valueOf(currDocNO));                    // add to the array list
            }
            catch (NumberFormatException ex){
                ex.printStackTrace();
            }
        }
        // order the result from the best to the worst (reverse order of the priority queue)
        Collections.reverse(rankedResults);

        endTime = System.currentTimeMillis();           // end time of hash map ordering
        printTime("*** Ranked results (results priority queue) in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

        return rankedResults;
    }

    /**
     * Function to retrieve length for each posting lists passed as parameter.
     *
     * @param terms    array of the terms
     * @return  an array of posting lists (ArrayList of posting). the array has length equal to the number of terms,
     *          and the i-th position in the array contains the posting list of the i-th term in the processedQuery
     */
    private static int[] retrieveLengthAllPostingLists(String[] terms)
    {
        int[] lengthPostingList = new int[terms.length];

        for (int i = 0; i < terms.length; i++)
        {
            lengthPostingList[i] = dictionary.getTermToTermStat().get(terms[i]).getDf();
        }

        return lengthPostingList;
    }

    /**
     * function to update the essential posting lists by the threshold passed as parameter
     *
     * @param sumTUBList        array of double to contain the sum of the term upper bound of the previous position
     * @param threshold         the current threshold passes as parameter
     * @param currentFirstEPL   the current position of the first Essential PL (index of beginning, to be updated)
     * @return  an integer that indicate the index of the first new essential posting list. The new essential posting
     *          list will be the posting lists with index between the returned integer and the last posting lists
     */
    private static int updateEssentialPositngLists (double[] sumTUBList, double threshold, int currentFirstEPL)
    {
        for (int i = (currentFirstEPL + 1); i < sumTUBList.length; i++)
        {
            if (sumTUBList[i] > threshold)
                return i - 1;
        }

        return sumTUBList.length-1;
    }

    /**
     * Function to initialize the SkipList objects related to the query terms. Used in case of only skipping enabled.
     *
     * @param orderedQueryTerm  the term of the query ordered by Term Upper Bound
     * @param postingLists      the Posting Lists of the query term
     * @return      the arraylist of SkipList related to the query term
     */
    private static SkipList[] SetAllSkipList(String[] orderedQueryTerm, ArrayList<Posting>[] postingLists)
    {
        SkipList[] tempSkipListArray;   // array of skipList, the skip list at i-th position is related to the i-th term of the query
        DictionaryElem tempDictElem;    // temp dictionary elem

        if ( (orderedQueryTerm.length == 0) || (postingLists.length == 0) )   // control check
            return null;

        tempSkipListArray = new SkipList[orderedQueryTerm.length];      // set array

        for (int i = 0; i < postingLists.length; i++)       // scan all query terms
        {
            //printDebug("SetAllSkipList -- SkipList related to the term: " + orderedQueryTerm[i] + " in position: " + i);
            tempDictElem = dictionary.getTermStat(orderedQueryTerm[i]); // get dictionary elem associated to term
            // create the skipList related to query's term
            tempSkipListArray[i] = new SkipList(tempDictElem.getSkipOffset(), tempDictElem.getSkipArrLen(),postingLists[i], tempDictElem.getDf());
        }

        return tempSkipListArray;
    }

    /**
     * Function to initialize the SkipList objects related to the query terms.
     * Used in case of both skipping and compression are enabled.
     *
     * @param processedQuery    the query of the user
     * @param orderedQueryTerm  the query term ordered by their Term Upper Bound (used in case of Max Score enabled)
     * @param maxScore          if 'true' Max Score enabled, if 'false' Max Score is not enabled
     * @return      the array of SkipList related to the query term
     */
    private static SkipList[] skipListInitCompAndSkip(ArrayList<String> processedQuery, String[] orderedQueryTerm, boolean maxScore)
    {
        SkipList[] tempSkipListArray;   // array of skipList, the skip list at i-th position is related to the i-th term of the query
        DictionaryElem tempDictElem;    // temp dictionary elem

        if (maxScore)       // Max Score enabled
        {
            int queryLen = orderedQueryTerm.length;
            if (queryLen == 0)   // control check
                return null;

            tempSkipListArray = new SkipList[queryLen];     // set array

            for (int i = 0; i < queryLen; i++)              // scan all query terms
            {
                if (dictionary.getTermToTermStat().containsKey(orderedQueryTerm[i]))  // check if the term is in the dictionary
                {
                    //printDebug("SetAllSkipList -- SkipList related to the term: " + orderedQueryTerm[i] + " in position: " + i);
                    tempDictElem = dictionary.getTermStat(orderedQueryTerm[i]); // get dictionary elem associated to term
                    // create the skipList related to query's term
                    tempSkipListArray[i] = new SkipList(tempDictElem.getSkipOffset(), tempDictElem.getSkipArrLen(),null, tempDictElem.getDf());
                }
            }
        }
        else                // simple DAAT
        {
            if (processedQuery.isEmpty())   // control check
                return null;

            tempSkipListArray = new SkipList[processedQuery.size()];      // set array

            for (int i = 0; i < processedQuery.size(); i++)       // scan all query terms
            {
                if (dictionary.getTermToTermStat().containsKey(processedQuery.get(i)))  // check if the term is in the dictionary
                {
                    //printDebug("SetAllSkipList -- SkipList related to the term: " + orderedQueryTerm[i] + " in position: " + i);
                    tempDictElem = dictionary.getTermStat(processedQuery.get(i)); // get dictionary elem associated to term
                    // create the skipList related to query's term
                    tempSkipListArray[i] = new SkipList(tempDictElem.getSkipOffset(), tempDictElem.getSkipArrLen(),null, tempDictElem.getDf());
                }
            }
        }

        return tempSkipListArray;
    }

    /**
     * Function that given the posting lists of each term in a given query returns an ordered list of the DocIDs
     * present in the all posting lists
     *
     * @param postingLists  the posting lists of each term in the query
     * @return  an ordered ArrayList of the DocIDs in the posting lists
     */
    private static ArrayList<Integer> DIDOrderedListOfQuery(ArrayList<Posting>[] postingLists, boolean isConjunctive)
    {
        PriorityQueue<Integer> tempPQMax = new PriorityQueue<>(Comparator.reverseOrder());  // priority queue from which to take the DID which will become the new max
        int[] postingListsIndex = getPostingListsIndex(postingLists);   // contain the current position index for the posting list of each term in the query
        ArrayList<Integer> orderedList = new ArrayList<>(); // the final ordered list of the DIDs
        TreeSet<Integer> tempList = new TreeSet<>();        // treeSet to temporarily contain the DIDs in an orderly manner
        int plsLen = postingLists.length;                   // take the number of posting list
        int[] lengthPostingList = new int[plsLen];          // array containing the length of the posting lists
        boolean allPostListScanned = false;     // indicate if all posting list are fully scanned
        boolean postingEnded = false;           // used only in conjunctive query (for optimization), indicates that one posting list is fully scanned
        int max = 0;                            // indicates the current max DID taken
        int currentDocID = 0;                   // var to contain the current DocID
        long startTime,endTime;                 // variables to calculate the execution time

        startTime = System.currentTimeMillis();             // start time to take th DocID list
        // take the len of each PL and the first DocID from PLs
        for (int i = 0; i < plsLen; i++)
        {
            lengthPostingList[i] = postingLists[i].size();          // add the len
            tempPQMax.add(postingLists[i].get(0).getDocId());       // add the DID into tempList
        }
        max = tempPQMax.peek(); // Get the maximum DID
        tempPQMax.clear();

        int currPLIndex = 0;                    // var for the index of the current PL at the current iteration
        if(isConjunctive)   // conjunctive
        {   // start - if - conjunctive -
            while(true)  // scan all posting list and insert DIDs into orderedList
            {   // start - main while - conj -
                allPostListScanned = true;      // set var, if not reset the while end
                // for each posting list I take all values less than max and add them to the hash map.
                for (int i = 0; i < plsLen; i++)
                {
                    currPLIndex = postingListsIndex[i];         // take the current index for the current posting list

                    if (currPLIndex >= lengthPostingList[i])    // posting list completely visited
                        postingEnded = true;    // set var
                    else                // posting list not completely visited
                    {
                        allPostListScanned = false;  // there is at least one posting not seen yet -> set var
                        currentDocID = postingLists[i].get(currPLIndex).getDocId();    // take current DID
                        // scan the current posting list until a DID greater than max (avoiding overflow)
                        while((currentDocID <= max) && (currPLIndex < lengthPostingList[i]))
                        {
                            tempList.add(currentDocID);             // take the DID
                            currPLIndex = ++postingListsIndex[i];   // update and take the current index for the current posting list
                            if (currPLIndex < lengthPostingList[i])
                                currentDocID = postingLists[i].get(currPLIndex).getDocId();    // take current DID
                        }
                    }
                }
                // in tempList there are all DID lower than max (and max) in ascending  order
                orderedList.addAll(tempList);       // add ordered DID in orderedList
                tempList.clear();                   // clear tempList

                // if the query is conjunctive and at least one posting list is fully scanned
                if (postingEnded && allPostListScanned)
                    break;      // exit to while. Stop the collection ofs DID, I have all the DIDs I need
                else    // no posting list has been fully visited, same number of PLs to be visited as at the beginning
                {
                    // take the new max (from the next DID of each posting list)
                    for (int i = 0; i < plsLen; i++)
                    {
                        if (postingListsIndex[i] < lengthPostingList[i])        // posting list not fully visited
                            tempPQMax.add(postingLists[i].get(postingListsIndex[i]).getDocId());    // add the DID into tempList
                    }
                    if (!tempPQMax.isEmpty())
                    {
                        max = tempPQMax.peek(); // Get the maximum DID
                        tempPQMax.clear();
                    }
                    else
                        break;  // exit to while
                }
            }   // end - main while - conj -
        }   // end - if - conjunctive -
        else    // disjunctive
        {   // start - else - disjunctive -
            while(!allPostListScanned)  // scan all posting list and insert DIDs into orderedList
            {   // start - main while - disj -
                allPostListScanned = true;      // set var, if not reset the while end
                // for each posting list I take all values less than max and add them to the hash map.
                for (int i = 0; i < plsLen; i++)
                {
                    currPLIndex = postingListsIndex[i];     // take the current index for the current posting list

                    if (currPLIndex < lengthPostingList[i])       // posting list not completely visited
                    {
                        allPostListScanned = false;  // there is at least one posting not seen yet -> set var
                        currentDocID = postingLists[i].get(currPLIndex).getDocId();    // take current DID
                        // scan the current posting list until a DID greater than max (avoiding overflow)
                        while((currentDocID <= max) && (currPLIndex < lengthPostingList[i]))
                        {
                            tempList.add(currentDocID);             // take the DID
                            currPLIndex = ++postingListsIndex[i];   // update and take the current index for the current posting list
                            if (currPLIndex < lengthPostingList[i])
                                currentDocID = postingLists[i].get(currPLIndex).getDocId();    // take current DID
                        }
                    }
                }
                // in tempList there are all DID lower than max (and max) in ascending  order
                orderedList.addAll(tempList);       // add ordered DID in orderedList
                tempList.clear();                   // clear tempList

                int numPLRemaining = 0;         // number of posting list that are not fully scanned
                // take the new max (from the next DID of each posting list)
                for (int i = 0; i < plsLen; i++)
                {
                    if (postingListsIndex[i] < lengthPostingList[i])        // posting list not fully visited
                    {
                        tempPQMax.add(postingLists[i].get(postingListsIndex[i]).getDocId()); // add the DID into tempList
                        ++numPLRemaining;
                    }
                }

                // check if thee is only one posting list not fully scanned
                if (numPLRemaining == 1)
                {
                    int index = 0;
                    // take the index of the not fully scanned posting lists
                    for (int i = 0; i < plsLen; i++)
                    {
                        if (postingListsIndex[i] < lengthPostingList[i])     // posting list not fully visited
                        {
                            index = i;  // take the index of the not fully scanned posting list
                            break;      // find the searched PL
                        }
                    }
                    // take all remaining DID in the posting list
                    while(postingListsIndex[index] < lengthPostingList[index])
                    {
                        orderedList.add(postingLists[index].get(postingListsIndex[index]).getDocId());  // add DID into orederdList
                        ++postingListsIndex[index];     // increment the related index
                    }
                    break;      // exit to while
                }
                else if ( (numPLRemaining != 0) && (!tempPQMax.isEmpty()) )      // there are more than one posting lists not fully scanned (and not is the last iteration of the while with tempPQMax empty)
                {
                    max = tempPQMax.peek(); // Get the maximum DID
                    tempPQMax.clear();
                }
                else    // tempList.size() = 0 -> is the last iteration of the while with tempList empty
                    break;      // exit to while
            }   // end - main while - disj -
        }   // end - else - disjunctive -
        endTime = System.currentTimeMillis();           // end time of DocID list ordering
        // shows query execution time
        printTime("*** ORDERED DID LIST (no PQ V.2.2) in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
        //printDebug("Ordered List (no PQ V.2.2) of DocID dim: " + orderedList.size());     // print orderedList
        return orderedList;
    }

    /**
     * Function to create an array of indexes for posting lists
     *
     * @param postingLists  the posting lists of each term in the query
     * @return  an array that contains the index for the current posting (position) for each posting lists of the term
     *          in the query
     */
    private static int[] getPostingListsIndex (ArrayList<Posting>[] postingLists)
    {
        int[] postingListsIndex = new int[postingLists.length];

        // set the index to 0 for each posting lists of the term in the query
        for (int i = 0; i < postingLists.length; i++)
        {
            postingListsIndex[i] = 0;       // set index to 0
        }

        return postingListsIndex;
    }

    /**
     * Function to create an array of term upper bound for posting lists(terms)
     *
     * @param processedQuery  ArrayList of the processed terms of the query
     * @return  an array that contains the term upper bound for each posting lists of the term in the query
     */
    private static double[] getPostingListsTermUpperBound (ArrayList<String> processedQuery)
    {
        double[] postingListsTermUB = new double[processedQuery.size()];

        // set the index to 0 for each posting lists of the term in the query
        for (int i = 0; i < processedQuery.size(); i++)
        {
            postingListsTermUB[i] = TermDocUpperBound.getTermUpperBound(processedQuery.get(i));       // set index to 0
            //printDebug("Pos: " + i + " term: " + processedQuery.get(i) + " term upper bound: " + postingListsTermUB[i]);
        }

        return postingListsTermUB;
    }

    /**
     * Function for boolean search, used in max score without skipping.
     *
     * @param tempList  list on which to do the binary search
     * @param targetDID the searched DID
     * @return      '-1' if the searched DID is not in the tempList
     *              'currentPos' the position of the searched DID in the tempList
     */
    private static int booleanSearch(List<Posting> tempList, int targetDID)
    {
        int startPos = 0;               // initial startPosition
        int endPos = tempList.size()-1; // initial endPosition
        int currentPos = 0;
        int currDID = 0;

        while (true)
        {
            if (startPos > endPos)
                return -1;          // not found

            currentPos = (startPos + endPos)/2;
            currDID = tempList.get(currentPos).getDocId();

            if (currDID == targetDID)
                return currentPos;
            else if (currDID > targetDID)
            {
                endPos = currentPos - 1;
            }
            else
            {
                startPos = currentPos + 1;
            }
        }
    }

    /**
     * Function to shows the user the ranked results (DocID) of the query executed
     *
     * @param rankedResults the results returned by the query
     */
    public static void printQueryResults(ArrayList<Integer> rankedResults)
    {
        if (!rankedResults.isEmpty())      // there are results
        {
            printUI("Query results:");
            for (int i = 0; i < rankedResults.size(); i++)
                printUI((i + 1) + " - " + rankedResults.get(i));
        }
        else                                // there aren't results
            printUI("No results found for this query.");
    }

    // -------- end: utilities function --------


    // -------- start: utilities for priority queue --------

    /**
     * Class to define orderedDIDBlock to insert in the priority queue.
     * These blocks and the queue will be used when compression and skipping are both enabled.
     * The compressed PLs will be decompressed block by block and the did and tf will be put in the priority queue and
     * then taken in order during DAAT.
     */
    private static class orderedDIDBlock
    {
        int DocID;                  // DocID
        int termFreq;               // the term frequency related the PL (of index indexPL) and the DocID
        int indexPL;                // the index of the term(PL) in the query

        // constructor with parameters
        public orderedDIDBlock(int DocID, int termFreq, int indexPL)
        {
            this.DocID = DocID;
            this.termFreq = termFreq;
            this.indexPL = indexPL;
        }

        public int getDocID() {
            return DocID;
        }
        public int getTermFreq() {
            return termFreq;
        }
        public int getIndexPL() {
            return indexPL;
        }

        @Override
        public String toString()
        {
            return "PB{" +
                    "\nDocID = '" + DocID + "'" +
                    ", term freq = '" + termFreq + "'" +
                    ", indexPL = '" + indexPL + "' }";
        }
    }

    /**
     * Class to compare the block, allows the order of the priority queue.
     */
    private static class CompareOrdDIDBlock implements Comparator<QueryProcessor.orderedDIDBlock>
    {
        @Override
        public int compare(QueryProcessor.orderedDIDBlock odb1, QueryProcessor.orderedDIDBlock odb2)
        {
            int scoreComparison = Double.compare(odb1.getDocID(), odb2.getDocID());     // comparing terms
            // if the term upper bound are equal, compare by position
            if (scoreComparison == 0) {
                // return order by both term upper bound and position (the TUB of the two blocks is equal)
                return Integer.compare(odb1.getIndexPL(), odb2.getIndexPL());
            }
            return scoreComparison;     // return order only by score (the score of the two blocks is different)
        }
    }

    /**
     * class to define termUpperBoundPostingList. The priority queue contains instances of termUpperBoundPostingList
     * representing the TUB for the term of the query.
     */
    private static class TermUpperBoundBlock
    {
        int termPosition;           // the position of the term in the query
        double termUpperBound;      // the value of the TermUpperBound

        // constructor with parameters
        public TermUpperBoundBlock(int termPosition, double termUpperBound)
        {
            this.termPosition = termPosition;
            this.termUpperBound = termUpperBound;
        }

        // get methods
        public int getTermPosition() {
            return termPosition;
        }

        public double getTermUpperBound() {
            return termUpperBound;
        }

        @Override
        public String toString() {
            return "PB{" +
                    "Query position = '" + termPosition + '\'' +
                    ", term upper bound =" + termUpperBound +
                    '}';
        }
    }

    /**
     * Class to compare the block, allows the order of the priority queue
     * The order is ascending order according to term upper bound. In case of equal term upper bound values will be made
     * order according to query position of the term. Term upper bound and position are sorted in ascending order.
     */
    private static class CompareTUBTerm implements Comparator<QueryProcessor.TermUpperBoundBlock> {
        @Override
        public int compare(QueryProcessor.TermUpperBoundBlock tub1, QueryProcessor.TermUpperBoundBlock tub2) {
            // comparing terms
            int scoreComparison = Double.compare(tub1.getTermUpperBound(), tub2.getTermUpperBound());
            // if the term upper bound are equal, compare by position
            if (scoreComparison == 0) {
                // return order by both term upper bound and position (the TUB of the two blocks is equal)
                return Integer.compare(tub1.getTermPosition(), tub2.getTermPosition());
            }

            return scoreComparison;     // return order only by score (the score of the two blocks is different)
        }
    }
    // -------- end: utilities for priority queue --------

    // -------- start: utilities for priority queue for the results --------
    /**
     * class to define ResultBlock. The priority queue contains instances of ResultBlock representing a result of the
     * score calculation for the research.
     */
    private static class ResultBlock
    {
        int DocID;          // DocID
        double score;       // the score of the scoring function for the doc identified by DocID

        // constructor with parameters
        public ResultBlock(int DocID, double score) {
            this.DocID = DocID;
            this.score = score;
        }

        // get methods
        public int getDID() {
            return DocID;
        }

        public double getScore() {
            return score;
        }

        @Override
        public String toString() {
            return "PB{" +
                    "DocID = '" + DocID + '\'' +
                    ", score =" + score +
                    '}';
        }
    }

    /**
     * Class to compare the block, allows the order of the priority queue
     * The order is ascending order according to doc score in case of equal score will be made order according to DocID.
     * DocIDs sorted in descending order. In this way in case of a tie, documents with the smallest DocID will be
     * considered better than those with the largest DocID.
     */
    private static class CompareTerm implements Comparator<QueryProcessor.ResultBlock>
    {
        @Override
        public int compare(QueryProcessor.ResultBlock rpb1, QueryProcessor.ResultBlock rpb2) {
            // comparing terms
            int scoreComparison = Double.compare(rpb1.getScore(), rpb2.getScore());
            // if the score are equal, compare by DocID
            if (scoreComparison == 0) {
                // return order by both score and DocID (the score of the two blocks is equal)
                return -Integer.compare(rpb1.getDID(), rpb2.getDID());  // '-' to obtain descending order, without it have ascending order
            }

            return scoreComparison;     // return order only by score (the score of the two blocks is different)
        }
    }
    // -------- end: utilities for priority queue for the results --------

    // -------- start: function to read collection of query --------
    /**
     * This function allows an automatic test of the resolution of preset queries saved in a file on disk.
     * This function will fetch and execute (from the specified file) a number of queries passed as a parameter.
     * These queries will be executed in both conjunctive and disjunctive modes, and a series of statistics will be
     * collected about them (total duration of the test, fastest and slowest conjunctive query, fastest and slowest
     * disjunctive query), which will be displayed at the end of the test.
     *
     * @param numQueries    the number of queries to be read from the file and executed
     * @param pathTest      string identified the path for the test queries
     * @param numTest       the number of test to do
     */
    public static void readQueryFromCollection(int numQueries, String pathTest, int numTest)
    {
        ArrayList<Integer> rankedResults;   // ArrayList that contain the ranked results of query
        int queryCount = 0;             // indicates how many queries have been made
        long startTime, endTime;        // variables to calculate the execution time
        long startTimeTest, endTimeTest;// variables to calculate the execution time of all test
        long fasterQueryCon = 100000;   // indicates the execution time for the fastest query in the collection (conjunctive)
        String quidFastCon = "";        // indicates the query ID of the fastest query in the collection (conjunctive)
        long slowerQueryCon = 0;        // indicates the execution time for the slowest query in the collection (conjunctive)
        String quidSlowCon = "";        // indicates the query ID of the slowest query in the collection (conjunctive)
        long fasterQueryDis = 100000;   // indicates the execution time for the fastest query in the collection (disjunctive)
        String quidFastDis = "";        // indicates the query ID of the fastest query in the collection (disjunctive)
        long slowerQueryDis = 0;        // indicates the execution time for the slowest query in the collection (disjunctive)
        String quidSlowDis = "";        // indicates the query ID of the slowest query in the collection (disjunctive)
        long avgExTimeCon = 0;          // indicate the average execution time for the queries in the collection (conjunctive)
        long avgExTimeDis = 0;          // indicate the average execution time for the queries in the collection (disjunctive)

        // control check for queries
        try {
            if (!queryStartControl())
                return;                 // there aren't all files needed for execute a query, function it's terminated
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        // read term upper bound if is needed
        if (Flags.isDynamicPruningEnabled() && TermDocUpperBound.termUpperBoundTableIsEmpty())
        {
            if (TermDocUpperBound.termUpperBoundFileExist())     // the file already exist
                TermDocUpperBound.readTermUpperBoundTableFromDisk();
            else                                                // the file not exist
                TermDocUpperBound.calculateTermsUpperBound(false);   // calculate term upper bound for each term of dictionary
        }

        // check if the collection file exists
        File file = new File(pathTest);
        if(!file.exists())
        {
            printError("The selected queries collection file does not exist.");
            return;
        }
        printUIMag(" Start query test... from: " + pathTest);         // control print
        String record;          // string to contain the queries and their result

        printUIMag("--------------------------------------------------------------------------------");
        startTimeTest = System.currentTimeMillis();         // start time of all test
        for (int i = 0; i < numTest; i++)
        {   // -- START - for - number of test -
            try
            {
                InputStream tarArchiveInputStream = new GzipCompressorInputStream(new FileInputStream(file));
                BufferedReader buffer_collection = new BufferedReader(new InputStreamReader(tarArchiveInputStream, StandardCharsets.UTF_8));

                printUIMag("-- Start query test number : " + i + " ---------------------------------------------");
                // scan all queries in the collection
                while (((record = buffer_collection.readLine()) != null) && (queryCount < numQueries))
                {   // -- START - while - i-th test -
                    if (record.isBlank())
                        continue;       // empty string or composed by whitespace characters or malformed

                    String[] queryProc = record.split("\t", 2);  // preprocess the query to obtain the result DocNO
                    String qid = queryProc[0];      // get the DocNO of the best result for the query

                    // print of the query and result obtained by search engine
                    printUIMag("---- Query number: " + queryCount + " -------------------------------------------- QueryID: " + qid + " ----");
                    printUIMag("The query is: " + queryProc[1]);
                    printUIMag("---- disjunctive mode ----");

                    startTime = System.currentTimeMillis();         // start time of execute query
                    rankedResults = queryManager(queryProc[1], false, 5);    // run the query in disjunctive mode
                    printQueryResults(rankedResults);
                    endTime = System.currentTimeMillis();           // end time of execute query
                    // shows query execution time
                    printTime("\nQuery (disjunctive mode) executes in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

                    // does queries collection statistics
                    if ((endTime - startTime) < fasterQueryDis) {
                        fasterQueryDis = (endTime - startTime);     // update faster time
                        quidFastDis = qid;                         // update quid
                    }
                    if ((endTime - startTime) > slowerQueryDis) {
                        slowerQueryDis = (endTime - startTime);     // update slower time
                        quidSlowDis = qid;                         //update quid
                    }
                    avgExTimeDis += (endTime - startTime);          // update avg execution time

                    printUIMag("---- conjunctive mode ----");
                    startTime = System.currentTimeMillis();         // start time of execute query
                    rankedResults = queryManager(queryProc[1], true, 5);    // run the query in conjunctive mode
                    printQueryResults(rankedResults);
                    endTime = System.currentTimeMillis();           // end time of execute query
                    // shows query execution time
                    printTime("\nQuery (conjunctive mode) executes in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
                    printUIMag("--------------------------------------------------------------------------------");

                    // does queries collection statistics
                    if ((endTime - startTime) < fasterQueryCon) {
                        fasterQueryCon = (endTime - startTime);     // update faster time
                        quidFastCon = qid;                         // update quid
                    }
                    if ((endTime - startTime) > slowerQueryCon) {
                        slowerQueryCon = (endTime - startTime);     // update slower time
                        quidSlowCon = qid;                         // update quid
                    }
                    avgExTimeCon += (endTime - startTime);          // update avg execution time

                    queryCount++;       // update counter
                }   // -- END - while - i-th test -

                queryCount = 0;                 // reset counter of the query
                buffer_collection.close();      // close buffer reader
                tarArchiveInputStream.close();  // close tar archive
            } catch (IOException e) {
                e.printStackTrace();
            }
        }   // -- END - for - number of test -
        endTimeTest = System.currentTimeMillis();         // start time of all test
        // print queries collection statistics
        printUIMag(" End query test... Executed: " + numTest + " times, from: " + pathTest);         // control print
        printTime(" All test executed in: "  + (endTimeTest - startTimeTest) + " ms (" + formatTime(startTimeTest, endTimeTest) + ")");
        printUIMag("--------------------------------------------------------------------------------");
        printTime("The fastest query (conjunctive mode) executes in " + fasterQueryCon + " ms and its QUID is " + quidFastCon);
        printTime("The slowest query (conjunctive mode) executes in " + slowerQueryCon + " ms and its QUID is " + quidSlowCon);
        printTime("The average queries execution time (conjunctive mode) is " + avgExTimeCon / ((long) numQueries * numTest) + " ms");

        printTime("\nThe fastest query (disjunctive mode) executes in " + fasterQueryDis + " ms and its QUID is " + quidFastDis);
        printTime("The slowest query (disjunctive mode) executes in " + slowerQueryDis + " ms and its QUID is " + quidSlowDis);
        printTime("The average queries execution time (disjunctive mode) is " + avgExTimeDis / ((long) numQueries * numTest) + " ms");
        printUIMag("--------------------------------------------------------------------------------");
    }

    /**
     * Function that automatically analyses the test query terms of a collection passed as a parameter.
     *
     * @param pathTest      string identified the path for the test queries
     */
    public static void readAndAnalyzeQueryFromCollection(String pathTest)
    {
        int queryCount = 0;             // indicates how many queries have been made
        int minTermBefore = 100000;     // min query len before preprocessing
        String quidMinBefore = "";      // indicates the query ID of the shortest query before preprocessing in the collection
        int minTermAfter = 100000;      // min query len after preprocessing
        String quidMinAfter = "";       // indicates the query ID of the shortest query after preprocessing in the collection
        float avgTermBefore = 0;        // avg query len before preprocessing
        float avgTermAfter = 0;         // avg query len after preprocessing
        int maxTermBefore = 0;          // max query len before preprocessing
        String quidMaxBefore = "";      // indicates the query ID of the longest query before preprocessing in the collection
        int maxTermAfter = 0;           // max query len after preprocessing
        String quidMaxAfter = "";       // indicates the query ID of the longest query after preprocessing in the collection
        long startTime, endTime;        // variables to calculate the execution time

        // control check for the collection file
        File collFile = new File(pathTest);               // blocks.txt
        if(!collFile.exists())
        {
            printError("The selected queries collection file does not exist.");
            return;
        }

        // -- control for file into disk
        if (!FileSystem.areThereAllMergedFiles() || !Flags.isThereFlagsFile() || !CollectionStatistics.isThereStatsFile())
        {
            printError("Error: missing required files.");
            return;
        }
        readFlagsFromDisk();                                // read flags from disk

        printUIMag("Term query analysis from: " + pathTest);         // control print
        String record;          // string to contain the queries and their result

        startTime = System.currentTimeMillis();         // start time of all test
        try
        {
            InputStream tarArchiveInputStream = new GzipCompressorInputStream(new FileInputStream(collFile));
            BufferedReader buffer_collection = new BufferedReader(new InputStreamReader(tarArchiveInputStream, StandardCharsets.UTF_8));

            // scan all queries in the collection
            while ((record = buffer_collection.readLine()) != null)
            {   // -- START - while - i-th test -
                if (record.isBlank())
                    continue;       // empty string or composed by whitespace characters or malformed

                String[] queryProc = record.split("\t", 2);  // preprocess the query to obtain the queryDID
                String qid = queryProc[0];      // get the DocNO of the best result for the query
                ArrayList<String> tokenList;    // for contain the preprocessed query
                int termNumBefore = 0;          // the number of the term before preproc for the current query
                int termNumAfter = 0;           // the number of the term after preproc for the current query

                // -- part only for clean the text
                tokenList = TextProcessor.removeAndCleanText(queryProc[1]); // clean
                termNumBefore = tokenList.size();           // take the size

                // queries collection statistics
                if (termNumBefore < minTermBefore)
                {
                    minTermBefore = termNumBefore;          // update faster time
                    quidMinBefore = qid;                    // update quid
                }
                if (termNumBefore > maxTermBefore)
                {
                    maxTermBefore = termNumBefore;          // update faster time
                    quidMaxBefore = qid;                    // update quid
                }
                avgTermBefore += termNumBefore;

                // -- part of preprocessing
                if (Flags.isSwsEnabled())   // Check if the filtering flag is enabled
                {
                    tokenList = TextProcessor.preprocessText(queryProc[1]); // preprocess the text according flags
                    termNumAfter = tokenList.size();        // take the size

                    if (termNumAfter < minTermAfter)
                    {
                        minTermAfter = termNumAfter;            // update faster time
                        quidMaxAfter = qid;                     // update quid
                    }
                    if (termNumAfter > maxTermAfter)
                    {
                        maxTermAfter = termNumAfter;            // update faster time
                        quidMaxAfter = qid;                     // update quid
                    }
                }
                avgTermAfter += termNumAfter;

                queryCount++;       // update counter
            }   // -- END - while - i-th test -

            buffer_collection.close();      // close buffer reader
            tarArchiveInputStream.close();  // close tar archive
        } catch (IOException e) {
            e.printStackTrace();
        }

        endTime = System.currentTimeMillis();         // start time of all test
        // print queries collection statistics
        printTime("Whole analysis executed in: "  + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
        printUIMag("--------------------------------------------------------------------------------");
        printTime("Analysis of the query before preprocessing (only the clean).");
        printTime("The shortest query has " + minTermBefore + " terms and its QUID is " + quidMinBefore);
        printTime("The longest query has " + maxTermBefore + " terms and its QUID is " + quidMaxBefore);
        printTime("The average number of terms before preprocessing is " + avgTermBefore / queryCount );

        if (Flags.isSwsEnabled())   // Check if the filtering flag is enabled
        {
            printTime("\nAnalysis of the query after preprocessing (with stopwords removal and stemming).");
            printTime("The shortest query has " + minTermAfter + " terms and its QUID is " + quidMinAfter);
            printTime("The longest query has " + maxTermAfter + " terms and its QUID is " + quidMaxAfter);
            printTime("The average number of terms after preprocessing is " + avgTermAfter / queryCount );
        }
        printUIMag("--------------------------------------------------------------------------------");
    }

    /**
     * Function allowing a user to enter a query that the system will execute a specified number of times, showing at
     * the end the times of the several executions and the average times.
     *
     * @param numTest       the number of test to do
     * @param sc    scanner to get the choice of the user inserted via keyboard
     */
    public static void testQuery(int numTest, Scanner sc)
    {
        ArrayList<Integer> rankedResults;   // ArrayList that contain the ranked results of query
        long startTime, endTime;        // variables to calculate the execution time
        long startTimeTest, endTimeTest;// variables to calculate the execution time of all test
        long fasterQueryCon = 100000;   // indicates the fastest execution time for the query (conjunctive)
        long slowerQueryCon = 0;        // indicates the slowest execution time for the query (conjunctive)
        long fasterQueryDis = 100000;   // indicates the fastest execution time for the query (disjunctive)
        long slowerQueryDis = 0;        // indicates the slowest execution time for the query (disjunctive)
        long avgExTimeCon = 0;          // indicate the average execution time for the queries(conjunctive)
        long avgExTimeDis = 0;          // indicate the average execution time for the queries (disjunctive)

        // control check for queries
        try {
            if (!queryStartControl())
                return;                 // there aren't all files needed for execute a query, function it's terminated
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        // read term upper bound if is needed
        if (Flags.isDynamicPruningEnabled() && TermDocUpperBound.termUpperBoundTableIsEmpty())
        {
            if (TermDocUpperBound.termUpperBoundFileExist())     // the file already exist
                TermDocUpperBound.readTermUpperBoundTableFromDisk();
            else                                                // the file not exist
                TermDocUpperBound.calculateTermsUpperBound(false);   // calculate term upper bound for each term of dictionary
        }

        // take the query from user
        printUI("Insert query: \n");
        String query = sc.nextLine();           // take user's query
        // control check of the query
        if (query == null || query.isEmpty())
        {
            printError("Error: the query is empty.");
            return;
        }

        // start test
        printUIMag("--------------------------------------------------------------------------------");
        startTimeTest = System.currentTimeMillis();         // start time of all test
        for (int i = 0; i < numTest; i++)
        {   // -- START - for - number of test -
            // print of the query and result obtained by search engine
            printUIMag("-- Start query test number : " + i + " ---------------------------------------------");
            printUIMag("---- disjunctive mode ----");

            startTime = System.currentTimeMillis();         // start time of execute query
            rankedResults = queryManager(query, false, 5);    // run the query in disjunctive mode
            printQueryResults(rankedResults);
            endTime = System.currentTimeMillis();           // end time of execute query
            // shows query execution time
            printTime("\nQuery (disjunctive mode) executes in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

            // does queries collection statistics
            if ((endTime - startTime) < fasterQueryDis)
                fasterQueryDis = (endTime - startTime);     // update faster time

            if ((endTime - startTime) > slowerQueryDis)
                slowerQueryDis = (endTime - startTime);     // update slower time

            avgExTimeDis += (endTime - startTime);          // update avg execution time

            printUIMag("---- conjunctive mode ----");
            startTime = System.currentTimeMillis();         // start time of execute query
            rankedResults = queryManager(query, true, 5);    // run the query in conjunctive mode
            printQueryResults(rankedResults);
            endTime = System.currentTimeMillis();           // end time of execute query
            // shows query execution time
            printTime("\nQuery (conjunctive mode) executes in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
            printUIMag("--------------------------------------------------------------------------------");

            // does queries collection statistics
            if ((endTime - startTime) < fasterQueryCon)
                fasterQueryCon = (endTime - startTime);     // update faster time

            if ((endTime - startTime) > slowerQueryCon)
                slowerQueryCon = (endTime - startTime);     // update slower time

            avgExTimeCon += (endTime - startTime);          // update avg execution time
        }   // -- END - for - number of test -

        endTimeTest = System.currentTimeMillis();         // start time of all test
        // print queries collection statistics
        printUIMag(" End query test... Executed: " + numTest + " times.");         // control print
        printTime(" All test executed in: "  + (endTimeTest - startTimeTest) + " ms (" + formatTime(startTimeTest, endTimeTest) + ")");
        printUIMag("--------------------------------------------------------------------------------");
        printTime("The fastest query execution (conjunctive mode) is " + fasterQueryCon + " ms");
        printTime("The slowest query execution (conjunctive mode) is " + slowerQueryCon + " ms");
        printTime("The average query execution time (conjunctive mode) is " + avgExTimeCon / numTest + " ms");

        printTime("\nThe fastest query execution (disjunctive mode) is " + fasterQueryDis + " ms");
        printTime("The slowest query execution (disjunctive mode) is " + slowerQueryDis + " ms");
        printTime("The average query execution time (disjunctive mode) is " + avgExTimeDis / numTest + " ms");
        printUIMag("--------------------------------------------------------------------------------");
    }
    // -------- end: function to read collection of query --------
}

/*
 * NOTE:
 * 0 - list of conditions (of the if):
 *      partialScore = 0 -> Document Upper bound < threshold -> also with the maximum score in the nonessential posting
 *                          lists the document score will be less than the minimum threshold to be among the top docs.
 *      firstEssPostListIndex = 0 -> all posting list have already been scanned
 *      resetScore = true -> conjunctive case, the DID must be in all posting lists
 *
 * 1 -
 */