package it.unipi.dii.aide.mircv;
import it.unipi.dii.aide.mircv.compression.Unary;
import it.unipi.dii.aide.mircv.compression.VariableBytes;
import it.unipi.dii.aide.mircv.data_structures.*;
import it.unipi.dii.aide.mircv.data_structures.Dictionary;
import it.unipi.dii.aide.mircv.utils.FileSystem;
import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.apache.commons.compress.compressors.gzip.GzipCompressorInputStream;

import java.io.*;
import java.nio.CharBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.StandardCharsets;
import java.util.*;

import static it.unipi.dii.aide.mircv.data_structures.CollectionStatistics.readCollectionStatsFromDisk;
import static it.unipi.dii.aide.mircv.data_structures.DataStructureHandler.readCompressedPostingListFromDisk;
import static it.unipi.dii.aide.mircv.data_structures.DataStructureHandler.readPostingListFromDisk;
import static it.unipi.dii.aide.mircv.data_structures.DataStructureHandler.readAndUncompressCompressedAndSkippedPLFromDisk;
import static it.unipi.dii.aide.mircv.data_structures.DataStructureHandler.readCompTFBlockFromDisk;
import static it.unipi.dii.aide.mircv.data_structures.DataStructureHandler.readCompDIDBlockFromDisk;

import static it.unipi.dii.aide.mircv.data_structures.DocumentElement.DOCELEM_SIZE;
import static it.unipi.dii.aide.mircv.data_structures.Flags.readFlagsFromDisk;
//import static it.unipi.dii.aide.mircv.data_structures.PartialIndexBuilder.dictionaryBlockOffsets;
import static it.unipi.dii.aide.mircv.utils.Constants.*;
import static it.unipi.dii.aide.mircv.utils.Constants.printDebug;
import static it.unipi.dii.aide.mircv.utils.Logger.collStats_logger;
import static java.lang.Math.min;

/**
 * Class to manage and execute query
 */
public final class QueryProcessor
{
    // indicate whether order all or only first "numberOfResults" results from hash table. TEST VARIABLE
    private static boolean orderAllHashMap = false;
    public static HashMap<Integer, DocumentElement> documentTable = new HashMap<>();    // hash table DocID to related DocElement
    public static HashMap<Integer, Double> termFreqWeightTable = new HashMap<>();    // hash table DocID to related DocElement
    static it.unipi.dii.aide.mircv.data_structures.Dictionary dictionary = new Dictionary();    // dictionary in memory
    private static ArrayList<String> termNotInCollection = new ArrayList<>();   // ArrayList that contain the term that are in the query but not in the collection
    static PriorityQueue<QueryProcessor.ResultBlock> resPQ;     // priority queue for the result of scoring function for the best numberOfResults docs

    // var for when both the compression and the skipping are enabled
    // priority queue for the ordered block of DID and Term freq to execute DAAT when the compression and skipping are enabled
    static PriorityQueue<QueryProcessor.orderedDIDBlock> ordDIDPQ = new PriorityQueue<>(new CompareOrdDIDBlock());
    static ArrayList<Posting>[] skipAndCompPLs;  // contains all the posting lists for each term of the query (case of compression and skipping enabled)
    static PriorityQueue<Integer> pqDID = new PriorityQueue<Integer>();    // priority queue for the DID when use max score when compression and skipping are enabled
    // ---------------------------------------------- start: functions -------------------------------------------------

    // ---------------------------------------- start: set and get functions -------------------------------------------

    /**
     *
     * @param pl
     * @param index
     */
    public static void setPLInSkipAndCompPLs(ArrayList<Posting> pl, int index)
    {
        if ((skipAndCompPLs != null) || (index >= skipAndCompPLs.length))
            return;

        skipAndCompPLs[index] = pl;  // set new uncompressed block of the posting list
    }

    /**
     *
     * @param index
     * @return
     */
    public static ArrayList<Posting> getPLInSkipAndCompPLs(int index)
    {
        if ((skipAndCompPLs != null) || (index >= skipAndCompPLs.length))
            return null;

        return skipAndCompPLs[index];
    }
    // ----------------------------------------- end: set and get functions --------------------------------------------

    /**
     * fuction to manage the query request. Prepare and execute the query and return the results.
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
            //printDebug("Query before processed: " + query);
            processedQuery = TextProcessor.preprocessText(query); // Preprocessing of document text
            //printDebug("Query after processed: " + processedQuery);

            // check if query is empty
            if (processedQuery.isEmpty() || (processedQuery.size() == 1 && processedQuery.get(0).equals("")))
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
     * function that checks whether all files and resources required to execute the query are available
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

        readFlagsFromDisk();
        readCollectionStatsFromDisk();
        //CollectionStatistics.printCollectionStatistics();

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

        if (termFreqWeightTable.isEmpty())
            readTFWeightFromDisk();
        else
            printDebug("termFreqWeightTable is already loaded.");

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

        if (dynamicPrun)        // dynamic pruning 'true' -> use Max Score
        {
            if (skipping && compression)
                DAATAlgMAXSCORESkipAndComp(processedQuery,scoringFunc,isConjunctive,numberOfResults);   // apply DAAT + MaxScore (comp+skipping)  to calculate the score of the Docs
            else if (skipping)
                DAATAlgMAXSCORESkipping(processedQuery,scoringFunc,isConjunctive,numberOfResults);  // apply DAAT + MaxScore (skipping) to calculate the score of the Docs
            else
                DAATAlgMAXSCORE(processedQuery,scoringFunc,isConjunctive,numberOfResults);  // apply DAAT + MaxScore (no skipping) to calculate the score of the Docs
        }
        else                    // dynamic pruning 'false' -> use DAAT alg
        {
            if (compression && skipping)
                DAATAlgCompSkip(processedQuery,scoringFunc,isConjunctive,numberOfResults);  // apply DAAT with compression and skipping
            else
                DAATAlgorithm(processedQuery,scoringFunc,isConjunctive,numberOfResults);        // apply DAAT to calculate the score of the Docs
        }
    }

    // ------------------------------------------- START - Execution Query alg -----------------------------------------
    /**
     * function for apply the Document at a Time algorithm
     *
     * @param scoringFunc       indicates the preference for scoring. if false use TFIDF, if true use BM25.
     * @param processedQuery    array list for containing the query term
     * @param isConjunctive     indicates whether the query is conjunctive or disjunctive type (default is disjunctive)
     * @param numberOfResults   indicated the max number of result to return to user
     */
    private static void DAATAlgorithm(ArrayList<String> processedQuery, boolean scoringFunc , boolean isConjunctive, int numberOfResults) throws FileNotFoundException
    {
        resPQ = new PriorityQueue<>(numberOfResults, new CompareTerm());    // length equal to the number of results to be returned to the user
        String[] terms = new String[processedQuery.size()];                 // array containing the terms of the query
        ArrayList<Integer> ordListDID;      // ordered list of the DocID present in the all posting lists of the term present in the query
        ArrayList<Posting>[] postingLists;  // contains all the posting lists for each term of the query
        double[] IDFweight;                 // array containing the IDF weight for each posting list
        int[] lengthPostingList;            // array containing the length of the posting lists
        int[] postingListsIndex;            // contain the current position index for the posting list of each term in the query
        Posting currentP;                   // support var
        int currentDID = 0;                 // DID of the current doc processed in algorithm
        double partialScore = 0;            // var that contain partial score
        int docScoreCalc = 0;               // indicates the number of documents whose score was calculated (0 to number of results requested by the user)
        boolean resetScore = false;         // used only in conjunctive case. indicates that the score must be set to 0 (the current Doc there aren't all the term of the query)
        int pLNotEmpty = 0;                 // contains the number of posting lists related to the query terms that aren't empty
        long startTime,endTime;             // variables to calculate the execution time

        startTime = System.currentTimeMillis();         // start time for retrieve all posting lists of the query
        postingLists = retrieveAllPostListsFromQuery(processedQuery);   // take all posting lists of query terms
        endTime = System.currentTimeMillis();           // end time for retrieve all posting lists of the query
        // shows query execution time
        printTime("\n*** Retrieved all posting lists in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

        pLNotEmpty = NumberOfPostingListNotEmpty(postingLists);     // take the number of list that are not empty
        // check the number of posting lists not empty and perform the best choice
        if (pLNotEmpty == 0)    // all terms in the query aren't in the dictionary or empty query
            return;     // exit
        else if (pLNotEmpty == 1)   // there is only 1 postingList (query with one term or query with more term but only one in dictionary)
        {
            if ((postingLists.length != 1) && (isConjunctive))  // case of query conjunctive and more term but only one in dictionary
                return;     // exit

            DAATOnePostingList(processedQuery, postingLists, scoringFunc, numberOfResults);   // execute DAAT algorithm
            return;     // exit
        }
        // more postingLists not empty
        ordListDID = DIDOrderedListOfQuery(postingLists, isConjunctive);        // take ordered list of DocID
        postingListsIndex = getPostingListsIndex(postingLists);     // get the index initialized
        // NEW VERSION
        for (int i = 0; i < processedQuery.size(); i++)
            terms[i] = processedQuery.get(i);
        lengthPostingList = retrieveLengthAllPostingLists(terms);   // take the length of each posting list
        IDFweight = calculateIDFWeight(lengthPostingList);          // calculate the IDF weight

        startTime = System.currentTimeMillis();           // start time of DAAT
        // scan all Doc retrieved and calculate score (TFIDF or BM25)
        for (Integer integer : ordListDID)
        {
            currentDID = integer;       // update the DID, document of which to calculate the score
            partialScore = 0;           // reset var
            resetScore = false;         // set to false

            // default case is query Disjunctive
            // take all values and calculating the scores in the posting related to currentDID
            for (int j = 0; j < postingLists.length; j++) {
                // check if the posting lists of j-th isn't at the end AND if the j-th term of the query is present in the doc identify by currentDID
                if ( (postingLists[j] != null) && (postingListsIndex[j] < postingLists[j].size()) && (postingLists[j].get(postingListsIndex[j]).getDocId() == currentDID)) {
                    currentP = postingLists[j].get(postingListsIndex[j]);              // take posting
                    postingListsIndex[j]++;                         // update index of current value
                    //System.out.println("DAAT, prescoring -- df = " + DataStructureHandler.postingListLengthFromTerm(processedQuery.get(j)));

                    // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                    if (scoringFunc)
                        partialScore += ScoringBM25(currentDID, currentP.getTermFreq(), IDFweight[j]);     // use BM25
                    else
                        partialScore += ScoringTFIDF(currentP.getTermFreq(), IDFweight[j]);               // use TFIDF
                    //printDebug("DAAT: posting del termine: " + processedQuery.get(j) + " della posting: " + j + " in pos: " + (postingListsIndex[j]-1) + " ha DID: " + currentDID + " and partialScore: " + partialScore);
                } else if (isConjunctive) {
                    // must take only the document in which there are all term (DID that compare in all posting lists of the terms)
                    resetScore = true;       // reset the partial score
                    //printDebug("Sono in conjuntive, azzero lo score. Posting list numero: " + j + " in pos: " + postingListsIndex[j]);
                    // if all postings in one posting lists have already been seen the next documents in the posting lists cannot contain all the terms in the query
                    if ((postingLists[j] == null) || (postingListsIndex[j] >= postingLists[j].size())) {
                        //printDebug("Query conjunctive, posting list numero: " + j + " finita. Si è in pos: " + postingListsIndex[j] + " su dimensione: " + postingLists[j].size());
                        endTime = System.currentTimeMillis();           // end time of DAAT
                        // shows query execution time
                        printTime("*** DAAT V.0.6 execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
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
                    docScoreCalc++;         // increment result in priority queue counter
                }
                else if (resPQ.peek().getScore() < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                {
                    // substitution of the block
                    //printDebug("Old block : DID = " + resPQ.peek().getDID()+ " score: " + resPQ.peek().getScore());
                    resPQ.poll();       // remove the first element
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                    //printDebug("New block : DID = " + currentDID+ " score: " + partialScore);
                }
            }
        }
        endTime = System.currentTimeMillis();           // end time of DAAT
        // shows DAAT execution time
        printTime("*** DAAT V.0.6 execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
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
        resPQ = new PriorityQueue<>(numberOfResults, new CompareTerm());    // length equal to the number of results to be returned to the user
        String[] terms = new String[processedQuery.size()];                 // array containing the terms of the query
        SkipList[] skipListArray;           // array of the Skip List reference related to the term of query
        int pLNotEmpty = 0;                 // contains the number of posting lists related to the query terms that aren't empty
        double[] IDFweight;                 // array containing the IDF weight for each posting list
        int[] lengthPostingList;            // array containing the length of the posting lists
        int[] blockIndex;                   // contain the index of the next block of skipList for each PL
        int[] elmOfBlockDone;               // contain the counter of the element belong to the same PL block
        int[] totalElemOfBlok;              // contain the total number of element of the same block in the PQ
        int previousDID = 0;                // DID of the previous doc processed in algorithm
        int currentDID = 0;                 // DID of the current doc processed in algorithm
        int currIndexPL = 0;                // indicates the index of the PL from which the current element taken from the PQ originates
        orderedDIDBlock currElem;           // the current (at each iteration) element poll from ordDIDPQ
        int docScoreCalc = 0;               // indicates the number of documents whose score was calculated (0 to number of results requested by the user)
        double partialScore = 0;            // var that contain partial score
        int countSameDID = 0;
        long startTime,endTime;             // variables to calculate the execution time

        // check the number of the term that have PL (which are in the dictionary)
        pLNotEmpty = numberOfQueryTermsInDictionary(processedQuery);
        // create the skip List reference related to the term of query
        skipListArray = skipListInitCompAndSkip(processedQuery, null, false);

        // check the number of posting lists not empty and perform the best choice
        if (pLNotEmpty == 0)        // all terms in the query aren't in the dictionary or empty query
            return;     // exit
        else if (pLNotEmpty == 1)   // there is only 1 postingList (query with one term or query with more term but only one in dictionary)
        {
            if ((processedQuery.size() != 1) && (isConjunctive))  // case of query conjunctive and more term but only one in dictionary
                return;     // exit

            // The PL is only one -> read and decompress the whole PL and use the classic optimization method
            startTime = System.currentTimeMillis();         // start time for retrieve the posting lists of the query
            ArrayList<Posting>[] postingLists = retrieveAllUncompPL(processedQuery, skipListArray); // get the uncompress PL
            endTime = System.currentTimeMillis();           // end time for retrieve the posting lists of the query
            // shows query execution time
            printTime("*** DAAT (comp+skipping) retrieved PL (case 1 PL) in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

            DAATOnePostingList(processedQuery, postingLists, scoringFunc, numberOfResults);   // execute DAAT algorithm
            return;     // exit
        }

        // -- more postingLists not empty --
        // 0) take the first block of the PL
        startTime = System.currentTimeMillis();         // start time for retrieve first block of PLs of the query
        retrieveFirstCompBlockOfPLFromQuery(processedQuery, null, skipListArray, false); // retrieve the first block of each PL and put value in the PQ
        endTime = System.currentTimeMillis();           // end time for retrieve first block of PLs of the query
        // shows query execution time
        printTime("*** DAAT (comp+skipping) retrieved first block of each PLs in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

        // 1) setting of all is need for the DAAT algorithm
        for (int i = 0; i < processedQuery.size(); i++)             // take the term of the query
            terms[i] = processedQuery.get(i);
        lengthPostingList = retrieveLengthAllPostingLists(terms);   // take the length of each posting list
        IDFweight = calculateIDFWeight(lengthPostingList);          // calculate the IDF weight
        // set the index block for each PL
        blockIndex = new int[processedQuery.size()];
        Arrays.fill(blockIndex, 1);                 // the first block (with index 0) has already been taken
        // set the counter and max number of elem of the same block in the PQ -> used to manage the load of next PL block from disk
        elmOfBlockDone = new int[processedQuery.size()];
        Arrays.fill(elmOfBlockDone, 0);             // at the beginning the elem taken by PQ are equal to 0
        totalElemOfBlok = new int[processedQuery.size()];
        // the max number of a block (a block has len equal SKIP_POINTERS_THRESHOLD if is not the last block and less
        // if is the last block but in this case there is no need to load the next block from memory
        Arrays.fill(totalElemOfBlok, SKIP_POINTERS_THRESHOLD);

        // 2) start DAAT
        startTime = System.currentTimeMillis();           // start time of DAAT (comp + skipping)
        previousDID = ordDIDPQ.peek().getDocID();
        int count = 0;
        while (!ordDIDPQ.isEmpty())
        {
            currElem = ordDIDPQ.poll();             // get the current element
            currentDID = currElem.getDocID();       // get current DID
            currIndexPL = currElem.getIndexPL();    // get the current PL index

            if (currentDID == previousDID)      // sum the current partial score with the previous
            {
                // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                if (scoringFunc)
                    partialScore += ScoringBM25(currentDID, currElem.getTermFreq(), IDFweight[currIndexPL]);     // use BM25
                else
                    partialScore += ScoringTFIDF(currElem.getTermFreq(), IDFweight[currIndexPL]);               // use TFIDF

                countSameDID++;     // update counter
            }
            else            // other DID, the previous partial score is final, save it.
            {
                if (isConjunctive && (countSameDID != processedQuery.size()))   // conjunctive and the current DID doesn't contain all the term in the query
                {
                    countSameDID = 1;   // reset counter
                }
                else
                {
                    countSameDID = 1;   // reset counter
                    // save score
                    if (partialScore != 0)
                    {
                        // insert without control into priority queue (is not full) or insert all results (orderAllHashMap = true)
                        if ((docScoreCalc < numberOfResults) || orderAllHashMap)
                        {
                            resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                            docScoreCalc++;         // increment result in priority queue counter
                        }
                        else if (resPQ.peek().getScore() < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                        {
                            // substitution of the block
                            //printDebug("Old block : DID = " + resPQ.peek().getDID()+ " score: " + resPQ.peek().getScore());
                            resPQ.poll();       // remove the first element
                            resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                            //printDebug("New block : DID = " + currentDID+ " score: " + partialScore);
                        }
                    }
                }

                // calculate SCORE (TFIDF or BM25) for this term and currentDID
                if (scoringFunc)
                    partialScore = ScoringBM25(currentDID, currElem.getTermFreq(), IDFweight[currIndexPL]);     // use BM25
                else
                    partialScore = ScoringTFIDF(currElem.getTermFreq(), IDFweight[currIndexPL]);               // use TFIDF

                previousDID = currentDID;       // update previous DID
            }

            //printDebug("PRIMA if -- counter elem: " + elmOfBlockDone[currIndexPL] + " of index: " + currIndexPL);
            // update the counter -> when a counter is equal to maximum value of the block load the next block of the related PL
            if (totalElemOfBlok[currIndexPL] == ++elmOfBlockDone[currIndexPL])  // all elem of a block have been taken
            {
                // load the new block and put elem into PQ
                assert skipListArray != null;
                //printDebug("All the elems for '" + processedQuery.get(currIndexPL) + "' are taken. Len of PQ: " + ordDIDPQ.size());
                retrieveCompBlockOfPL(processedQuery.get(currIndexPL), skipListArray, blockIndex[currIndexPL], currIndexPL, false);
                blockIndex[currIndexPL]++;          // increment the counter of block
                //printDebug("Now the len of PQ: " + ordDIDPQ.size() + " now the bloc index is: " + blockIndex[currIndexPL]);
                elmOfBlockDone[currIndexPL] = 0;    // reset the counter of elem
            }
            //printDebug("DOPO if -- counter elem: " + elmOfBlockDone[currIndexPL] + " of index: " + currIndexPL);
            count++;    // update iteration counter
        }
        if (isConjunctive && (countSameDID != processedQuery.size()))   // conjunctive and the current DID doesn't contain all the term in the query
        {
            partialScore = 0;   // reset counter
        }
        // save the last partial score
        if ((partialScore != 0))
        {
            // insert without control into priority queue (is not full) or insert all results (orderAllHashMap = true)
            if ((docScoreCalc < numberOfResults) || orderAllHashMap)
            {
                resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                docScoreCalc++;         // increment result in priority queue counter
            }
            else if (resPQ.peek().getScore() < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
            {
                // substitution of the block
                //printDebug("Old block : DID = " + resPQ.peek().getDID()+ " score: " + resPQ.peek().getScore());
                resPQ.poll();       // remove the first element
                resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                //printDebug("New block : DID = " + currentDID+ " score: " + partialScore);
            }
        }
        endTime = System.currentTimeMillis();           // end time of DAAT (comp + skipping)
        // shows DAAT (comp + skipping) execution time
        printTime("*** DAAT (comp+skipping) execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
        //printDebug("The total iteration od DAAT alg. is: " + count);
    }

    /**
     * function for apply the Document at a Time algorithm in the special case of only one posting list not empty
     *
     * @param postingLists      array list for containing the posting lists of the query's term
     * @param scoringFunc       indicates the preference for scoring. If false use TFIDF, if true use BM25.
     * @param processedQuery    array list for containing the query term
     * @param numberOfResults   indicated the max number of result to return to user
     */
    private static void DAATOnePostingList (ArrayList<String> processedQuery, ArrayList<Posting>[] postingLists, boolean scoringFunc, int numberOfResults)
    {
        String[] terms = new String[processedQuery.size()]; // array containing the terms of the query
        double partialScore = 0;    // var that contain partial score
        long startTime,endTime;     // variables to calculate the execution time
        int docScoreCalc = 0;       // indicates the number of documents whose score was calculated (0 to number of results requested by the user)
        double[] IDFweight;         // array containing the IDF weight for each posting list
        int[] lengthPostingList;    // array containing the

        // initialize
        for (int i = 0; i < processedQuery.size(); i++)
            terms[i] = processedQuery.get(i);
        lengthPostingList = retrieveLengthAllPostingLists(terms);   // take the length of each posting list
        IDFweight = calculateIDFWeight(lengthPostingList);          // calculate the IDF weight
        //printDebug("IDFweight length: " + IDFweight.length + " value: " + IDFweight[0]);

        startTime = System.currentTimeMillis();         // start time of DAAT
        int index = 0;
        for (int i = 0; i < postingLists.length; i++)   // take the index of the not fully scanned posting lists
        {
            // (term that there isn't in collection -> posting list == null) OR (posting list completely visited)
            if (postingLists[i] == null)
                continue;           // go to next posting list
            index = i;          // take the index of the not fully scanned posting list
        }

        // optimization -> there is one posting list -> the DID are already sort
        for (Posting p : postingLists[index])
        {
            partialScore = 0;           // reset var
            // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
            if (scoringFunc)
                partialScore += ScoringBM25(p.getDocId(),p.getTermFreq(), IDFweight[index]  );     // use BM25
            else
                partialScore += ScoringTFIDF(p.getTermFreq(), IDFweight[index]);               // use TFIDF

            // save score
            if (partialScore != 0)
            {
                // insert without control into priority queue (is not full) or insert all results (orderAllHashMap = true)
                if ((docScoreCalc < numberOfResults) || orderAllHashMap)
                {
                    resPQ.add(new QueryProcessor.ResultBlock(p.getDocId(), partialScore));     // add to priority queue
                    docScoreCalc++;         // increment result in priority queue counter
                }
                else if (resPQ.peek().getScore() < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                {
                    // substitution of the block
                    //printDebug("Old block : DID = " + resPQ.peek().getDID()+ " score: " + resPQ.peek().getScore());
                    resPQ.poll();       // remove the first element
                    resPQ.add(new QueryProcessor.ResultBlock(p.getDocId(), partialScore));     // add to priority queue
                    //printDebug("New block : DID = " + currentDID+ " score: " + partialScore);
                }
            }
        }
        endTime = System.currentTimeMillis();           // end time of DAAT
        // shows DAAT execution time
        printTime("*** DAAT V.0.5 (only 1 postingList) execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
    }

    /**
     * Function for apply the Document at a Time algorithm with WAND algorithm as dynamic pruning algorithm.
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
        String[] terms = new String[processedQuery.size()];                 // array containing the terms of the query
        ArrayList<Posting>[] postingLists;  // contains all the posting lists for each term of the query
        ArrayList<Integer> ordListDID;      // ordered list of the DocID present in the all posting lists of the term present in the query
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
        long startTime,endTime;             // variables to calculate the execution time

        startTime = System.currentTimeMillis();         // start time for retrieve all posting lists of the query
        postingLists = retrieveAllPostListsFromQuery(processedQuery);   // take all posting lists of query terms
        endTime = System.currentTimeMillis();           // end time for retrieve all posting lists of the query
        // shows query execution time
        printTime("\n*** Retrieved all posting lists in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

        pLNotEmpty = NumberOfPostingListNotEmpty(postingLists);     // take the number of list that are not empty
        // check the number of posting lists not empty and perform the best choice
        if (pLNotEmpty == 0)    // all terms in the query aren't in the dictionary or empty query
            return;     // exit
        else if (pLNotEmpty == 1)   // there is only 1 postingList (query with one term or query with more term but only one in dictionary)
        {
            if ((postingLists.length != 1) && (isConjunctive))  // case of query conjunctive and more term but only one in dictionary
                return;     // exit

            DAATOnePostingList(processedQuery, postingLists, scoringFunc, numberOfResults);   // execute DAAT algorithm
            return;     // exit
        }

        // more postingLists not empty, use WAND algorithm
        ordListDID = DIDOrderedListOfQuery(postingLists, isConjunctive);    // take ordered list of DocID
        postingListsIndex = getPostingListsIndex(postingLists);             // get the index initialized
        postingListsIndexWAND = getPostingListsIndex(postingLists);         // get the WAND index initialized
        termUpperBoundList = getPostingListsTermUpperBound(processedQuery); // get the term upper bound for each term(postinglist)
        // initialize array for improvement TFIDF and BM25 scoring
        for (int i = 0; i < processedQuery.size(); i++)                     // get query terms
            terms[i] = processedQuery.get(i);
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

    // ---- NEW VERSION OF DAAT WITH MAXSCORE -- START ----
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
        String[] orderedQueryTerm = new String[processedQuery.size()];      // contains ...
        double[] termUpperBoundList  = new double[processedQuery.size()];   // contains all the term upper bound for each term of the query
        double[] sumTUBList  = new double[processedQuery.size()];           // array containing the sum of TUB, the value at the i-th position is the sum of TUBs from position 0 to (i-1)
        ArrayList<Posting>[] postingLists;  // contains all the posting lists for each term of the query
        ArrayList<Integer> ordListDID;      // ordered list of the DocID present in the all posting lists of the term present in the query
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
        long startTime,endTime;             // variables to calculate the execution time

        startTime = System.currentTimeMillis();         // start time for retrieve all posting lists of the query
        postingLists = retrieveAllPostingListsMaxScore(processedQuery,orderedQueryTerm,termUpperBoundList,sumTUBList);   // take all posting lists of query terms
        endTime = System.currentTimeMillis();           // end time for retrieve all posting lists of the query
        // shows query execution time
        printTime("\n*** Retrieved all posting lists in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

        pLNotEmpty = NumberOfPostingListNotEmpty(postingLists);     // take the number of list that are not empty
        // check the number of posting lists not empty and perform the best choice
        if (pLNotEmpty == 0)    // all terms in the query aren't in the dictionary or empty query
            return;     // exit
        else if (pLNotEmpty == 1)   // there is only 1 postingList (query with one term or query with more term but only one in dictionary)
        {
            if ((postingLists.length != 1) && (isConjunctive))  // case of query conjunctive and more term but only one in dictionary
                return;     // exit

            DAATOnePostingList(processedQuery, postingLists, scoringFunc, numberOfResults);   // execute DAAT algorithm
            return;     // exit
        }

        // more postingLists not empty, use MAXSCORE algorithm
        ordListDID = DIDOrderedListOfQuery(postingLists, isConjunctive);        // take ordered list of DocID
        postingListsIndex = getPostingListsIndex(postingLists);                 // get the index initialized
        lengthPostingList = retrieveLengthAllPostingLists(orderedQueryTerm);    // take the length of each posting list
        IDFweight = calculateIDFWeight(lengthPostingList);                      // calculate the IDF weight
        // control print
        //printDebug("orderedQueryTerm -> " + Arrays.toString(orderedQueryTerm));
        //printDebug("termUpperBoundList -> " + Arrays.toString(termUpperBoundList));
        //printDebug("sumTUBList -> " + Arrays.toString(sumTUBList));
        //printDebug("lengthPostingList -> " + Arrays.toString(lengthPostingList));
        /*
        for (int i = 0; i < skipListArray.length; i++)
        {
            printDebug("Control print skipping list of term: " + orderedQueryTerm[i]);
            skipListArray[i].testReadAllSkip();
        }//*/

        //printDebug("START DAAT + MAXSCORE\n");
        startTime = System.currentTimeMillis();           // start time of DAAT + MAX SCORE
        // MaxScore algorithm - scan all Doc retrieved and calculate score (TFIDF or BM25)
        for (Integer currentDID : ordListDID)
        {   // -- start - for 0: DID --
            //printDebug("START cycle with DID: " + currentDID);
            partialScore = 0;           // reset var
            resetScore = false;         // set to false

            //printDebug("-- START EP - with first EP: " + firstEssPostListIndex + " of DID: " + currentDID);
            // scan the essential posting lists, default case is query Disjunctive
            for (int j = firstEssPostListIndex; j < postingLists.length; j++)
            {   // -- start - for 0.1: EPL --
                //printDebug("-- EP - search DID: " + currentDID + " in term: '" + orderedQueryTerm[j] + "' of posting: " + j + " in pos: " + (postingListsIndex[j]));
                // check if the posting lists of j-th isn't at the end AND if the j-th term of the query is present in the doc identify by currentDID
                if ( (postingLists[j] != null) && (postingListsIndex[j] < postingLists[j].size()) && (postingLists[j].get(postingListsIndex[j]).getDocId() == currentDID))
                {
                    currentP = postingLists[j].get(postingListsIndex[j]);              // take posting
                    postingListsIndex[j]++;                         // update index of current value
                    // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                    if (scoringFunc)
                        partialScore += ScoringBM25(currentDID,currentP.getTermFreq(), IDFweight[j]);     // use BM25
                    else
                        partialScore += ScoringTFIDF(currentP.getTermFreq(), IDFweight[j]);               // use TFIDF
                    //printDebug("---- EP - find DID: " + currentDID + " and partialScore: " + partialScore + " in term: '" + orderedQueryTerm[j] + "' of posting: " + j + " in pos: " + (postingListsIndex[j]-1));
                }
                else if (isConjunctive)     // must take only the document in which there are all term (DID that compare in all posting lists of the terms)
                {
                    resetScore = true;       // reset the partial score
                    // if all postings in one posting lists have already been seen the next documents in the posting lists cannot contain all the terms in the query
                    if ((postingLists[j] == null) || (postingListsIndex[j] >= postingLists[j].size())) {
                        //printDebug("Query conjunctive, posting list numero: " + j + " finita. Si è in pos: " + postingListsIndex[j] + " su dimensione: " + postingLists[j].size());
                        endTime = System.currentTimeMillis();           // end time of DAAT
                        // shows query execution time
                        printTime("*** DAAT + MAX SCORE V.1 execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
                        return;             // exit from function
                    }
                }
            }   // -- end - for 0.1: EPL --
            //printDebug("-- END EP -> DID : " + currentDID + " partialScore: " + partialScore + " and threshold: " + threshold + " reset score: " + resetScore);

            // Conditions under which analysis of nonessential posting lists can be skipped -- SEE NOTE 0 --
            if ( (partialScore == 0) || resetScore )
                continue;       // go to next iteration, the current doc can't be among the top result

            // scan non essential posting lists
            if (firstEssPostListIndex != 0)
            {   // -- start - if: NoEPL --
                currentDocUpperBound = partialScore + sumTUBList[firstEssPostListIndex];    // calculate the current DUB
                //printDebug("-- START NonEP - DocUB: " + currentDocUpperBound + " and threshold: " + threshold);
                // check if the doc has no zero possibility to have a score greater than threshold
                if (currentDocUpperBound <= threshold)
                    continue;                           // go to next iteration with next DID

                // update the score
                for (int i = 0; i < firstEssPostListIndex; i++)
                {
                    if (isConjunctive)
                        resetScore = true;       // reset the partial score

                    if ((postingLists[i] != null) && (postingListsIndex[i] < postingLists[i].size()))
                    {
                        //printDebug("---- NonEP HOP -> search the DID " + currentDID + " for the term: '" + orderedQueryTerm[i] + "' in the posting list: " + i + " in pos: " + postingListsIndex[i] + " and max len: " + postingLists[i].size() + " and did is: " + postingLists[i].get(postingListsIndex[i]).getDocId());
                        /*
                        // -- START -- NEW PART --
                        if (postingLists[i].get(postingListsIndex[i]).getDocId() < currentDID)
                        {
                            int pos = skipListArray[i].nextGEQ(currentDID);
                            printDebug("Uses skipping -> postingListSize: " + postingLists[i].size() + " search DID: " + currentDID +
                                    " and nextGEQ return position " + pos + " that have DID: " + (pos < postingLists[i].size() ? postingLists[i].get(pos).getDocId() : " out of bound"));
                        }
                        //*/
                        // -- END -- NEW PART --

                        // check first position
                        if (postingLists[i].get(postingListsIndex[i]).getDocId() > currentDID)
                        {
                            currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                            //printDebug("------ NonEP - First position > currentDID. update DUB " + currentDocUpperBound + " partialScore: " + partialScore + " in PL: " + i + " in pos: " + postingListsIndex[i]);
                            continue;       // go to the next step
                        }
                        //printDebug("**** pos+range: " + (postingListsIndex[i] + range) + " ha DID: " + postingLists[i].get(postingListsIndex[i]+range).getDocId());
                        while (((postingListsIndex[i] + range) < postingLists[i].size()) && (postingLists[i].get(postingListsIndex[i]+range).getDocId() <= currentDID))
                        {
                            //printDebug("-------- NonEP - search(RANGE) DID: " + currentDID + " in PL: " + i + " in pos: " + postingListsIndex[i] + " with did value: " + postingLists[i].get(postingListsIndex[i]).getDocId());
                            postingListsIndex[i] += range;         // update index of current value
                        }
                        //printDebug("------ NonEP - Use booleanSearch for the index of currentDID.");
                        // check the size of the posting list
                        if ((postingListsIndex[i] + range) < postingLists[i].size())
                            tempList = postingLists[i].subList(postingListsIndex[i],postingListsIndex[i]+range);
                        else
                            tempList = postingLists[i].subList(postingListsIndex[i],postingLists[i].size()-1);

                        newIndex = booleanSearch(tempList, currentDID);        // get the new index
                        if (newIndex == -1)     // the searched DID there isn't in the posting list
                        {
                            postingListsIndex[i]++;                             // update the index with the new value
                            currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                            //printDebug("------ NonEP - CurrentDID there isn't in this posting list. update DUB " + currentDocUpperBound + " partialScore: " + partialScore + " in PL: " + i + " in pos: " + postingListsIndex[i]);
                        }
                        else                    // the searched DID there is in the posting list
                        {
                            postingListsIndex[i] += newIndex;       // update the index with the new value
                            // find the searched DID - update the partialScore and the currentDocUpperBound
                            currentP = postingLists[i].get(postingListsIndex[i]);              // take posting
                            //printDebug("------ NonEP -> FIND the DID " + currentDID + " old DUB: " + currentDocUpperBound + " and old partialScore: " + partialScore);
                            currentDocUpperBound -= partialScore;               // update currentDocUpperBound
                            // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                            if (scoringFunc)
                                partialScore += ScoringBM25(currentDID,currentP.getTermFreq(), IDFweight[i]);   // use BM25
                            else
                                partialScore += ScoringTFIDF(currentP.getTermFreq(), IDFweight[i]);             // use TFIDF
                            currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                            currentDocUpperBound += partialScore;               // update currentDocUpperBound

                            resetScore = false;         // reset the partial score
                            //printDebug("------ NonEP -> FIND the DID " + currentDID + " in posting: " + i + " in pos: " + postingListsIndex[i] + " with valueDID: " + postingLists[i].get(postingListsIndex[i]).getDocId() + " update DUB: " + currentDocUpperBound + " and partialScore: " + partialScore);
                            postingListsIndex[i]++;     // update the index of the postin list
                        }
                    } // if all postings in one posting lists have already been seen the next documents in the posting lists cannot contain all the terms in the query
                    else if (isConjunctive)
                    {
                        //printDebug("Query conjunctive, posting list numero: " + j + " finita. Si è in pos: " + postingListsIndex[j] + " su dimensione: " + postingLists[j].size());
                        endTime = System.currentTimeMillis();           // end time of DAAT
                        // shows query execution time
                        printTime("*** DAAT + MAX SCORE V.1 execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
                        return;             // exit from function
                    }
                    //printDebug("-- END NonEP - DocUB: " + currentDocUpperBound + " with score: " + partialScore + " and threshold: " + threshold);
                    // check if the doc has no zero possibility to have a score greater than threshold
                    if ((currentDocUpperBound <= threshold) || resetScore)
                        break;
                }
            }   // -- end - if: NoEPL --
            // save score
            if (resetScore)
                continue;       // go to the next iteration (next Doc)
            // insert without control into priority queue (is not full) or insert all results (orderAllHashMap = true)
            if (docScoreCalc < numberOfResults)
            {
                resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                docScoreCalc++;                         // increment result in priority queue counter
                if (docScoreCalc == numberOfResults)
                    threshold = resPQ.peek().getScore();    // update threshold
                //printDebug("-- SCORING - Add result: " + partialScore);
            }
            else if (threshold < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
            {
                // substitution of the block
                resPQ.poll();                           // remove the first element
                resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                threshold = resPQ.peek().getScore();    // update threshold
                // calculate new essential posting lists and update firstEssPostListIndex
                firstEssPostListIndex = updateEssentialPositngLists(sumTUBList, threshold);
                //printDebug("-- **** New threshold: " + threshold + " new first essential posting list: " + firstEssPostListIndex);
            }
        }   // -- end - for: DID --

        endTime = System.currentTimeMillis();           // end time of DAAT
        // shows DAAT execution time
        printTime("*** DAAT + MAX SCORE V.1 execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
    }

    // ---- NEW VERSION OF DAAT WITH MAXSCORE -- START ----
    /**
     * Function for apply the Document at a Time algorithm with max Score algorithm as dynamic pruning algorithm.
     * This function apply the Max Score in the case of skipping enabled and compression disabled.
     *
     * @param processedQuery    array list for containing the query term
     * @param scoringFunc       indicates the preference for scoring. if false use TFIDF, if true use BM25.
     * @param isConjunctive     indicates whether the query is conjunctive or disjunctive type (default is disjunctive)
     * @param numberOfResults   indicated the max number of result to return to user
     * @throws FileNotFoundException
     */
    private static void DAATAlgMAXSCORESkipping(ArrayList<String> processedQuery, boolean scoringFunc , boolean isConjunctive, int numberOfResults) throws FileNotFoundException
    {
        resPQ = new PriorityQueue<>(numberOfResults, new CompareTerm());    // length equal to the number of results to be returned to the user
        String[] orderedQueryTerm = new String[processedQuery.size()];      // contains ...
        double[] termUpperBoundList  = new double[processedQuery.size()];   // contains all the term upper bound for each term of the query
        double[] sumTUBList  = new double[processedQuery.size()];           // array containing the sum of TUB, the value at the i-th position is the sum of TUBs from position 0 to (i-1)
        ArrayList<Posting>[] postingLists;  // contains all the posting lists for each term of the query
        ArrayList<Integer> ordListDID;      // ordered list of the DocID present in the all posting lists of the term present in the query
        double[] IDFweight;                 // array containing the IDF weight for each posting list
        int[] lengthPostingList;            // array containing the length of the posting lists
        int[] postingListsIndex;            // contain the current position index for the posting list of each term in the query
        SkipList[] skipListArray;           // array of the Skip List reference related to the term of query
        Posting currentP;                   // support var
        int firstEssPostListIndex = 0;      // indicates the index of the first (current) essential posting list
        double threshold = 0;               // var that contain the current threshold for MaxScore (is the minimum score value to be in the current best result)
        double currentDocUpperBound = 0;    // current document upper bound (used in max score algorithm for early stopping)
        double partialScore = 0;            // var that contain partial score
        int docScoreCalc = 0;               // indicates the number of documents whose score was calculated (0 to number of results requested by the user)
        boolean resetScore = false;         // used only in conjunctive case. indicates that the score must be set to 0 (the current Doc there aren't all the term of the query)
        int pLNotEmpty = 0;                 // contains the number of posting lists related to the query terms that aren't empty
        int postListCurrDID = 0;            // the DID in the current position of the posting list (used in no-essential PL)
        long startTime,endTime;             // variables to calculate the execution time

        startTime = System.currentTimeMillis();         // start time for retrieve all posting lists of the query
        postingLists = retrieveAllPostingListsMaxScore(processedQuery,orderedQueryTerm,termUpperBoundList,sumTUBList);   // take all posting lists of query terms
        endTime = System.currentTimeMillis();           // end time for retrieve all posting lists of the query
        // shows query execution time
        printTime("\n*** Retrieved all posting lists in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

        pLNotEmpty = NumberOfPostingListNotEmpty(postingLists);     // take the number of list that are not empty
        // check the number of posting lists not empty and perform the best choice
        if (pLNotEmpty == 0)    // all terms in the query aren't in the dictionary or empty query
            return;     // exit
        else if (pLNotEmpty == 1)   // there is only 1 postingList (query with one term or query with more term but only one in dictionary)
        {
            if ((postingLists.length != 1) && (isConjunctive))  // case of query conjunctive and more term but only one in dictionary
                return;     // exit

            DAATOnePostingList(processedQuery, postingLists, scoringFunc, numberOfResults);   // execute DAAT algorithm
            return;     // exit
        }

        // more postingLists not empty, use MAXSCORE algorithm
        ordListDID = DIDOrderedListOfQuery(postingLists, isConjunctive);        // take ordered list of DocID
        postingListsIndex = getPostingListsIndex(postingLists);                 // get the index initialized
        lengthPostingList = retrieveLengthAllPostingLists(orderedQueryTerm);    // take the length of each posting list
        IDFweight = calculateIDFWeight(lengthPostingList);                      // calculate the IDF weight
        skipListArray = SetAllSkipList(orderedQueryTerm, postingLists);         // create the skip List reference related to the term of query
        // control print
        //printDebug("orderedQueryTerm -> " + Arrays.toString(orderedQueryTerm));
        //printDebug("termUpperBoundList -> " + Arrays.toString(termUpperBoundList));
        //printDebug("sumTUBList -> " + Arrays.toString(sumTUBList));
        //printDebug("lengthPostingList -> " + Arrays.toString(lengthPostingList));
        /*
        for (int i = 0; i < skipListArray.length; i++)
        {
            printDebug("Control print skipping list of term: " + orderedQueryTerm[i]);
            skipListArray[i].testReadAllSkip();
        }//*/

        //printDebug("START DAAT + MAXSCORE\n");
        startTime = System.currentTimeMillis();           // start time of DAAT + MAX SCORE
        // MaxScore algorithm - scan all Doc retrieved and calculate score (TFIDF or BM25)
        for (Integer currentDID : ordListDID)
        {   // -- start - for 0: DID --
            //printDebug("START cycle with DID: " + currentDID);
            partialScore = 0;           // reset var
            resetScore = false;         // set to false

            //printDebug("-- START EP - with first EP: " + firstEssPostListIndex + " of DID: " + currentDID);
            // scan the essential posting lists, default case is query Disjunctive
            for (int j = firstEssPostListIndex; j < postingLists.length; j++)
            {   // -- start - for 0.1: EPL --
                //printDebug("-- EP - search DID: " + currentDID + " in term: '" + orderedQueryTerm[j] + "' of posting: " + j + " in pos: " + (postingListsIndex[j]));
                // check if the posting lists of j-th isn't at the end AND if the j-th term of the query is present in the doc identify by currentDID
                if ( (postingLists[j] != null) && (postingListsIndex[j] < postingLists[j].size()) && (postingLists[j].get(postingListsIndex[j]).getDocId() == currentDID))
                {
                    currentP = postingLists[j].get(postingListsIndex[j]);              // take posting
                    postingListsIndex[j]++;                         // update index of current value
                    // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                    if (scoringFunc)
                        partialScore += ScoringBM25(currentDID,currentP.getTermFreq(), IDFweight[j]);     // use BM25
                    else
                        partialScore += ScoringTFIDF(currentP.getTermFreq(), IDFweight[j]);               // use TFIDF
                    //printDebug("---- EP - find DID: " + currentDID + " and partialScore: " + partialScore + " in term: '" + orderedQueryTerm[j] + "' of posting: " + j + " in pos: " + (postingListsIndex[j]-1));
                }
                else if (isConjunctive)     // must take only the document in which there are all term (DID that compare in all posting lists of the terms)
                {
                    resetScore = true;       // reset the partial score
                    // if all postings in one posting lists have already been seen the next documents in the posting lists cannot contain all the terms in the query
                    if ((postingLists[j] == null) || (postingListsIndex[j] >= postingLists[j].size())) {
                        //printDebug("Query conjunctive, posting list numero: " + j + " finita. Si è in pos: " + postingListsIndex[j] + " su dimensione: " + postingLists[j].size());
                        endTime = System.currentTimeMillis();           // end time of DAAT
                        // shows query execution time
                        printTime("*** DAAT + MAX SCORE V.2 (skipping) execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
                        return;             // exit from function
                    }
                }
            }   // -- end - for 0.1: EPL --
            //printDebug("-- END EP -> DID : " + currentDID + " partialScore: " + partialScore + " and threshold: " + threshold + " reset score: " + resetScore);

            // Conditions under which analysis of nonessential posting lists can be skipped -- SEE NOTE 0 --
            if ( (partialScore == 0) || resetScore )
                continue;       // go to next iteration, the current doc can't be among the top result

            // scan non essential posting lists
            if (firstEssPostListIndex != 0)
            {   // -- start - if: NoEPL --
                currentDocUpperBound = partialScore + sumTUBList[firstEssPostListIndex];    // calculate the current DUB
                //printDebug("-- START NonEP - DocUB: " + currentDocUpperBound + " and threshold: " + threshold);
                // check if the doc has no zero possibility to have a score greater than threshold
                if (currentDocUpperBound <= threshold)
                    continue;                           // go to next iteration with next DID

                // update the score
                for (int i = 0; i < firstEssPostListIndex; i++)
                {
                    if (isConjunctive)
                        resetScore = true;       // reset the partial score

                    if ((postingLists[i] != null) && (postingListsIndex[i] < postingLists[i].size()))
                    {
                        //printDebug("---- NonEP HOP -> search the DID " + currentDID + " for the term: '" + orderedQueryTerm[i] + "' in the posting list: " + i + " in pos: " + postingListsIndex[i] + " and max len: " + postingLists[i].size() + " and did is: " + postingLists[i].get(postingListsIndex[i]).getDocId());
                        postListCurrDID = postingLists[i].get(postingListsIndex[i]).getDocId();

                        // check first position
                        if (postListCurrDID == currentDID)
                        {
                            // find the searched DID - update the partialScore and the currentDocUpperBound
                            currentP = postingLists[i].get(postingListsIndex[i]);              // take posting
                            //printDebug("------ NonEP -> FIND the DID " + currentDID + " old DUB: " + currentDocUpperBound + " and old partialScore: " + partialScore);
                            currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                            currentDocUpperBound -= partialScore;               // update currentDocUpperBound
                            // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                            if (scoringFunc)
                                partialScore += ScoringBM25(currentDID, currentP.getTermFreq(), IDFweight[i]);   // use BM25
                            else
                                partialScore += ScoringTFIDF(currentP.getTermFreq(), IDFweight[i]);             // use TFIDF
                            currentDocUpperBound += partialScore;               // update currentDocUpperBound

                            resetScore = false;         // reset the partial score
                            //printDebug("------ NonEP -> FIND the DID " + currentDID + " in posting: " + i + " in pos: " + postingListsIndex[i] + " with valueDID: " + postingLists[i].get(postingListsIndex[i]).getDocId() + " update DUB: " + currentDocUpperBound + " and partialScore: " + partialScore);
                            postingListsIndex[i]++;     // update the index of the postin list
                            //if (currentDocUpperBound <= threshold)
                            //    break;
                            //else
                            continue;
                        }
                        else if (postListCurrDID > currentDID)      // the searched DID is not in this posting list
                        {
                            currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                            //printDebug("------ NonEP - First position > currentDID. update DUB " + currentDocUpperBound + " partialScore: " + partialScore + " in PL: " + i + " in pos: " + postingListsIndex[i]);
                            continue;       // go to the next step
                        }
                        else    // postListCurrDID < currentDID, search in the posting listw
                        {
                            assert skipListArray != null;
                            postingListsIndex[i] = skipListArray[i].nextGEQ(currentDID, postingListsIndex[i]);
                            //printDebug("Uses skipping -> postingListSize: " + postingLists[i].size() + " search DID: " + currentDID + " and nextGEQ return position " + pos + " that have DID: " + (pos < postingLists[i].size() ? postingLists[i].get(pos).getDocId() : " out of bound"));
                        }

                        // check the index returned by nextGEQ
                        if (postingListsIndex[i] >= postingLists[i].size())  // check for out of bound in case of reaching the end of the list
                        {
                            //printDebug("NextGEQ return index: " + postingListsIndex[i] + " greater than size: " + postingLists[i].size());
                            continue;
                        }

                        postListCurrDID = postingLists[i].get(postingListsIndex[i]).getDocId(); // take the did

                        if (postListCurrDID != currentDID)
                        {
                            // should be always greater than currentDID
                            currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                            //printDebug("------ NonEP - CurrentDID there isn't in this posting list. update DUB " + currentDocUpperBound + " partialScore: " + partialScore + " in PL: " + i + " in pos: " + postingListsIndex[i]);
                        }
                        else                    // the searched DID there is in the posting list
                        {
                            // find the searched DID - update the partialScore and the currentDocUpperBound
                            currentP = postingLists[i].get(postingListsIndex[i]);              // take posting
                            //printDebug("------ NonEP -> FIND the DID " + currentDID + " old DUB: " + currentDocUpperBound + " and old partialScore: " + partialScore);
                            currentDocUpperBound -= partialScore;               // update currentDocUpperBound
                            currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                            // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                            if (scoringFunc)
                                partialScore += ScoringBM25(currentDID, currentP.getTermFreq(), IDFweight[i]);   // use BM25
                            else
                                partialScore += ScoringTFIDF(currentP.getTermFreq(), IDFweight[i]);             // use TFIDF
                            currentDocUpperBound += partialScore;               // update currentDocUpperBound

                            resetScore = false;         // reset the partial score
                            //printDebug("------ NonEP -> FIND the DID " + currentDID + " in posting: " + i + " in pos: " + postingListsIndex[i] + " with valueDID: " + postingLists[i].get(postingListsIndex[i]).getDocId() + " update DUB: " + currentDocUpperBound + " and partialScore: " + partialScore);
                            postingListsIndex[i]++;     // update the index of the postin list
                        }
                    } // if all postings in one posting lists have already been seen the next documents in the posting lists cannot contain all the terms in the query
                    else if (isConjunctive)
                    {
                        //printDebug("Query conjunctive, posting list numero: " + j + " finita. Si è in pos: " + postingListsIndex[j] + " su dimensione: " + postingLists[j].size());
                        endTime = System.currentTimeMillis();           // end time of DAAT
                        // shows query execution time
                        printTime("*** DAAT + MAX SCORE V.2 (skipping) execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
                        return;             // exit from function
                    }
                    //printDebug("-- END NonEP - DocUB: " + currentDocUpperBound + " with score: " + partialScore + " and threshold: " + threshold);
                    // check if the doc has no zero possibility to have a score greater than threshold
                    if ((currentDocUpperBound <= threshold) || resetScore)
                        break;
                }
            }   // -- end - if: NoEPL --
            // save score
            if (resetScore)
                continue;       // go to the next iteration (next Doc)
            // insert without control into priority queue (is not full) or insert all results (orderAllHashMap = true)
            if (docScoreCalc < numberOfResults)
            {
                resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                docScoreCalc++;                         // increment result in priority queue counter
                if (docScoreCalc == numberOfResults)
                    threshold = resPQ.peek().getScore();    // update threshold
                //printDebug("-- SCORING - Add result: " + partialScore);
            }
            else if (threshold < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
            {
                // substitution of the block
                resPQ.poll();                           // remove the first element
                resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                threshold = resPQ.peek().getScore();    // update threshold
                // calculate new essential posting lists and update firstEssPostListIndex
                firstEssPostListIndex = updateEssentialPositngLists(sumTUBList, threshold);
                //printDebug("-- **** New threshold: " + threshold + " new first essential posting list: " + firstEssPostListIndex);
            }
        }   // -- end - for: DID --

        endTime = System.currentTimeMillis();           // end time of DAAT
        // shows DAAT execution time
        printTime("*** DAAT + MAX SCORE V.2 (skipping) execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
    }

    /**
     * Function for apply the Document at a Time algorithm with max Score algorithm as dynamic pruning algorithm.
     * This function apply the Max Score in the case of both compression and skipping are enabled.
     *
     * @param processedQuery    array list for containing the query term
     * @param scoringFunc       indicates the preference for scoring. If false use TFIDF, if true use BM25.
     * @param isConjunctive     indicates whether the query is conjunctive or disjunctive type (default is disjunctive)
     * @param numberOfResults   indicated the max number of result to return to user
     * @throws FileNotFoundException
     */
    private static void DAATAlgMAXSCORESkipAndComp(ArrayList<String> processedQuery, boolean scoringFunc , boolean isConjunctive, int numberOfResults) throws FileNotFoundException
    {
        resPQ = new PriorityQueue<>(numberOfResults, new CompareTerm());    // length equal to the number of results to be returned to the user
        String[] orderedQueryTerm = new String[processedQuery.size()];      // contains the query term ordered by TermUpperBound (ascending order)
        double[] termUpperBoundList  = new double[processedQuery.size()];   // contains all the term upper bound for each term of the query
        double[] sumTUBList  = new double[processedQuery.size()];           // array containing the sum of TUB, the value at the i-th position is the sum of TUBs from position 0 to (i-1)
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
        long startTime,endTime;             // variables to calculate the execution time

        // check the number of the term that have PL (which are in the dictionary)
        pLNotEmpty = numberOfQueryTermsInDictionary(processedQuery);
        // check the number of posting lists not empty and perform the best choice
        if (pLNotEmpty == 0)        // all terms in the query aren't in the dictionary or empty query
            return;     // exit
        else if (pLNotEmpty == 1)   // there is only 1 postingList (query with one term or query with more term but only one in dictionary)
        {
            if ((processedQuery.size() != 1) && (isConjunctive))  // case of query conjunctive and more term but only one in dictionary
                return;     // exit

            // create the skip List reference related to the term of query
            skipListArray = skipListInitCompAndSkip(processedQuery, null, false);
            // The PL is only one -> read and decompress the whole PL and use the classic optimization method
            startTime = System.currentTimeMillis();         // start time for retrieve the posting lists of the query
            ArrayList<Posting>[] postingLists = retrieveAllUncompPL(processedQuery, skipListArray); // get the uncompress PL
            endTime = System.currentTimeMillis();           // end time for retrieve the posting lists of the query
            // shows query execution time
            printTime("*** MAX SCORE (comp+skipping) retrieved PL (case 1 PL) in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

            DAATOnePostingList(processedQuery, postingLists, scoringFunc, numberOfResults);   // execute DAAT algorithm
            return;     // exit
        }

        // -- more postingLists not empty --
        // 0) take the first block of the PL
        skipAndCompPLs = new ArrayList[processedQuery.size()];
        setAllUtilsListMAxScoreCompAndSkipping(processedQuery, orderedQueryTerm, termUpperBoundList, sumTUBList);   // calculate utilities arrays
        // create the skip List reference related to the term of query
        skipListArray = skipListInitCompAndSkip(processedQuery, orderedQueryTerm, true);

        startTime = System.currentTimeMillis();         // start time for retrieve first block of PLs of the query
        retrieveFirstCompBlockOfPLFromQuery(processedQuery,orderedQueryTerm, skipListArray, true); // retrieve the first block of each PL and put value in the PQ
        endTime = System.currentTimeMillis();           // end time for retrieve first block of PLs of the query
        printTime("*** MAX SCORE (comp+skipping) retrieved first block of each PLs in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

        blockIndex = new int[processedQuery.size()];    // set the index block for each PL
        Arrays.fill(blockIndex, 1);                 // the first block (with index 0) has already been taken
        postingListsIndex = getPostingListsIndex(skipAndCompPLs);               // get the index initialized
        lengthPostingList = retrieveLengthAllPostingLists(orderedQueryTerm);    // take the length of each posting list
        IDFweight = calculateIDFWeight(lengthPostingList);                      // calculate the IDF weight
        // control print
        /*
        printDebug("orderedQueryTerm -> " + Arrays.toString(orderedQueryTerm));
        printDebug("termUpperBoundList -> " + Arrays.toString(termUpperBoundList));
        printDebug("sumTUBList -> " + Arrays.toString(sumTUBList));
        printDebug("lengthPostingList -> " + Arrays.toString(lengthPostingList));
        */
        /*
        for (int i = 0; i < skipListArray.length; i++)
        {
            printDebug("Control print skipping list of term: " + orderedQueryTerm[i]);
            //skipListArray[i].testReadAllSkip();       // ++++++++ usare solo con posting list non compressa
        }
        //*/

        startTime = System.currentTimeMillis();           // start time of DAAT + MAX SCORE (comp + skipping)
        /*
        // control check for conjunctive case
        if (!isConjunctive)
        {   // -- start - IF disjunctive
            // MaxScore algorithm - scan all Doc retrieved and calculate score (TFIDF or BM25)
            while (!pqDID.isEmpty())
            {   // -- start - while 0: DID --
                currentDID = pqDID.poll();  // take the current DID
                //printDebug("START cycle with DID: " + currentDID);
                partialScore = 0;           // reset var

                //printDebug("-- START EP - with first EP: " + firstEssPostListIndex + " of DID: " + currentDID);
                // scan the essential posting lists, default case is query Disjunctive
                for (int j = firstEssPostListIndex; j < skipAndCompPLs.length; j++)
                {   // -- start - for 0.1: EPL --
                    //printDebug("-- EP - search DID: " + currentDID + " in term: '" + orderedQueryTerm[j] + "' of posting: " + j + " in pos: " + (postingListsIndex[j]));
                    if (skipAndCompPLs[j] != null)
                    {
                        // check if it needs to upload another block of the posting list
                        if (postingListsIndex[j] >= skipAndCompPLs[j].size())
                        {
                            // load the new block
                            if (retrieveCompBlockOfPL(orderedQueryTerm[j], skipListArray, blockIndex[j], j, true))
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

                        // check if the posting lists of j-th isn't at the end AND if the j-th term of the query is present in the doc identify by currentDID
                        if ( (postingListsIndex[j] < skipAndCompPLs[j].size()) && (skipAndCompPLs[j].get(postingListsIndex[j]).getDocId() == currentDID))
                        {
                            currentP = skipAndCompPLs[j].get(postingListsIndex[j]);              // take posting
                            postingListsIndex[j]++;                         // update index of current value
                            // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                            if (scoringFunc)
                                partialScore += ScoringBM25(currentDID,currentP.getTermFreq(), IDFweight[j]);     // use BM25
                            else
                                partialScore += ScoringTFIDF(currentP.getTermFreq(), IDFweight[j]);               // use TFIDF
                            //printDebug("---- EP - find DID: " + currentDID + " and partialScore: " + partialScore + " in term: '" + orderedQueryTerm[j] + "' of posting: " + j + " in pos: " + (postingListsIndex[j]-1));
                        }
                    }
                }   // -- end - for 0.1: EPL --

                //printDebug("-- END EP -> DID : " + currentDID + " partialScore: " + partialScore + " and threshold: " + threshold + " reset score: " + resetScore);
                // Conditions under which analysis of nonessential posting lists can be skipped -- SEE NOTE 0 --
                if ( partialScore == 0 )
                    continue;       // go to next iteration, the current doc can't be among the top result

                // scan non essential posting lists
                if (firstEssPostListIndex != 0)
                {   // -- start - if: NoEPL --
                    currentDocUpperBound = partialScore + sumTUBList[firstEssPostListIndex];    // calculate the current DUB
                    //printDebug("-- START NonEP - DocUB: " + currentDocUpperBound + " and threshold: " + threshold);
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
                                if (retrieveCompBlockOfPL(orderedQueryTerm[i], skipListArray, blockIndex[i], i, true))
                                {
                                    blockIndex[i]++;            // increment the counter of block
                                    postingListsIndex[i] = 0;   // reset the index of the posting list
                                }
                                else    // was the last block, the PL is over
                                {
                                    skipAndCompPLs[i] = null;   // set to null the posting list
                                    continue;
                                }
                            }

                            //printDebug("---- NonEP HOP -> search the DID " + currentDID + " for the term: '" + orderedQueryTerm[i] + "' in the posting list: " + i + " in pos: " + postingListsIndex[i] + " and max len: " + postingLists[i].size() + " and did is: " + postingLists[i].get(postingListsIndex[i]).getDocId());
                            postListCurrDID = skipAndCompPLs[i].get(postingListsIndex[i]).getDocId();

                            // check first position
                            if (postListCurrDID == currentDID)
                            {
                                // find the searched DID - update the partialScore and the currentDocUpperBound
                                currentP = skipAndCompPLs[i].get(postingListsIndex[i]);              // take posting
                                //printDebug("------ NonEP -> FIND the DID " + currentDID + " old DUB: " + currentDocUpperBound + " and old partialScore: " + partialScore);
                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                currentDocUpperBound -= partialScore;               // update currentDocUpperBound
                                // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                                if (scoringFunc)
                                    partialScore += ScoringBM25(currentDID, currentP.getTermFreq(), IDFweight[i]);   // use BM25
                                else
                                    partialScore += ScoringTFIDF(currentP.getTermFreq(), IDFweight[i]);             // use TFIDF
                                currentDocUpperBound += partialScore;               // update currentDocUpperBound

                                //printDebug("------ NonEP -> FIND the DID " + currentDID + " in posting: " + i + " in pos: " + postingListsIndex[i] + " with valueDID: " + postingLists[i].get(postingListsIndex[i]).getDocId() + " update DUB: " + currentDocUpperBound + " and partialScore: " + partialScore);
                                postingListsIndex[i]++;     // update the index of the postin list
                                continue;
                            }
                            else if (postListCurrDID > currentDID)      // the searched DID is not in this posting list
                            {
                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                //printDebug("------ NonEP - First position > currentDID. update DUB " + currentDocUpperBound + " partialScore: " + partialScore + " in PL: " + i + " in pos: " + postingListsIndex[i]);
                                continue;       // go to the next step
                            }
                            else    // postListCurrDID < currentDID, search in the posting listw
                            {
                                assert skipListArray != null;
                                postingListsIndex[i] = skipListArray[i].nextGEQCompSkip(currentDID, postingListsIndex[i], i,orderedQueryTerm[i]);
                                //printDebug("Uses skipping -> postingListSize: " + postingLists[i].size() + " search DID: " + currentDID + " and nextGEQ return position " + pos + " that have DID: " + (pos < postingLists[i].size() ? postingLists[i].get(pos).getDocId() : " out of bound"));
                            }
                            // check the index returned by nextGEQ
                            if (postingListsIndex[i] >= skipAndCompPLs[i].size())  // check for out of bound in case of reaching the end of the list
                            {
                                //printDebug("NextGEQ return index: " + postingListsIndex[i] + " greater than size: " + postingLists[i].size());
                                continue;
                            }
                            postListCurrDID = skipAndCompPLs[i].get(postingListsIndex[i]).getDocId(); // take the did
                            // check if the current target has been found or not
                            if (postListCurrDID != currentDID)
                            {
                                // should be always greater than currentDID
                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                //printDebug("------ NonEP - CurrentDID there isn't in this posting list. update DUB " + currentDocUpperBound + " partialScore: " + partialScore + " in PL: " + i + " in pos: " + postingListsIndex[i]);
                            }
                            else                    // the searched DID there is in the posting list
                            {
                                // find the searched DID - update the partialScore and the currentDocUpperBound
                                currentP = skipAndCompPLs[i].get(postingListsIndex[i]);              // take posting
                                //printDebug("------ NonEP -> FIND the DID " + currentDID + " old DUB: " + currentDocUpperBound + " and old partialScore: " + partialScore);
                                currentDocUpperBound -= partialScore;               // update currentDocUpperBound
                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                                if (scoringFunc)
                                    partialScore += ScoringBM25(currentDID, currentP.getTermFreq(), IDFweight[i]);   // use BM25
                                else
                                    partialScore += ScoringTFIDF(currentP.getTermFreq(), IDFweight[i]);             // use TFIDF
                                currentDocUpperBound += partialScore;               // update currentDocUpperBound

                                //printDebug("------ NonEP -> FIND the DID " + currentDID + " in posting: " + i + " in pos: " + postingListsIndex[i] + " with valueDID: " + postingLists[i].get(postingListsIndex[i]).getDocId() + " update DUB: " + currentDocUpperBound + " and partialScore: " + partialScore);
                                postingListsIndex[i]++;     // update the index of the postin list
                            }
                        } // if all postings in one posting lists have already been seen the next documents in the posting lists cannot contain all the terms in the query
                        //printDebug("-- END NonEP - DocUB: " + currentDocUpperBound + " with score: " + partialScore + " and threshold: " + threshold);
                        // check if the doc has no zero possibility to have a score greater than threshold
                        if (currentDocUpperBound <= threshold)
                            break;
                    }   // -- start - for: scan NoEPLs --
                }   // -- end - if: NoEPL --
                // save score
                // insert without control into priority queue (is not full) or insert all results (orderAllHashMap = true)
                if (docScoreCalc < numberOfResults)
                {
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                    docScoreCalc++;                         // increment result in priority queue counter
                    if (docScoreCalc == numberOfResults)
                        threshold = resPQ.peek().getScore();    // update threshold
                    //printDebug("-- SCORING - Add result: " + partialScore);
                }
                else if (threshold < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                {
                    // substitution of the block
                    resPQ.poll();                           // remove the first element
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                    threshold = resPQ.peek().getScore();    // update threshold
                    // calculate new essential posting lists and update firstEssPostListIndex
                    firstEssPostListIndex = updateEssentialPositngLists(sumTUBList, threshold);
                    //printDebug("-- **** New threshold: " + threshold + " new first essential posting list: " + firstEssPostListIndex);
                }
            }   // -- end - while 0: DID --
        }   // -- end - IF disjunctive
        else
        {   // -- start - ELSE conjunctive
            // MaxScore algorithm - scan all Doc retrieved and calculate score (TFIDF or BM25)
            while (!pqDID.isEmpty())
            {   // -- start - while 0: DID --
                currentDID = pqDID.poll();  // take the current DID
                //printDebug("START cycle with DID: " + currentDID);
                partialScore = 0;           // reset var
                resetScore = false;         // set to false

                //printDebug("-- START EP - with first EP: " + firstEssPostListIndex + " of DID: " + currentDID);
                // scan the essential posting lists, default case is query Disjunctive
                for (int j = firstEssPostListIndex; j < skipAndCompPLs.length; j++)
                {   // -- start - for 0.1: EPL --
                    //printDebug("-- EP - search DID: " + currentDID + " in term: '" + orderedQueryTerm[j] + "' of posting: " + j + " in pos: " + (postingListsIndex[j]));
                    if (skipAndCompPLs[j] != null)
                    {   // -- start - if 0.1.1: check PL null --
                        // check if it needs to upload another block of the posting list
                        if (postingListsIndex[j] >= skipAndCompPLs[j].size())
                        {
                            // load the new block
                            if (retrieveCompBlockOfPL(orderedQueryTerm[j], skipListArray, blockIndex[j], j, true))
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

                        // check if the posting lists of j-th isn't at the end AND if the j-th term of the query is present in the doc identify by currentDID
                        if ( (postingListsIndex[j] < skipAndCompPLs[j].size()) && (skipAndCompPLs[j].get(postingListsIndex[j]).getDocId() == currentDID))
                        {
                            currentP = skipAndCompPLs[j].get(postingListsIndex[j]);              // take posting
                            postingListsIndex[j]++;                         // update index of current value
                            // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                            if (scoringFunc)
                                partialScore += ScoringBM25(currentDID,currentP.getTermFreq(), IDFweight[j]);     // use BM25
                            else
                                partialScore += ScoringTFIDF(currentP.getTermFreq(), IDFweight[j]);               // use TFIDF
                            //printDebug("---- EP - find DID: " + currentDID + " and partialScore: " + partialScore + " in term: '" + orderedQueryTerm[j] + "' of posting: " + j + " in pos: " + (postingListsIndex[j]-1));
                        }
                        else if (isConjunctive)
                            resetScore = true;       // reset the partial score

                    }   // -- end - if 0.1.1: check PL null --
                    else if (isConjunctive)     // must take only the document in which there are all term (DID that compare in all posting lists of the terms)
                    {
                        //printDebug("Query conjunctive, posting list numero: " + j + " finita. Si è in pos: " + postingListsIndex[j] + " su dimensione: " + postingLists[j].size());
                        endTime = System.currentTimeMillis();           // end time of DAAT
                        // shows query execution time
                        printTime("*** MAX SCORE (comp+skipping) execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
                        return;             // exit from function
                    }
                }   // -- end - for 0.1: EPL --

                //printDebug("-- END EP -> DID : " + currentDID + " partialScore: " + partialScore + " and threshold: " + threshold + " reset score: " + resetScore);
                // Conditions under which analysis of nonessential posting lists can be skipped -- SEE NOTE 0 --
                if ( (partialScore == 0) || resetScore )
                    continue;       // go to next iteration, the current doc can't be among the top result

                // scan non essential posting lists
                if (firstEssPostListIndex != 0)
                {   // -- start - if: NoEPL --
                    currentDocUpperBound = partialScore + sumTUBList[firstEssPostListIndex];    // calculate the current DUB
                    //printDebug("-- START NonEP - DocUB: " + currentDocUpperBound + " and threshold: " + threshold);
                    // check if the doc has no zero possibility to have a score greater than threshold
                    if (currentDocUpperBound <= threshold)
                        continue;                           // go to next iteration with next DID

                    // update the score
                    for (int i = 0; i < firstEssPostListIndex; i++)
                    {   // -- start - for: scan NoEPLs --
                        if (isConjunctive)
                            resetScore = true;       // reset the partial score

                        if (skipAndCompPLs[i] != null)
                        {
                            // check if it needs to upload another block of the posting list
                            if (postingListsIndex[i] >= skipAndCompPLs[i].size())
                            {
                                // load the new block
                                if (retrieveCompBlockOfPL(orderedQueryTerm[i], skipListArray, blockIndex[i], i, true))
                                {
                                    blockIndex[i]++;            // increment the counter of block
                                    postingListsIndex[i] = 0;   // reset the index of the posting list
                                }
                                else    // was the last block, the PL is over
                                {
                                    skipAndCompPLs[i] = null;   // set to null the posting list
                                    continue;
                                }
                            }

                            //printDebug("---- NonEP HOP -> search the DID " + currentDID + " for the term: '" + orderedQueryTerm[i] + "' in the posting list: " + i + " in pos: " + postingListsIndex[i] + " and max len: " + postingLists[i].size() + " and did is: " + postingLists[i].get(postingListsIndex[i]).getDocId());
                            postListCurrDID = skipAndCompPLs[i].get(postingListsIndex[i]).getDocId();

                            // check first position
                            if (postListCurrDID == currentDID)
                            {
                                // find the searched DID - update the partialScore and the currentDocUpperBound
                                currentP = skipAndCompPLs[i].get(postingListsIndex[i]);              // take posting
                                //printDebug("------ NonEP -> FIND the DID " + currentDID + " old DUB: " + currentDocUpperBound + " and old partialScore: " + partialScore);
                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                currentDocUpperBound -= partialScore;               // update currentDocUpperBound
                                // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                                if (scoringFunc)
                                    partialScore += ScoringBM25(currentDID, currentP.getTermFreq(), IDFweight[i]);   // use BM25
                                else
                                    partialScore += ScoringTFIDF(currentP.getTermFreq(), IDFweight[i]);             // use TFIDF
                                currentDocUpperBound += partialScore;               // update currentDocUpperBound

                                resetScore = false;         // reset the partial score
                                //printDebug("------ NonEP -> FIND the DID " + currentDID + " in posting: " + i + " in pos: " + postingListsIndex[i] + " with valueDID: " + postingLists[i].get(postingListsIndex[i]).getDocId() + " update DUB: " + currentDocUpperBound + " and partialScore: " + partialScore);
                                postingListsIndex[i]++;     // update the index of the postin list
                                continue;
                            }
                            else if (postListCurrDID > currentDID)      // the searched DID is not in this posting list
                            {
                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                //printDebug("------ NonEP - First position > currentDID. update DUB " + currentDocUpperBound + " partialScore: " + partialScore + " in PL: " + i + " in pos: " + postingListsIndex[i]);
                                continue;       // go to the next step
                            }
                            else    // postListCurrDID < currentDID, search in the posting listw
                            {
                                assert skipListArray != null;
                                postingListsIndex[i] = skipListArray[i].nextGEQCompSkip(currentDID, postingListsIndex[i], i,orderedQueryTerm[i]);
                                //printDebug("Uses skipping -> postingListSize: " + postingLists[i].size() + " search DID: " + currentDID + " and nextGEQ return position " + pos + " that have DID: " + (pos < postingLists[i].size() ? postingLists[i].get(pos).getDocId() : " out of bound"));
                            }
                            // check the index returned by nextGEQ
                            if (postingListsIndex[i] >= skipAndCompPLs[i].size())  // check for out of bound in case of reaching the end of the list
                            {
                                //printDebug("NextGEQ return index: " + postingListsIndex[i] + " greater than size: " + postingLists[i].size());
                                continue;
                            }
                            postListCurrDID = skipAndCompPLs[i].get(postingListsIndex[i]).getDocId(); // take the did
                            // check if the current target has been found or not
                            if (postListCurrDID != currentDID)
                            {
                                // should be always greater than currentDID
                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                //printDebug("------ NonEP - CurrentDID there isn't in this posting list. update DUB " + currentDocUpperBound + " partialScore: " + partialScore + " in PL: " + i + " in pos: " + postingListsIndex[i]);
                            }
                            else                    // the searched DID there is in the posting list
                            {
                                // find the searched DID - update the partialScore and the currentDocUpperBound
                                currentP = skipAndCompPLs[i].get(postingListsIndex[i]);              // take posting
                                //printDebug("------ NonEP -> FIND the DID " + currentDID + " old DUB: " + currentDocUpperBound + " and old partialScore: " + partialScore);
                                currentDocUpperBound -= partialScore;               // update currentDocUpperBound
                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                                if (scoringFunc)
                                    partialScore += ScoringBM25(currentDID, currentP.getTermFreq(), IDFweight[i]);   // use BM25
                                else
                                    partialScore += ScoringTFIDF(currentP.getTermFreq(), IDFweight[i]);             // use TFIDF
                                currentDocUpperBound += partialScore;               // update currentDocUpperBound

                                resetScore = false;         // reset the partial score
                                //printDebug("------ NonEP -> FIND the DID " + currentDID + " in posting: " + i + " in pos: " + postingListsIndex[i] + " with valueDID: " + postingLists[i].get(postingListsIndex[i]).getDocId() + " update DUB: " + currentDocUpperBound + " and partialScore: " + partialScore);
                                postingListsIndex[i]++;     // update the index of the postin list
                            }
                        } // if all postings in one posting lists have already been seen the next documents in the posting lists cannot contain all the terms in the query
                        else if (isConjunctive)
                        {
                            //printDebug("Query conjunctive, posting list numero: " + j + " finita. Si è in pos: " + postingListsIndex[j] + " su dimensione: " + postingLists[j].size());
                            endTime = System.currentTimeMillis();           // end time of DAAT
                            // shows query execution time
                            printTime("*** DAAT + MAX SCORE (comp+skipping) execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
                            return;             // exit from function
                        }
                        //printDebug("-- END NonEP - DocUB: " + currentDocUpperBound + " with score: " + partialScore + " and threshold: " + threshold);
                        // check if the doc has no zero possibility to have a score greater than threshold
                        if ((currentDocUpperBound <= threshold) || resetScore)
                            break;
                    }   // -- start - for: scan NoEPLs --
                }   // -- end - if: NoEPL --
                // save score
                if (resetScore)
                    continue;       // go to the next iteration (next Doc)
                // insert without control into priority queue (is not full) or insert all results (orderAllHashMap = true)
                if (docScoreCalc < numberOfResults)
                {
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                    docScoreCalc++;                         // increment result in priority queue counter
                    if (docScoreCalc == numberOfResults)
                        threshold = resPQ.peek().getScore();    // update threshold
                    //printDebug("-- SCORING - Add result: " + partialScore);
                }
                else if (threshold < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                {
                    // substitution of the block
                    resPQ.poll();                           // remove the first element
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                    threshold = resPQ.peek().getScore();    // update threshold
                    // calculate new essential posting lists and update firstEssPostListIndex
                    firstEssPostListIndex = updateEssentialPositngLists(sumTUBList, threshold);
                    //printDebug("-- **** New threshold: " + threshold + " new first essential posting list: " + firstEssPostListIndex);
                }
            }   // -- end - while 0: DID --
        }   // -- end - ELSE conjunctive
        */
        ///*
        // MaxScore algorithm - scan all Doc retrieved and calculate score (TFIDF or BM25)
            while (!pqDID.isEmpty())
            {   // -- start - while 0: DID --
                currentDID = pqDID.poll();  // take the current DID
                //printDebug("START cycle with DID: " + currentDID);
                partialScore = 0;           // reset var
                resetScore = false;         // set to false

                //printDebug("-- START EP - with first EP: " + firstEssPostListIndex + " of DID: " + currentDID);
                // scan the essential posting lists, default case is query Disjunctive
                for (int j = firstEssPostListIndex; j < skipAndCompPLs.length; j++)
                {   // -- start - for 0.1: EPL --
                    //printDebug("-- EP - search DID: " + currentDID + " in term: '" + orderedQueryTerm[j] + "' of posting: " + j + " in pos: " + (postingListsIndex[j]));
                    if (skipAndCompPLs[j] != null)
                    {   // -- start - if 0.1.1: check PL null --
                        // check if it needs to upload another block of the posting list
                        if (postingListsIndex[j] >= skipAndCompPLs[j].size())
                        {
                            // load the new block
                            if (retrieveCompBlockOfPL(orderedQueryTerm[j], skipListArray, blockIndex[j], j, true))
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

                        // check if the posting lists of j-th isn't at the end AND if the j-th term of the query is present in the doc identify by currentDID
                        if ( (postingListsIndex[j] < skipAndCompPLs[j].size()) && (skipAndCompPLs[j].get(postingListsIndex[j]).getDocId() == currentDID))
                        {
                            currentP = skipAndCompPLs[j].get(postingListsIndex[j]);              // take posting
                            postingListsIndex[j]++;                         // update index of current value
                            // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                            if (scoringFunc)
                                partialScore += ScoringBM25(currentDID,currentP.getTermFreq(), IDFweight[j]);     // use BM25
                            else
                                partialScore += ScoringTFIDF(currentP.getTermFreq(), IDFweight[j]);               // use TFIDF
                            //printDebug("---- EP - find DID: " + currentDID + " and partialScore: " + partialScore + " in term: '" + orderedQueryTerm[j] + "' of posting: " + j + " in pos: " + (postingListsIndex[j]-1));
                        }
                        else if (isConjunctive)
                            resetScore = true;       // reset the partial score

                    }   // -- end - if 0.1.1: check PL null --
                    else if (isConjunctive)     // must take only the document in which there are all term (DID that compare in all posting lists of the terms)
                    {
                        //printDebug("Query conjunctive, posting list numero: " + j + " finita. Si è in pos: " + postingListsIndex[j] + " su dimensione: " + postingLists[j].size());
                        endTime = System.currentTimeMillis();           // end time of DAAT
                        // shows query execution time
                        printTime("*** MAX SCORE (comp+skipping) execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
                        return;             // exit from function
                    }
                }   // -- end - for 0.1: EPL --

                //printDebug("-- END EP -> DID : " + currentDID + " partialScore: " + partialScore + " and threshold: " + threshold + " reset score: " + resetScore);
                // Conditions under which analysis of nonessential posting lists can be skipped -- SEE NOTE 0 --
                if ( (partialScore == 0) || resetScore )
                    continue;       // go to next iteration, the current doc can't be among the top result

                // scan non essential posting lists
                if (firstEssPostListIndex != 0)
                {   // -- start - if: NoEPL --
                    currentDocUpperBound = partialScore + sumTUBList[firstEssPostListIndex];    // calculate the current DUB
                    //printDebug("-- START NonEP - DocUB: " + currentDocUpperBound + " and threshold: " + threshold);
                    // check if the doc has no zero possibility to have a score greater than threshold
                    if (currentDocUpperBound <= threshold)
                        continue;                           // go to next iteration with next DID

                    // update the score
                    for (int i = 0; i < firstEssPostListIndex; i++)
                    {   // -- start - for: scan NoEPLs --
                        if (isConjunctive)
                            resetScore = true;       // reset the partial score

                        if (skipAndCompPLs[i] != null)
                        {
                            // check if it needs to upload another block of the posting list
                            if (postingListsIndex[i] >= skipAndCompPLs[i].size())
                            {
                                // load the new block
                                if (retrieveCompBlockOfPL(orderedQueryTerm[i], skipListArray, blockIndex[i], i, true))
                                {
                                    blockIndex[i]++;            // increment the counter of block
                                    postingListsIndex[i] = 0;   // reset the index of the posting list
                                }
                                else    // was the last block, the PL is over
                                {
                                    skipAndCompPLs[i] = null;   // set to null the posting list
                                    continue;
                                }
                            }

                            //printDebug("---- NonEP HOP -> search the DID " + currentDID + " for the term: '" + orderedQueryTerm[i] + "' in the posting list: " + i + " in pos: " + postingListsIndex[i] + " and max len: " + postingLists[i].size() + " and did is: " + postingLists[i].get(postingListsIndex[i]).getDocId());
                            postListCurrDID = skipAndCompPLs[i].get(postingListsIndex[i]).getDocId();

                            // check first position
                            if (postListCurrDID == currentDID)
                            {
                                // find the searched DID - update the partialScore and the currentDocUpperBound
                                currentP = skipAndCompPLs[i].get(postingListsIndex[i]);              // take posting
                                //printDebug("------ NonEP -> FIND the DID " + currentDID + " old DUB: " + currentDocUpperBound + " and old partialScore: " + partialScore);
                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                currentDocUpperBound -= partialScore;               // update currentDocUpperBound
                                // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                                if (scoringFunc)
                                    partialScore += ScoringBM25(currentDID, currentP.getTermFreq(), IDFweight[i]);   // use BM25
                                else
                                    partialScore += ScoringTFIDF(currentP.getTermFreq(), IDFweight[i]);             // use TFIDF
                                currentDocUpperBound += partialScore;               // update currentDocUpperBound

                                resetScore = false;         // reset the partial score
                                //printDebug("------ NonEP -> FIND the DID " + currentDID + " in posting: " + i + " in pos: " + postingListsIndex[i] + " with valueDID: " + postingLists[i].get(postingListsIndex[i]).getDocId() + " update DUB: " + currentDocUpperBound + " and partialScore: " + partialScore);
                                postingListsIndex[i]++;     // update the index of the postin list
                                continue;
                            }
                            else if (postListCurrDID > currentDID)      // the searched DID is not in this posting list
                            {
                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                //printDebug("------ NonEP - First position > currentDID. update DUB " + currentDocUpperBound + " partialScore: " + partialScore + " in PL: " + i + " in pos: " + postingListsIndex[i]);
                                continue;       // go to the next step
                            }
                            else    // postListCurrDID < currentDID, search in the posting listw
                            {
                                //printDebug("Use skipping");
                                assert skipListArray != null;
                                postingListsIndex[i] = skipListArray[i].nextGEQCompSkip(currentDID, postingListsIndex[i], i,orderedQueryTerm[i]);
                                //printDebug("Uses skipping -> postingListSize: " + postingLists[i].size() + " search DID: " + currentDID + " and nextGEQ return position " + pos + " that have DID: " + (pos < postingLists[i].size() ? postingLists[i].get(pos).getDocId() : " out of bound"));
                            }
                            // check the index returned by nextGEQ
                            if (postingListsIndex[i] >= skipAndCompPLs[i].size())  // check for out of bound in case of reaching the end of the list
                            {
                                //printDebug("NextGEQ return index: " + postingListsIndex[i] + " greater than size: " + postingLists[i].size());
                                continue;
                            }
                            postListCurrDID = skipAndCompPLs[i].get(postingListsIndex[i]).getDocId(); // take the did
                            // check if the current target has been found or not
                            if (postListCurrDID != currentDID)
                            {
                                // should be always greater than currentDID
                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                //printDebug("------ NonEP - CurrentDID there isn't in this posting list. update DUB " + currentDocUpperBound + " partialScore: " + partialScore + " in PL: " + i + " in pos: " + postingListsIndex[i]);
                            }
                            else                    // the searched DID there is in the posting list
                            {
                                // find the searched DID - update the partialScore and the currentDocUpperBound
                                currentP = skipAndCompPLs[i].get(postingListsIndex[i]);              // take posting
                                //printDebug("------ NonEP -> FIND the DID " + currentDID + " old DUB: " + currentDocUpperBound + " and old partialScore: " + partialScore);
                                currentDocUpperBound -= partialScore;               // update currentDocUpperBound
                                currentDocUpperBound -= termUpperBoundList[i];      // update currentDocUpperBound
                                // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                                if (scoringFunc)
                                    partialScore += ScoringBM25(currentDID, currentP.getTermFreq(), IDFweight[i]);   // use BM25
                                else
                                    partialScore += ScoringTFIDF(currentP.getTermFreq(), IDFweight[i]);             // use TFIDF
                                currentDocUpperBound += partialScore;               // update currentDocUpperBound

                                resetScore = false;         // reset the partial score
                                //printDebug("------ NonEP -> FIND the DID " + currentDID + " in posting: " + i + " in pos: " + postingListsIndex[i] + " with valueDID: " + postingLists[i].get(postingListsIndex[i]).getDocId() + " update DUB: " + currentDocUpperBound + " and partialScore: " + partialScore);
                                postingListsIndex[i]++;     // update the index of the postin list
                            }
                        } // if all postings in one posting lists have already been seen the next documents in the posting lists cannot contain all the terms in the query
                        else if (isConjunctive)
                        {
                            //printDebug("Query conjunctive, posting list numero: " + j + " finita. Si è in pos: " + postingListsIndex[j] + " su dimensione: " + postingLists[j].size());
                            endTime = System.currentTimeMillis();           // end time of DAAT
                            // shows query execution time
                            printTime("*** DAAT + MAX SCORE (comp+skipping) execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
                            return;             // exit from function
                        }
                        //printDebug("-- END NonEP - DocUB: " + currentDocUpperBound + " with score: " + partialScore + " and threshold: " + threshold);
                        // check if the doc has no zero possibility to have a score greater than threshold
                        if ((currentDocUpperBound <= threshold) || resetScore)
                            break;
                    }   // -- start - for: scan NoEPLs --
                }   // -- end - if: NoEPL --
                // save score
                if (resetScore)
                    continue;       // go to the next iteration (next Doc)
                // insert without control into priority queue (is not full) or insert all results (orderAllHashMap = true)
                if (docScoreCalc < numberOfResults)
                {
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                    docScoreCalc++;                         // increment result in priority queue counter
                    if (docScoreCalc == numberOfResults)
                        threshold = resPQ.peek().getScore();    // update threshold
                    //printDebug("-- SCORING - Add result: " + partialScore);
                }
                else if (threshold < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                {
                    // substitution of the block
                    resPQ.poll();                           // remove the first element
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                    threshold = resPQ.peek().getScore();    // update threshold
                    // calculate new essential posting lists and update firstEssPostListIndex
                    firstEssPostListIndex = updateEssentialPositngLists(sumTUBList, threshold);
                    //printDebug("-- **** New threshold: " + threshold + " new first essential posting list: " + firstEssPostListIndex);
                }
            }   // -- end - while 0: DID --
        //*/
        endTime = System.currentTimeMillis();           // end time of DAAT
        // shows DAAT execution time
        printTime("*** MAX SCORE (comp+skipping) execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
    }
    ///* ---- NEW VERSION OF DAAT -- END ----

    // -------------------------------------------- END - Execution Query alg ------------------------------------------

    /**
     * Function to calculate the number of the query terms that are not in the dictionary.
     *
     * @param processedQuery    array list for containing the query term
     * @return      return the number of query terms that are in the dictionary.
     *              (i.e. the number of query terms that have a posting list not null).
     */
    private static int numberOfQueryTermsInDictionary(ArrayList<String> processedQuery)
    {
        int counter = 0;    // counter of the query terms in the dictionary

        for (String term : processedQuery)      // scan all query terms
        {
            // there is a posting list for the query term == the term is in the collection
            if (dictionary.getTermToTermStat().containsKey(term))
                counter++;
        }

        return counter;
    }

    /**
     * Function that given the posting lists of each term in a given query returns the number of posting lists that
     * aren't empty
     *
     * @param postingLists  the posting lists of each term in the query
     * @return the number of posting lists that aren't empty
     */
    private static int NumberOfPostingListNotEmpty (ArrayList<Posting>[] postingLists)
    {
        int postingListNotEmpty = 0;    // contains the number of posting lists related to the query terms that aren't empty

        // control check for empty posting lists (the terms are not present in the document collection)
        if (postingLists.length == 0)
        {
            printUI("The terms in query there aren't in collection.");
            return 0;     // exit to function
        }
        // scan all posting lists
        for (int i = 0; i < postingLists.length; i++)
        {
            if (postingLists[i] == null)    // term that there isn't in collection -> posting list == null
                continue;                   // go to next posting list
            postingListNotEmpty++;          // increment counter
        }

        return postingListNotEmpty;
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

        processedQuery = new ArrayList<>();
        processedQuery.add(term);           // insert the term

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
        for (Posting p : postingLists[0])
        {
            partialScore = 0;               // reset var
            // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
            if (scoringFunc)
                partialScore += ScoringBM25(p.getDocId(), p.getTermFreq(), IDFWeight); // use BM25
            else
                partialScore += ScoringTFIDF(p.getTermFreq(), IDFWeight);              // use TFIDF

            // save score if is the new current term upper bound
            if ( partialScore > maxScore)
            {
                maxScore = partialScore;
            }
            // compute term freq statistics
            if (computeStats)
                CollectionStatistics.addTFOccToTermFreqTable(p.getTermFreq());
        }
        //printDebug("Term upper bound for term: " + term + " is: " + maxScore);  // control print
        return maxScore;        // return term upper bound
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

        // control to avoid log and division to 0
        if (termFreq == 0)
            return (double) 0;

        TFweight = termFreqWeightTable.get(termFreq);   // calculate TF weight, new version
        scoreTFIDF = TFweight * IDFweight;          // calculate TFIDF weight from Tf and IDF weight values
        //printDebug("ScoringTFIDF - TFweight = " + TFweight + " IDFweight = " + IDFweight + " scoreTFIDF = " + scoreTFIDF);

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
                termFreqWeightTable.put(i,TFweight);    // put into hash table
                //printDebug("termFreqWeightTable ( " + i + " , " + TFweight + " )");
                // store
                buffer.putDouble(TFweight);             // write termFreqWeight
            }
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
        printDebug("The maxTF is: " + maxTF + " and the size of termFreqWeightTable is: " + termFreqWeightTable.size());
        endTime = System.currentTimeMillis();
        printTime("TermFreqWeight calculated and stored in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
    }

    /**
     * Function to read from disk the precomputed TermFreqWeight.
     */
    public static void readTFWeightFromDisk()
    {
        long startTime, endTime;
        int tf = 1;

        printLoad("Loading all useful values for termFreqWeight from disk...");

        startTime = System.currentTimeMillis();
        if (!termFreqWeightTable.isEmpty())
            termFreqWeightTable.clear();

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
                termFreqWeightTable.put(tf, buffer.getDouble());
                tf++;
            }
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
        endTime = System.currentTimeMillis();
        printTime("TermFreqWeight loaded in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
        //printDebug("The maxTF is: " + CollectionStatistics.getMaxTermFreq() + " and the size of termFreqWeightTable is: " + termFreqWeightTable.size());
        //printDebug("TermFreqWeight hashtable: " + termFreqWeightTable);
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

        // control to avoid log and division to 0
        if (termFreq == 0)
            return (double) 0;

        denominator = documentTable.get(DocID).getDenomPartBM25() + termFreq;
        scoreBM25 = (termFreq / denominator) * IDFweight;      // calculate TFIDF weight from Tf and IDF weight values
        //printDebug("ScoringBM25 - docLen = " + docLen + " denominator = " + denominator + " IDFweight = " + IDFweight + " scoreBM25 = " + scoreBM25);

        return scoreBM25;
    }

    // ------------------------ end: scoring function ------------------------

    // ------------------------ start: function to retrieve PL of the term in query ------------------------

    /**
     * Function to retrieve all the posting lists for each term of the query passed as parameter.
     * This function is used when compression and skipping are not enabled or only one of them is enabled.
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
                // there is a posting list for the query term == the term is in the collection
                if (dictionary.getTermToTermStat().containsKey(term))
                {
                    //printDebug("DAAT: retrieve posting list of  " + term);
                    DictionaryElem de = dictionary.getTermToTermStat().get(term);

                    if(Flags.isCompressionEnabled() && !Flags.considerSkippingBytes())  // if the compression is enabled and the skipping is not enabled
                        postingLists[iterator] = readCompressedPostingListFromDisk(de.getOffsetDocId(),de.getOffsetTermFreq(), de.getTermFreqSize(), de.getDocIdSize(), de.getDf(), docIdChannel, termFreqChannel); //read compressed posting list
                    else    // take the postingList of term
                        postingLists[iterator] = readPostingListFromDisk(de.getOffsetDocId(),de.getOffsetTermFreq(),de.getDf(),docIdChannel,termFreqChannel);
                }
                else        // there isn't a posting list for the query term == the term isn't in the collection
                {
                    termNotInCollection.add(term);  // update array list of the term not in collection
                }
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
        // array of arrayList (posting list) that contain all the posting lists for each term in the query
        ArrayList<Posting>[] postingLists = new ArrayList[processedQuery.size()];
        // priority queue for ordering the posting list according to term upper bound
        PriorityQueue<QueryProcessor.TermUpperBoundBlock> pq = new PriorityQueue<>(processedQuery.size(), new CompareTUBTerm());
        TermUpperBoundBlock tempTUBblock;   // the current element of the pq taken at each iteration
        double sumTUB = 0;                  // used to calculate the value for the 'sumTUBList'
        int iterator = 0;                   // iterator for saving posting lists term in correct position

        // control check
        if ( (orderedQueryTerm.length != processedQuery.size()) || (termUpperBoundList.length != processedQuery.size()) || (sumTUBList.length != processedQuery.size()))
        {
            printError("Error in retrieveAllPostingListsMaxScore: wrong length in orderedQueryTerm or in termUpperBoundList or in sumTUBList.");
            return postingLists;
        }

        // retrieve the term upper bound for each posting lists and put into PQ
        for (int i = 0; i < processedQuery.size(); i++)
            pq.add(new QueryProcessor.TermUpperBoundBlock(i, TermDocUpperBound.getTermUpperBound(processedQuery.get(i))));     // add to priority queue

        // extract the ordered posting lists and related terms and insert them in the array of the term
        for (int i = 0; i < processedQuery.size(); i++)
        {
            tempTUBblock = pq.poll();    // get block
            assert tempTUBblock != null;
            orderedQueryTerm[i] = processedQuery.get(tempTUBblock.getTermPosition());   // get the term
            termUpperBoundList[i] = tempTUBblock.getTermUpperBound();                   // get term upper bound
            sumTUBList[i] = sumTUB;                 // get the term upper bound sum of the previous posting lists
            sumTUB += termUpperBoundList[i];
            //printDebug("Ordered query term -> Position: " + i + " term: '" + orderedQueryTerm[i] + "' with TUB: " + termUpperBoundList[i]); // control print
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

                // there is a posting list for the query term == the term is in the collection
                if (dictionary.getTermToTermStat().containsKey(term))
                {
                    if(Flags.isCompressionEnabled())
                        postingLists[iterator] = readCompressedPostingListFromDisk(de.getOffsetDocId(),de.getOffsetTermFreq(), de.getTermFreqSize(), de.getDocIdSize(), de.getDf(), docIdChannel, termFreqChannel); //read compressed posting list
                    else    // take the postingList of term
                        postingLists[iterator] = readPostingListFromDisk(de.getOffsetDocId(),de.getOffsetTermFreq(),de.getDf(),docIdChannel,termFreqChannel);
                }
                else        // there isn't a posting list for the query term == the term isn't in the collection
                {
                    termNotInCollection.add(term);  // update array list of the term not in collection
                }
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
        // priority queue for ordering the posting list according to term upper bound
        PriorityQueue<QueryProcessor.TermUpperBoundBlock> pq = new PriorityQueue<>(processedQuery.size(), new CompareTUBTerm());
        TermUpperBoundBlock tempTUBblock;   // the current element of the pq taken at each iteration
        double sumTUB = 0;                  // used to calculate the value for the 'sumTUBList'

        // control check
        if ( (orderedQueryTerm.length != processedQuery.size()) || (termUpperBoundList.length != processedQuery.size()) || (sumTUBList.length != processedQuery.size()))
            printError("Error in retrieveAllPostingListsMaxScore: wrong length in orderedQueryTerm or in termUpperBoundList or in sumTUBList.");

        // retrieve the term upper bound for each posting lists and put into PQ
        for (int i = 0; i < processedQuery.size(); i++)
            pq.add(new QueryProcessor.TermUpperBoundBlock(i, TermDocUpperBound.getTermUpperBound(processedQuery.get(i))));     // add to priority queue

        // extract the ordered posting lists and related terms and insert them in the array of the term
        for (int i = 0; i < processedQuery.size(); i++)
        {
            tempTUBblock = pq.poll();               // get block
            assert tempTUBblock != null;
            orderedQueryTerm[i] = processedQuery.get(tempTUBblock.getTermPosition());   // get the term
            termUpperBoundList[i] = tempTUBblock.getTermUpperBound();                   // get term upper bound
            sumTUBList[i] = sumTUB;                 // get the term upper bound sum of the previous posting lists
            sumTUB += termUpperBoundList[i];        // update sumTUB
            //printDebug("Ordered query term -> Position: " + i + " term: '" + orderedQueryTerm[i] + "' with TUB: " + termUpperBoundList[i]); // control print
        }
    }

    /**
     * Function to read from disk and decompress the first compressed block of the posting list for each query term.
     * If is for DAAT the term freq e DID read are put in a priority queue (parameter maxScore = false).
     * If is for Max Score the term freq and DID are used for creating posting lists, and DIDs are also used to form the
     * priority queue for algorithm execution (parameter maxScore = true).
     *
     * @param processedQuery    ArrayList  of the processed terms of the query (used only in the case of simple DAAT)
     * @param orderedQueryTerm  array of string to contain the query terms sorted by term upper bound. (used only in the case of Max Score)
     * @param slArr             the array of SkipInfo to initialize with object related to the query term
     * @param maxScore          if 'true' is the case of max score enabled, if 'false' max score is not enabled use simple DAAT
     * @return  an array of posting lists (ArrayList of posting). the array has length equal to the number of terms,
     *          and the i-th position in the array contains the posting list of the i-th term in the processedQuery
     */
    private static void retrieveFirstCompBlockOfPLFromQuery(ArrayList<String> processedQuery, String[] orderedQueryTerm, SkipList[] slArr, boolean maxScore)
    {
        byte[] tf;                      // array for the compressed TermFreq list
        byte[] docids;                  // array for the compressed TermFreq list
        int iterator = 0;               // iterator for saving posting lists term in correct position
        int currDID = 0;                // the current value of uncompressed DID (case max score)

        try(
                // open complete files to read the postingList
                RandomAccessFile docidFile = new RandomAccessFile(DOCID_FILE, "rw");
                RandomAccessFile termfreqFile = new RandomAccessFile(TERMFREQ_FILE, "rw");
                // FileChannel
                FileChannel docIdChannel = docidFile.getChannel();
                FileChannel termFreqChannel = termfreqFile.getChannel()
        ) {
            // check if is the case of DAAT or Max Score
            if (maxScore)   // Max Score case
            {
                // take posting list for each term in query ordered by Term Upper Bound
                for (String term : orderedQueryTerm)
                {
                    // there is a posting list for the query term == the term is in the collection
                    if (dictionary.getTermToTermStat().containsKey(term))
                    {
                        //printDebug("MAX SCORE (comp+skipping): retrieve posting list of  " + term);
                        DictionaryElem de = dictionary.getTermToTermStat().get(term);

                        // get the compressed block
                        tf = readCompTFBlockFromDisk(slArr[iterator], 0, de.getOffsetTermFreq(), de.getTermFreqSize(), de.getSkipArrLen(), termFreqChannel);
                        docids = readCompDIDBlockFromDisk(slArr[iterator], 0, de.getOffsetDocId(), de.getDocIdSize(), de.getSkipArrLen(), docIdChannel);

                        int numTFComp = min(SKIP_POINTERS_THRESHOLD, de.getDf());
                        // decompress the block
                        ArrayList<Integer> uncompressedTf = Unary.integersDecompression(tf, numTFComp);  // decompress term freq
                        ArrayList<Integer> uncompressedDocid = VariableBytes.integersDecompression(docids,true);    // decompress DocID

                        // add the block to the related PL
                        skipAndCompPLs[iterator] = new ArrayList<>();    // decompressed posting list
                        for (int i = 0; i < numTFComp; i++)
                        {
                            currDID = uncompressedDocid.get(i);     // get DID
                            // add the posting to the posting list
                            skipAndCompPLs[iterator].add(new Posting(currDID, uncompressedTf.get(i)));
                            // add the block to PQ (pqDID)
                            if (!pqDID.contains(currDID))
                                pqDID.add(currDID);                 // add to priority queue
                        }
                        slArr[iterator].setCurrPostList(skipAndCompPLs[iterator]);  // set the current uncompressed block of PL in SkipList instance
                    }
                    else        // there isn't a posting list for the query term == the term isn't in the collection
                    {
                        termNotInCollection.add(term);  // update array list of the term not in collection
                    }
                    iterator++;                 // update iterator
                }
            }
            else            // DAAT case
            {
                // take posting list for each term in query ordered by entry in the query
                for (String term : processedQuery)
                {
                    // there is a posting list for the query term == the term is in the collection
                    if (dictionary.getTermToTermStat().containsKey(term))
                    {
                        //printDebug("DAAT: retrieve posting list of  " + term);
                        DictionaryElem de = dictionary.getTermToTermStat().get(term);

                        // get the compressed block
                        tf = readCompTFBlockFromDisk(slArr[iterator], 0,de.getOffsetTermFreq(), de.getTermFreqSize(), de.getSkipArrLen(), termFreqChannel);
                        docids = readCompDIDBlockFromDisk(slArr[iterator], 0, de.getOffsetDocId(), de.getDocIdSize(), de.getSkipArrLen(), docIdChannel);

                        int numTFComp = min(SKIP_POINTERS_THRESHOLD, de.getDf());
                        // decompress the block
                        ArrayList<Integer> uncompressedTf = Unary.integersDecompression(tf, numTFComp);  // decompress term freq
                        ArrayList<Integer> uncompressedDocid = VariableBytes.integersDecompression(docids,true);    // decompress DocID

                        // add the block to PQ (ordDIDPQ)
                        for (int i = 0; i < numTFComp; i++)
                            ordDIDPQ.add(new QueryProcessor.orderedDIDBlock(uncompressedDocid.get(i),uncompressedTf.get(i),iterator));     // add to priority queue
                    }
                    else        // there isn't a posting list for the query term == the term isn't in the collection
                    {
                        termNotInCollection.add(term);  // update array list of the term not in collection
                    }
                    iterator++;                 // update iterator
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Function to retrieve one compressed block of a posting list and uncompress it.
     *
     * @param term          term related to the posting list
     * @param slArr         the array of the SkipList of the query terms
     * @param blockIndex    the position of the block to load from disk
     * @param indexPL       the index of the posting list
     * @param maxScore      if 'true' is the case of max score enabled, if 'false' max score is not enabled use simple DAAT
     * @return              'true' if a block has been read, decompressed, and added, 'false' otherwhise
     */
    private static boolean retrieveCompBlockOfPL(String term, SkipList[] slArr, int blockIndex, int indexPL, boolean maxScore)
    {
        byte[] tf;                      // array for the compressed TermFreq list
        byte[] docids;                  // array for the compressed TermFreq list
        int currDID = 0;                // the current value of uncompressed DID (case max score)

        // control check
        if (indexPL >= slArr.length)
        {
            printError("retrieveCompBlockOfPL -> error in the parameter. indexPL: " + indexPL + " and blockIndex: " + blockIndex);
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
            // there is a posting list for the query term == the term is in the collection
            if (dictionary.getTermToTermStat().containsKey(term))
            {
                DictionaryElem de = dictionary.getTermToTermStat().get(term);
                // check if is the case of DAAT or Max Score
                if (maxScore)   // Max Score case
                {
                    //printDebug("MAXSCORE: retrieve posting list of  " + term);
                    // get the compressed block (the control check of the block index is in this function)
                    tf = readCompTFBlockFromDisk(slArr[indexPL], blockIndex,de.getOffsetTermFreq(), de.getTermFreqSize(), de.getSkipArrLen(), termFreqChannel);
                    docids = readCompDIDBlockFromDisk(slArr[indexPL], blockIndex, de.getOffsetDocId(), de.getDocIdSize(), de.getSkipArrLen(), docIdChannel);
                    int numTFComp = min(SKIP_POINTERS_THRESHOLD, (de.getDf() - (SKIP_POINTERS_THRESHOLD * blockIndex)));

                    if ( (tf == null) || (docids == null) )
                        return false;
                    // decompress the block
                    ArrayList<Integer> uncompressedTf = Unary.integersDecompression(tf, numTFComp);  // decompress term freq
                    ArrayList<Integer> uncompressedDocid = VariableBytes.integersDecompression(docids,true);    // decompress DocID
                    //printDebug("Request Block: " + blockIndex + " related to the term: " + term + " with skipArr len: " + de.getSkipArrLen());
                    //printDebug("uncompressedTf len: " + uncompressedTf.size() + " uncompressedDocid len: " + uncompressedDocid.size());
                    // add the block to the related PL
                    skipAndCompPLs[indexPL] = new ArrayList<>();    // decompressed posting list
                    for (int i = 0; i < numTFComp; i++)
                    {
                        currDID = uncompressedDocid.get(i);     // get DID
                        // add the posting to the posting list
                        skipAndCompPLs[indexPL].add(new Posting(currDID, uncompressedTf.get(i)));
                        // add the block to PQ (pqDID)
                        if (!pqDID.contains(currDID))
                            pqDID.add(currDID);                 // add to priority queue
                    }
                    slArr[indexPL].setCurrPostList(skipAndCompPLs[indexPL]);  // set the current uncompressed block of PL in SkipList instance
                }
                else            // DAAT case
                {
                    //printDebug("DAAT: retrieve posting list of  " + term);
                    // get the compressed block (the control check of the block index is in this function)
                    tf = readCompTFBlockFromDisk(slArr[indexPL], blockIndex,de.getOffsetTermFreq(), de.getTermFreqSize(), de.getSkipArrLen(), termFreqChannel);
                    docids = readCompDIDBlockFromDisk(slArr[indexPL], blockIndex, de.getOffsetDocId(), de.getDocIdSize(), de.getSkipArrLen(), docIdChannel);

                    int numTFComp = min(SKIP_POINTERS_THRESHOLD, (de.getDf() - (SKIP_POINTERS_THRESHOLD * blockIndex)));
                    if ( (tf == null) || (docids == null) )
                        return false;
                    // decompress the block
                    ArrayList<Integer> uncompressedTf = Unary.integersDecompression(tf, numTFComp);  // decompress term freq
                    ArrayList<Integer> uncompressedDocid = VariableBytes.integersDecompression(docids,true);    // decompress DocID
                    //printDebug("Request Block: " + blockIndex + " related to the term: " + term + " with skipArr len: " + de.getSkipArrLen());
                    //printDebug("uncompressedTf len: " + uncompressedTf.size() + " uncompressedDocid len: " + uncompressedDocid.size());
                    // add the block to PQ (ordDIDPQ)
                    for (int i = 0; i < numTFComp; i++)
                        ordDIDPQ.add(new QueryProcessor.orderedDIDBlock(uncompressedDocid.get(i),uncompressedTf.get(i),indexPL));     // add to priority queue
                }
            }
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
                // there is a posting list for the query term == the term is in the collection
                if (dictionary.getTermToTermStat().containsKey(term))
                {
                    //printDebug("DAAT: retrieve posting list of  " + term);
                    DictionaryElem de = dictionary.getTermToTermStat().get(term);
                    // read adn uncompress the whole posting list related to term
                    postingLists[iterator] = readAndUncompressCompressedAndSkippedPLFromDisk(slArr[iterator], de.getOffsetDocId(), de.getOffsetTermFreq(), de.getTermFreqSize(), de.getDocIdSize(), de.getSkipArrLen(), de.getDf(), docIdChannel, termFreqChannel);
                }
                else        // there isn't a posting list for the query term == the term isn't in the collection
                {
                    termNotInCollection.add(term);  // update array list of the term not in collection
                }
                iterator++;                 // update iterator
            }
            return postingLists;

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
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
        ArrayList<Integer> rankedResults = new ArrayList<>();   // array list to contain the top "numResults" docs
        String currDocNO;           // indicates the DocNO of the current document in the result (top 'numResults' docs)
        long startTime, endTime;    // variables to calculate the execution time

        if (numResults <= 0)        // control check
            return rankedResults;

        startTime = System.currentTimeMillis();         // start time of hash map ordering

        QueryProcessor.ResultBlock currentResPQ;        // var that contain the resultBlock extract from pq in the current iteration
        ArrayList<Integer> results = new ArrayList<Integer>();       // Create an ArrayList object

        while(!resPQ.isEmpty())                         // control if the priority queue for results is empty
        {
            //printDebug("Taken: " + resPQ.peek().getDID() + " with score: " + resPQ.peek().getScore());
            currentResPQ = resPQ.poll();                                    // take the lowest element (score and DID)
            currDocNO = documentTable.get(currentResPQ.getDID()).getDocno();// take the DocNo related to the DID
            try{
                results.add(Integer.valueOf(currDocNO));                    // add to the array list
            }
            catch (NumberFormatException ex){
                ex.printStackTrace();
            }
        }
        // order the result from the best to the worst (reverse order of the priority queue)
        rankedResults = new ArrayList<Integer>(results);     // Create an ArrayList object
        Collections.reverse(rankedResults);
        //printDebug("orderedResults " + rankedResults);

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
            assert terms[i] != null;
            if (dictionary.getTermToTermStat().get(terms[i]) != null)
                lengthPostingList[i] = dictionary.getTermToTermStat().get(terms[i]).getDf();
            else
                lengthPostingList[i] = 1;
            //printDebug("The length for the posting list of the term: '" + terms[i] + "' is : " + lengthPostingList[i]);
        }

        return lengthPostingList;
    }

    /**
     * function to update the essential posting lists by the threshold passed as parameter
     *
     * @param sumTUBList    array of double to contain the sum of the term upper bound of the previous position
     * @param threshold     the current threshold passes as parameter
     * @return  an integer that indicate the index of the first new essential posting list. The new essential posting
     *          list will be the posting lists with index between the returned integer and the last posting lists
     */
    private static int updateEssentialPositngLists (double[] sumTUBList, double threshold)
    {
        for (int i = 1; i < sumTUBList.length; i++)
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
            if (postingLists[i] == null)    // control check
            {
                //printDebug("SetAllSkipList -- term: " + orderedQueryTerm[i] + " not in collection.");
                tempSkipListArray[i] = null;
                continue;
            }
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
            if (orderedQueryTerm.length == 0)   // control check
                return null;

            tempSkipListArray = new SkipList[orderedQueryTerm.length];  // set array

            for (int i = 0; i < orderedQueryTerm.length; i++)       // scan all query terms
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
    private static ArrayList<Integer> DIDOrderedListOfQuery(ArrayList<Posting>[] postingLists, boolean isConjunctive) throws FileNotFoundException
    {
        LinkedHashMap<Integer, Integer> hashDocID = new LinkedHashMap<>();  //hashmap to get all DocID without copies
        ArrayList<Integer> tempList = new ArrayList<>();    // create arrayList to contain temporarily the get DID
        ArrayList<Integer> orderedList = new ArrayList<>(); // the final ordered list of the DIDs
        boolean allPostListScanned = false;     // indicate if all posting list are fully scanned
        boolean postingEnded = false;           // used only in conjunctive query (for optimization), indicates that one posting list is fully scanned
        int max = 0;                            // indicates the current max DID taken
        int currentDocID = 0;                   // var to contain the current DocID
        long startTime,endTime;                 // variables to calculate the execution time

        int[] postingListsIndex = getPostingListsIndex(postingLists);   // contain the current position index for the posting list of each term in the query

        // NEW VERSION -- hash map V.2.2 -- start ------------------------------------------------------------------
        startTime = System.currentTimeMillis();             // start time to take th DocID list
        // take first DocID from posting lists
        for (int i = 0; i < postingLists.length; i++)
        {
            if (postingLists[i] == null)    // term that there isn't in collection -> posting list == null
                continue;                   // go to next posting list
            tempList.add(postingLists[i].get(0).getDocId());        // add the DID into tempList
        }

        Collections.sort(tempList);             // order the list of the first DID
        max = tempList.get(tempList.size()-1);  // take the maximum DID (the last one in the order tempList) and set MAX var
        tempList.clear();                       // clear the tempList
        //printDebug("First max = " + max);

        // scan all posting list and insert DIDs into orderedList
        while(!allPostListScanned) // scan all posting lists and insert DID into priority queue
        {
            allPostListScanned = true;      // set var, if not reset the while end

            // for each posting list I take all values less than max and add them to the hash map.
            for (int i = 0; i < postingLists.length; i++)
            {
                // (term that there isn't in collection -> posting list == null) OR (posting list completely visited)
                if ((postingLists[i] == null) || (postingListsIndex[i] >= postingLists[i].size()))
                {
                    postingEnded = true;    // set var
                    continue;               // go to next posting list
                }

                allPostListScanned = false;  // there is at least one posting not seen yet -> set var
                currentDocID = postingLists[i].get(postingListsIndex[i]).getDocId();    // take current DID
                // scan the current posting list until a DID greater than max (avoiding overflow)
                while((currentDocID <= max) && (postingListsIndex[i] < postingLists[i].size()))
                {
                    if (!hashDocID.containsKey(currentDocID))
                        hashDocID.put(currentDocID,1);  // put DocID in hashtable

                    // update index and currentDID
                    postingListsIndex[i]++;     // increment the related index
                    if (postingListsIndex[i] < postingLists[i].size())
                        currentDocID = postingLists[i].get(postingListsIndex[i]).getDocId();    // take current DID
                }
            }
            //printDebug("------- Step for ------");

            // in hashDocID there are all DID lower than max
            for (Map.Entry<Integer, Integer> entry : hashDocID.entrySet()) {
                tempList.add(entry.getKey());        // insert into orderedList the DIDs taken before
            }
            hashDocID.clear();                  // clear HashMap
            Collections.sort(tempList);         // order the list of DocID
            orderedList.addAll(tempList);       // add ordered DID in orderedList
            tempList.clear();                   // clear templist

            // if the query is conjunctive and at least one posting list is fully scanned
            if (isConjunctive && postingEnded)
                break;      // exit to while. Stop the collection od DID, I have all the DIDs I need

            // take the new max (from the next DID of each posting list)
            for (int i = 0; i < postingLists.length; i++)
            {
                // (term that there isn't in collection -> posting list == null) OR (posting list completely visited)
                if ((postingLists[i] == null) || (postingListsIndex[i] >= postingLists[i].size()))
                    continue;           // go to next posting list

                tempList.add(postingLists[i].get(postingListsIndex[i]).getDocId());        // add the DID into tempList
            }

            // check if thee is only one posting list not fully scanned
            if (tempList.size() == 1)
            {
                //printDebug("Remain only one list");
                if (isConjunctive)
                    break;      // exit to while. Stop the collection od DID, I have all the DIDs I need

                int index = 0;

                for (int i = 0; i < postingLists.length; i++)   // take the index of the not fully scanned posting lists
                {
                    // (term that there isn't in collection -> posting list == null) OR (posting list completely visited)
                    if ((postingLists[i] == null) || (postingListsIndex[i] >= postingLists[i].size()))
                        continue;           // go to next posting list
                    index = i;          // take the index of the not fully scanned posting list
                }
                //printDebug("Remain only the list: " + index + " - total postinglists: " + postingLists.length + " are in the position (of posting list): " + postingListsIndex[index] + " of " + postingLists[index].size());

                while(postingListsIndex[index] < postingLists[index].size())    // take all remaining DID in the posting list
                {
                    orderedList.add(postingLists[index].get(postingListsIndex[index]).getDocId());  // add DID into orederdList
                    postingListsIndex[index]++;     // increment the related index
                }

                allPostListScanned = true;      // set var, if not reset the while end
            }
            else if (tempList.size() != 0)      // there are more than one posting lists not fully scanned (and not is the last iteration of the while with tempList empty)
            {
                Collections.sort(tempList);             // order the list of the first DID
                max = tempList.get(tempList.size()-1);  // take the maximum DID (the last one in the order tempList) and set MAX var
                tempList.clear();                       // clear the tempList
                //printDebug("New max = " + max);
            }
            else    // tempList.size() = 0 -> is the last iteration of the while with tempList empty
            {
                allPostListScanned = true;      // set var, if not reset the while end
            }
        }

        endTime = System.currentTimeMillis();           // end time of DocID list ordering
        // shows query execution time
        printTime("*** ORDERED DID LIST (no PQ V.2.2) in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
        //printDebug("Ordered List (no PQ V.2.2) of DocID dim: " + orderedList.size());     // print orderedList
        // NEW VERSION -- hash map V.2.2 -- end ------------------------------------------------------------------*/
        return orderedList;
    }

    /**
     * function to create an array of indexes for posting lists
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
     */
    private static int booleanSearch(List<Posting> tempList, int targetDID)
    {
        int startPos = 0;
        int endPos = tempList.size()-1;
        int currentPos = 0;

        while (true)
        {
            if (startPos > endPos)
                return -1;          // not found

            currentPos = (startPos + endPos)/2;

            if (tempList.get(currentPos).getDocId() == targetDID)
                return currentPos;
            else if (tempList.get(currentPos).getDocId() > targetDID)
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
        if (rankedResults.size() != 0)      // there are results
        {
            printUI("Query results:");
            for (int i = 0; i < rankedResults.size(); i++)
                printUI((i + 1) + " - " + rankedResults.get(i));
        }
        else                                // there aren't results
            printUI("No results found for this query.");
    }

    /**
     * Function to return the dictionary loaded in memory
     */
    public static HashMap<String, DictionaryElem> getDictionary()
    {
        return dictionary.getTermToTermStat();
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
     * This function allows an automatic test of the resolution of preset queries saved in a file on disk
     * (in this case msmmarco-test2020-queries). This function will fetch and execute (from the specified file) a number
     * of queries passed as a parameter. These queries will be executed in both conjunctive and disjunctive modes,
     * and a series of statistics will be collected about them (total duration of the test, fastest and slowest
     * conjunctive query, fastest and slowest disjunctive query), which will be displayed at the end of the test.
     *
     * @param numQueries    the number of queries to be read from the file and executed
     * @param pathTest      string identified the path for the test queries
     */
    public static void readQueryFromCollection(int numQueries, String pathTest)
    {
        ArrayList<Integer> rankedResults;   // ArrayList that contain the ranked results of query
        int queryCount = 0;             // indicates how many queries have been made
        long startTime, endTime;        // variables to calculate the execution time
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

        // read term upper bound
        if(TermDocUpperBound.termUpperBoundTableIsEmpty())
        {
            if(TermDocUpperBound.termUpperBoundFileExist())     // the file already exist
                TermDocUpperBound.readTermUpperBoundTableFromDisk();
            else                                                // the file not exist
                TermDocUpperBound.calculateTermsUpperBound(false);   // calculate term upper bound for each term of dictionary
        }

        printUIMag(" Start query test... from: " + pathTest);         // control print
        printUIMag("--------------------------------------------------------------------------------");
        File file = new File(pathTest);
        try (
                InputStream tarArchiveInputStream = new GzipCompressorInputStream(new FileInputStream(file));
        ) {
            BufferedReader buffer_collection;
            buffer_collection = new BufferedReader(new InputStreamReader(tarArchiveInputStream, StandardCharsets.UTF_8));
            String record;          // string to contain the queries and their result

            // scan all queries in the collection
            while (((record = buffer_collection.readLine()) != null) && (queryCount < numQueries))
            {
                if (record.isBlank())
                    continue;       // empty string or composed by whitespace characters or malformed

                String[] queryProc = record.split("\t", 2);  // preprocess the query to obtain the result DocNO
                String qid = queryProc[0];      // get the DocNO of the best result for the query

                // print of the query and result obtained by search engine
                printUIMag("---- Query number: " + queryCount + " -------------------------------------------- QueryID: " + qid + " ----");
                printUIMag("The query is: " + queryProc[1]);
                printUIMag("---- disjunctive mode ----");

                startTime = System.currentTimeMillis();         // start time of execute query
                rankedResults = queryManager(queryProc[1],false,5);    // run the query in disjunctive mode
                printQueryResults(rankedResults);
                endTime = System.currentTimeMillis();           // end time of execute query
                // shows query execution time
                printTime("\nQuery (disjunctive mode) executes in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

                // does queries collection statistics
                if ((endTime - startTime) < fasterQueryDis)
                {
                    fasterQueryDis = (endTime - startTime);     // update faster time
                    quidFastDis = qid;                         // update quid
                }
                if ((endTime - startTime) > slowerQueryDis)
                {
                    slowerQueryDis = (endTime - startTime);     // update slower time
                    quidSlowDis = qid;                         //update quid
                }
                avgExTimeDis += (endTime - startTime);          // update avg execution time

                printUIMag("---- conjunctive mode ----");
                startTime = System.currentTimeMillis();         // start time of execute query
                rankedResults = queryManager(queryProc[1],true,5);    // run the query in conjunctive mode
                printQueryResults(rankedResults);
                endTime = System.currentTimeMillis();           // end time of execute query
                // shows query execution time
                printTime("\nQuery (conjunctive mode) executes in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
                printUIMag("--------------------------------------------------------------------------------");

                // does queries collection statistics
                if ((endTime - startTime) < fasterQueryCon)
                {
                    fasterQueryCon = (endTime - startTime);     // update faster time
                    quidFastCon = qid;                         // update quid
                }
                if ((endTime - startTime) > slowerQueryCon)
                {
                    slowerQueryCon = (endTime - startTime);     // update slower time
                    quidSlowCon = qid;                         // update quid
                }
                avgExTimeCon += (endTime - startTime);          // update avg execution time

                queryCount++;       // update counter
            }

            // print queries collection statistics
            printUIMag(" End query test... from: " + pathTest);         // control print
            printUIMag("--------------------------------------------------------------------------------");
            printTime("The fastest query (conjunctive mode) executes in " + fasterQueryCon + " ms and its QUID is " + quidFastCon);
            printTime("The slowest query (conjunctive mode) executes in " + slowerQueryCon + " ms and its QUID is " + quidSlowCon);
            printTime("The average queries execution time (conjunctive mode) is " + avgExTimeCon/numQueries + " ms");

            printTime("\nThe fastest query (disjunctive mode) executes in " + fasterQueryDis + " ms and its QUID is " + quidFastDis);
            printTime("The slowest query (disjunctive mode) executes in " + slowerQueryDis + " ms and its QUID is " + quidSlowDis);
            printTime("The average queries execution time (disjunctive mode) is " + avgExTimeDis/numQueries + " ms");
            printUIMag("--------------------------------------------------------------------------------");
        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }
    // -------- end: function to read collection of query --------
}

/*
 * NOTE:
 * 0 - list of conditions (of the if):
 *      partialScore = 0 -> also with the maximum score in the nonessential posting lists the document score will be
 *                          less than the minimum threshold to be among the top docs.
 *      firstEssPostListIndex = 0 -> all posting list have already been scanned
 *      resetScore = true -> conjunctive case, the DID must be in all posting lists
 *
 * 1 -
 */
