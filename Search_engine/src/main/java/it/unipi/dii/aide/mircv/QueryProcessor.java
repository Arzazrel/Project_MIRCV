package it.unipi.dii.aide.mircv;
import it.unipi.dii.aide.mircv.compression.Unary;
import it.unipi.dii.aide.mircv.data_structures.*;
import it.unipi.dii.aide.mircv.data_structures.Dictionary;
import it.unipi.dii.aide.mircv.utils.FileSystem;
import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.apache.commons.compress.compressors.gzip.GzipCompressorInputStream;

import java.io.*;
import java.nio.channels.FileChannel;
import java.nio.charset.StandardCharsets;
import java.util.*;

import static it.unipi.dii.aide.mircv.data_structures.CollectionStatistics.readCollectionStatsFromDisk;
import static it.unipi.dii.aide.mircv.data_structures.DataStructureHandler.readCompressedPostingListFromDisk;
import static it.unipi.dii.aide.mircv.data_structures.DataStructureHandler.readPostingListFromDisk;
import static it.unipi.dii.aide.mircv.data_structures.Flags.readFlagsFromDisk;
//import static it.unipi.dii.aide.mircv.data_structures.PartialIndexBuilder.dictionaryBlockOffsets;
import static it.unipi.dii.aide.mircv.utils.Constants.*;
import static it.unipi.dii.aide.mircv.utils.Constants.printDebug;

/**
 * Class to manage and execute query
 */
public final class QueryProcessor {

    // indicate whether order all or only first "numberOfResults" results from hash table. TEST VARIABLE
    private static boolean orderAllHashMap = false;
    public static HashMap<Integer, DocumentElement> documentTable = new HashMap<>();    // hash table DocID to related DocElement
    static it.unipi.dii.aide.mircv.data_structures.Dictionary dictionary = new Dictionary();    // dictionary in memory
    private static ArrayList<String> termNotInCollection = new ArrayList<>();   // ArrayList that contain the term that are in the query but not in the collection
    // variable to BM25
    private static double k = 1.2;      // typical values between 1,2 and 2
    private static double b = 0.75;     // typical values around 0.75
    private static double avgDocLen;     // average document length
    static PriorityQueue<QueryProcessor.PostingBlock> pq;       // priority queue for ordering the Docs in posting lists by DocID (ascending order)
    static PriorityQueue<QueryProcessor.ResultBlock> resPQ;     // priority queue for the result of scoring function for the best numberOfResults docs

    /**
     * fuction to manage the query request. Prepare and execute the query and return the results.
     *
     * @param query             is the query of the users (in words)
     * @param isConjunctive     indicates whether the query is of conjunctive type
     * @param isDisjunctive     indicates whether the query is of disjunctive type
     * @param numberOfResults   the number of results to be returned by the query
     * @return  an ArrayList of integer that representing an ordered list of DocIDs
     */
    public static ArrayList<Integer> queryManager(String query, boolean isConjunctive, boolean isDisjunctive, int numberOfResults)
    {
        ArrayList<Integer> rankedResults = new ArrayList<>();   // ArrayList that contain the ranked results of query
        ArrayList<String> processedQuery;                       // array list for containing the query term

        // take user's choices that affecting the query execution
        boolean scoringFunc = Flags.isScoringEnabled();      // take user's choice about using scoring function
        //boolean scoringFunc = true;      // to debug use, to activate BM25

        //printDebug("User choice for scoring is: " + scoringFunc);
        if (scoringFunc)
            avgDocLen = CollectionStatistics.getTotDocLen() / CollectionStatistics.getNDocs();  // set average doc len

        try{
            // processed the query to obtain the term
            printDebug("Query before processed: " + query);
            processedQuery = TextProcessor.preprocessText(query); // Preprocessing of document text
            printDebug("Query after processed: " + processedQuery);
            //printDebug(dictionary.getTermStat("0000").toString());

            // check if query is empty
            if (processedQuery.isEmpty() || (processedQuery.size() == 1 && processedQuery.get(0).equals("")))
            {
                printError("Error: query is empty, please retry.");     // mex of error
                return rankedResults;
            }

            // control for correct form
            if ( (isConjunctive && isDisjunctive) || !(isConjunctive || isDisjunctive))     // query is Conjunctive or Disjunctive cannot be both or neither
            {
                printError("Error: query is Conjunctive or Disjunctive cannot be both or neither.");  // mex of error
                return rankedResults;
            }

            DAATAlgorithm(processedQuery,scoringFunc,isConjunctive, isDisjunctive,numberOfResults);        // apply DAAT to calculate the score of the Docs
            //resPQ.clear();
            //DAATAlgorithmV1(processedQuery,scoringFunc,isConjunctive, isDisjunctive,numberOfResults);        // apply DAAT V.1 to calculate the score of the Docs

            rankedResults = getRankedResults(numberOfResults);          // get ranked results

            // check if array list is empty
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
    public static boolean queryStartControl() throws IOException {
        // -- control for file into disk
        if (!FileSystem.areThereAllMergedFiles() ||
                !Flags.isThereFlagsFile() ||
                !CollectionStatistics.isThereStatsFile()) {
            printError("Error: missing required files.");
            return false;
        }

        readFlagsFromDisk();
        readCollectionStatsFromDisk();

        // -- control for structures in memory - if not load them from disk
        if (!dictionary.dictionaryIsSet())
        {
            long startTime = System.currentTimeMillis();
            dictionary.readDictionaryFromDisk();
            long endTime = System.currentTimeMillis();
            printTime( "Dictionary loaded in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
        }
        if(documentTable.isEmpty())
        {
            long startTime = System.currentTimeMillis();
            DataStructureHandler.readDocumentTableFromDisk(false);
            long endTime = System.currentTimeMillis();
            printTime("Document Table loaded in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
        }

        return true;
    }

    /**
     * function for apply the Document at a Time algorithm
     *
     * @param scoringFunc       indicates the preference for scoring. if false use TFIDF, if true use BM25.
     * @param processedQuery    array list for containing the query term
     * @param isConjunctive     indicates whether the query is of conjunctive type
     * @param isDisjunctive     indicates whether the query is of disjunctive type
     * @param numberOfResults   indicated the max number of result to return to user
     */
    private static void DAATAlgorithm(ArrayList<String> processedQuery, boolean scoringFunc , boolean isConjunctive, boolean isDisjunctive, int numberOfResults) throws FileNotFoundException
    {
        resPQ = new PriorityQueue<>(numberOfResults, new CompareTerm());    // length equal to the number of results to be returned to the user
        int docScoreCalc = 0;               // indicates the number of documents whose score was calculated (0 to number of results requested by the user)
        boolean resetScore = false;         // used only in conjunctive case. indicates that the score must be set to 0 (the current Doc there aren't all the term of the query)
        int pLNotEmpty = 0;                 // contains the number of posting lists related to the query terms that aren't empty
        int df = 0;                         // contains df of the term (used in score function)
        double min = 0;

        ArrayList<Integer> ordListDID;      // ordered list of the DocID present in the all posting lists of the term present in the query
        ArrayList<Posting>[] postingLists;  // contains all the posting lists for each term of the query
        Posting currentP;                   // support var
        int[] postingListsIndex;            // contain the current position index for the posting list of each term in the query
        int currentDID = 0;                 // DID of the current doc processed in algorithm
        double partialScore = 0;            // var that contain partial score
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
        postingListsIndex = getPostingListsIndex(postingLists);     // get the index initialized   NEW VERSION

        startTime = System.currentTimeMillis();           // start time of DAAT
        // scan all Doc retrieved and calculate score (TFIDF or BM25)
        for (Integer integer : ordListDID) {

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
                    String term = processedQuery.get(j);
                    assert term != null;
                    df = dictionary.getTermToTermStat().get(term).getDf();
                    if (scoringFunc)
                        partialScore += ScoringBM25(currentDID,currentP.getTermFreq(), df);     // use BM25
                    else
                        partialScore += ScoringTFIDF(currentP.getTermFreq(), df);               // use TFIDF
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
                        printTime("*** DAAT V.0.5 execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
                        return;             // exit from function
                    }
                }
            }
            //printDebug("---- for hop ----");

            // save score
            if ((partialScore != 0) && !resetScore)
            {
                // insert without control into priority queue (is not full) or insert all results (orderAllHashMap = true)
                if ((docScoreCalc < numberOfResults) || orderAllHashMap)
                {
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                    //min = partialScore;
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
                /*else if (min < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                {
                    // substitution of the block
                    //printDebug("Old block : DID = " + resPQ.peek().getDID()+ " score: " + resPQ.peek().getScore());
                    resPQ.poll();       // remove the first element
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                    //printDebug("New block : DID = " + currentDID+ " score: " + partialScore);
                    min = partialScore;
                }*/

            }
        }
        endTime = System.currentTimeMillis();           // end time of DAAT
        // shows DAAT execution time
        printTime("*** DAAT V.0.5 execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
    }

    /**
     * function that given the posting lists of each term in a given query returns the number of posting lists that
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
     * function for apply the Document at a Time algorithm in the special case of only one posting list not empty
     *
     * @param postingLists      array list for containing the posting lists of the query's term
     * @param scoringFunc       indicates the preference for scoring. If false use TFIDF, if true use BM25.
     * @param processedQuery    array list for containing the query term
     * @param numberOfResults   indicated the max number of result to return to user
     */
    private static void DAATOnePostingList (ArrayList<String> processedQuery, ArrayList<Posting>[] postingLists, boolean scoringFunc, int numberOfResults)
    {
        double partialScore = 0;    // var that contain partial score
        long startTime,endTime;     // variables to calculate the execution time
        int docScoreCalc = 0;       // indicates the number of documents whose score was calculated (0 to number of results requested by the user)
        int df = 0;                 // contains df of the term (used in score function)

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
            String term = processedQuery.get(index);        // take the related term
            assert term != null;
            df = dictionary.getTermToTermStat().get(term).getDf();              // retrieve df for the term
            if (scoringFunc)
                partialScore += ScoringBM25(p.getDocId(), p.getTermFreq(), df); // use BM25
            else
                partialScore += ScoringTFIDF(p.getTermFreq(), df);              // use TFIDF

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

    ///* ---- NEW VERSION OF DAAT -- START ----
    private static void DAATAlgorithmV1(ArrayList<String> processedQuery, boolean scoringFunc , boolean isConjunctive, boolean isDisjunctive, int numberOfResults) throws FileNotFoundException
    {
        resPQ = new PriorityQueue<>(numberOfResults, new CompareTerm());    // length equal to the number of results to be returned to the user
        ArrayList<Integer> tempList = new ArrayList<>();    // create arrayList to contain temporarily the get DID
        ArrayList<Posting>[] postingLists;  // contains all the posting lists for each term of the query
        int[] postingListsIndex ;           // contain the current position index for the posting list of each term in the query
        double partialScore = 0;            // var that contain partial score
        int currentDID = 0;                 // DID of the current doc processed in algorithm
        int df = 0;                         // contains df of the term (used in score function)
        int min = 0;                        // indicates the current min DID taken
        int docScoreCalc = 0;               // indicates the number of documents whose score was calculated (0 to number of results requested by the user)
        boolean resetScore = false;         // used only in conjunctive case. indicates that the score must be set to 0 (the current Doc there aren't all the term of the query)
        boolean allPostListScanned = false; // indicate if all posting list are fully scanned
        boolean postingEnded = false;       // used only in conjunctive query (for optimization), indicates that one posting list is fully scanned
        long startTime,endTime;             // variables to calculate the execution time

        startTime = System.currentTimeMillis();         // start time for retrieve all posting lists of the query
        postingLists = retrieveAllPostListsFromQuery(processedQuery);   // take all posting lists of query terms
        endTime = System.currentTimeMillis();           // end time for retrieve all posting lists of the query
        // shows query execution time
        printTime("\n*** Retrieved all posting lists in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

        // -- DAAT V.1.5 -- start --
        startTime = System.currentTimeMillis();         // start time of DAAT

        int postingListNotEmpty = 0;    // contains the number of posting lists related to the query terms that aren't empty
        // control check for empty posting lists (the terms are not present in the document collection)
        if (postingLists.length == 0)
        {
            printUI("The terms in query there aren't in collection.");
            return;     // exit to function
        }
        // scan all posting lists
        for (int i = 0; i < postingLists.length; i++)
        {
            if (postingLists[i] == null)    // term that there isn't in collection -> posting list == null
                continue;                   // go to next posting list
            postingListNotEmpty++;          // increment counter
        }

        if ((postingListNotEmpty == 1) && ( !isConjunctive || (postingListNotEmpty == postingLists.length)))   // only one posting list to scan
        {
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
                String term = processedQuery.get(index);        // take the related term
                assert term != null;
                df = dictionary.getTermToTermStat().get(term).getDf();              // retrieve df for the term
                if (scoringFunc)
                    partialScore += ScoringBM25(p.getDocId(), p.getTermFreq(), df); // use BM25
                else
                    partialScore += ScoringTFIDF(p.getTermFreq(), df);              // use TFIDF

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
            printTime("*** DAAT V.1.5 execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
            // -- VDAAT V.1.5 -- end --
            return;
        }
        else if ( isConjunctive && (postingListNotEmpty != postingLists.length))    // case of conjunctive query and at least one query term is not in dictionary
        {
            endTime = System.currentTimeMillis();           // end time of DAAT
            // shows DAAT execution time
            printTime("*** DAAT V.1.5 execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
            // -- VDAAT V.1.5 -- end --
            return;
        }

        if (!isConjunctive)         // disjunctive case
        {
                //printDebug("More than one list.");
                postingListsIndex = getPostingListsIndex(postingLists); // initialize the index for the posting lists

                for (int i = 0; i < postingLists.length; i++)           // take first DocID from posting lists
                {
                    if (postingLists[i] == null)    // term that there isn't in Dictionary -> posting list == null
                        continue;                   // go to next posting list

                    tempList.add(postingLists[i].get(0).getDocId());        // add the DID into tempList
                }
                //printDebug("First tempList = " + tempList);
                Collections.sort(tempList);             // order the list of the first DID
                //printDebug("First ordered tempList = " + tempList);
                min = tempList.get(0);  // take the minimum DID (the first one in the order tempList) and set MIN var
                tempList.clear();                       // clear the tempList
                //printDebug("First min = " + min);

                // scan all posting lists and perform DAAT and calculate the score for each Doc
                while(!allPostListScanned)
                {   // start - while 0 (execution DAAT algoritm for more posting lists)
                    allPostListScanned = true;      // set var, if not reset the while end
                    partialScore = 0;           // reset var
                    // for each posting list I calculate the score for the Doc with DID equal to the current MIN
                    for (int i = 0; i < postingLists.length; i++)
                    {   // start - for 1 (scan posting list to calculate score)

                        // (term that there isn't in collection -> posting list == null) OR (posting list completely visited)
                        if ((postingLists[i] == null) || (postingListsIndex[i] >= postingLists[i].size()))
                        {
                            //printDebug("Setto postingEnded");
                            continue;                   // go to next posting list
                        }
                        allPostListScanned = false;  // there is at least one posting not seen yet -> set var

                        currentDID = postingLists[i].get(postingListsIndex[i]).getDocId();    // take current DID
                        if (currentDID == min)    // Same DID of the min, is the DID for this iteration of DAAT
                        {
                            // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                            String term = processedQuery.get(i);        // take the related term
                            assert term != null;
                            df = dictionary.getTermToTermStat().get(term).getDf();  // retrieve df for the term
                            if (scoringFunc)
                                partialScore += ScoringBM25(currentDID, postingLists[i].get(postingListsIndex[i]).getTermFreq(), df);   // use BM25
                            else
                                partialScore += ScoringTFIDF(postingLists[i].get(postingListsIndex[i]).getTermFreq(), df);              // use TFIDF

                            postingListsIndex[i]++;     // increment the related index
                        }
                    }   // end - for 1 (scan posting list to calculate score)
                    //printDebug("------- Step for ------");

                    // save score
                    if (partialScore != 0)
                    {
                        // insert without control into priority queue (is not full) or insert all results (orderAllHashMap = true)
                        if ((docScoreCalc < numberOfResults) || orderAllHashMap)
                        {
                            resPQ.add(new QueryProcessor.ResultBlock(min, partialScore));     // add to priority queue
                            docScoreCalc++;         // increment result in priority queue counter
                        }
                        else if (resPQ.peek().getScore() < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                        {
                            // substitution of the block
                            //printDebug("Old block : DID = " + resPQ.peek().getDID()+ " score: " + resPQ.peek().getScore());
                            resPQ.poll();       // remove the first element
                            resPQ.add(new QueryProcessor.ResultBlock(min, partialScore));     // add to priority queue
                            //printDebug("New block : DID = " + currentDID+ " score: " + partialScore);
                        }
                    }

                    // take the new min (from the next DID of each posting list)
                    for (int i = 0; i < postingLists.length; i++)
                    {
                        // (term that there isn't in collection -> posting list == null) OR (posting list completely visited)
                        if ((postingLists[i] == null) || (postingListsIndex[i] >= postingLists[i].size()))
                            continue;           // go to next posting list

                        tempList.add(postingLists[i].get(postingListsIndex[i]).getDocId());        // add the DID into tempList
                    }

                    // check if thee is only one posting list not fully scanned
                    if (tempList.size() == 1)
                    {   // start - if 0.1 (remain one posting list)
                        //printDebug("Remain only one list");
                        int index = 0;
                        for (int i = 0; i < postingLists.length; i++)   // take the index of the not fully scanned posting lists
                        {
                            // (term that there isn't in collection -> posting list == null) OR (posting list completely visited)
                            if ((postingLists[i] == null) || (postingListsIndex[i] >= postingLists[i].size()))
                                continue;           // go to next posting list
                            index = i;          // take the index of the not fully scanned posting list
                        }
                        //printDebug("Remain only the list: " + index + " - total postinglists: " + postingLists.length + " are in the position (of posting list): " + postingListsIndex[index] + " of " + postingLists[index].size());
                        while(postingListsIndex[index] < postingLists[index].size())    // take all remaining posting in the posting list
                        {   // start - while 0.1 (calculate score when remain only one posting list)
                            partialScore = 0;           // reset var
                            // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                            String term = processedQuery.get(index);        // take the related term
                            assert term != null;
                            df = dictionary.getTermToTermStat().get(term).getDf();  // retrieve df for the term
                            currentDID = postingLists[index].get(postingListsIndex[index]).getDocId();  // take current DID
                            if (scoringFunc)
                                partialScore += ScoringBM25(currentDID, postingLists[index].get(postingListsIndex[index]).getTermFreq(), df);   // use BM25
                            else
                                partialScore += ScoringTFIDF(postingLists[index].get(postingListsIndex[index]).getTermFreq(), df);  // use TFIDF

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

                            postingListsIndex[index]++;     // increment the related index
                        }   // end - while 0.1 (calculate score when remain only one posting list)

                        allPostListScanned = true;      // set var, if not reset the while end
                    }   // end - if 0.1 (remain one posting list)
                    else if (tempList.size() != 0)      // there are more than one posting lists not fully scanned (and not is the last iteration of the while with tempList empty)
                    {
                        Collections.sort(tempList);             // order the list of the first DID
                        min = tempList.get(0);  // take the maximum DID (the last one in the order tempList) and set MAX var
                        tempList.clear();                       // clear the tempList
                        //printDebug("New min = " + min);
                    }
                    else    // tempList.size() = 0 -> is the last iteration of the while with tempList empty
                    {
                        allPostListScanned = true;      // set var, if not reset the while end
                        //printDebug("End");
                    }
                }   // end - while 0 (execution DAAT algoritm for more posting lists)
        }
        else if (postingLists.length == postingListNotEmpty)    // cunjunctive case
        {
                //printDebug("More than one list.");
                postingListsIndex = getPostingListsIndex(postingLists); // initialize the index for the posting lists

                for (int i = 0; i < postingLists.length; i++)           // take first DocID from posting lists
                {
                    if (postingLists[i] == null)    // term that there isn't in Dictionary -> posting list == null
                            return;
                    tempList.add(postingLists[i].get(0).getDocId());        // add the DID into tempList
                }
                //printDebug("First tempList = " + tempList);
                Collections.sort(tempList);             // order the list of the first DID
                //printDebug("First ordered tempList = " + tempList);
                min = tempList.get(0);  // take the minimum DID (the first one in the order tempList) and set MIN var
                tempList.clear();                       // clear the tempList
                //printDebug("First min = " + min);

                // scan all posting lists and perform DAAT and calculate the score for each Doc
                while(!allPostListScanned)
                {   // start - while 0 (execution DAAT algoritm for more posting lists)
                    allPostListScanned = true;      // set var, if not reset the while end
                    partialScore = 0;           // reset var
                    resetScore = false;          // reset var
                    // for each posting list I calculate the score for the Doc with DID equal to the current MIN
                    for (int i = 0; i < postingLists.length; i++)
                    {   // start - for 1 (scan posting list to calculate score)
                        // (term that there isn't in collection -> posting list == null) OR (posting list completely visited)
                        if ((postingLists[i] == null) || (postingListsIndex[i] >= postingLists[i].size()))
                        {
                            //printDebug("Setto postingEnded");
                            postingEnded = true;        // set var
                            resetScore = true;          // reset the partial score
                            continue;                   // go to next posting list
                        }
                        allPostListScanned = false;  // there is at least one posting not seen yet -> set var

                        currentDID = postingLists[i].get(postingListsIndex[i]).getDocId();    // take current DID
                        if (currentDID == min)    // Same DID of the min, is the DID for this iteration of DAAT
                        {
                            // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                            String term = processedQuery.get(i);        // take the related term
                            assert term != null;
                            df = dictionary.getTermToTermStat().get(term).getDf();  // retrieve df for the term
                            if (scoringFunc)
                                partialScore += ScoringBM25(currentDID, postingLists[i].get(postingListsIndex[i]).getTermFreq(), df);   // use BM25
                            else
                                partialScore += ScoringTFIDF(postingLists[i].get(postingListsIndex[i]).getTermFreq(), df);              // use TFIDF

                            postingListsIndex[i]++;     // increment the related index
                        }
                        else  // at least one posting list doesn't have this doc hence this doc does not contain at least one term from the query
                        {
                            resetScore = true;          // reset the partial score
                        }
                    }   // end - for 1 (scan posting list to calculate score)
                    //printDebug("------- Step for ------");

                    // save score
                    if ((partialScore != 0) && !resetScore)
                    {
                        // insert without control into priority queue (is not full) or insert all results (orderAllHashMap = true)
                        if ((docScoreCalc < numberOfResults) || orderAllHashMap)
                        {
                            resPQ.add(new QueryProcessor.ResultBlock(min, partialScore));     // add to priority queue
                            docScoreCalc++;         // increment result in priority queue counter
                        }
                        else if (resPQ.peek().getScore() < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                        {
                            // substitution of the block
                            //printDebug("Old block : DID = " + resPQ.peek().getDID()+ " score: " + resPQ.peek().getScore());
                            resPQ.poll();       // remove the first element
                            resPQ.add(new QueryProcessor.ResultBlock(min, partialScore));     // add to priority queue
                            //printDebug("New block : DID = " + currentDID+ " score: " + partialScore);
                        }
                    }

                    // if the query is conjunctive and at least one posting list is fully scanned
                    if (postingEnded)
                    {
                        //printDebug("query conjunctive e postingEnded settato esco");
                        break;      // exit to while. Stop posting list scanning and score calculating. I have the score of all the DIDs I need
                    }

                    // take the new min (from the next DID of each posting list)
                    for (int i = 0; i < postingLists.length; i++)
                    {
                        // (term that there isn't in collection -> posting list == null) OR (posting list completely visited)
                        if ((postingLists[i] == null) || (postingListsIndex[i] >= postingLists[i].size()))
                            continue;           // go to next posting list

                        tempList.add(postingLists[i].get(postingListsIndex[i]).getDocId());        // add the DID into tempList
                    }

                    // check if thee is only one posting list not fully scanned
                    if (tempList.size() == 1)
                    {   // start - if 0.1 (remain one posting list)
                        break;      // exit to while. Stop the collection od DID, I have all the DIDs I need
                    }   // end - if 0.1 (remain one posting list)
                    else if (tempList.size() != 0)      // there are more than one posting lists not fully scanned (and not is the last iteration of the while with tempList empty)
                    {
                        Collections.sort(tempList);             // order the list of the first DID
                        min = tempList.get(0);  // take the maximum DID (the last one in the order tempList) and set MAX var
                        tempList.clear();                       // clear the tempList
                        //printDebug("New min = " + min);
                    }
                    else    // tempList.size() = 0 -> is the last iteration of the while with tempList empty
                    {
                        allPostListScanned = true;      // set var, if not reset the while end
                        //printDebug("End");
                    }
                }   // end - while 0 (execution DAAT algoritm for more posting lists)
        }

        endTime = System.currentTimeMillis();           // end time of DAAT
        // shows DAAT execution time
        printTime("*** DAAT V.1.5 execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
        // -- VDAAT V.1.5 -- end --

        /* -- DAAT V.1 -- start --
        startTime = System.currentTimeMillis();         // start time of DAAT
        if (postingLists.length == 1)       // there is only one posting (query with one term)
        {   // start - if 0
            //printDebug("Only one list.");
            if (postingLists[0] == null)    // term that there isn't in collection -> posting list == null
                return;
            // optimization -> there is one posting list -> the DID are already sort
            for (Posting p : postingLists[0])
            {
                partialScore = 0;           // reset var
                // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                String term = processedQuery.get(0);        // take the related term
                assert term != null;
                df = dictionary.getTermToTermStat().get(term).getDf();  // retrieve df for the term
                if (scoringFunc)
                    partialScore += ScoringBM25(p.getDocId(), p.getTermFreq(), df);   // use BM25
                else
                    partialScore += ScoringTFIDF(p.getTermFreq(), df);              // use TFIDF

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
        }   // end - if 0
        else    // there are more posting lists (query composed by more term)
        {   // start - else 0 (more postinglists)
            //printDebug("More than one list.");
            postingListsIndex = getPostingListsIndex(postingLists); // initialize the index for the posting lists

            for (int i = 0; i < postingLists.length; i++)           // take first DocID from posting lists
            {
                if (postingLists[i] == null)    // term that there isn't in Dictionary -> posting list == null
                {
                    if (isConjunctive)   // conjunctive query case, if on posting lists is empty it's impossible that any doc can have all the query terms
                        return;
                    continue;                   // go to next posting list
                }
                tempList.add(postingLists[i].get(0).getDocId());        // add the DID into tempList
            }
            //printDebug("First tempList = " + tempList);
            Collections.sort(tempList);             // order the list of the first DID
            //printDebug("First ordered tempList = " + tempList);
            min = tempList.get(0);  // take the minimum DID (the first one in the order tempList) and set MIN var
            tempList.clear();                       // clear the tempList
            //printDebug("First min = " + min);

            // scan all posting lists and perform DAAT and calculate the score for each Doc
            while(!allPostListScanned)
            {   // start - while 0 (execution DAAT algoritm for more posting lists)
                allPostListScanned = true;      // set var, if not reset the while end
                partialScore = 0;           // reset var
                resetScore = false;          // reset var
                // for each posting list I calculate the score for the Doc with DID equal to the current MIN
                for (int i = 0; i < postingLists.length; i++)
                {   // start - for 1 (scan posting list to calculate score)

                    // (term that there isn't in collection -> posting list == null) OR (posting list completely visited)
                    if ((postingLists[i] == null) || (postingListsIndex[i] >= postingLists[i].size()))
                    {
                        //printDebug("Setto postingEnded");
                        postingEnded = true;        // set var
                        if (isConjunctive)
                            resetScore = true;      // reset the partial score
                        continue;                   // go to next posting list
                    }

                    allPostListScanned = false;  // there is at least one posting not seen yet -> set var

                    currentDID = postingLists[i].get(postingListsIndex[i]).getDocId();    // take current DID
                    if (currentDID == min)    // Same DID of the min, is the DID for this iteration of DAAT
                    {
                        // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                        String term = processedQuery.get(i);        // take the related term
                        assert term != null;
                        df = dictionary.getTermToTermStat().get(term).getDf();  // retrieve df for the term
                        if (scoringFunc)
                            partialScore += ScoringBM25(currentDID, postingLists[i].get(postingListsIndex[i]).getTermFreq(), df);   // use BM25
                        else
                            partialScore += ScoringTFIDF(postingLists[i].get(postingListsIndex[i]).getTermFreq(), df);              // use TFIDF

                        postingListsIndex[i]++;     // increment the related index
                    }
                    else if(isConjunctive)  // at least one posting list doesn't have this doc hence this doc does not contain at least one term from the query
                    {
                        resetScore = true;          // reset the partial score
                    }
                }   // end - for 1 (scan posting list to calculate score)
                //printDebug("------- Step for ------");

                // save score
                if ((partialScore != 0) && !resetScore)
                {
                    // insert without control into priority queue (is not full) or insert all results (orderAllHashMap = true)
                    if ((docScoreCalc < numberOfResults) || orderAllHashMap)
                    {
                        resPQ.add(new QueryProcessor.ResultBlock(min, partialScore));     // add to priority queue
                        docScoreCalc++;         // increment result in priority queue counter
                    }
                    else if (resPQ.peek().getScore() < partialScore)    // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                    {
                        // substitution of the block
                        //printDebug("Old block : DID = " + resPQ.peek().getDID()+ " score: " + resPQ.peek().getScore());
                        resPQ.poll();       // remove the first element
                        resPQ.add(new QueryProcessor.ResultBlock(min, partialScore));     // add to priority queue
                        //printDebug("New block : DID = " + currentDID+ " score: " + partialScore);
                    }
                }

                // if the query is conjunctive and at least one posting list is fully scanned
                if (isConjunctive && postingEnded)
                {
                    //printDebug("query conjunctive e postingEnded settato esco");
                    break;      // exit to while. Stop posting list scanning and score calculating. I have the score of all the DIDs I need
                }

                // take the new min (from the next DID of each posting list)
                for (int i = 0; i < postingLists.length; i++)
                {
                    // (term that there isn't in collection -> posting list == null) OR (posting list completely visited)
                    if ((postingLists[i] == null) || (postingListsIndex[i] >= postingLists[i].size()))
                        continue;           // go to next posting list

                    tempList.add(postingLists[i].get(postingListsIndex[i]).getDocId());        // add the DID into tempList
                }

                // check if thee is only one posting list not fully scanned
                if (tempList.size() == 1)
                {   // start - if 0.1 (remain one posting list)
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
                    while(postingListsIndex[index] < postingLists[index].size())    // take all remaining posting in the posting list
                    {   // start - while 0.1 (calculate score when remain only one posting list)
                        partialScore = 0;           // reset var
                        // calculate SCORE (TFIDF or BM25) for this term and currentDID and sum to partial score
                        String term = processedQuery.get(index);        // take the related term
                        assert term != null;
                        df = dictionary.getTermToTermStat().get(term).getDf();  // retrieve df for the term
                        currentDID = postingLists[index].get(postingListsIndex[index]).getDocId();  // take current DID
                        if (scoringFunc)
                            partialScore += ScoringBM25(currentDID, postingLists[index].get(postingListsIndex[index]).getTermFreq(), df);   // use BM25
                        else
                            partialScore += ScoringTFIDF(postingLists[index].get(postingListsIndex[index]).getTermFreq(), df);  // use TFIDF

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

                        postingListsIndex[index]++;     // increment the related index
                    }   // end - while 0.1 (calculate score when remain only one posting list)

                    allPostListScanned = true;      // set var, if not reset the while end
                }   // end - if 0.1 (remain one posting list)
                else if (tempList.size() != 0)      // there are more than one posting lists not fully scanned (and not is the last iteration of the while with tempList empty)
                {
                    Collections.sort(tempList);             // order the list of the first DID
                    min = tempList.get(0);  // take the maximum DID (the last one in the order tempList) and set MAX var
                    tempList.clear();                       // clear the tempList
                    //printDebug("New min = " + min);
                }
                else    // tempList.size() = 0 -> is the last iteration of the while with tempList empty
                {
                    allPostListScanned = true;      // set var, if not reset the while end
                    //printDebug("End");
                }
            }   // end - while 0 (execution DAAT algoritm for more posting lists)
        }   // start - else 0 (more postinglists)

        endTime = System.currentTimeMillis();           // end time of DAAT
        // shows DAAT execution time
        printTime("*** DAAT V.1 execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
        // -- DAAT V.1 -- end -- */
    }
    ///* ---- NEW VERSION OF DAAT -- END ----

    

    // -------- start: scoring function --------

    /**
     * function to calculate TFIDF for one term and one document
     *
     * @param termFreq          term frequency of the term in the document
     * @param postListLength    number of documents in which the term occurs
     * @return  the TFIDF score for one term and one document. The total score for a document will be the sum of the
     *          result of this function for each term that is both in the document and in the query
     */
    private static Double ScoringTFIDF(int termFreq, int postListLength)
    {
        double TFweight, IDFweight, scoreTFIDF;     // variables to calculate the TFIDF score value

        // control to avoid log and division to 0
        if (termFreq == 0 || postListLength == 0)
            return (double) 0;

        TFweight = (1 + Math.log10(termFreq));      // calculate TF weight
        IDFweight = Math.log10(((double) CollectionStatistics.getNDocs() / postListLength));    // calculate IDF weight
        scoreTFIDF = TFweight * IDFweight;          // calculate TFIDF weight from Tf and IDF weight values

        //printDebug("ScoringTFIDF - TFweight = " + TFweight + " IDFweight = " + IDFweight + " scoreTFIDF = " + scoreTFIDF);

        return scoreTFIDF;
    }

    /**
     * function to calculate BM25 for one term and one document
     *
     * @param DocID             DocID of the document processed
     * @param termFreq          term frequency of the term in the document
     * @param postListLength    number of documents in which the term occurs
     * @return  the BM25 score for one term and one document. The total score for a document will be the sum of the
     *          result of this function for each term that is both in the document and in the query
     */
    private static Double ScoringBM25(int DocID, int termFreq, int postListLength)
    {
        double docLen, denominator, IDFweight, scoreBM25;     // variables to calculate the BM25 score value

        // control to avoid log and division to 0
        if (termFreq == 0 || postListLength == 0)
            return (double) 0;

        docLen = documentTable.get(DocID).getDoclength();   // get doc length

        denominator = k * ((1 - b) + b * (docLen / avgDocLen)) + termFreq;                      // calculate TF weight
        IDFweight = Math.log10(((double) CollectionStatistics.getNDocs() / postListLength));    // calculate IDF weight
        scoreBM25 = (termFreq / denominator) * IDFweight;      // calculate TFIDF weight from Tf and IDF weight values

        //printDebug("ScoringBM25 - docLen = " + docLen + " denominator = " + denominator + " IDFweight = " + IDFweight + " scoreBM25 = " + scoreBM25);

        return scoreBM25;
    }

    // -------- end: scoring function --------

    // -------- start: utilities function --------

    /**
     * function to elaborate all docs and related scores to obtain the ranked list of results
     *
     * @param numResults    number of result(DocID) to return to the user
     * @return  an ordered ArrayList that represent the top numResults results for the query
     */
    private static ArrayList<Integer> getRankedResults(int numResults)
    {
        ArrayList<Integer> rankedResults = new ArrayList<>();   // array list to contain the top "numResults" docs
        long startTime, endTime;            // variables to calculate the execution time

        //control check
        if (numResults <= 0)
            return rankedResults;

        startTime = System.currentTimeMillis();         // start time of hash map ordering

        QueryProcessor.ResultBlock currentResPQ;        // var that contain the resultBlock extract from pq in the current iteration
        ArrayList<Integer> results = new ArrayList<Integer>();       // Create an ArrayList object

        while(!resPQ.isEmpty())                         // control if the priority queue for results is empty
        {
            currentResPQ = resPQ.poll();                            // take lowest element (score and DID)
            results.add(currentResPQ.getDID());                     // add to the array list
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
     * function to retrieve all the posting lists for each term of the query passed as parameter
     *
     * @param processedQuery    ArrayList of the processed ter of the query
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
     * function that given the posting lists of each term in a given query returns an ordered list of the DocIDs
     * present in the all posting lists
     *
     * @param postingLists  the posting lists of each term in the query
     * @return  an ordered ArrayList of the DocIDs in the posting lists
     */
    private static ArrayList<Integer> DIDOrderedListOfQuery(ArrayList<Posting>[] postingLists, boolean isConjunctive) throws FileNotFoundException {

        // ordered list of the DocID present in the all posting lists of the term present in the query
        ArrayList<Integer> orderedList = new ArrayList<>();
        LinkedHashMap<Integer, Integer> hashDocID = new LinkedHashMap<>();  //hashmap to get all DocID without copies
        long startTime,endTime;                         // variables to calculate the execution time
        int currentDocID = 0;                           // var to contain the current DocID

        //printDebug("The isConjunctive passed as parameter has this value: " + isConjunctive);
        /* print posting lists
        for (int i = 0; i < postingLists.length; i++)
        {
            if (postingLists[i] == null)    // term that there isn't in collection -> posting list == null
                continue;                   // go to next posting list

            printDebug("PostingLists: " + postingLists[i]);
        }
        //*/

        /* print how many DID there are in more than one posting list
        int count = 0;
        for (int i = 0; i < postingLists.length; i++)
        {
            if (postingLists[i] == null)    // term that there isn't in collection -> posting list == null
                continue;                   // go to next posting list
            // scan all DocID in the i-th posting list
            for (Posting p : postingLists[i])
            {
                currentDocID = p.getDocId();            // take DocID in the current posting
                //System.out.println(ANSI_YELLOW + "\n*** Posting: " + i + " - DID: " + currentDocID);
                if (i == 0)
                    hashDocID.put(currentDocID,1);      // add DocID
                    // control check for duplicate DocID, do only after first posting list
                else if (hashDocID.containsKey(currentDocID))
                {
                    printDebug("Termine a comune: " + currentDocID);
                    count++;
                }
            }
        }
        printDebug("Termini presenti sia nella prima lista che nella seconda: " + count);
        hashDocID.clear();
        //*/

        ///* NEW VERSION -- hash map V.2.2 -- start ------------------------------------------------------------------
        int[] postingListsIndex = getPostingListsIndex(postingLists);   // contain the current position index for the posting list of each term in the query
        ArrayList<Integer> tempList = new ArrayList<>();    // create arrayList to contain temporarily the get DID
        boolean allPostListScanned = false; // indicate if all posting list are fully scanned
        boolean postingEnded = false;       // used only in conjunctive query (for optimization), indicates that one posting list is fully scanned
        int max = 0;                        // indicates the current max DID taken

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
        //System.out.println("Ordered List (no PQ V.2.2): " + orderedList);     // print orderedList
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
     * function to shows the user the ranked results (DocID) of the query executed
     *
     * @param rankedResults the results returned by the query
     */
    public static void printQueryResults(ArrayList<Integer> rankedResults)
    {
        if (rankedResults.size() != 0)      // there are results
        {
            System.out.println(ANSI_CYAN + "Query results:" + ANSI_RESET);
            for (int i = 0; i < rankedResults.size(); i++)
                printUI((i + 1) + " - " + rankedResults.get(i));
        }
        else                                // there aren't results
            printUI("No results found for this query.");
    }
    // -------- end: utilities function --------


    // -------- start: utilities for priority queue --------
    /**
     * class to define PostingBlock. The priority queue contains instances of PostingBlock
     */
    private static class PostingBlock {
        int DocID;                  // DocID
        int indexOfPostingList;     // reference to the posting list (index in the array of posting lists of the query) containing DcoID

        // constructor with parameters
        public PostingBlock(int DocID, int indexOfPostingList) {
            this.DocID = DocID;
            this.indexOfPostingList = indexOfPostingList;
        }

        public int getDID() {
            return DocID;
        }

        public int getIndex() {
            return indexOfPostingList;
        }

        @Override
        public String toString() {
            return "PB{" +
                    "DocID = '" + DocID + '\'' +
                    ", index of pl =" + indexOfPostingList +
                    '}';
        }
    }
    /**
     * class to compare the block, allows the order of the priority queue
     */
    private static class ComparePostingBlockTerm implements Comparator<QueryProcessor.PostingBlock> {
        @Override
        public int compare(QueryProcessor.PostingBlock pb1, QueryProcessor.PostingBlock pb2) {
            // comparing terms
            int DocIDComparison = Integer.compare(pb1.getDID(), pb2.getDID());
            // if the DocID are equal, compare by block number
            if (DocIDComparison == 0) {
                return Integer.compare(pb1.getIndex(), pb2.getIndex());
            }

            return DocIDComparison;
        }
    }

    /**
     * class to define PostingBlock. The priority queue contains instances of PostingBlock
     */
    private static class DIDBlock {
        int DocID;                  // DocID

        // constructor with parameters
        public DIDBlock(int DocID) {
            this.DocID = DocID;
        }

        public int getDID() {
            return DocID;
        }

        @Override
        public String toString() {
            return "PB{" + "DocID = '" + DocID + '\'' + '}';
        }
    }
    /**
     * class to compare the block, allows the order of the priority queue
     */
    private static class CompareDIDBlock implements Comparator<QueryProcessor.DIDBlock> {
        @Override
        public int compare(QueryProcessor.DIDBlock pb1, QueryProcessor.DIDBlock pb2) {
            // comparing terms
            return Integer.compare(pb1.getDID(), pb2.getDID());
        }
    }
    // -------- end: utilities for priority queue --------

    // -------- start: utilities for priority queue for the results --------
    /**
     * class to define ResultBlock. The priority queue contains instances of ResultBlock representing a result of the score calculation for the research
     */
    private static class ResultBlock {
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
     * class to compare the block, allows the order of the priority queue
     * The order is ascending order according to doc score in case of equal score will be made order according to DocID.
     * DocIDs sorted in descending order. In this way in case of a tie, documents with the smallest DocID will be
     * considered better than those with the largest DocID.
     */
    private static class CompareTerm implements Comparator<QueryProcessor.ResultBlock> {
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
     * @param numQueries the number of queries to be read from the file and executed
     */
    public static void readQueryFromCollection(int numQueries)
    {
        ArrayList<Integer> rankedResults;   // ArrayList that contain the ranked results of query
        int queryCount = 0;             // indicates how many queries have been made
        long startTime, endTime;        // variables to calculate the execution time
        long fasterQueryCon = 100000;   // indicates the execution time for the fastest query in the collection (conjunctive)
        String quidFastCon = "";        // indicates the QUID of the fastest query in the collection (conjunctive)
        long slowerQueryCon = 0;        // indicates the execution time for the slowest query in the collection (conjunctive)
        String quidSlowCon = "";        // indicates the QUID of the slowest query in the collection (conjunctive)
        long fasterQueryDis = 100000;   // indicates the execution time for the fastest query in the collection (disjunctive)
        String quidFastDis = "";        // indicates the QUID of the fastest query in the collection (disjunctive)
        long slowerQueryDis = 0;        // indicates the execution time for the slowest query in the collection (disjunctive)
        String quidSlowDis = "";        // indicates the QUID of the slowest query in the collection (disjunctive)
        long avgExTimeCon = 0;          // indicate the average execution time for the queries in the collection (conjunctive)
        long avgExTimeDis = 0;          // indicate the average execution time for the queries in the collection (disjunctive)

        // control check for queries
        try {
            if (!queryStartControl()) {
                return;                 // there aren't all files needed for execute a query, function it's terminated
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        printDebug(" Start query test...");         // control print
        File file = new File(QUERIES_COLLECTION_PATH); // file that contain the queries
        try (
                FileReader fr = new FileReader(file);  //reads the file
        ) {
            BufferedReader buffer_collection;
            buffer_collection = new BufferedReader(fr);
            String record;          // string to contain the queries and their result

            // scan all queries in the collection
            while (((record = buffer_collection.readLine()) != null) && (queryCount < numQueries)) {

                if (record.isBlank()) {
                    continue;       // empty string or composed by whitespace characters or malformed
                }
                String[] queryProc = record.split("\t", 2);  // preprocess the query to obtain the result DocNO
                String quid = queryProc[0];      // get the DocNO of the best result for the query

                // print of the query and result obtained by search engine
                printDebug("---- Query number: " + queryCount + " -------------------------------------------- QUID: " + quid + " ----");
                printDebug("---- disjunctive mode ----");

                startTime = System.currentTimeMillis();         // start time of execute query
                rankedResults = queryManager(queryProc[1],false,true,5);    // run the query in disjunctive mode
                printQueryResults(rankedResults);
                endTime = System.currentTimeMillis();           // end time of execute query
                // shows query execution time
                printTime("\nQuery (disjunctive mode) executes in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");

                // does queries collection statistics
                if ((endTime - startTime) < fasterQueryDis)
                {
                    fasterQueryDis = (endTime - startTime);     // update faster time
                    quidFastDis = quid;                         // update quid
                }
                if ((endTime - startTime) > slowerQueryDis)
                {
                    slowerQueryDis = (endTime - startTime);     // update slower time
                    quidSlowDis = quid;                         //update quid
                }
                avgExTimeDis += (endTime - startTime);          // update avg execution time

                printDebug("---- conjunctive mode ----");
                startTime = System.currentTimeMillis();         // start time of execute query
                rankedResults = queryManager(queryProc[1],true,false,5);    // run the query in conjunctive mode
                printQueryResults(rankedResults);
                endTime = System.currentTimeMillis();           // end time of execute query
                // shows query execution time
                printTime("\nQuery (conjunctive mode) executes in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
                printDebug("--------------------------------------------------------------------------------");

                // does queries collection statistics
                if ((endTime - startTime) < fasterQueryCon)
                {
                    fasterQueryCon = (endTime - startTime);     // update faster time
                    quidFastCon = quid;                         // update quid
                }
                if ((endTime - startTime) > slowerQueryCon)
                {
                    slowerQueryCon = (endTime - startTime);     // update slower time
                    quidSlowCon = quid;                         // update quid
                }
                avgExTimeCon += (endTime - startTime);          // update avg execution time

                queryCount++;       // update counter
            }

            // print queries collection statistics
            printTime("The fastet query (conjunctive mode) executes in " + fasterQueryCon + " ms and its QUID is " + quidFastCon);
            printTime("The slowest query (conjunctive mode) executes in " + slowerQueryCon + " ms and its QUID is " + quidSlowCon);
            printTime("The average queries execution time (conjunctive mode) is " + avgExTimeCon/numQueries + " ms");

            printTime("\nThe fastet query (disjunctive mode) executes in " + fasterQueryDis + " ms and its QUID is " + quidFastDis);
            printTime("The slowest query (disjunctive mode) executes in " + slowerQueryDis + " ms and its QUID is " + quidSlowDis);
            printTime("The average queries execution time (disjunctive mode) is " + avgExTimeDis/numQueries + " ms");

            printDebug(" End query test...");         // control print
            printDebug("--------------------------------------------------------------------------------");
        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }
    // -------- end: function to read collection of query --------
}

/*
 * NOTE:
 * 0 -
 */
