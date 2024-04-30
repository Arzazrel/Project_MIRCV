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
    // HashMap for containing the DocID and document score related. DID -> doc score
    private static final HashMap<Integer, Double> tableDAAT = new HashMap<>();
    // HashMap for containing the top "numberOfResults" DocID of the document related to score value. doc score -> ArrayList of DID
    private static final HashMap<Double, ArrayList<Integer>> scoreToDocID = new HashMap<>();
    // HashMap for containing the score values for which "numberOfResults" docs have already been found. Score -> true or false
    private static final HashMap<Double, Boolean> scoreWithMaxDoc = new HashMap<>();
    public static HashMap<Integer, DocumentElement> documentTable = new HashMap<>();    // hash table DocID to related DocElement
    static it.unipi.dii.aide.mircv.data_structures.Dictionary dictionary = new Dictionary();    // dictionary in memory
    private static ArrayList<String> termNotInCollection = new ArrayList<>();   // ArrayList that contain the term that are in the query but not in the collection
    // variable to BM25
    private static double k = 1.2;      // typical values between 1,2 and 2
    private static double b = 0.75;     // typical values around 0.75
    private static double avgDocLen;     // average document length

    //static PriorityQueue<QueryProcessor.PostingBlock> pq;
    static PriorityQueue<QueryProcessor.ResultBlock> resPQ;     // priority queue for the result of scoring function for the best numberOfResults docs
    private static boolean scoringFunc;

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

        printDebug("User choice for scoring is: " + scoringFunc);
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

            DAATAlgorithm(processedQuery,scoringFunc,isConjunctive, isDisjunctive,numberOfResults);        // apply DAAT, result in tableDAAT

            rankedResults = getRankedResults(numberOfResults);          // get ranked results
            tableDAAT.clear();                                          // clear HashMap
            scoreToDocID.clear();                                       // clear HashMap
            scoreWithMaxDoc.clear();                                    // clear HashMa
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
    private static void DAATAlgorithm(ArrayList<String> processedQuery, boolean scoringFunc , boolean isConjunctive, boolean isDisjunctive, int numberOfResults) throws FileNotFoundException {

        ArrayList<Integer> ordListDID;          // ordered list of the DocID present in the all posting lists of the term present in the query
        ArrayList<Posting>[] postingLists;      // contains all the posting lists for each term of the query
        Posting currentP;                       // support var
        int currentDID = 0;                     // DID of the current doc processed in algorithm
        int[] postingListsIndex;                // contain the current position index for the posting list of each term in the query
        double partialScore = 0;                // var that contain partial score
        long startTime,endTime;                 // variables to calculate the execution time
        // new
        resPQ = new PriorityQueue<>(numberOfResults, new CompareTerm());    // length equal to the number of results to be returned to the user
        int docScoreCalc = 0;                   // indicates the number of documents whose score was calculated (0 to number of results requested by the user)

        startTime = System.currentTimeMillis();           // end time of hash map ordering
        postingLists = retrieveAllPostListsFromQuery(processedQuery);   // take all posting lists of query terms
        endTime = System.currentTimeMillis();           // end time of hash map ordering
        // shows query execution time
        System.out.println(ANSI_YELLOW + "\n*** Retrieved all posting lists in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")" + ANSI_RESET);

        // control check for empty posting lists (the terms are not present in the document collection)
        if (postingLists.length == 0)
        {
            printUI("The term in query there aren't in collection.");
            return;     // exit to function
        }

        ordListDID = DIDOrderedListOfQuery(postingLists);           // take ordered list of DocID
        postingListsIndex = getPostingListsIndex(postingLists);     // get the index initialized   NEW VERSION

        startTime = System.currentTimeMillis();           // end time of hash map ordering

        // scan all Doc retrieved and calculate score (TFIDF or BM25)
        for (Integer integer : ordListDID) {

            currentDID = integer;     // update the DID, document of which to calculate the score
            partialScore = 0;                   // reset var

            // default case is query Disjunctive
            // take all values and calculating the scores in the posting related to currentDID
            for (int j = 0; j < postingLists.length; j++) {
                // check if the posting lists of j-th isn't at the end AND if the j-th term of the query is present in the doc identify by currentDID
                if ( (postingLists[j] != null) && (postingListsIndex[j] < postingLists[j].size()) && (postingLists[j].get(postingListsIndex[j]).getDocId() == currentDID)) {
                    currentP = postingLists[j].get(postingListsIndex[j]);              // take posting
                    postingListsIndex[j]++;                         // update index of current value

                    //System.out.println("DAAT, prescoring -- df = " + DataStructureHandler.postingListLengthFromTerm(processedQuery.get(j)));

                    // calculate TFIDF for this term and currentDID and sum to partial score
                    String term = processedQuery.get(j);
                    assert term != null;
                    int df = dictionary.getTermToTermStat().get(term).getDf();
                    if (scoringFunc)
                        partialScore += ScoringBM25(currentDID,currentP.getTermFreq(), df);     // use BM25
                    else
                        partialScore += ScoringTFIDF(currentP.getTermFreq(), df);               // use TFIDF

//                    printDebug("DAAT: posting del termine: " + processedQuery.get(j) + " in array pos: " + j + " ha DID: " + currentDID + " and partialScore: " + partialScore);
                } else if (isConjunctive) {
                    // must take only the document in which there are all term (DID that compare in all posting lists of the terms)
                    partialScore = 0;       // reset the partial score
                    // if all postings in one posting lists have already been seen the next documents in the posting lists cannot contain all the terms in the query
                    if ((postingListsIndex[j] >= postingLists[j].size()) || (postingLists[j] == null)) {
                        endTime = System.currentTimeMillis();           // end time of hash map ordering
                        // shows query execution time
                        System.out.println(ANSI_YELLOW + "\n*** DAAT execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")" + ANSI_RESET);
                        return;             // exit from function
                    } else
                        break;              // exit from the for and go to next Document
                }
            }

            // save score
            if (partialScore != 0) {
                /*
                // -- old version: not priority queue for the result --
                //tableDAAT.put(currentDID,partialScore);     // add DID and related score to HashMap   OLD VERSION
                if (!scoreWithMaxDoc.containsKey(partialScore)) {
                    tableDAAT.put(currentDID, partialScore);     // add DID and related score to HashMap     NEW VERSION
                    addToScoreToDocID(partialScore, currentDID, numberOfResults); // add DID to the related DID in hashmap
                }
                //printDebug("Final TFIDF scoring for DID = " + currentDID + " is: " + tableDAAT.get(currentDID));
                // -- END -- old version: not priority queue for the result --
                */

                // -- START -- new version: priority queue for the result --
                if (docScoreCalc < numberOfResults)     // insert without control into priority queue
                {
                    resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                    docScoreCalc++;         // increment
                }
                else        // number of user-requested results achieved, check whether the current doc is within the best docs to return (score greater than the first item in the priority queue)
                {
                    if (resPQ.peek().getScore() < partialScore) // substitution of the block
                    {
                        //printDebug("Old block : DID = " + resPQ.peek().getDID()+ " score: " + resPQ.peek().getScore());
                        resPQ.poll();       // remove the first element
                        resPQ.add(new QueryProcessor.ResultBlock(currentDID, partialScore));     // add to priority queue
                        //printDebug("New block : DID = " + currentDID+ " score: " + partialScore);
                    }
                }
                // -- END -- new version: priority queue for the result --
            }
        }

        endTime = System.currentTimeMillis();           // end time of hash map ordering
        // shows query execution time
        System.out.println(ANSI_YELLOW + "\n*** DAAT execute in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")" + ANSI_RESET);
    }

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
        ArrayList<Double> orderedList = new ArrayList<>();      // contain scores of all docs
        long startTime, endTime;            // variables to calculate the execution time

        //control check
        if (numResults < 0 /*|| tableDAAT.isEmpty()*/)
            return rankedResults;

        // ---- START -- new version with priority queue ----
        startTime = System.currentTimeMillis();         // start time of hash map ordering

        QueryProcessor.ResultBlock currentresPQ;     // var that contain the resultBlock extract from pq in the current iteration
        ArrayList<Integer> results = new ArrayList<Integer>();  // Create an ArrayList object

        while(!resPQ.isEmpty())
        {
            currentresPQ = resPQ.poll();                            // take lowest element (score and DID)
            results.add(currentresPQ.getDID());                     // add to the array list
        }
        // order the result from the best to the worst (reverse order of the priority queue)
        rankedResults = new ArrayList<Integer>(results);     // Create an ArrayList object
        Collections.reverse(rankedResults);
        printDebug("orderedResults " + rankedResults);

        endTime = System.currentTimeMillis();           // end time of hash map ordering
        System.out.println(ANSI_YELLOW + "\n*** Ranked results (results priority queue) in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")" + ANSI_RESET);
        // ---- END -- new version with priority queue ----

        /*
        // take ranked list of DocID
        for (Map.Entry<Integer, Double> entry : tableDAAT.entrySet()) {
            orderedList.add(entry.getValue());
        }

        System.out.println("\n*** REMOVE DUPLICATES orderedLISt - before size: " + orderedList.size());
        startTime = System.currentTimeMillis();         // start time of hash map ordering
        Set<Double> set = new HashSet<>(orderedList);
        orderedList.clear();
        orderedList.addAll(set);
        endTime = System.currentTimeMillis();           // end time of hash map ordering
        System.out.println(ANSI_YELLOW + "\n*** REMOVE DUPLICATE orderedList in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")" + ANSI_RESET);
        System.out.println("\n*** REMOVE DUPLICATES orderedLISt - after size: " + orderedList.size());

        startTime = System.currentTimeMillis();         // start time of hash map ordering
        orderedList.sort(Collections.reverseOrder());
        endTime = System.currentTimeMillis();           // end time of hash map ordering
        System.out.println(ANSI_YELLOW + "\n*** ORDER orderedList in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")" + ANSI_RESET);

        // true in testing phase -> order and show all results (required long time)
        if (orderAllHashMap)
        {
            startTime = System.currentTimeMillis();         // start time of hash map ordering
            for (double num : orderedList) {
                for (Map.Entry<Integer, Double> entry : tableDAAT.entrySet()) {
                    if (entry.getValue() == num && !rankedResults.contains(entry.getKey())) {
                        rankedResults.add(entry.getKey());
                    }
                }
            }

            printDebug("Total ranked results: " + rankedResults);

            // if the ranked results are more than numResults, cut the last results
            if (rankedResults.size() > numResults)
            {
                List<Integer> ord = rankedResults.subList(0,numResults);    // retrieve only the first numResults DocID
                rankedResults = new ArrayList<>(ord);
                System.out.println("Cut ranked results: " + rankedResults);
            }
            endTime = System.currentTimeMillis();           // end time of hash map ordering
            // shows query execution time
            printTime("\n*** TOTAL HashMap ordered in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")");
        }
        else
        {
            int iterator = 1;                   // iterator to stop the ordering
            ArrayList<Integer> retrievDocIDs;   // list to get the arraylist of DocID related to a score

            startTime = System.currentTimeMillis();         // start time of hash map ordering
            //remove duplicate
            for (double num : orderedList)
            {
                if (scoreToDocID.containsKey(num))
                {
                    // take DocID of the document that have num as score
                    retrievDocIDs = scoreToDocID.get(num);
                    scoreToDocID.remove(num);
                    //System.out.println("\n*** retrieveDocIDs size: " + retrievDocIDs.size() + "\nretrieveDocIDs array list: " + retrievDocIDs);
                    // scan all DocID retrieved
                    for (Integer i : retrievDocIDs)
                    {
                        rankedResults.add(i);
                        if (iterator < numResults)
                            iterator++;
                        else
                        {
                            endTime = System.currentTimeMillis();           // end time of hash map ordering
                            // shows query execution time
                            System.out.println(ANSI_YELLOW + "\n*** PARTIAL HashMap ordered in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")" + ANSI_RESET);
                            return rankedResults;
                        }
                    }
                }
                else
                {
                    System.out.println(ANSI_RED + "ERROR in scoreToDocID hashMap there isn't a Doc for the score: " + num + ANSI_RESET);
                }
            }
        }
        */

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
//                printDebug("DAAT: retrieve posting list of  " + term);

                DictionaryElem de = dictionary.getTermToTermStat().get(term);

                // there is a posting list for the query term == the term is in the collection
                if (dictionary.getTermToTermStat().containsKey(term))
                {
                    if(Flags.isCompressionEnabled())
                        postingLists[iterator] = readCompressedPostingListFromDisk(de.getOffsetDocId(),de.getOffsetTermFreq(), de.getTermFreqSize(), de.getDocIdSize(), de.getDf(), docIdChannel, termFreqChannel); //read compressed posting list
                    else // take the postingList of term
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
    private static ArrayList<Integer> DIDOrderedListOfQuery(ArrayList<Posting>[] postingLists) throws FileNotFoundException {

        // ordered list of the DocID present in the all posting lists of the term present in the query
        ArrayList<Integer> orderedList = new ArrayList<>();
        LinkedHashMap<Integer, Integer> hashDocID = new LinkedHashMap<>();  //hashmap to get all DocID without copies
        long startTime,endTime;                         // variables to calculate the execution time
        int currentDocID = 0;                           // var to contain the current DocID

        ///* OLD VERSION -- start

        startTime = System.currentTimeMillis();         // start time to take th DocID list
        // scan all posting lists passed as parameters
        for (int i = 0; i < postingLists.length; i++)
        {
            if (postingLists[i] == null)    // term that there isn't in collection -> posting list == null
                continue;                   // go to next posting list
            // scan all DocID in the i-th posting list
            for (Posting p : postingLists[i])
            {
                currentDocID = p.getDocId();            // take DocID in the current posting
                if (i == 0)
                    hashDocID.put(currentDocID,1);      // add DocID
                // control check for duplicate DocID, do only after first posting list
                else if (!hashDocID.containsKey(currentDocID))
                {
                    hashDocID.put(currentDocID,1);      // add DocID
                }
            }
        }
        for (Map.Entry<Integer, Integer> entry : hashDocID.entrySet()) {
            orderedList.add(entry.getKey());
        }
        endTime = System.currentTimeMillis();          // end time to take th DocID list
        System.out.println(ANSI_YELLOW + "\n*** TAKE DID LIST (no PQ) in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")" + ANSI_RESET);

        startTime = System.currentTimeMillis();         // start time of DocID list ordering
        Collections.sort(orderedList);          // order the list of DocID
        endTime = System.currentTimeMillis();           // end time of DocID list ordering
        // shows query execution time
        System.out.println(ANSI_YELLOW + "\n*** ORDERED DID LIST (no PQ) in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")" + ANSI_RESET);

        //printDebug("Ordered List of DocID for the query:  " + orderedList);     // print orderedList

        System.out.println("Ordered List (no PQ) of DocID dim: " + orderedList.size());     // print orderedList
        hashDocID.clear();          // clear linkHashMap
        // OLD VERSION - end */

        /* NEW VERSION (priority queue) - start
        // clear the previus work without PQ
        orderedList = new ArrayList<>();
        hashDocID.clear();

        // create PQ
        //PriorityQueue<QueryProcessor.PostingBlock> pq = new PriorityQueue<>(postingLists.length, new CompareTerm());
        pq = new PriorityQueue<>(postingLists.length, new CompareTerm());
        int[] postingListsIndex = getPostingListsIndex(postingLists); // contain the current position index for the posting list of each term in the query
        QueryProcessor.PostingBlock currentPostingBlock;     // var that contain the PostBlock extract from pq in the current iteration

        startTime = System.currentTimeMillis();         // start time to take th DocID list
        // take first DocID from posting lists
        for (int i = 0; i < postingLists.length; i++)
        {
            if (postingLists[i] == null)    // term that there isn't in collection -> posting list == null
                continue;                   // go to next posting list
            pq.add(new QueryProcessor.PostingBlock(postingLists[i].get(0).getDocId(), i));     // add to the priority queue a TermBlock element (term + its blocks number)
        }

        while(!pq.isEmpty()) //
        {
            //qSystem.out.println("PQ:\n" + pq);               // print priority queue
            currentPostingBlock = pq.poll();                // take lowest element (DID and index)
            currentDocID = currentPostingBlock.getDID();
            hashDocID.put(currentDocID,1);  // put DocID in hashtable
            // scann all current position in posting list and update indexes
            for (int i = 0; i < postingLists.length; i++)
            {
                // check if the DocID in the posting list is the same in currentPostingBlock and check if there is another posting in the posting lists
                if ( (postingLists[i] != null) && (postingListsIndex[i] < postingLists[i].size()) && (currentDocID == postingLists[i].get(postingListsIndex[i]).getDocId()) )
                {
                    // check if is the posting list of the currentPostingBlock
                    if ( i != currentPostingBlock.getIndex())
                    {
                        pq.poll();                  // remove one element in pq
                    }
                    postingListsIndex[i]++;     // update index
                    // check if there is another posting in the posting lists
                    if (postingListsIndex[i] < postingLists[i].size())
                        pq.add(new QueryProcessor.PostingBlock(postingLists[i].get(postingListsIndex[i]).getDocId(), i));  // insert new posting in pq
                }
            }
        }
        // pass from hashMap to ArrayList
        for (Map.Entry<Integer, Integer> entry : hashDocID.entrySet()) {
            orderedList.add(entry.getKey());
        }
        endTime = System.currentTimeMillis();           // end time of DocID list ordering
        // shows query execution time
        System.out.println(ANSI_YELLOW + "\n*** TAKE AND ORDERED DID LIST (PQ) in " + (endTime - startTime) + " ms (" + formatTime(startTime, endTime) + ")" + ANSI_RESET);

        System.out.println("Ordered List (PQ) of DocID dim: " + orderedList.size());     // print orderedList
        // NEW VERSION (priority queue) - end */

        return orderedList;
    }

    /**
     * function to put in hashMap the document related its score
     *
     * @param score             score of the document
     * @param DocID             DocID of the document
     * @param numberOfResults   maximum number of document to keep for each score
     */
    private static void addToScoreToDocID (double score, int DocID,int numberOfResults)
    {
        ArrayList<Integer> values;

        if (scoreToDocID.containsKey(score))        // contains key, add DocID to the arrayList
        {
            if (scoreToDocID.get(score).size() >= numberOfResults)      // SEE NOTE 0
            {
                scoreWithMaxDoc.put(score,true);                        // SEE NOTE 1
                return;
            }
            // get value from HashMap
            ArrayList<Integer> oldValues = scoreToDocID.get(score);
            //System.out.println("*** addToScoreToDocID: old value = " + scoreToDocID.get(score) + " for score: " + score);
            values = new ArrayList<>(oldValues.subList(0,oldValues.size()));
            values.add(DocID);                  // add current DID
            scoreToDocID.put(score,values);     // update to hashMap
            //System.out.println("*** addToScoreToDocID: new value = " + scoreToDocID.get(score) + " for score: " + score);
        }
        else        // First DocID create an arrayList with one element
        {
            values = new ArrayList<>();
            values.add(DocID);
            scoreToDocID.put(score,values);     // add to hashMap
            //System.out.println("*** addToScoreToDocID: new value = " + scoreToDocID.get(score) + " for score: " + score);
        }
    }

    // new version to substitute remove with get in the posting lists

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
    /*
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
    */
    /**
     * class to compare the block, allows the order of the priority queue
     */
    /*
    private static class CompareTerm implements Comparator<QueryProcessor.PostingBlock> {
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
    */
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

                startTime = System.currentTimeMillis();         // start time of execute query
                queryManager(queryProc[1],true,false,5);    // run the query in conjunctive mode
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
            printTime("\nThe fastet query (conjunctive mode) executes in " + fasterQueryCon + " ms and its QUID is " + quidFastCon);
            printTime("\nThe slowest query (conjunctive mode) executes in " + slowerQueryCon + " ms and its QUID is " + quidSlowCon);
            printTime("\nThe average queries execution time (conjunctive mode) is " + avgExTimeCon/numQueries + " ms");

            printTime("\nThe fastet query (disjunctive mode) executes in " + fasterQueryDis + " ms and its QUID is " + quidFastDis);
            printTime("\nThe slowest query (disjunctive mode) executes in " + slowerQueryDis + " ms and its QUID is " + quidSlowDis);
            printTime("\nThe average queries execution time (disjunctive mode) is " + avgExTimeDis/numQueries + " ms");

            printDebug(" End query test...");         // control print
        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }
    // -------- end: function to read collection of query --------
}

/*
 * NOTE:
 * 0 - The idea behind the reasoning is that if I have to return to the user the "numberOfResults" documents with the
 *     best result, if for each result I collect at least the first "numberOfResults" documents with that result and
 *     don't register the others at the end I will have given the same best "numberOfResults" results as if I had
 *     registered all the documents.
 * 1 - the maximum required number of documents (numberOfResults) have been scored, advise don't take any more
 *     documents with this score
 */
