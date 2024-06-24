package it.unipi.dii.aide.mircv.data_structures;

import it.unipi.dii.aide.mircv.TextProcessor;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import static it.unipi.dii.aide.mircv.utils.Constants.*;
import static it.unipi.dii.aide.mircv.data_structures.DataStructureHandler.*;

import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.apache.commons.compress.compressors.gzip.GzipCompressorInputStream;

/**
 * Class to make the partial inverted index, dictionary and document table.
 */
public final class PartialIndexBuilder
{
    // Data structures initialization
    static HashMap<Integer, DocumentElement> documentTable = new HashMap<>();     // hash table DocID to related DocElement
    static Dictionary dictionary = new Dictionary();                              // dictionary in memory
    static HashMap<String, ArrayList<Posting>> invertedIndex = new HashMap<>();   // hash table Term to related Posting list
    static ArrayList<Long> dictionaryBlockOffsets = new ArrayList<>();            // Offsets of the dictionary blocks

    /**
     * Implements the SPIMI algorithm for indexing large collections.
     */
    public static void SPIMIalgorithm()
    {
        long memoryAvailable = (long) (Runtime.getRuntime().maxMemory() * MEMORY_THRESHOLD);    // amount of memory which can be used
        int docCounter = 1;         // counter for DocID
        int termCounter = 0;        // counter for TermID
        int totDocLen = 0;          // variable for the sum of the lengths of all documents
        double avgDocLen;           // average document length (used in BM25 scoring function)
        int emptyDocs = 0;          // number of empty docs in the collection
        int minlenDoc = 100000;     // len of the shortest doc in the collection
        int maxLenDoc = 0;          // len of the longest doc in the collection
        int maxTermFreq = 0;        // max termFreq in the collection
        int tempCurrTF = 0;         //
        int malformedDocs = 0;      //

        printDebug("The memory available for each block is: " + memoryAvailable + " bites ( " + (double)((double)memoryAvailable/(double)1073741824) + " GB)");

        File file = new File(COLLECTION_PATH);
        try (
            final TarArchiveInputStream tarArchiveInputStream = new TarArchiveInputStream(new GzipCompressorInputStream(new FileInputStream(file)));
        ) {
            TarArchiveEntry tarArchiveEntry = tarArchiveInputStream.getNextTarEntry();
            BufferedReader buffer_collection;
            if(tarArchiveEntry == null)     // empty collection
                return;
            buffer_collection = new BufferedReader(new InputStreamReader(tarArchiveInputStream, StandardCharsets.UTF_8));
            String record;                  // string to contain the document

            // scan all documents in the collection
            while ((record = buffer_collection.readLine()) != null)
            {   // -- start - while 0 - scan all docs --
                int separator = record.indexOf("\t");       // check for malformed line, no \t
                if (record.isBlank() || separator == -1)    // empty string or composed by whitespace characters or malformed
                {
                    malformedDocs++;
                    continue;       // malformed doc -> go to the next doc
                }

                ArrayList<String> preprocessed = TextProcessor.preprocessText(record); // Preprocessing of document text
                String docno = preprocessed.remove(0);      // get the DocNO of the current document

                // check if document is empty
                if (preprocessed.isEmpty() || (preprocessed.size() == 1 && preprocessed.get(0).equals("")))
                {
                    emptyDocs++;
                    continue;              // skip to next while iteration (next document)
                }

                // to collect collection statistics
                if (preprocessed.size() < minlenDoc)
                {
                    minlenDoc = preprocessed.size();
                }
                if (preprocessed.size() > maxLenDoc)
                {
                    maxLenDoc = preprocessed.size();
                }

                DocumentElement de = new DocumentElement(docno, docCounter, preprocessed.size());   // create new Document element
                documentTable.put(docCounter, de);      // add current Document into Document Table in memory
                totDocLen += preprocessed.size();       // update total doc length, add current document length (value will be stored in collection statistics)
                // scan all term in the current document
                for (String term : preprocessed)
                {   // -- start - for 0 - scan all term in current doc --
                    // control check if the length of the current term is greater than the maximum allowed
                    if(term.length() > TERM_DIM)
                        term = term.substring(0,TERM_DIM);      // truncate term

                    // control check if the current term has already been found or is the first time
                    if (!dictionary.getTermToTermStat().containsKey(term))
                        termCounter++;                          // update TermID counter

                    assert !term.equals("");
                    DictionaryElem dictElem = dictionary.getOrCreateTerm(term,termCounter);     // Dictionary build

                    // check if the term is already find in this doc or it is the first time
                    if(addTerm(term, docCounter, 0))
                        dictElem.addDf(1);  // update document frequency (number of docs in which there is the term)
                    dictElem.addCf(1);  // update collection frequency (number of occurrences of the term in the collection)

                    N_POSTINGS++;       // update number of partial postings to save in the file

                    tempCurrTF = invertedIndex.get(term).get(invertedIndex.get(term).size() - 1).getTermFreq();
                    if (tempCurrTF > maxTermFreq)
                        maxTermFreq = tempCurrTF;
                }   // -- end - for 0 - scan all term in current doc --
                docCounter++;       // update DocID counter

                if(Runtime.getRuntime().totalMemory() > memoryAvailable)
                {
                    printDebug("********** Memory full **********");
                    storeIndexAndDictionaryIntoDisk();  //store index and dictionary to disk
                    storeDocumentTableIntoDisk();       // store document table one document at a time for each block
                    freeMemory();   // delete information in document table, dictionary,and inverted index
                    System.gc();    // effort JVM
                    printDebug("********** Free memory **********");
                    N_POSTINGS = 0;     // new partial index, reset number of postings in the block
                }
            }   // -- end - while 0 - scan all docs --
            DataStructureHandler.storeBlockOffsetsIntoDisk();   // store into file all the blocks offset

            printDebug("Malformed docs: " + malformedDocs);
            printDebug("The number of empty docs is: " + emptyDocs + " the shortest doc have len of: " + minlenDoc + " the longest doc have len of: " + maxLenDoc + " the max TermFreq is: " + maxTermFreq);
            // calculate and store collection statistics values
            avgDocLen = (double) totDocLen / docCounter;    // set average doc len

            CollectionStatistics.setNDocs(docCounter);          // set total number of Document in the collection
            CollectionStatistics.setTotDocLen(totDocLen);       // set the sum of the all document length in the collection
            CollectionStatistics.setAvgDocLen(avgDocLen);       // set the Doc average len
            CollectionStatistics.setEmptyDocs(emptyDocs);       // set the number of empty docs in the collection
            CollectionStatistics.setMinLenDoc(minlenDoc);       // set the len of the shortest doc in the collection
            CollectionStatistics.setMaxLenDoc(maxLenDoc);       // set the len of the longest doc in the collection
            CollectionStatistics.setMaxTermFreq(maxTermFreq);   // set the max termFreq in the collection
            CollectionStatistics.storeCollectionStatsIntoDisk();    // store collection statistics into disk

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /***
     * Add the current term to the inverted index
     *
     * @param term  the term passed as parameter
     * @param docId the docID of the current document
     * @param tf    the
     * @return false if the term has been already encountered in the current document,
     *         true if the term has been encountered for the first time in the current document or if the term was for
     *              the first time encountered
     * ***/
    private static boolean addTerm(String term, int docId, int tf)
    {
        // Initialize term frequency to 1 if tf is not provided (tf = 0 during index construction)
        int termFreq = (tf != 0) ? tf : 1;

        // Get or create the PostingList associated with the term
        if(!invertedIndex.containsKey(term))
            invertedIndex.put(term, new ArrayList<>());     // inverted index doesn't contain the term -> adds it

        int size = invertedIndex.get(term).size();          // take the size of posting lists related to the term

        // Check if the posting list is empty or if the last posting is for a different document
        if (invertedIndex.get(term).isEmpty() || invertedIndex.get(term).get(size - 1).getDocId() != docId)
        {
            invertedIndex.get(term).add(new Posting(docId, termFreq));  // Add a new posting for the current doc

            // Print term frequency and term frequency in the current posting (only during index construction)
            if (tf != 0)
                printDebug("SPIMI(add term): term: " + term + "TF: " + tf + " TERMFREQ: " + termFreq);

            return true;    // it's a new doc for this term -> Increment df
        }
        else    // term
        {
            invertedIndex.get(term).get(size - 1).addTermFreq(1); // Increment the term frequency for the current doc
            return false;   // this term has already been found in this doc -> no need to increment df
        }
    }

    // method to free memory by deleting the information in document table, dictionary,and inverted index
    private static void freeMemory()
    {
        documentTable.clear();
        dictionary.getTermToTermStat().clear();
        invertedIndex.clear();
    }

}
