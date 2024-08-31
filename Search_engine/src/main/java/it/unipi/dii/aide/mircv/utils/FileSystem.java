package it.unipi.dii.aide.mircv.utils;

import org.apache.commons.io.FileUtils;
import java.io.*;
import java.nio.channels.FileChannel;
import java.util.ArrayList;

import it.unipi.dii.aide.mircv.data_structures.SkipInfo;
import static it.unipi.dii.aide.mircv.utils.Constants.*;

public final class FileSystem
{
    private FileSystem() {
        throw new UnsupportedOperationException();
    }

    /**
     * Function to delete all the file except "stopwords.txt", "collection.tsv", and "msmarco-test2020-queries.tsv"
     * that are in resources.
     */
    public static void file_cleaner()
    {
        try {
            // delete folders
            File partial_folder = new File(PARTIAL_FOLDER);
            if(partial_folder.exists())     // check if exist the folder
                FileUtils.cleanDirectory(partial_folder);       // delete files in partial folder
            else
                partial_folder.mkdir();     // create the empty folder

            File merged_folder = new File(MERGED_FOLDER);
            if(merged_folder.exists())
                FileUtils.cleanDirectory(merged_folder);        // delete files in merged folder
            else
                merged_folder.mkdir();      // create the empty folder

            File upperBound_folder = new File(UPPERBOUND_FOLDER);
            if(upperBound_folder.exists())
                FileUtils.cleanDirectory(upperBound_folder);    // delete files in upperBound folder
            else
                upperBound_folder.mkdir();  // create the empty folder

            // delete files
            File flags = new File(FLAGS_FILE);
            if(flags.exists())
                FileUtils.delete(flags);                    // delete flags file if exist

            File collectionStatistics = new File(STATS_FILE);
            if(collectionStatistics.exists())
                FileUtils.delete(collectionStatistics);     // delete collection statistics file if exist

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * function to delete all the temporary files in the partial folder
     */
    public static void delete_tempFiles()
    {
        File partial_directory = new File(PARTIAL_FOLDER);
        try {
            FileUtils.cleanDirectory(partial_directory);    // delete files in partial folder
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * function to delete all the merged files in the merged folder
     */
    public static void delete_mergedFiles()
    {
        File dict = new File(DICTIONARY_FILE);
        File docid = new File(DOCID_FILE);
        File termfreq = new File(TERMFREQ_FILE);

        if(dict.exists() && docid.exists() && termfreq.exists()) // if the files exist delete them
        {
            try {
                FileUtils.delete(dict);         // delete dictionary file
                FileUtils.delete(docid);        // delete docid file
                FileUtils.delete(termfreq);     // delete termfreq file
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    /**
     * Function that check if there are all .txt files in "/resources/merged" folder
     * The file that controls are: dictionary.txt, docId.txt, documentTable.txt, termFreq.txt
     *
     * @return  true -> there are all merged files into disk
     *          false -> there aren't all merged files into disk
     */
    public static boolean areThereAllMergedFiles()
    {
        // define all file
        File docTable = new File(DOCTABLE_FILE);        // documentTable.txt
        File dict = new File(DICTIONARY_FILE);          // dictionary.txt"
        File docDID = new File(DOCID_FILE);             // docId.txt
        File docTF = new File(TERMFREQ_FILE);           // termFreq.txt

        return docTable.exists() && dict.exists() && docDID.exists() && docTF.exists();
    }

    /**
     * Function to show the size of all file saved into disk. It can be useful for analyze and compare the size of the
     * files (and consequently of the inverted index, dictionary, doc table, skip info adn Term Upper Bound) of different version and with
     * different techniques enabled (skipping, compression and so on)
     */
    public static void printFileSize()
    {
        printUI("File sizes, if any, will be shown (divided by category).");

        printUI("Partial files:");      // size of the partial files
        File block = new File(BLOCKOFFSETS_FILE);               // blocks.txt
        if(block.exists())
            printSize(formatSize("blocks", block.length()));

        File partDict = new File(PARTIAL_DICTIONARY_FILE);      // partial_dictionary.txt"
        if(partDict.exists())
            printSize(formatSize("partial_dictionary", partDict.length()));

        File partDocDID = new File(PARTIAL_DOCID_FILE);         // partial_docId.txt
        if(partDocDID.exists())
            printSize(formatSize("partial_docId", partDocDID.length()));

        File partDocTF = new File(PARTIAL_TERMFREQ_FILE);       // partial_termFreq.txt
        if(partDocTF.exists())
            printSize(formatSize("partial_termFreq", partDocTF.length()));

        printUI("Merged files:");       // size of the merged files
        File docTable = new File(DOCTABLE_FILE);        // documentTable.txt
        if(docTable.exists())
            printSize(formatSize("documentTable", docTable.length()));

        File dict = new File(DICTIONARY_FILE);          // dictionary.txt"
        if(dict.exists())
            printSize(formatSize("dictionary", dict.length()));

        File docDID = new File(DOCID_FILE);             // docId.txt
        if(docDID.exists())
            printSize(formatSize("docId", docDID.length()));

        File docTF = new File(TERMFREQ_FILE);           // termFreq.txt
        if(docTF.exists())
            printSize(formatSize("termFreq", docTF.length()));

        File skip = new File(SKIP_FILE);                // skipInfo
        if(skip.exists())
            printSize(formatSize("skipInfo", skip.length()));

        printUI("TUB files:");                      // size of term upper bound file
        File termUB = new File(TERMUPPERBOUND_FILE);
        if(termUB.exists())
            printSize(formatSize("termsUpperBound", termUB.length()));

        File termFreqWeight = new File(TERMFREQWEIGHT_FILE);    // size of the precomputed termFreqWeight file
        if(termUB.exists())
            printSize(formatSize("termFreqWeight", termFreqWeight.length()));

        printUI("Generic files:");      // size of the other files (collection statistic and flags files)
        File flags = new File(FLAGS_FILE);
        if(flags.exists())
            printSize(formatSize("flags", flags.length()));

        File collectionStatistics = new File(STATS_FILE);       // size of the collection statistics
        if(collectionStatistics.exists())
            printSize(formatSize("collectionStatistics", collectionStatistics.length()));
    }

    // ------------------------------- start: method not used (for future implementation) ------------------------------
    /**
     * Function to save docids or tf posting list into file (in order to compare before and after compression)
     *
     * @param postings      array list of the DID
     * @param tempFileName  name of the file
     * @throws FileNotFoundException
     */
    public static void saveDocsInFile(ArrayList<Integer> postings, String tempFileName) throws FileNotFoundException
    {
        // Create a file
        File outputf = new File(tempFileName);

        try (PrintWriter outputWriter = new PrintWriter(outputf))
        {
            for (int i = 0; i < postings.size(); i++)
            {
                printDebug("posting" + i + ": " + postings.get(i));
                outputWriter.print(postings.get(i));
                outputWriter.println(); // Add a newline character
            }
        }
    }

    /**
     *
     * @param si
     * @param tempFileName
     */
    public static void saveDocsInFileSkipInfo(SkipInfo si, String tempFileName)
    {
        // ----------- debug file ---------------
        File outputf = new File(tempFileName);

        try(PrintWriter outputWriter = new PrintWriter(outputf);)
        {
            outputWriter.print(si.toString());
            outputWriter.println();
        } catch (FileNotFoundException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     *
     * @param data
     * @param fileName
     */
    public static void saveStructureToFile(ArrayList<String> data, String fileName)
    {
        try (FileWriter writer = new FileWriter(DEBUG_FOLDER + fileName, false))
        {
            for (String line : data) {
                writer.write(line);
                writer.write(System.lineSeparator());
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void appendStringToFile(String data, String fileName)
    {
        try (FileWriter writer = new FileWriter(DEBUG_FOLDER + fileName, true))
        {
            writer.write(data);
            writer.write(System.lineSeparator());
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    // ------------------------------- end: method not used (for future implementation) ------------------------------
}
