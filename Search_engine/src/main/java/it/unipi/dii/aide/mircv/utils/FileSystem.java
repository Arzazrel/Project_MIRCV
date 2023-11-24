package it.unipi.dii.aide.mircv.utils;

import java.io.*;

import it.unipi.dii.aide.mircv.data_structures.SkipInfo;
import org.apache.commons.io.FileUtils;

import java.util.ArrayList;

import static it.unipi.dii.aide.mircv.utils.Constants.*;


public final class FileSystem {

    private FileSystem() {
        throw new UnsupportedOperationException();
    }

    /**
     * function to delete all the file except "stopwords.txt", "collection.tsv", and "msmarco-test2020-queries.tsv"
     * that are in resources.
     */
    public static void file_cleaner() {

        try {
            File partial_folder = new File(PARTIAL_FOLDER);
            FileUtils.cleanDirectory(partial_folder);
            File merged_folder = new File(MERGED_FOLDER);
            FileUtils.cleanDirectory(merged_folder);

            File flags = new File(FLAGS_FILE);
            if(flags.exists())
                FileUtils.delete(flags);

        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public static void delete_tempFiles() {

        File partial_directory = new File(PARTIAL_FOLDER);
        try {
            FileUtils.cleanDirectory(partial_directory);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public static void delete_mergedFiles() {

        File dict = new File(DICTIONARY_FILE);
        File docid = new File(DOCID_FILE);
        File termfreq = new File(TERMFREQ_FILE);
        if(dict.exists() && docid.exists() && termfreq.exists()) {
            try {
                FileUtils.delete(dict);
                FileUtils.delete(docid);
                FileUtils.delete(termfreq);
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

    // function to save docids  or tf posting list into file (in order to compare before and after compression)
    public static void saveDocsInFile(ArrayList<Integer> postings, String tempFileName) throws FileNotFoundException {
        // Create a file
        File outputf = new File(tempFileName);

        try (PrintWriter outputWriter = new PrintWriter(outputf)) {
            for (int i = 0; i < postings.size(); i++) {
                printDebug("posting" + i + ": " + postings.get(i));
                outputWriter.print(postings.get(i));
                outputWriter.println(); // Add a newline character
            }
        }
    }

    public static void saveDocsInFileSkipInfo(SkipInfo si, String tempFileName) {
        // ----------- debug file ---------------
        File outputf = new File(tempFileName);

        try(PrintWriter outputWriter = new PrintWriter(outputf);) {

            outputWriter.print(si.toString());
            outputWriter.println();
        } catch (FileNotFoundException e) {
            throw new RuntimeException(e);
        }
    }


//    public static void saveStructureToFile(ArrayList<String> data, String fileName) {
//        try (FileWriter writer = new FileWriter(DEBUG_FOLDER + fileName, false)) {
//            for (String line : data) {
//                writer.write(line);
//                writer.write(System.lineSeparator());
//            }
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
//    }
//
//    public static void appendStringToFile(String data, String fileName) {
//        try (FileWriter writer = new FileWriter(DEBUG_FOLDER + fileName, true)) {
//            writer.write(data);
//            writer.write(System.lineSeparator());
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
//    }

}
