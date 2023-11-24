package it.unipi.dii.aide.mircv;

import it.unipi.dii.aide.mircv.data_structures.Flags;
import org.tartarus.snowball.ext.PorterStemmer;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class TextProcessor {

    // Global variable to store stopwords
    public static List<String> globalStopwords;

    // Main function to execute text preprocessing
    public static ArrayList<String> preprocessText(String input) throws IOException {
        ArrayList<String> tokenList;

        /* Preprocess the Text */
        input = cleanText(input); // Clean the text by removing URLs, HTML tags, punctuation, etc.
        tokenList = tokenizeText(input); // Tokenize the cleaned text into individual words

        /* Remove stop-words and perform stemming */
        if (Flags.isSwsEnabled()) { // Check if the filtering flag is enabled
            tokenList = removeStopwords(tokenList); // Remove common stopwords
            tokenList = applyStemming(tokenList); // Perform stemming on the remaining words
        }

        return tokenList; // Return the preprocessed tokens
    }

    // Function to remove unicode characters from a string
    private static String removeNonASCIIChars(String input) {
        String cleanedInput;
        byte[] inputBytes = input.getBytes(StandardCharsets.UTF_8);
        cleanedInput = new String(inputBytes, StandardCharsets.UTF_8);

        // Define a pattern to match any Unicode characters outside the ASCII range
        Pattern nonASCIICharsPattern = Pattern.compile("[^\\x00-\\x7F]",
                Pattern.UNICODE_CASE | Pattern.CANON_EQ | Pattern.CASE_INSENSITIVE);

        Matcher nonASCIICharsMatcher = nonASCIICharsPattern.matcher(cleanedInput);
        cleanedInput = nonASCIICharsMatcher.replaceAll(" "); // Replace non-ASCII characters with a space

        return cleanedInput;
    }

    // Function to clean the input text
    private static String cleanText(String input) {
        /* Remove URLs */
        input = input.replaceAll("https?://\\S+\\s?", " "); // Replace URLs with spaces

        /* Convert to lowercase */
        input = input.toLowerCase(); // Convert text to lowercase

        /* Remove HTML tags */
        input = input.replaceAll("<[^>]*>", ""); // Remove HTML tags

        /* Remove punctuation */
        input = input.replaceAll("\\p{Punct}", " "); // Replace punctuation with spaces

        /* Remove non-ASCII characters */
        input = removeNonASCIIChars(input); // Remove non-ASCII characters

        /* Remove extra whitespaces */
        input = input.replaceAll("\\s+", " "); // Replace multiple spaces with a single space

        return input;
    }

    // Function to tokenize a string into individual words
    private static ArrayList<String> tokenizeText(String input) {
        return new ArrayList<>(Arrays.asList(input.toLowerCase().split(" "))); // Tokenize by splitting on spaces
    }

    // Function to apply stemming on tokens
    private static ArrayList<String> applyStemming(ArrayList<String> tokens) {
        PorterStemmer stemmer = new PorterStemmer();
        for (int i = 0; i < tokens.size(); i++) {
            stemmer.setCurrent(tokens.get(i)); // Set the current word for stemming
            stemmer.stem(); // Perform stemming
            tokens.set(i, stemmer.getCurrent()); // Replace the word with its stemmed form
        }
        return tokens; // Return the tokens after stemming
    }

    // Function to remove stopwords from tokens
    private static ArrayList<String> removeStopwords(ArrayList<String> tokens) throws IOException {
        globalStopwords = Files.readAllLines(Paths.get("src/main/resources/stopwords.txt")); // Read stopwords from a file
        tokens.removeAll(globalStopwords); // Remove stopwords from the token list
        return tokens; // Return tokens after stopwords removal
    }
}
