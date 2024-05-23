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

public class TextProcessor
{
    public static List<String> globalStopwords;         // list of strings to store stopwords

    /**
     * Function to preprocessing the query of the user (in word).
     * According to the flags of the user return the term of the query from the query in word.
     *
     * @param input             is the query of the users (in words)
     * @return  an ArrayList of string that contains the term of the query after preprocessing
     */
    public static ArrayList<String> preprocessText(String input) throws IOException
    {
        ArrayList<String> tokenList;

        // Preprocess the Text
        input = cleanText(input);           // Clean the text by removing URLs, HTML tags, punctuation, etc...
        tokenList = tokenizeText(input);    // Tokenize the cleaned text into individual words

        // Remove stop-words and perform stemming
        if (Flags.isSwsEnabled())   // Check if the filtering flag is enabled
        {
            tokenList = removeStopwords(tokenList);     // Remove common stopwords
            tokenList = applyStemming(tokenList);       // Perform stemming on the remaining words
        }

        return tokenList;   // Return the preprocessed tokens
    }

    /**
     * Function to remove unicode characters from a string.
     *
     * @param input             is the query of the users (in words)
     * @return  a string that contains the string without non ASCII characters
     */
    private static String removeNonASCIIChars(String input)
    {
        String cleanedInput;        // to contain the cleaned input
        byte[] inputBytes = input.getBytes(StandardCharsets.UTF_8);     // convert string into array of bytes
        cleanedInput = new String(inputBytes, StandardCharsets.UTF_8);  // encode in UTF-8 the string

        // Define a pattern to match any Unicode characters outside the ASCII range
        Pattern nonASCIICharsPattern = Pattern.compile("[^\\x00-\\x7F]",
                Pattern.UNICODE_CASE | Pattern.CANON_EQ | Pattern.CASE_INSENSITIVE);

        Matcher nonASCIICharsMatcher = nonASCIICharsPattern.matcher(cleanedInput);  // find all the non-ASCII characters that are present in input (if there are)
        cleanedInput = nonASCIICharsMatcher.replaceAll(" "); // Replace non-ASCII characters with a space

        return cleanedInput;
    }

    /**
     * Function to clean the input text.
     *
     * @param input             is the query of the users (in words)
     * @return  a string that contains the string without not permitted or useless characters
     */
    private static String cleanText(String input) {

        input = input.replaceAll("https?://\\S+\\s?", " "); // Replace URLs with spaces

        input = input.toLowerCase();        // Convert text to lowercase

        input = input.replaceAll("<[^>]*>", "");        // Remove HTML tags

        input = input.replaceAll("\\p{Punct}", " ");    // Replace punctuation with spaces

        input = removeNonASCIIChars(input); // Remove non-ASCII characters

        input = input.replaceAll("\\s+", " ");  // Replace multiple(extra) spaces with a single space

        return input;
    }

    /**
     * Function to tokenize a string into individual words.
     *
     * @param input             is the query of the users (in words)
     * @return  an arraylist of string that contains the term(string) of the query
     */
    private static ArrayList<String> tokenizeText(String input) {
        return new ArrayList<>(Arrays.asList(input.toLowerCase().split(" "))); // Tokenize by splitting on spaces
    }

    /**
     * Function to apply stemming on tokens.
     *
     * @param tokens             is the query of the users (in token)
     * @return  an arraylist of string that contains the token of the query after the stemming
     */
    private static ArrayList<String> applyStemming(ArrayList<String> tokens)
    {
        PorterStemmer stemmer = new PorterStemmer();    // Initialize the porterStemmer to perform stemming

        for (int i = 0; i < tokens.size(); i++)     // Scan all tokens
        {
            stemmer.setCurrent(tokens.get(i));      // Set the current word for stemming
            stemmer.stem();                         // Perform stemming
            tokens.set(i, stemmer.getCurrent());    // Replace the word with its stemmed form
        }
        return tokens;      // Return the tokens after stemming
    }

    /**
     * Function to remove stopwords from tokens.
     *
     * @param tokens             is the query of the users (in token)
     * @return  an arraylist of string that contains the token of the query after removing the stopwords
     */
    private static ArrayList<String> removeStopwords(ArrayList<String> tokens) throws IOException
    {
        globalStopwords = Files.readAllLines(Paths.get("src/main/resources/stopwords.txt")); // Read stopwords from a file
        tokens.removeAll(globalStopwords); // Remove stopwords from the token list

        return tokens; // Return tokens after stopwords removal
    }
}
