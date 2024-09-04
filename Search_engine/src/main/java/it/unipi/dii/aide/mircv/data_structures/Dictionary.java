
package it.unipi.dii.aide.mircv.data_structures;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.CharBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.StandardCharsets;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.stream.Collectors;

import static it.unipi.dii.aide.mircv.data_structures.DictionaryElem.getDictElemSize;
import static it.unipi.dii.aide.mircv.data_structures.DocumentElement.DOCELEM_SIZE;
import static it.unipi.dii.aide.mircv.utils.Constants.*;
/**
 * Dictionary class
 */
public class Dictionary
{
    // HashMap between a term and its statistics, contained in a DictionaryElem.
    private HashMap<String, DictionaryElem> termToTermStat;

    // constructor without parameters
    public Dictionary() {
        this.termToTermStat = new HashMap<>();
    }

    /**
     * Function which returns, if present, the DictionaryElem associated with the term passed as a parameter.
     * Otherwise, it creates a new DictionaryElem associated with the term, inserts it in the HashMap and returns it.
     *
     * @param term          the term (string) searched or to be created
     */
    public DictionaryElem getOrCreateTerm(String term)
    {
        return termToTermStat.computeIfAbsent(term, t -> new DictionaryElem(term));
    }

    /**
     * Function that return the DictionaryElem associated to the term passed as parameter
     *
     * @param term  the term (string) passed as parameter of which to obtain the information
     * @return the dictionary element related to the term passed as parameter
     */
    public DictionaryElem getTermStat(String term) {
        return termToTermStat.get(term);
    }

    /**
     * Function that return the Dictionary (hash map)
     *
     * @return the whole hash map containing the dictionary
     */
    public HashMap<String, DictionaryElem> getTermToTermStat() {
        return termToTermStat;
    }

    /**
     * Function to sort in lexicographic order the term in the dictionary
     */
    public void sort()
    {
        termToTermStat = getTermToTermStat().entrySet().stream()
                .sorted(Map.Entry.comparingByKey())
                .collect(Collectors.toMap(
                        Map.Entry::getKey,
                        Map.Entry::getValue,
                        (a, b) -> { throw new AssertionError(); },
                        LinkedHashMap::new
                ));
    }

    /**
     * Function that return if the dictionary in memory is set or not
     */
    public boolean dictionaryIsSet()
    {
        return !termToTermStat.isEmpty();  // the hash map in dictionary is empty, the dictionary isn't set
    }

    /**
     *  Function to read whole Dictionary from disk.
     */
    public void readDictionaryFromDisk()
    {
        // Define maximum chunk size (default = 80% of the total memory available)
        final long CHUNK_SIZE = (long) (Runtime.getRuntime().maxMemory() * MEMORY_THRESHOLD);
        int dictElemSize = getDictElemSize();   // size of each element of the dictionary
        long dictSize = 0;          // the size of the DocumentTable
        long position = 0;          // current position
        MappedByteBuffer buffer;    // the buffer for the reading

        printLoad("Loading dictionary from disk..."); // control print
        try (
                FileChannel channel = new RandomAccessFile(DICTIONARY_FILE, "rw").getChannel()
        ) {
            dictSize = channel.size();                          // size of the dictionary saved into disk

            while (position < dictSize)
            {
                long remaining = dictSize - position;               // how much is left to read
                long chunkSize = Math.min(remaining, CHUNK_SIZE);   // get the dimension for the current chunk to read
                // map the current chunk in memory
                buffer = channel.map(FileChannel.MapMode.READ_ONLY, position, chunkSize);

                int offset = 0;             // reset the offset for the current chunk
                while (offset < chunkSize)
                {
                    if (buffer.remaining() < dictElemSize)
                        break;              // exit if there are not enough bytes for another element

                    DictionaryElem dictElem = new DictionaryElem();
                    // Decoding the characters for the term
                    byte[] termBytes = new byte[TERM_DIM];
                    buffer.get(termBytes);                          // read the term
                    String term = new String(termBytes, StandardCharsets.UTF_8).trim();
                    if (term.isEmpty())                             // control check for term
                        break;

                    dictElem.setTerm(term);                         // read term
                    dictElem.setDf(buffer.getInt());                // read and set Df
                    dictElem.setCf(buffer.getInt());                // read and set Cf
                    dictElem.setOffsetTermFreq(buffer.getLong());   // read and set offset Tf
                    dictElem.setOffsetDocId(buffer.getLong());      // read and set offset DID
                    if(Flags.isCompressionEnabled())        // check if the compression flag is enabled
                    {
                        dictElem.setTermFreqSize(buffer.getInt());  // read dimension in byte of compressed DocID of the PL
                        dictElem.setDocIdSize(buffer.getInt());     // read dimension in byte of compressed termFreq of the PL
                    }
                    if(Flags.considerSkippingBytes())       // if skipping is enabled
                    {
                        dictElem.setSkipOffset(buffer.getLong());   // read offset of the skip element
                        dictElem.setSkipArrLen(buffer.getInt());    // read len of the skip array (equal to the number of skipping block)
                    }

                    termToTermStat.put(term, dictElem);             // add DictionaryElem into memory
                    offset += dictElemSize;                         // update offset
                }
                position += offset;                 // update the position (the byte read)
            }
            //printLoad("Vocabulary size: " + termToTermStat.size());
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}

