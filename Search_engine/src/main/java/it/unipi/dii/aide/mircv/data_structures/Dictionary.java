
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
    public DictionaryElem getOrCreateTerm(String term, int termCounter)
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
     * Function to read whole Dictionary from disk
     */
    public void readDictionaryFromDisk()
    {
        long position = 0;          // indicate the position where read at each iteration
        MappedByteBuffer buffer;    // get first term of the block

        printLoad("Loading dictionary from disk...");       // control print
        try (
                FileChannel channel = new RandomAccessFile(DICTIONARY_FILE, "rw").getChannel()
        ) {
            long len = channel.size();          // size of the dictionary saved into disk

            // scan all Dictionary Element saved into disk
            while(position < len)
            {
                buffer = channel.map(FileChannel.MapMode.READ_ONLY, position, getDictElemSize());// read one DictionaryElem
                position += getDictElemSize();                      // update read position

                DictionaryElem dictElem = new DictionaryElem();     // create new DictionaryElem

                CharBuffer.allocate(TERM_DIM);              //allocate a charbuffer of the dimension reserved to docno
                CharBuffer charBuffer = StandardCharsets.UTF_8.decode(buffer);

                if(charBuffer.toString().split("\0").length == 0)       // control check of the term size
                    continue;

                String term = charBuffer.toString().split("\0")[0];     // split using end string character

                /*
                if(term.equals("epstein"))      // debug print
                    printDebug("TERM: " + term);
                //*/

                dictElem.setTerm(term);                         // read term
                buffer.position(TERM_DIM);                      // skip docno
                dictElem.setDf(buffer.getInt());                // read and set Df
                dictElem.setCf(buffer.getInt());                // read and set Cf
                dictElem.setOffsetTermFreq(buffer.getLong());   // read and set offset Tf
                dictElem.setOffsetDocId(buffer.getLong());      // read and set offset DID
                if(Flags.isCompressionEnabled())    // check if the compression flag is enabled
                {
                    dictElem.setTermFreqSize(buffer.getInt());  //
                    dictElem.setDocIdSize(buffer.getInt());     //
                }
                dictElem.setSkipOffset(buffer.getLong());       //
                dictElem.setSkipOffset(buffer.getInt());        //
                termToTermStat.put(term, dictElem);             // add DictionaryElem into memory
            }

            printLoad("vocabulary size: " + termToTermStat.size());
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}

