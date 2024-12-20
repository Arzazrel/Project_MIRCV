package it.unipi.dii.aide.mircv.data_structures;

import java.io.IOException;
import java.nio.CharBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.StandardCharsets;

import static it.unipi.dii.aide.mircv.utils.Constants.*;
import static it.unipi.dii.aide.mircv.utils.Logger.dict_logger;

/**
 *  Stores unique terms and their statistics.
 */
public class DictionaryElem
{
    private String term;        // 20 byte
    private int df;             // document frequency, number of documents in which there is the term
    private int cf;             // collection frequency, number of occurrences of the term in the collection
    private long offsetTermFreq;// starting point of the posting list of the term in the term freq file (set in partial(SPIMI) and total index)
    private long offsetDocId;   // starting point of the posting list of the term in the docid file (set in partial(SPIMI) and total index)
    // compression
    private int docIdSize;      // dimension in byte of compressed DocID of the posting list
    private int termFreqSize;   // dimension in byte of compressed termFreq of the posting list
    // skipping
    private long skipOffset;    // offset of the skip element
    private int skipArrLen;     // len of the skip array (equal to the number of skipping block)

    /**
     * Constructor without argument.
     */
    DictionaryElem()
    {
        this.df = 0;
        this.cf = 0;
        this.term = "";
        this.docIdSize = 0;
        this.termFreqSize = 0;
        this.skipOffset = -1;
        this.skipArrLen = -1;
    }

    /**
        Constructor with TermID parameter. Called when the first occurrence of a term is found.
        Is the first occurrence found of the term in the collection, will be in at least one document and present
        once in the collection for these set df and cf to 1.
     */
    public DictionaryElem(String term) {
        this.term = term;
        this.df = 0;
        this.cf = 0;
        this.docIdSize = 0;
        this.termFreqSize = 0;
        this.skipOffset = 0;
        this.skipArrLen = 0;
    }

    // add the quantity passed as a parameter to the current Df
    public void addDf(int n) { setDf(getDf() + n); }

    // add the quantity passed as a parameter to the current Cf
    public void addCf(int n) { setCf(getCf() + n); }

    // -------- start get methods --------
    public void setDf(int df) { this.df = df; }

    public void setCf(int cf) { this.cf = cf; }

    public void setTerm(String term) { this.term = term; }

    public void setOffsetTermFreq(long offsetTermFreq) {
        this.offsetTermFreq = offsetTermFreq;
    }

    public void setOffsetDocId(long offsetDocId) {
        this.offsetDocId = offsetDocId;
    }

    public void setDocIdSize(int docIdSize) {
        this.docIdSize = docIdSize;
    }

    public void setTermFreqSize(int termFreqSize) {
        this.termFreqSize = termFreqSize;
    }

    public void setSkipOffset(long skipOffset) { this.skipOffset = skipOffset; }

    public void setSkipArrLen(int skipArrLen) { this.skipArrLen = skipArrLen; }

    // -------- start set methods --------

    public int getDf() { return df; }

    public int getCf() { return cf; }

    public String getTerm() {
        return term;
    }

    public long getOffsetTermFreq() {
        return offsetTermFreq;
    }

    public long getOffsetDocId() {
        return offsetDocId;
    }

    public int getDocIdSize() {
        return docIdSize;
    }

    public int getTermFreqSize() {
        return termFreqSize;
    }

    public long getSkipOffset() { return skipOffset; }

    public int getSkipArrLen() { return skipArrLen; }

    /**
     * Function that returns the size of a DictionaryElem.
     * The size change in according to the value of the flags for skipping and compression.
     *
     * @return  the size of dictionaryElem as int
     */
    static int getDictElemSize()
    {
        // if compression case, need to store 2 more integers (dimension of compressed DocID and Term Frequency values)
        int DICT_ELEM_SIZE = TERM_DIM + 2*Integer.BYTES + 2*Long.BYTES;
        return DICT_ELEM_SIZE + ((Flags.isCompressionEnabled()) ? 2*Integer.BYTES : 0)
                + ((Flags.considerSkippingBytes()) ? (Long.BYTES + Integer.BYTES) : 0);
    }

    @Override
    public String toString()
    {
        return "DictionaryElem{" +
                "term='" + term + '\'' +
                ", df=" + df +
                ", cf=" + cf +
                ", offsetTermFreq=" + offsetTermFreq +
                ", offsetDocId=" + offsetDocId +
                ", docIdSize=" + docIdSize +
                ", termFreqSize=" + termFreqSize +
                ", offsetSkip=" + skipOffset +
                ", skipSize=" + skipArrLen +
                '}';
    }

    /**
     * function to store one dictionary elem into disk
     *
     * @param channel   indicate the file where to write
     */
    void storeDictionaryElemIntoDisk(FileChannel channel)
    {
        try {
            MappedByteBuffer buffer;
            buffer = channel.map(FileChannel.MapMode.READ_WRITE, channel.size(), getDictElemSize());

            if(buffer == null)      // Buffer not created
                return;

            CharBuffer charBuffer = CharBuffer.allocate(TERM_DIM);  // allocate bytes for docno
            //put every char into charbuffer
            for(int i = 0; i < term.length(); i++)
                charBuffer.put(i, term.charAt(i));

            // write docno, docid and doclength into document file
            buffer.put(StandardCharsets.UTF_8.encode(charBuffer));      // write term
            buffer.putInt(df);                                          // write Df
            buffer.putInt(cf);                                          // write Cf
            buffer.putLong(offsetTermFreq);                             // write offset Tf
            buffer.putLong(offsetDocId);
            // if in merge phase, need to store also the size of DocID and Term Frequency compressed values
            if (Flags.isCompressionEnabled())       // if compression is enabled
            {
                buffer.putInt(termFreqSize);            // write dimension in byte of compressed DocID of the PL
                buffer.putInt(docIdSize);               // write dimension in byte of compressed termFreq of the PL
            }
            // write offset DID
            if(Flags.considerSkippingBytes())       // if skipping is enabled
            {
                buffer.putLong(skipOffset);             // write offset of the skip element
                buffer.putInt(skipArrLen);              // write len of the skip array (equal to the number of skipping block)
            }

            if(debug)
                dict_logger.logInfo(this.toString());   // save log info

            PARTIAL_DICTIONARY_OFFSET += getDictElemSize();       // update offset

            if(term.equals(""))
                printError("Empty term.");

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * function to read one Dictionary Element from disk
     *
     * @param start     offset of the document reading from document file
     * @param channel   indicate the file from which to read
     * @return a DictionaryElem with the value read from disk
     */
    public void readDictionaryElemFromDisk(long start, FileChannel channel)
    {
        try {
            MappedByteBuffer buffer = channel.map(FileChannel.MapMode.READ_ONLY, start,  getDictElemSize()); // get first term of the block
            CharBuffer.allocate(TERM_DIM);              //allocate a charbuffer of the dimension reserved to docno
            CharBuffer charBuffer = StandardCharsets.UTF_8.decode(buffer);

            if(!(charBuffer.toString().split("\0").length == 0))
                term = charBuffer.toString().split("\0")[0];   // read term

            buffer.position(TERM_DIM);          // skip over term position
            df = buffer.getInt();               // read Df
            cf = buffer.getInt();               // read Cf
            offsetTermFreq = buffer.getLong();  // read offset Tf
            offsetDocId = buffer.getLong();     // read offset DID
            if (Flags.isCompressionEnabled())   // if compression is enabled
            {
                termFreqSize = buffer.getInt(); // read dimension in byte of compressed DocID of the PL
                docIdSize = buffer.getInt();    // read dimension in byte of compressed termFreq of the PL
            }
            if(Flags.considerSkippingBytes())   // if skipping is enabled
            {
                skipOffset = buffer.getLong();  // read offset of the skip element
                skipArrLen = buffer.getInt();   // read len of the skip array (equal to the number of skipping block)
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
