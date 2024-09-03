package it.unipi.dii.aide.mircv.data_structures;

import java.io.IOException;
import java.nio.CharBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.StandardCharsets;
import java.util.Objects;

import static it.unipi.dii.aide.mircv.utils.Constants.*;

/**
 * class to
 */
public class DocumentElement
{
    // -- static variable
    public static final int DOCNO_DIM = 10;         // Length of docno (in bytes), maximum 10 Bytes
    public final static int DOCELEM_SIZE = INT_BYTES + DOCNO_DIM + INT_BYTES + DOUBLE_BYTES;   // Size in bytes of docid, docno, and doclength

    // -- not static variable (single document element)
    private int docid;      // docID of the document
    private String docno;   // docNo of the document
    private int doclength;  // length (number of word) of the document

    // the denominator in BM25 is = k * ((1 - b) + b * (docLen / avgDocLen)) + termFreq;
    private double denomPartBM25;   // first part of denominator, saved for optimization purpose = k * ((1 - b) + b * (docLen / avgDocLen))

    /**
     * Constructor with parameter.
     *
     * @param docno     the docNo of the document
     * @param docid     the docID of the document
     * @param doclength the length (number of words) of the document
     */
    public DocumentElement(String docno, int docid, int doclength)
    {
        this.docno = docno;
        this.doclength = doclength;
        this.docid = docid;
        this.denomPartBM25 = 0;
    }

    /**
     * Constructor with all parameters.
     *
     * @param docno         the docNo of the document
     * @param docid         the docID of the document
     * @param doclength     the length (number of words) of the document
     * @param denomPartBM25 // first part of denominator in BM25
     */
    public DocumentElement(String docno, int docid, int doclength, double denomPartBM25)
    {
        this.docno = docno;
        this.doclength = doclength;
        this.docid = docid;
        this.denomPartBM25 = denomPartBM25;
    }

    /**
     * Constructor without parameters.
     */
    public DocumentElement() {
        this.docno = "";
        this.docid = 0;
        this.doclength = 0;
        this.denomPartBM25 = 0;
    }

    @Override
    public String toString()
    {
        return "DocumentElement\n{" +
                "\n - docno = " + docno +
                "\n - docid = " + docid +
                "\n - doclength = " + doclength +
                "\n - denomPartBM25 = " + denomPartBM25 +
                "\n}";
    }

    // Override of the equals method to compare two instances of DocumentElement
    @Override
    public boolean equals(Object obj)
    {
        if (this == obj)    // check if they are the same object
            return true;

        if (obj == null || getClass() != obj.getClass())    // obj is null or not a Person
            return false;

        DocumentElement obj2 = (DocumentElement) obj;
        return ( docno.equals(obj2.docno) && (docid == obj2.docid) && (doclength == obj2.doclength) && (denomPartBM25 == obj2.denomPartBM25) );
    }

    // Override of the hashCode method
    @Override
    public int hashCode() {
        return Objects.hash(docno, docid, doclength, denomPartBM25);    // create hash code
    }

    // ---- start method get and set ----

    public int getDoclength() {
        return doclength;
    }

    public void setDoclength(int doclength) {
        this.doclength = doclength;
    }

    public String getDocno() {
        return docno;
    }

    public void setDocno(String docno) {
        this.docno = docno;
    }

    public int getDocid() {
        return docid;
    }

    public void setDocid(int docid) {
        this.docid = docid;
    }

    public double getDenomPartBM25() {
        return denomPartBM25;
    }

    public void setDenomPartBM25(double denomPartBM25) {
        this.denomPartBM25 = denomPartBM25;
    }


    /**
     * Function to read one Document Element from disk and set the values read into the variables of this instance
     *
     * @param start     offset of the document reading from document file
     * @param channel   indicate the file from which to read
     */
    public void readDocumentElementFromDisk(int start, FileChannel channel) throws IOException
    {
        MappedByteBuffer buffer = channel.map(FileChannel.MapMode.READ_ONLY, start, DOCELEM_SIZE);

        if(buffer == null)      // Buffer not created
            return;

        CharBuffer.allocate(DOCNO_DIM); //allocate a charbuffer of the dimension reserved to docno
        CharBuffer charBuffer = StandardCharsets.UTF_8.decode(buffer);

        if(charBuffer.toString().split("\0").length == 0)   // check that the docNo is not empty
            return;

        docno =  charBuffer.toString().split("\0")[0];      //split using end string character
        buffer.position(DOCNO_DIM);                               // skip docno
        docid = buffer.getInt();                                  // read the docID
        doclength = buffer.getInt();                              // read the length of the document
        denomPartBM25 = buffer.getDouble();                       // read the part of denominator of BM25
    }

    /**
     * Function to read one Document Element from disk and set the values read into the variables of this instance.
     * Version with MappedByteBuffer not for each document element, version much faster version than the one with the
     * buffer for each element.
     *
     * @param start     offset of the document reading from document file
     * @param buffer    indicate the buffer (MappedByteBuffer) from which to read
     */
    public void readDocumentElementFromDisk(int start, MappedByteBuffer buffer)
    {
        buffer.position(start);

        // check if there is not enough data in the buffer to read a complete DocumentElement
        if (buffer.remaining() < DOCELEM_SIZE)
            return;

        // Decoding the characters for the docno
        byte[] docnoBytes = new byte[DOCNO_DIM];
        buffer.get(docnoBytes);                     // read the docno
        String docnoStr = new String(docnoBytes, StandardCharsets.UTF_8).trim();
        if (docnoStr.isEmpty())                     // control check for docno
            return;

        this.docno = docnoStr;                      // take the docno
        this.docid = buffer.getInt();               // read the docID
        this.doclength = buffer.getInt();           // read the length of the document
        this.denomPartBM25 = buffer.getDouble();    // read the part of denominator of BM25
    }

}
