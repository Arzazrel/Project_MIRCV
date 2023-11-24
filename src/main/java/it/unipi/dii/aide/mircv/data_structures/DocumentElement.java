package it.unipi.dii.aide.mircv.data_structures;

import java.io.IOException;
import java.nio.CharBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.StandardCharsets;

import static it.unipi.dii.aide.mircv.utils.Constants.*;

public class DocumentElement {

    public static final int DOCNO_DIM = 10;                     // Length of docno (in bytes)
    public final static int DOCELEM_SIZE = 4 + DOCNO_DIM + 4;   // Size in bytes of docid, docno, and doclength

    //single document element with correspondent docid and relative length
    private int docid;
    private String docno;
    private int doclength;


    public DocumentElement(String docno, int docid, int doclength) {
        this.docno = docno;
        this.doclength = doclength;
        this.docid = docid;
    }

    public DocumentElement() {
        this.docno = "";
        this.docid = 0;
        this.doclength = 0;
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


    /**
     * function to read one Document Element from disk
     *
     * @param start     offset of the document reading from document file
     * @param channel   indicate the file from which to read
     * @return a DocumentElement with the value read from disk
     */
    public void readDocumentElementFromDisk(int start, FileChannel channel) throws IOException {

        MappedByteBuffer buffer = channel.map(FileChannel.MapMode.READ_ONLY, start, DOCELEM_SIZE);

        if(buffer == null)      // Buffer not created
            return;

        CharBuffer.allocate(DOCNO_DIM); //allocate a charbuffer of the dimension reserved to docno
        CharBuffer charBuffer = StandardCharsets.UTF_8.decode(buffer);

        if(charBuffer.toString().split("\0").length == 0)
            return;

        docno =  charBuffer.toString().split("\0")[0]; //split using end string character
        buffer.position(DOCNO_DIM);             //skip docno
        docid = buffer.getInt();
        doclength = buffer.getInt();

    }
}
