package it.unipi.dii.aide.mircv.data_structures;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import static it.unipi.dii.aide.mircv.utils.Constants.*;

public class SkipInfo
{
    public static final int SKIPPING_INFO_SIZE = 2 * LONG_BYTES + INT_BYTES;    // the size
    private int maxDocId;       // the maximum DocID in the skipping block
    private long docIdOffset;   // offset of the first docID in the next skipping block
    private long freqOffset;    // offset of the first termFreq in the next skipping block

    /**
     * Constructor without parameters.
     */
    public SkipInfo()
    {
        this.maxDocId = 0;
        this.docIdOffset = 0;
        this.freqOffset = 0;
    }

    /**
     * Constructor with parameters.
     *
     * @param maxDocId      the maximum DocID in the skipping block
     * @param docIdOffset   the offset of the first docID in the next skipping block
     * @param freqOffset    offset of the first termFreq in the next skipping block
     */
    public SkipInfo(int maxDocId, long docIdOffset, long freqOffset)
    {
        this.maxDocId = maxDocId;
        this.docIdOffset = docIdOffset;
        this.freqOffset = freqOffset;
    }

    public long getMaxDocId() { return maxDocId; }

    public void setMaxDocId(int maxDocId) { this.maxDocId = maxDocId; }

    public long getDocIdOffset() { return docIdOffset; }

    public void setDocIdOffset(long docIdOffset) {
        this.docIdOffset = docIdOffset;
    }

    public long getFreqOffset() {
        return freqOffset;
    }

    public void setFreqOffset(long freqOffset) {
        this.freqOffset = freqOffset;
    }

    @Override
    public String toString() {
        return "SkipInfo{" +
                "maxDocId=" + maxDocId +
                ", docIdOffset=" + docIdOffset +
                ", freqOffset=" + freqOffset +
                '}';
    }

    /**
     * Function to store a skip Info instance into disk.
     *
     * @param skipFileChannel   // the channel for the file into disk
     * @throws IOException
     */
    public void storeSkipInfoToDisk(FileChannel skipFileChannel) throws IOException
    {
        ByteBuffer skipPointsBuffer = ByteBuffer.allocate(SKIPPING_INFO_SIZE);
        skipFileChannel.position(skipFileChannel.size());

        skipPointsBuffer.putInt(this.maxDocId);         // write maxDocID
        skipPointsBuffer.putLong(this.docIdOffset);     // write docID offset
        skipPointsBuffer.putLong(this.freqOffset);      // write freq offset

        skipPointsBuffer = ByteBuffer.wrap(skipPointsBuffer.array());   // wrap in buffer

        while(skipPointsBuffer.hasRemaining())
            skipFileChannel.write(skipPointsBuffer);    // write in the file (SKIP_FILE = skipInfo)
    }

    /**
     * Function to read a skip Info instance from disk.
     *
     * @param start             //
     * @param skipFileChannel   // the channel for the file into disk
     * @throws IOException
     */
    public void readSkipInfoFromDisk(long start, FileChannel skipFileChannel) throws IOException
    {
        ByteBuffer skipPointsBuffer = ByteBuffer.allocate(SKIPPING_INFO_SIZE);

        skipFileChannel.position(start);                    // set the position for the read start

        while (skipPointsBuffer.hasRemaining())             // read from the file (skipFile)
            skipFileChannel.read(skipPointsBuffer);

        skipPointsBuffer.rewind();                          // reset to 0 the current position of the buffer

        this.setMaxDocId(skipPointsBuffer.getInt());        // read MaxDocID
        this.setDocIdOffset(skipPointsBuffer.getLong());    // read docID offset
        this.setFreqOffset(skipPointsBuffer.getLong());     // read term freq offset
    }
}
