package it.unipi.dii.aide.mircv.data_structures;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;

public class SkipInfo {

    public static final int SKIPPING_INFO_SIZE = 3 * Long.BYTES;    // the size

    private long maxDocId;      // the maximum DocID in the skipping block
    private long docIdOffset;   // offset of the first docID in the next skipping block
    private long freqOffset;    // offset of the first termFreq in the next skipping block

    /**
     *
     *
     * @param maxDocId
     * @param docIdOffset
     * @param freqOffset
     */
    public SkipInfo(long maxDocId, long docIdOffset, long freqOffset)
    {
        this.maxDocId = maxDocId;
        this.docIdOffset = docIdOffset;
        this.freqOffset = freqOffset;
    }

    public long getMaxDocId() { return maxDocId; }

    public void setMaxDocId(long maxDocId) { this.maxDocId = maxDocId; }

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
     *
     * @param skipFileChannel
     * @throws IOException
     */
    public void storeSkipInfoToDisk(FileChannel skipFileChannel) throws IOException
    {
        ByteBuffer skipPointsBuffer = ByteBuffer.allocate(SKIPPING_INFO_SIZE);
        skipFileChannel.position(skipFileChannel.size());

        skipPointsBuffer.putLong(this.maxDocId);        // write maxDocID
        skipPointsBuffer.putLong(this.docIdOffset);     // write docID offset
        skipPointsBuffer.putLong(this.freqOffset);      // write freq offset

        skipPointsBuffer = ByteBuffer.wrap(skipPointsBuffer.array());   // wrap in buffer

        while(skipPointsBuffer.hasRemaining())
            skipFileChannel.write(skipPointsBuffer);    // write in the file (SKIP_FILE = skipInfo)
    }

    /**
     *
     * @param start
     * @param skipFileChannel
     * @throws IOException
     */
    public void readSkipInfoFromDisk(long start, FileChannel skipFileChannel) throws IOException
    {
        ByteBuffer skipPointsBuffer = ByteBuffer.allocate(SKIPPING_INFO_SIZE);

        skipFileChannel.position(start);

        while (skipPointsBuffer.hasRemaining())
            skipFileChannel.read(skipPointsBuffer);

        skipPointsBuffer.rewind();
        this.setMaxDocId(skipPointsBuffer.getLong());
        this.setDocIdOffset(skipPointsBuffer.getLong());
        this.setFreqOffset(skipPointsBuffer.getLong());
    }

}
