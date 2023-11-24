package it.unipi.dii.aide.mircv.data_structures;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;

public class SkipInfo {

    public static final int SKIPPING_INFO_SIZE = 3 * Long.BYTES;

    private long maxDocId;
    private long docIdOffset;
    private long freqOffset;

    public SkipInfo(long maxDocId, long docIdOffset, long freqOffset) {
        this.maxDocId = maxDocId;
        this.docIdOffset = docIdOffset;
        this.freqOffset = freqOffset;
    }

    public long getMaxDocId() {
        return maxDocId;
    }

    public void setMaxDocId(long maxDocId) {
        this.maxDocId = maxDocId;
    }

    public long getDocIdOffset() {
        return docIdOffset;
    }

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

    public void storeSkipInfoToDisk(FileChannel skipFileChannel) throws IOException {

        ByteBuffer skipPointsBuffer = ByteBuffer.allocate(SKIPPING_INFO_SIZE);
        skipFileChannel.position(skipFileChannel.size());

        skipPointsBuffer.putLong(this.maxDocId);
        skipPointsBuffer.putLong(this.docIdOffset);
        skipPointsBuffer.putLong(this.freqOffset);

        skipPointsBuffer = ByteBuffer.wrap(skipPointsBuffer.array());

        while(skipPointsBuffer.hasRemaining())
            skipFileChannel.write(skipPointsBuffer);
    }

    public void readSkipInfoFromDisk(long start, FileChannel skipFileChannel) throws IOException {
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
