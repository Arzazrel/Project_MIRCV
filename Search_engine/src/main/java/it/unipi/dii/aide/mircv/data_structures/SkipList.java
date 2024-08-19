package it.unipi.dii.aide.mircv.data_structures;

import it.unipi.dii.aide.mircv.QueryProcessor;
import it.unipi.dii.aide.mircv.compression.Unary;
import it.unipi.dii.aide.mircv.compression.VariableBytes;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.CharBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;

import static it.unipi.dii.aide.mircv.data_structures.DataStructureHandler.readCompDIDBlockFromDisk;
import static it.unipi.dii.aide.mircv.data_structures.DataStructureHandler.readCompTFBlockFromDisk;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static it.unipi.dii.aide.mircv.utils.Constants.*;

/**
 * Class to use skipping and implement the DocID search with skipping in te posting list.
 */
public class SkipList
{
    private ArrayList<SkipInfo> points;             // contains the array ok skipInfo points
    private int skipArrLen;                         // len of the skipping points
    private int skipInterval;                       // number of postings in each skipping block
    private int pointsIndex;                        // the index of the current position for the skip array
    private Iterator<SkipInfo> skipElemIterator;    // iterator for the skipping array (points)
    private ArrayList<Posting> currPostList;        // contains the related posting list
    private int postListIndex;                      // the index of the current position for the posting list
    private int totalPostListLen;                   // the len of the total posting list (used when both skip and compression are enabled)

    /**
     * Constructor without parameters.
     */
    public SkipList()
    {
        this.points = null;
        this.skipArrLen = 0;
        this.skipInterval = 0;
        this.pointsIndex = 0;
        this.skipElemIterator = null;
        this.currPostList = null;
        this.postListIndex = 0;
        this.totalPostListLen = 0;
    }

    /**
     * Constructor with parameter that taken a posting list and skip parameter prepare all is need to use the skipping.
     *
     * @param skipOffset    offset of the beginning of the skipInfo Array in the skipInfo file in the disk
     * @param skipArrLen    len of the skip array (equal to the number of skipInfo block)
     * @param postList      whole posting list
     * @param lenPL         the len of the total posting list
     */
    public SkipList(long skipOffset,int skipArrLen, ArrayList<Posting> postList, int lenPL)
    {
        // set the indexes, they will be updated at each iteration next
        pointsIndex = 0;
        postListIndex = 0;

        totalPostListLen = lenPL;                                       // set the len of the whole PL
        currPostList = postList;                                        // get the whole posting list
        skipInterval = SKIP_POINTERS_THRESHOLD;                         // calculate the skip interval
        this.skipArrLen = skipArrLen;                                   // take skipArrLen

        if (Flags.considerSkippingBytes())
            points = getSkipArrayFromDisk(skipOffset, skipArrLen);      // get the array of skipInfo
        else
        {
            points = null;
            //printDebug("---- This posting list has size = " + currPostList.size() + " that is too small for activate the skipping (skipping threshold = " + SKIP_POINTERS_THRESHOLD);
        }

        //if ((currPostList!=null) && (currPostList.size() >= SKIP_POINTERS_THRESHOLD))
        //if (currPostList != null)
        if (Flags.considerSkippingBytes())
        {
            skipElemIterator = points.iterator();                       // initialize the iterator
            if (skipElemIterator.hasNext())                             // set the iterator
                skipElemIterator.next();
                //printDebug("Constructor first skipping block: " + skipElemIterator.next().getMaxDocId());                              // go in first position with iterator
        }
        //printDebug("-- SkipList constructor -- This skipList element:\n" + this);    // debug print
    }

    @Override
    public String toString()
    {
        return "SkipList\n{" +
                "\n Skipping info:" +
                "\n - points (size) = " + (points == null ? " null" : points.size()) +
                "\n - skipInterval = " + skipInterval +
                "\n - skipArrLen = " + skipArrLen +
                "\n Posting list (len) = " + (currPostList == null ? " null" : currPostList.size()) +
                "\n Whole posting list (len)" + totalPostListLen +
                "\n Indexes:" +
                "\n - pointsIndex = " + pointsIndex +
                "\n - postListIndex = " + postListIndex +
                "\n}";
    }

    //-------------------------------------------- start: get and set --------------------------------------------------

    /**
     * Function to return the Skip Info block with index passed as parameter.
     *
     * @param index     the index of the block in the SkipInfo array.
     * @return          return the Skip Info if there is or null otherwise
     */
    public SkipInfo getSkipBlockInfo(int index)
    {
        if ((points != null) && (index < points.size()))
            return points.get(index);
        else
            return null;
    }

    /**
     * Function to read all the SkipInfo block from disk.
     *
     * @param skipOffset    the offset indicates the start of the Skip Info array
     * @param skipArrLen    the len of the SkipInfo array
     * @return
     */
    private ArrayList<SkipInfo> getSkipArrayFromDisk (long skipOffset,int skipArrLen)
    {
        ArrayList<SkipInfo> tempSkipList = new ArrayList<>();   // temporary skip array
        SkipInfo tempSkipInfo;      // temp var for the Skip Info object read
        int index = 0;              // take the current position of the SkipInfo array
        long offset = skipOffset;   // var to indicates the offset of each SkipInfo element in the disk

        //printDebug("-- getSkipArrayFromDisk function(SkipList):");
        try (
                FileChannel channel = new RandomAccessFile(SKIP_FILE, "rw").getChannel()
        ) {
            // get all skipInfo related to the posting list
            while(index < skipArrLen)
            {
                tempSkipInfo = new SkipInfo();                      // create new temp object
                tempSkipInfo.readSkipInfoFromDisk(offset, channel); // read from disk
                tempSkipList.add(tempSkipInfo);                     // add the read skipInfo element into array
                index++;                                            // update index
                offset += SkipInfo.SKIPPING_INFO_SIZE;              // update offset, pointer to the beginning of next block
                //printDebug("---- Read the block in position : " + index + "\n------" + tempSkipInfo);// debug print
            }
        } catch (IOException e)
        {
            e.printStackTrace();
        }

        return tempSkipList;
    }

    public int getSkipArrLen() {  return this.skipArrLen;  }

    public void setCurrPostList(ArrayList<Posting> currPostList)
    {
        this.currPostList = currPostList;
    }

    //--------------------------------------------- end: get and set ---------------------------------------------------

    //--------------------------------- start: function for only skipping enabled --------------------------------------

    /**
     * Function to advance the iterator forward to the next posting with a document identifier greater than or equal to
     * the searched one.
     *
     * @param docID         the DocID to search
     * @param currentPos    the current position in the PL, helps to restrict the research area at each iteration
     * @return              the position in PL of the DocID passed as parameter or the smallest number greater than the
     *                      number searched
     */
    public int nextGEQ(int docID, int currentPos)
    {
        SkipInfo currentSI;     // to contain the current SkipInfo block
        int searchIndex;        // contains the index of the DocID searched or the greater one
        int startBlockPos;      //
        int endBlockPos;        //

        // check the next hop
        if(skipElemIterator == null)        // no blocks or finished
            return postListIndex;

        // check if the posting list length is enough for skipping
        if (currPostList.size() >= SKIP_POINTERS_THRESHOLD)
        {
            currentSI = points.get(pointsIndex);        // initialize the current SkipInfo
            //printDebug("++ IN nextGEQ - 0 iteration -> search DID: " + docID + " skipArrayPosition: " + pointsIndex + " wit maxDID: " + currentSI.getMaxDocId());
            // use skipping to find the searched DocID
            while (currentSI.getMaxDocId() < docID)
            {
                if (!skipElemIterator.hasNext())    // last block and the searched id there isn't in the posting list
                    return currPostList.size();         // return the outbound of 1 position to indicate that the searched DocID there isn't in the posting and it is ended

                currentSI = skipElemIterator.next();    // hop to next position
                pointsIndex++;                          // update the index
                //printDebug("++++ IN nextGEQ -> skipArrayPosition: " + pointsIndex + " with maxDID: " + currentSI.getMaxDocId());
            }

            // the searched DocID is in the current skip block (if there is)
            // the current skip block(if isn't the last block) start in: startPos = pointsIndex * skipInterval and end in: endPos = ((pointsIndex + 1) * skipInterval) - 1
            startBlockPos = pointsIndex * skipInterval;
            postListIndex = max(startBlockPos, currentPos);                         // SEE NOTE 0
            endBlockPos = min((currPostList.size()-1), ((pointsIndex + 1) * skipInterval) - 1);
            //printDebug("++ IN nextGEQ end iteration -> search: " + docID + " and startblock: " + startBlockPos + " postlistIndex: " + postListIndex + " -> effective maxDID: " + currPostList.get(endBlockPos).getDocId());

            searchIndex = booleanSearch(docID, endBlockPos);    // search the index of the searched DocID
            //printDebug("Used skipping -> found position: " + searchIndex + " with DocID: " + currPostList.get(searchIndex).getDocId());
        }
        else        // posting list too small, normal binary search
        {
            searchIndex = booleanSearch(docID, (currPostList.size()-1));    // search the index of the searched DocID
            //printDebug("Used booleanSearch -> found position: " + searchIndex + " with DocID: " + currPostList.get(searchIndex).getDocId());
        }

        return searchIndex;
    }

    /**
     * Function to do the binary search in the block.
     *
     * @param targetDID the searched DocID
     * @param maxPos    the position of the end of the block
     * @return      the index of the searched DID or the next greater one
     */
    private int booleanSearch(int targetDID, int maxPos)
    {
        int startPos = postListIndex;   // set the start position (the current position on the block to avoid the scan of all block is not needed)
        int endPos = maxPos;            // set the last position
        int currentPos;

        while (true)
        {
            if (startPos > endPos)      // end of the research, the searched DID is not in the list
            {
                //printDebug("++ in boolean search ++ did not found.");
                //printDebug("++++ startPos: " + startPos + " adn DID: " + currPostList.get(startPos).getDocId());
                //printDebug("++++ currentPos: " + currentPos + " adn DID: " + currPostList.get(currentPos).getDocId());
                postListIndex = startPos;     // update index for the posting list
                return startPos;  // not found (return the position with the first DocID greater than the searched one)
            }

            currentPos = (startPos + endPos)/2;

            if (currPostList.get(currentPos).getDocId() == targetDID)
            {
                postListIndex = currentPos;     // update index for the posting list
                return currentPos;              // search DocID found
            }
            else if (currPostList.get(currentPos).getDocId() > targetDID)
            {
                endPos = currentPos - 1;
            }
            else
            {
                startPos = currentPos + 1;
            }
        }
    }

    //---------------------------------- end: function for only skipping enabled ---------------------------------------

    //------------------------- start: function for both skipping and compression enabled ------------------------------
    /**
     * Function to advance the iterator forward to the next posting with a document identifier greater than or equal to
     * the searched one. In the case of both skipping and compression enabled.
     * In this case 'currPostList' contains one decompressed block of the posting list. At the creation of the SkipList
     * instance is loaded the first block of the PL. At each iteration of this block if there is a hop in other skip block
     * that block will be decompressed and put in 'currPostList' adn in the related PL in 'QueryProcessor'.
     *
     * @param docID         the DocID to search
     * @param currentPos    the current position in the PL, helps to restrict the research area at each iteration
     * @param plIndex       the index of the PL in the array of PLs (the index of the term in the ordered query term)
     * @param term          the term of the query related to this instance
     * @return              the position in PL of the DocID passed as parameter or the smallest number greater than the
     *                      number searched
     */
    public int nextGEQCompSkip(int docID, int currentPos, int plIndex, String term)
    {
        SkipInfo currentSI;     // to contain the current SkipInfo block
        int searchIndex;        // contains the index of the DocID searched or the greater one
        int endBlockPos;        // the last position of the current block for the binary search
        int startPointsIndex;   //

        // check the next hop
        if( (skipElemIterator == null) || (currPostList == null))   // no blocks or finished
            return postListIndex;

        //printDebug("Cerco DID: " + docID + " e posizione corrente: " + currentPos);
        startPointsIndex = pointsIndex; // set the current pointsIndex before eventual skip
        if (currentPos != SKIP_POINTERS_THRESHOLD)
            postListIndex = currentPos;     // set the current position in the current block of PL
        else
            postListIndex = 0;              // set the current position in the current block of PL
        currentSI = points.get(pointsIndex);        // initialize the current SkipInfo
        //printDebug("++ IN nextGEQ - 0 iteration -> search DID: " + docID + " skipArrayPosition: " + pointsIndex + " wit maxDID: " + currentSI.getMaxDocId());
        // use skipping to find the searched DocID
        while (currentSI.getMaxDocId() < docID)
        {
            postListIndex = 0;     // update the start position, it will be loaded a new block -> start position = 0
            if (!skipElemIterator.hasNext())    // last block and the searched id there isn't in the posting list
                return currPostList.size();         // return the outbound of 1 position to indicate that the searched DocID there isn't in the posting and it is ended

            currentSI = skipElemIterator.next();    // hop to next position
            pointsIndex++;                          // update the index
            //printDebug("++++ IN nextGEQ -> skipArrayPosition: " + pointsIndex + " with maxDID: " + currentSI.getMaxDocId());
        }

        // skip from a skipping block and another
        if (startPointsIndex != pointsIndex)
            readAndAddUncompBlockPL(term, plIndex);

        endBlockPos = min(skipInterval, ( totalPostListLen - (pointsIndex * skipInterval)));    // take the end position
        //printDebug("++ IN nextGEQ end iteration -> search: " + docID + " postlistIndex (startPos): " + postListIndex + " -> effective maxDID: " + currPostList.get(endBlockPos).getDocId());
        searchIndex = booleanSearch(docID, endBlockPos);    // search the index of the searched DocID

        return searchIndex;
    }

    /**
     * Function to load the new uncompressed block set it in 'QueryProcessor'
     *
     * @param term      the term of the query related to this instance
     * @param indexPL   the index of the PL in the array of PLs (the index of the term in the ordered query term)
     */
    private void readAndAddUncompBlockPL(String term, int indexPL)
    {
        byte[] tf;                      // array for the compressed TermFreq list
        byte[] docids;                  // array for the compressed TermFreq list

        try(
                // open complete files to read the postingList
                RandomAccessFile docidFile = new RandomAccessFile(DOCID_FILE, "rw");
                RandomAccessFile termfreqFile = new RandomAccessFile(TERMFREQ_FILE, "rw");
                // FileChannel
                FileChannel docIdChannel = docidFile.getChannel();
                FileChannel termFreqChannel = termfreqFile.getChannel()
        ) {
            DictionaryElem de = QueryProcessor.getDictionary().get(term);
            // get the compressed block
            tf = readCompTFBlockFromDisk(this, pointsIndex,de.getOffsetTermFreq(), de.getTermFreqSize(), de.getSkipArrLen(), termFreqChannel);
            docids = readCompDIDBlockFromDisk(this, pointsIndex, de.getOffsetDocId(), de.getDocIdSize(), de.getSkipArrLen(), docIdChannel);

            int numTFComp = min(SKIP_POINTERS_THRESHOLD, (de.getDf() - (SKIP_POINTERS_THRESHOLD * pointsIndex)));
            if ( (tf == null) || (docids == null) )     // control check
                return;
            // decompress the block
            ArrayList<Integer> uncompressedTf = Unary.integersDecompression(tf, numTFComp);  // decompress term freq
            ArrayList<Integer> uncompressedDocid = VariableBytes.integersDecompression(docids,true);    // decompress DocID
            //printDebug("readAndAddUncompBlockPL for term: " + "'" + term + "' -> Request Block:" + pointsIndex + " with skipArr len: " + de.getSkipArrLen());
            //printDebug("uncompressedTf len: " + uncompressedTf.size() + " uncompressedDocid len: " + uncompressedDocid.size());
            //printDebug("First DID uncompressed is: " + uncompressedDocid.get(0));
            // add the block to the related PL
            currPostList.clear();
            for (int i = 0; i < numTFComp; i++)
            {
                // add the posting to the posting list
                currPostList.add(new Posting(uncompressedDocid.get(i), uncompressedTf.get(i)));
                // SEE NOTE 1
            }
            QueryProcessor.setPLInSkipAndCompPLs(currPostList, indexPL);      // set the PL in 'QueryProcessor'
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
    //-------------------------- end: function for both skipping and compression enabled -------------------------------

    // ------------------------------------------- start: test functions -----------------------------------------------
    /**
     * Function to read the first did and compare the last DID saved with the real of all skipInfo of the SkipList.
     */
    public void testReadAllSkip()
    {
        int i = 0;
        int startBlockPos, endBlockPos;

        // check if the posting list length is enough for skipping
        if (currPostList.size() >= SKIP_POINTERS_THRESHOLD)
        {
            printDebug("testReadAllSkip:");
            while (i < points.size())           // scann all skipping block of the SkipList
            {
                startBlockPos = i * skipInterval;       // calculate the position in the PL of the first DID of the block
                endBlockPos = min(((i + 1) * skipInterval) - 1, currPostList.size() - 1);   // calculate the position in the PL of the last DID of the block
                printDebug("-- skipArrayPosition: " + i + " with first DID: " + currPostList.get(startBlockPos).getDocId() + " with maxDID: " + points.get(i).getMaxDocId() + " and real maxDID: " + currPostList.get(endBlockPos).getDocId());
                i++;
            }
        }
        else
            printDebug("The posting list is too small, the skipping is not used.");
    }

    // -------------------------------------------- end: test functions ------------------------------------------------
}

/*
 * -- NOTE 0 --
 * this way when several consecutive searches will target the same block you will not do the binary search each time on
 * the whole block but from time to time in a smaller portion.
 *
 * -- NOTE 1 --
 * there is no need to add found DIDs to the priority queue of DIDs to be scrolled through during the DAAT because if
 * present in other documents the DIDs will have already been added and if they are not present (since this method is
 * only invoked for nonessential posting lists) they can be skipped because they will not be able to reach the minimum
 * score to pass the top k results threshold.
 *
 */