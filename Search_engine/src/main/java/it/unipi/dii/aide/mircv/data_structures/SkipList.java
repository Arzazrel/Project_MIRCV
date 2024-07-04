package it.unipi.dii.aide.mircv.data_structures;

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
    }

    /**
     * Constructor with parameter that taken a posting list and skip parameter prepare all is need to use the skipping.
     *
     * @param skipOffset    offset of the beginning of the skipInfo Array in the skipInfo file in the disk
     * @param skipArrLen    len of the skip array (equal to the number of skipInfo block)
     * @param postList      whole posting list
     */
    public SkipList(long skipOffset,int skipArrLen, ArrayList<Posting> postList, int lenPL)
    {
        // set the indexes, they will be updated at each iteration next
        pointsIndex = 0;
        postListIndex = 0;

        currPostList = postList;                                        // get the whole posting list
        skipInterval = SKIP_POINTERS_THRESHOLD;                         // calculate the skip interval

        this.skipArrLen = skipArrLen;                                   // take skipArrLen e lo assegno alla var

        if (lenPL > SKIP_POINTERS_THRESHOLD)
            points = getSkipArrayFromDisk(skipOffset, skipArrLen);      // get the array of skipInfo
        else
        {
            points = null;
            //printDebug("---- This posting list has size = " + currPostList.size() + " that is too small for activate the skipping (skipping threshold = " + SKIP_POINTERS_THRESHOLD);
        }

        if ((currPostList!=null) && (currPostList.size() >= SKIP_POINTERS_THRESHOLD))
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
                "\n Posting list (size) = " + (currPostList == null ? " null" : currPostList.size()) +
                "\n Indexes:" +
                "\n - pointsIndex = " + pointsIndex +
                "\n - postListIndex = " + postListIndex +
                "\n}";
    }

    /**
     *
     * @param index
     * @return
     */
    public SkipInfo getSkipBlockInfo(int index)
    {
        if ((points != null) && (index < points.size()))
            return points.get(index);
        else
            return null;
    }

    /**
     *
     * @param skipOffset
     * @param skipArrLen
     * @return
     */
    private ArrayList<SkipInfo> getSkipArrayFromDisk (long skipOffset,int skipArrLen)
    {
        ArrayList<SkipInfo> tempSkipList = new ArrayList<>();   // temporary skip array
        int index = 0;          // take the current position of the SkipInfo array
        SkipInfo tempSkipInfo;  //
        long offset;            //

        //printDebug("-- getSkipArrayFromDisk function(SkipList):");
        try (
                FileChannel channel = new RandomAccessFile(SKIP_FILE, "rw").getChannel()
        ) {
            offset = skipOffset;

            // get all skipInfo related to the posting list
            while(index < skipArrLen)
            {
                tempSkipInfo = new SkipInfo();      // create new temp object
                tempSkipInfo.readSkipInfoFromDisk(offset, channel); // read from disk
                tempSkipList.add(tempSkipInfo);     // add the read skipInfo element into array
                index++;                            // update index
                offset += tempSkipInfo.SKIPPING_INFO_SIZE;
                //printDebug("---- Read the block in position : " + index + "\n------" + tempSkipInfo);// debug print
            }
        } catch (IOException e)
        {
            e.printStackTrace();
        }

        return tempSkipList;
    }

    /**
     * Function to advance the iterator forward to the next posting with a document identifier greater than or equal to
     * the searched one.
     *
     * @param docID
     * @param currentPos
     * @return
     */
    public int nextGEQ(int docID, int currentPos)
    {
        int searchIndex = -1;    //
        SkipInfo currentSI;
        int startBlockPos = 0;
        int endBlockPos = currPostList.size()-1;

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
        int startPos = postListIndex;
        int endPos = maxPos;
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


    // ------------------------------------------- start: test functions -----------------------------------------------
    /**
     *
     */
    public void testReadAllSkip()
    {
        int i = 0;
        int startBlockPos, endBlockPos;

        // check if the posting list length is enough for skipping
        if (currPostList.size() >= SKIP_POINTERS_THRESHOLD)
        {
            printDebug("testReadAllSkip:");
            while (i < points.size())
            {
                startBlockPos = i * skipInterval;
                endBlockPos = min(((i + 1) * skipInterval) - 1, currPostList.size() - 1);
                printDebug("-- skipArrayPosition: " + i + " with maxDID: " + points.get(i).getMaxDocId() + " and real maxDID: " + currPostList.get(endBlockPos).getDocId());
                i++;
            }
        }
        else
            printDebug("The posting list is too small, the skipping is not used.");
    }
}

/*
 * -- NOTE 0 --
 * this way when several consecutive searches will target the same block you will not do the binary search each time on
 * the whole block but from time to time in a smaller portion.
 *
 *
 *
 */