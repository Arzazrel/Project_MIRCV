package it.unipi.dii.aide.mircv.data_structures;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;

public class SkipList {

    public ArrayList<SkipInfo> points;

    public Iterator<Posting> postingIterator;
    public Iterator<SkipInfo> skipElemIterator;

    private final LinkedHashMap<Integer, SkipInfo> skipPoints;
    private Posting currPosting;

    private double maxTFIDF;
    private double maxBM25;

    public SkipList() {
        this.skipPoints = new LinkedHashMap<>();
        this.postingIterator = null;
        this.skipElemIterator = null;
    }


    //moves sequentially the iterator to the next posting
    public void next() {

        if(postingIterator.hasNext()) {
            currPosting = postingIterator.next();
            return;
        }

        if(skipElemIterator == null || !skipElemIterator.hasNext()) {      //no blocks or finished
            currPosting = null;
            return;
        }

        SkipInfo si = skipElemIterator.next();

//        postingIterator =

    }


    //returns the score of the current posting
    public double score() {
        return 0.0;
    }

    // advances the iterator forward to the next posting with a document identifier greater than or equal to
    // d ⇒ skipping
    public void nextGEQ(long docID) throws IOException {

    }

}
