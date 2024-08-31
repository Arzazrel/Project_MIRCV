package it.unipi.dii.aide.mircv.data_structures;

/**
 * Class defining a posting (unit element of a posting list).
 */
public class Posting
{
    private final int docId;    // DocID (recommended delta code compression)
    private int termFreq;       // frequency of the term in the document(recommended unary code compression)

    /**
     * Create a posting with a specified document ID. Constructor with parameter.
     *
     * @param docId The document ID associated with the posting.
     * @param termFreq the occurrence of term in the document associated with the posting.
     */
    public Posting(int docId, int termFreq) {
        this.docId = docId;
        this.termFreq = termFreq;
    }

    /**
     * Function that increment the term frequency of a term.
     *
     * @param n quantity to be added
     */
    public void addTermFreq(int n){
        this.termFreq += n;
    }

    // ---- start method get and set ----

    public int getDocId() { return docId; }

    public int getTermFreq() { return termFreq; }

    @Override
    public String toString() {
        return "Posting{" +
                "docId=" + docId +
                ", termFreq=" + termFreq +
                '}';
    }
}