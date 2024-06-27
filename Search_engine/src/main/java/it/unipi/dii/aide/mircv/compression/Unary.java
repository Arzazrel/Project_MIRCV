package it.unipi.dii.aide.mircv.compression;

import it.unipi.dii.aide.mircv.data_structures.Posting;

import java.io.IOException;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.ArrayList;

import static it.unipi.dii.aide.mircv.utils.Constants.printDebug;
import static it.unipi.dii.aide.mircv.utils.Constants.verbose;

/**
 * Class to the compression of positive integer using Unary code.
 * Unary code represents an integer x > 0 as U(x) = 1^(x-1)0 , i.e. a series of ones (of length x-1) followed by a 0.
 * Es. U(1) = 0 , U(2) = 10 , U(3) = 110 , and so on...
 * SEE NOTE 0 at the end of the class
 */
public class Unary
{
    /**
     * Function to compress the Term Frequency values using Unary compression code.
     *
     * @param termFreqToCompress             ArrayList of Term Frequency of the posting list
     * @return  a byte array with the Unary compression of the input Term Frequency values
     */
    public static byte[] integersCompression(ArrayList<Integer> termFreqToCompress)
    {
        // get the number of Bytes of the compressed form of the integer list passed as parameter.
        int numBytes = computeByteOfCompressedList(termFreqToCompress);
        byte[] compressedResult = new byte[numBytes];   // array of bytes for the compressed integer list
        int byteToWrite = 0;        // index for the bytes to write
        int bitPosition = 0;        // index for the current bit to write
        compressedResult[byteToWrite] = 0;  // set to 0

        for (int num : termFreqToCompress)      // scan all integer of the list
        {
            //System.out.println("num: " + num);
            for (int j = 0; j < num - 1; j++)   // insert (num -1) '1' and then '0'
            {
                //set correspondent bit to 1 using OR bit a bit and left shift. SEE NOTE 1
                compressedResult[byteToWrite] = (byte) (compressedResult[byteToWrite] | (1 << (7 - bitPosition)));
                bitPosition++;          // update the bit position to write
                // check if all bits in the byte have been written
                if (bitPosition == 8)
                {
                    byteToWrite++;      // go to next byte
                    compressedResult[byteToWrite] = 0;  // set the new byte to 0 in all bits
                    bitPosition = 0;     // reset bit value
                }
            }
            bitPosition++;      // write the 0 to the end of series of ones to represent the num
            // check if all bits in the byte have been written
            if (bitPosition == 8)
            {
                byteToWrite++;      // go to next byte
                compressedResult[byteToWrite] = 0;  // set the new byte to 0 in all bits
                bitPosition = 0;     // reset bit value
            }
        }
        /*
        for (byte b : compressedResult) {
           System.out.print(Integer.toBinaryString(b & 0xFF) + " ");
        }*/
        //unaryToInt(compressedResult);
        //integerArrayDecompression(compressedResult, 3);
        printDebug("Unary compression -- the list passed len: " + termFreqToCompress.size() + " int with size: " + (termFreqToCompress.size()*4) + " Bytes -> after compression the size is: " + numBytes + " Bytes.");
        printCompressedList(compressedResult);      // debug print
        return compressedResult;
    }

    /**
     * fuction to compress the Term Frequency values using Unary compression
     *
     * @param compressedArray           array containing the compressed Term Frequency values
     * @param totNum                    total number of integers to decompress
     * @return  an ArrayList containing the decompressed Term Frequency values
     */
    public static ArrayList<Integer> integersDecompression(byte[] compressedArray, int totNum)
    {
        ArrayList<Integer> decompressedList = new ArrayList<>();    // contains the list of the decompressed integer
        int currentBit = 7;         // index for the current bit to read
        int currentValue = 1;       // indicates the current value of the integer that is being decompressed
        int nIntegers = 0;          // number of integer decompressed

        // scan all compressed bytes of the list
        // take the current byte of the list
        for (byte currentByteValue : compressedArray)
        {
            while (true)
            {
                // Read current bit
                int bit = (currentByteValue >> currentBit) & 1; // take the value of the (8 - currentBit)-th bit

                if (bit == 1)   // If bit is 1, increment current value
                {
                    currentValue++;     // increment the current value (the decompressed value will be = 1 + # of 1
                } else            // If bit is 0, add current value to decompressed list (compression end of the current integer)
                {
                    decompressedList.add(currentValue);     // add the decompressed integer value into list

                    currentValue = 1;               // Reset current value
                    nIntegers++;                    // increment the number of decompressed integer
                    if (nIntegers == totNum)        // check if the current decompressed value is the last one
                    {
                        printDebug("Unary Decompression -- the compressed list len: " + totNum + " int with size: " + compressedArray.length + " Bytes -> after decompression the size is: " + (decompressedList.size()*4) + " Bytes.");
                        printDecompressedList(decompressedList);      // debug print
                        return decompressedList;        // end of the decompression
                    }
                }

                currentBit--;           // update the position of bit to read

                if (currentBit < 0)     // check if all byte read
                {
                    currentBit = 7;     // reset the position to read (first bit of the next byte)
                    break;              // go to the next byte to read
                }
            }
        }

        return decompressedList;
    }

    // ---- start: utilities functions ----

    /**
     * Function to compute the number of bytes needed for the compressed form (in unary code) of a list of integer passed
     * as parameter.
     *
     * @param termFreqToCompress    the list of integer to compress
     * @return  the number of bytes needed for the compressed form (in Unary code) of the list passed as parameter. This
     *          value is the len in bytes of the compressed integer list.
     */
    private static int computeByteOfCompressedList (ArrayList<Integer> termFreqToCompress)
    {
        int numBits = 0;    // indicates the number of bits of the compressed form of the integer list passed as parameter.
        int numBytes = 0;   // indicates the number of Bytes of the compressed form of the integer list passed as parameter.

        for (Integer freqTerm : termFreqToCompress)   // scan all integer int he list
            numBits += freqTerm;    // update the number of beets needed. In Unary #bit of compressed form = integer value

        numBytes = (int) Math.ceil((numBits / 8)) + (numBits % 8 == 0? 0 : 1);    // compute the number of byte needed

        return numBytes;
    }

    /**
     * Function to print the bit of the compressed list passed as parameter.
     *
     * @param byteArray     the list of compressed integer, in Unary code, to be printed.
     */
    private static void printCompressedList(byte[] byteArray)
    {
        for (byte b : byteArray)
            System.out.print(String.format("%8s", Integer.toBinaryString(b & 0xFF)).replace(' ', '0') + " ");
            //System.out.print(String.format("%8s", Integer.toBinaryString(b & 0xFF)).replace(' ', '0') + " ");
        System.out.println(" ");
    }

    /**
     * Function to print a list of integer of the decompressed list passed as parameter.
     *
     * @param intArray      the list of decompressed integer, in Unary code, to be printed.
     */
    private static void printDecompressedList(ArrayList<Integer> intArray)
    {
        for (int num : intArray)
            System.out.print(num + " ");
        System.out.println(" ");
    }

    // ---- end: utilities functions ----

}

/*
 * -- NOTE 0 --
 * The Unary code is very well for small integer. Below is shown a list of the top 20 values for TermFreq (from the
 * posting lists) with their relative occurrence in raw form and percentage of the collection used for this examination.
 *
 * Top 20 Term Frequency occurrence list (from the most common to the rarest):
 * -- pos 1 -> TermFreqValue: 1 , occurrence: 245012338 (row value) , occurrence: 73.98%
 * -- pos 2 -> TermFreqValue: 2 , occurrence: 57270188 (row value) , occurrence: 17.29%
 * -- pos 3 -> TermFreqValue: 3 , occurrence: 14645204 (row value) , occurrence: 4.42%
 * -- pos 4 -> TermFreqValue: 4 , occurrence: 7707034 (row value) , occurrence: 2.33%
 * -- pos 5 -> TermFreqValue: 5 , occurrence: 2805138 (row value) , occurrence: 0.85%
 * -- pos 6 -> TermFreqValue: 6 , occurrence: 1762460 (row value) , occurrence: 0.53%
 * -- pos 7 -> TermFreqValue: 7 , occurrence: 741176 (row value) , occurrence: 0.22%
 * -- pos 8 -> TermFreqValue: 8 , occurrence: 522791 (row value) , occurrence: 0.16%
 * -- pos 9 -> TermFreqValue: 9 , occurrence: 238715 (row value) , occurrence: 0.07%
 * -- pos 10 -> TermFreqValue: 10 , occurrence: 179306 (row value) , occurrence: 0.05414301501810548%
 * -- pos 11 -> TermFreqValue: 11 , occurrence: 88686 (row value) , occurrence: 0.026779513401089212%
 * -- pos 12 -> TermFreqValue: 12 , occurrence: 71628 (row value) , occurrence: 0.021628701101563022%
 * -- pos 13 -> TermFreqValue: 13 , occurrence: 36448 (row value) , occurrence: 0.01100579239612678%
 * -- pos 14 -> TermFreqValue: 14 , occurrence: 30443 (row value) , occurrence: 0.009192530122785545%
 * -- pos 15 -> TermFreqValue: 15 , occurrence: 17050 (row value) , occurrence: 0.00514839662955338%
 * -- pos 16 -> TermFreqValue: 16 , occurrence: 14018 (row value) , occurrence: 0.004232857709858023%
 * -- pos 17 -> TermFreqValue: 17 , occurrence: 7860 (row value) , occurrence: 0.002373395748286778%
 * -- pos 18 -> TermFreqValue: 18 , occurrence: 6522 (row value) , occurrence: 0.0019693749453341433%
 * -- pos 19 -> TermFreqValue: 19 , occurrence: 3778 (row value) , occurrence: 0.0011408001446599806%
 * -- pos 20 -> TermFreqValue: 20 , occurrence: 3168 (row value) , occurrence: 9.566053092331442E-4%
 *
 * An integer is on 4Bytes (32bits) as we can see above the most common term frequency values are the smaller one.
 * In particular the value '1' is the most common with more than 70% of the total occurrence.
 * In unary code one is the value with the smallest length in the compressed form and is the most common.
 * The second common value is the second smallest value in unary and so on, the unary code is perfect for the composition
 * of term frequency in this collection.
 * With this
 *
 * -- NOTE 1 --
 * Execute the or (bit-wise) between the value of the bit of the current byte and the left shift of 1 of (7-bitPosition)
 * position. 1 is equal 0000 0001 and the left shift of q from (7-bitPosition) is equal to write 1 in the bitPosition bit.
 * Than the or is done for ad the 1 in the target position without modifying the bit values in the other position.
 */