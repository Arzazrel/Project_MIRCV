package it.unipi.dii.aide.mircv.compression;

import java.util.ArrayList;

/**
 * Class to the compression of positive integer using VariableBytes code.
 * VariableBytes code represents an integer x > 0 as VB(x) in a variable number of bytes. In each byte the first bit is
 * a control bit, set to '1' if is the last byte of the number or set '0' otherwise (there are other encode bytes)
 *
 * Es. U(1) = 10000001 , U(5) = 10000101 , U(824) = 00000110 10111000
 * SEE NOTE 0 at the end of the class
 */
public class VariableBytes
{
    /**
     * Function to compress the DocID values using Variable Bytes compression.
     *
     * @param docIDsToCompress  ArrayList of DocIDs of the posting list to compress
     * @param dGaps             if 'true' use dGaps , if 'false' not use dGaps. When dGaps the first integer of the list
     *                          is compressed normally while the other integers will be compressed the difference with
     *                          the previous integer.
     * @return  a byte array with the Variable Bytes compression of the input DocID values
     */
    public static byte[] integersCompression(ArrayList<Integer> docIDsToCompress, boolean dGaps)
    {
        int nDocId = docIDsToCompress.size();   // get the number of integer to compress
        int previousDID = 0;        // value of previous DID, used when dGaps is enabled
        int byteSize = 0;           // indicates the number of bits for representing the current number
        int numToRoundUp = 6;       // the number to round up the number of bytes needed for each integer to compress
        int currentIndex;           // indicates the bytes index for the current integer to compress
        int currDocIDVal;           // the value of the current DocID
        int currentByte = 0;        // indicates the bytes, of compressed integer, saved
        int valueMask = 0x7F;       // mask to obtain the value bits from compressed byte (01111111)
        int currDiff = 0;
        byte b;

        // numberOfLeadingZeros(num) -> returns the total number of '0' bits preceding the first '1' bit
        // scan all integer to compress and compute the total number of bytes requested for the compression (number of bits / 7) one bit for each bytes will be for control
        if (dGaps)      // consider the gap between the raw integer
        {
            for (Integer iDsToCompress : docIDsToCompress)
            {
                byteSize += ((32 - Integer.numberOfLeadingZeros((iDsToCompress - previousDID))) + numToRoundUp) / 7;   //set byte size (removing leading zeros)
                previousDID = iDsToCompress;        // update previous DID
            }
        }
        else            // consider the raw integer
            for (Integer iDsToCompress : docIDsToCompress)
                byteSize += ((32 - Integer.numberOfLeadingZeros(iDsToCompress)) + numToRoundUp) / 7;   //set byte size (removing leading zeros)


        byte[] docIDCompressed = new byte[byteSize];    // initialize the array for the compressed list

        previousDID = 0;            // reset previous DID,
        for(int i = 0; i < nDocId; i++)     // scan all list of integer to compress
        {
            currentIndex = 0;                           // reset index
            currDocIDVal = docIDsToCompress.get(i);     // get the DocID value
            if (dGaps)      // consider the gap between the raw integer
            {
                currDiff = currDocIDVal - previousDID;
                byteSize = ((32 - Integer.numberOfLeadingZeros(currDiff) + numToRoundUp) / 7); // bytes size for the compressed form
                previousDID = currDocIDVal;                 // update previous DID
            }
            else            // consider the raw integer
                byteSize = ((32 - Integer.numberOfLeadingZeros(currDocIDVal) + numToRoundUp) / 7); // bytes size for the compressed form

            //for all the bytes of the compressed form of the current integer
            while (currentIndex < byteSize)
            {
                if (dGaps)      // consider the gap between the raw integer
                    b = (byte) ((currDiff >> ((byteSize - currentIndex - 1) * 7)) & valueMask); // move to right and AND bit to bit with 0111111
                else
                    b = (byte) ((currDocIDVal >> ((byteSize - currentIndex - 1) * 7)) & valueMask); // move to right and AND bit to bit with 0111111

                if (currentIndex != byteSize - 1)   // it is not the last byte
                    b = (byte)(b & valueMask);          // Set to '0' the most meaningful bit (AND 01111111 = 127)
                else                                // it is the last byte
                    b = (byte)(b | (1 << 7));           // Set to '1' the most meaningful bit (OR 10000000)

                docIDCompressed[currentByte++] = b;     // save the bytes in the list of compressed integer
                currentIndex++;                         // increment bytes index
            }
        }

        return docIDCompressed;     // return compressed list
    }

    /**
     * Function to decompress the DocID values using Variable Bytes compression.
     *
     * @param docIDsToDecompress    array containing the compressed DocID values
     * @param dGaps                 if 'true' use dGaps , if 'false' not use dGaps. When dGaps the first integer of the
     *                              list is decompressed normally while the other compressed integers will be the
     *                              difference with the previous integer. To obtain the raw value it's needed to sum
     *                              with the previous raw value.
     * @return  an ArrayList containing the decompressed DocID values
     */
    public static ArrayList<Integer> integersDecompression(byte[] docIDsToDecompress, boolean dGaps)
    {
        ArrayList<Integer> decompressedList = new ArrayList<>();  // list of decompressed integer
        int shift = 0;          // is the number of bits to shift in the current iteration to decompress the current compressed integer
        int currentIndex = 0;   // index for the list of compressed integer
        int num = 0;            // the current decompressed integer
        int valueMask = 0x7F;   // mask to obtain the value bits from compressed byte (01111111)
        int currentValue;       // the value of the current compressed byte

        //for all the bytes in docIDsToDecompress byte array
        while (currentIndex < docIDsToDecompress.length)
        {
            byte currentByte = docIDsToDecompress[currentIndex];    // get the current byte to decompress
            currentValue = currentByte & valueMask;                 // get last 7 bit from byte

            num = num << shift;         // update the position of the bits, the first time for each integer is 0 << 0
            num |= currentValue;        // Combine bits in the result (recompose the value of the bits)
            shift = 7;                  // update shift, the next byte read will be not the first, shift to 7 bits
            if ((currentByte & (-128)) == -128)     // (and 10000000) most meaningful bit is '1', so it's the last byte
            {
                if (dGaps && !decompressedList.isEmpty())      // consider the gap between the raw integer
                    num += decompressedList.get((decompressedList.size() - 1));   // add the previous value

                decompressedList.add(num);  // add the decompressed number
                shift = 0;                  // reset shift value
                num = 0;                    // reset value for decompressed integer
            }
            currentIndex++;                 // increment index for the list of compressed integer
        }

        return decompressedList;
    }

    /**
     * Function to compress the DocID values using Variable Bytes compression.
     * This method "updates" the number on which to start the difference every 'sizeBlock' positions.
     * That is, at each 'sizeBlock' elements saves the number in raw compressed form and for the next number will
     * be saved the difference with it.
     *
     * @param docIDsToCompress  ArrayList of DocIDs of the posting list to compress
     * @param sizeBlock         Size of the blocks for the DGaps. At each sizeBlock integer compressed
     * @return  a byte array with the Variable Bytes compression of the input DocID values
     */
    public static byte[] intCompDGapsBlock(ArrayList<Integer> docIDsToCompress, int sizeBlock)
    {
        int nDocId = docIDsToCompress.size();   // get the number of integer to compress
        int previousDID = 0;        // value of previous DID
        int byteSize = 0;           // indicates the number of bits for representing the current number
        int numToRoundUp = 6;       // the number to round up the number of bytes needed for each integer to compress
        int currentIndex;           // indicates the bytes index for the current integer to compress
        int currDocIDVal;           // the value of the current DocID
        int currentByte = 0;        // indicates the bytes, of compressed integer, saved
        int valueMask = 0x7F;       // mask to obtain the value bits from compressed byte (01111111)
        int currDiff = 0;
        int counter = 0;            // the counter for the iteration related to sizeBlock
        byte b;


        // numberOfLeadingZeros(num) -> returns the total number of '0' bits preceding the first '1' bit
        // scan all integer to compress and compute the total number of bytes requested for the compression (number of bits / 7) one bit for each bytes will be for control

        for (Integer iDsToCompress : docIDsToCompress)
        {
            counter++;      // increment counter
            if (counter == sizeBlock)       // update the number for the difference
            {
                previousDID = 0;        // reset the previous
                counter = 0;            // reset the counter
            }
            byteSize += ((32 - Integer.numberOfLeadingZeros((iDsToCompress - previousDID))) + numToRoundUp) / 7;   //set byte size (removing leading zeros)
            previousDID = iDsToCompress;        // update previous DID
        }

        byte[] docIDCompressed = new byte[byteSize];    // initialize the array for the compressed list

        counter = 0;                // reset the counter
        previousDID = 0;            // reset previous DID
        for(int i = 0; i < nDocId; i++)     // scan all list of integer to compress
        {
            currentIndex = 0;                           // reset index
            currDocIDVal = docIDsToCompress.get(i);     // get the DocID value
            counter++;                                  // increment counter
            if (counter == sizeBlock)                   // update the number for the difference
            {
                previousDID = 0;                        // reset the previous
                counter = 0;                            // reset the counter
            }

            currDiff = currDocIDVal - previousDID;
            byteSize = ((32 - Integer.numberOfLeadingZeros(currDiff) + numToRoundUp) / 7); // bytes size for the compressed form
            previousDID = currDocIDVal;                 // update previous DID

            //for all the bytes of the compressed form of the current integer
            while (currentIndex < byteSize)
            {
                b = (byte) ((currDiff >> ((byteSize - currentIndex - 1) * 7)) & valueMask); // move to right and AND bit to bit with 0111111

                if (currentIndex != byteSize - 1)   // it is not the last byte
                    b = (byte)(b & valueMask);          // Set to '0' the most meaningful bit (AND 01111111 = 127)
                else                                // it is the last byte
                    b = (byte)(b | (1 << 7));           // Set to '1' the most meaningful bit (OR 10000000)

                docIDCompressed[currentByte++] = b;     // save the bytes in the list of compressed integer
                currentIndex++;                         // increment bytes index
            }
        }

        return docIDCompressed;     // return compressed list
    }

    /**
     * Function to decompress the DocID values using Variable Bytes compression.
     * This method "updates" the number on which to start the difference every 'sizeBlock' positions.
     * That is, at each 'sizeBlock' elements saves the number in raw compressed form and for the next number will
     * be saved the difference with it.
     *
     * @param docIDsToDecompress    array containing the compressed DocID values
     * @param sizeBlock             Size of the blocks for the DGaps. At each sizeBlock integer compressed
     * @return  an ArrayList containing the decompressed DocID values
     */
    public static ArrayList<Integer> intDecompDGapsBlock(byte[] docIDsToDecompress, int sizeBlock)
    {
        ArrayList<Integer> decompressedList = new ArrayList<>();  // list of decompressed integer
        int shift = 0;          // is the number of bits to shift in the current iteration to decompress the current compressed integer
        int currentIndex = 0;   // index for the list of compressed integer
        int num = 0;            // the current decompressed integer
        int valueMask = 0x7F;   // mask to obtain the value bits from compressed byte (01111111)
        int currentValue;       // the value of the current compressed byte
        int counter = 0;        // the counter for the iteration related to sizeBlock
        int previousDID = 0;    // value of previous DID
        byte currentByte;

        while (currentIndex < docIDsToDecompress.length)        //for all the bytes in docIDsToDecompress byte array
        {
            currentByte = docIDsToDecompress[currentIndex];     // get the current byte to decompress
            currentValue = currentByte & valueMask;             // get last 7 bit from byte


            num = num << shift;         // update the position of the bits, the first time for each integer is 0 << 0
            num |= currentValue;        // Combine bits in the result (recompose the value of the bits)
            shift = 7;                  // update shift, the next byte read will be not the first, shift to 7 bits
            if ((currentByte & (-128)) == -128)     // (and 10000000) most meaningful bit is '1', so it's the last byte
            {
                counter++;                  // increment counter
                if (counter == sizeBlock)           // update the number for the difference
                {
                    previousDID = 0;                // reset the previous
                    counter = 0;                    // reset the counter
                }
                num += previousDID;         // add the previous value

                decompressedList.add(num);  // add the decompressed number
                previousDID = num;          // update previous DID
                shift = 0;                  // reset shift value
                num = 0;                    // reset value for decompressed integer
            }
            currentIndex++;                 // increment index for the list of compressed integer
        }

        return decompressedList;
    }

    // ---- start: utilities functions ----

    /**
     * Function to print the bit of the compressed list passed as parameter.
     *
     * @param byteArray     the list of compressed integer, in Unary code, to be printed.
     */
    public static void printCompressedList(byte[] byteArray)
    {
        for (byte b : byteArray)
            System.out.print(String.format("%8s", Integer.toBinaryString(b & 0xFF)).replace(' ', '0') + " ");

        System.out.println(" ");
    }

    /**
     * Function to print a list of integer of the decompressed list passed as parameter.
     *
     * @param intArray      the list of decompressed integer, in Unary code, to be printed.
     */
    public static void printDecompressedList(ArrayList<Integer> intArray)
    {
        for (int num : intArray)
            System.out.print(num + " ");
        System.out.println(" ");
    }
}

/*
 * -- NOTE 0 --
 * - Statistics of DocID gap in posting list:
 * -- Generic:
 * ---- the average gap between DID of the same posting list: 46747.79260327121
 * ---- the min avg gap between DID of the same posting list: 0.0
 * ---- the max avg gap between DID of the same posting list: 4026500.0
 * -- With partition in block of the posting lists (es. skipping or compression). Skipping enabled false and compression enabled false :
 * ---- the average gap between DID of the same block: 42112.16373854081
 * ---- the min avg gap between DID of the same block: 0.0
 * ---- the max avg gap between DID of the same block: 36616.62857142857
 *
 * in this collection the DocID value go from 1 to 8841702. the type of the DocID is int which is represented on 4 bytes
 * Integer on 4Byte in two's complement can have values from -2147483648(2^-31) to 2147483647(2^31 -1).
 * For the DocID in this collection in the worst case the number can be represented in at most 24 bits (in two's complement).
 */
