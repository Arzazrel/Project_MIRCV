package it.unipi.dii.aide.mircv.compression;

import java.sql.Array;
import java.util.ArrayList;

public class VariableBytes {

    //variable bytes representation of integer

    /**
     * fuction to compress the Term Frequency values using Variable Bytes compression
     *
     * @param docidsToCompress             ArrayList of DocIDs of the posting list to compress
     * @return  a byte array with the Variable Bytes compression of the input DocID values
     */
    public static byte[] integersCompression(ArrayList<Integer> docidsToCompress){

        int nDocId = docidsToCompress.size();
        int byteSize = 0;
        int currentIndex;
        int currentByte = 0;

        for(int i = 0; i < nDocId; i++)
            byteSize += ((32 - Integer.numberOfLeadingZeros(docidsToCompress.get(i))) + 6) / 7; //set byte size (removing leading zeros)

        byte[] docidCompressed = new byte[byteSize];

        for(int i = 0; i < nDocId; i++) {
            currentIndex = 0;
            byteSize = 0;

            int docidVal = docidsToCompress.get(i);

            byteSize = ((32 - Integer.numberOfLeadingZeros(docidVal) + 6) / 7); //size of current integer

            //for all the bytes of the current integer
            while (currentIndex < byteSize) {

                byte b = (byte) ((docidVal >> (currentIndex * 7)) & 0x7F); // move to right and AND bit to bit with 0111111

                if (currentIndex != byteSize - 1)
                    b = (byte)(b | (1 << 7)); // Set to 1 the most meaningful bit, if not the last byte

                docidCompressed[currentByte++] = b;
                currentIndex++;

            }
        }
        return docidCompressed;
    }

    /**
     * fuction to compress the Term Frequency values using Unary compression
     *
     * @param docidsToDecompress           array containing the compressed DocID values
     * @return  an ArrayList containing the decompressed DocID values
     */
    public static ArrayList<Integer> integersDecompression(byte[] docidsToDecompress) {
        ArrayList<Integer> result = new ArrayList<>();
        int shift = 0;
        int currentIndex = 0;
        int num = 0;

        //for all the bytes in docidsToDecompress byte array
        while (currentIndex < docidsToDecompress.length) {

            byte currentByte = docidsToDecompress[currentIndex];

            int value = currentByte & 0x7F; // get last 7 bit from byte

            num |= (value << shift); // Combine bits in the result (recompose the value of the bits)
            shift += 7;
            if ((currentByte & 0x80) == 0) {
                // Most meaningful bit is 0, so it's the last byte
                result.add(num);
                shift = 0;
                num = 0;
            }

            currentIndex++;
        }
        return result;
    }


}
