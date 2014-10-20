package ldpc;

import java.io.IOException;
import java.util.*;

/**
 * Contains Main method.
 * @author Sasank Chilamkurthy
 */
public class LDPC {

    public static void main(String[] args) throws IOException {
		//Let's do a simulation to check post FEC BER 
        BSC channel = new BSC(0.02);
        ldpcDecoder decoder = new ldpcDecoder("204.33.486.txt", 0.02);


        int numberOferrors = 0;
        boolean[] rec = new boolean[204];
        for(int i = 0; i< 1000; i++){
            
            boolean[] msg = new boolean[204]; 
            boolean[] sent = channel.send(msg);//send 0 vecor
            
            rec = decoder.Decode(sent);
            
            for(int j = 0; j<204; j++){
                if(rec[j] ==true){
                    numberOferrors = numberOferrors+1;
                    break;
                }
            }
        }
        
        System.out.println("Bit error rate after decoding = "+ ((double)numberOferrors)/204000);
        System.out.println("Compare this to ber before decoding = "+ 0.02);
        
    }
    
    public static boolean[][] toBoolean(int[][] H){
        boolean [][] Hb = new boolean[H.length][H[0].length];
        for(int i = 0; i< H.length; i++){
            for(int j = 0; j< H[0].length; j++){
                Hb[i][j] = (H[i][j]!=0);
            }
        }
        return Hb;
    }
    
    public static boolean[] toBoolean(int[] H){
        boolean [] Hb = new boolean[H.length];
        for(int i = 0; i< H.length; i++){
                Hb[i] = (H[i]!=0);
        }
        return Hb;
    }
    
    public static int[] toInt(boolean[] H){
        int[] Hi = new int[H.length];
        for(int i = 0; i< H.length; i++){
                Hi[i] = (H[i]==false)?0:1;
        }
        return Hi;
    }
}
