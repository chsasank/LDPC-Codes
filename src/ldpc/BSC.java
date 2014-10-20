package ldpc;
import java.math.*;

/** 
 * This is a model for binary symmetric channel. Error occurs with probability 'error_probability'
 * @author Sasank Chilamkurthy
 */
public class BSC {
    double error_probability;
    
    /**
     * Construct a BSC channel with error probability f
     * @param f
     */
    public BSC(double f){
        error_probability = f;
    }
    
    /**
     * Send boolean array x through BSC. 
     * Returns the received array
     * @param x
     * @return
     */
    public boolean[] send(boolean[] x){
        boolean [] x_temp = x;
        for(int i = 0; i<x.length;i++)
            if(Math.random()<error_probability){
                x_temp[i] = !x[i];
            }
        return x_temp;
    }
}
