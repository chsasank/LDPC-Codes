package ldpc;
import java.util.*;
import java.math.*;
import java.io.*;
/**
 * This is decoder interface.
 * First creates a graph with given matrix H.
 * Methods will be like decode(vector) or something.
 * Syndrome type decoding
 * Reference: McKay chapter on ldpc codes
 * @author Sasank Chilamkurthy
 */
public class ldpcDecoder {
    private final ArrayList <ArrayList<Integer>> checks;  //adjacencies for 'check' nodes
    private final ArrayList <ArrayList<Integer>> bits;    //adjacencies for 'bit' nodes
    private final int Max_iterations = 20;
    private boolean[] n_estimate;
    private Table<double[]> r;  //Messages from 'check' to 'bit' nodes
    private Table<double[]> q;  //Messages from 'bit' to 'check' nodes
    private double[] q0;
    private boolean[] z;        //syndrome
    private double f;           //BSC error prabability
    
    /**
     * Initialize decoder with parity check matrix H
     * @param H
     * @param bsc_error_probability
     */
    public ldpcDecoder(boolean[][] H, double bsc_error_probability){
        
        //H is m x n matrix
        int m = H.length;   //rows  
        int n = H[0].length;//columns
        //(m<n)
       
        checks = new ArrayList<>();
        bits = new ArrayList<>();
        r = new Table();
        q = new Table();
        n_estimate = new boolean[n];
        q0 = new double[n];
        z = new boolean[m];
        f = bsc_error_probability;
        
        //Filling adjacency list for check nodes.
        for(int i = 0; i< m; i++){
            ArrayList bitList = new ArrayList<>();
            for(int j = 0; j<n;j++){
                if(H[i][j] == true) bitList.add(j);
            }
            checks.add(bitList);
        }
        
        //Filling adjacency list for bit nodes.
        for(int j = 0;j<n;j++){
            ArrayList checkList = new ArrayList<>();
            for(int i = 0; i<m;i++){
                if(H[i][j] == true) checkList.add(i);
            }
            bits.add(checkList);
        }
    }
    
    /**
     * Initializes decoder with parity check matrix read from a specially formated file.
     * Files should be in format specified here: http://www.inference.phy.cam.ac.uk/mackay/codes/alist.html
     * @param filename
     * @param bsc_error_probability
     * @throws IOException
     */
    public ldpcDecoder(String filename, double bsc_error_probability) throws IOException{
        checks = new ArrayList<>();
        bits = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(filename)))
            {
                String sCurrentLine = br.readLine();
                Scanner scan;
                Scanner sizes = new Scanner(sCurrentLine);
                int n = sizes.nextInt();
                int m = sizes.nextInt();
                
                sCurrentLine = br.readLine();   //ignore max weights
                sCurrentLine = br.readLine();   //ignore individual row weights
                sCurrentLine = br.readLine();   //ignore individual column weights
                
                //Fill bits adajacency
                for(int i = 0; i< n; i++){
                    ArrayList checklist = new ArrayList();
                    sCurrentLine = br.readLine();
                    scan = new Scanner(sCurrentLine);
                    while(scan.hasNextInt()){
                        int t;
                        t = scan.nextInt();
                        if(t>=1) {
                            checklist.add(t-1 );
                        }
                    }
                    bits.add(checklist);
                }
                
                
                //Fill checks adjacency
                for(int i =0 ; i< m; i++){
                    ArrayList bitlist = new ArrayList();
                    sCurrentLine = br.readLine();
                    scan = new Scanner(sCurrentLine);
                    while(scan.hasNextInt()){
                        int t;
                        t = scan.nextInt();
                        if(t>=1){
                            bitlist.add(t -1 );
                        }
                    }
                    checks.add(bitlist);
                }
            }
            catch (IOException e) {
                    e.printStackTrace();
            }
        
        int m = checks.size();
        int n = bits.size();
        r = new Table();
        q = new Table();
        n_estimate = new boolean[n];
        q0 = new double[n];
        z = new boolean[m];
        f = bsc_error_probability;
        
    }

    public boolean[] Decode(boolean[] x){
        if(x.length != bits.size() ) System.out.println("input must be size "+ bits.size());
        
        //computing syndrome
        z = computeSyndrome(x);
        
        boolean[] zeros = new boolean[checks.size()];      
        if(isEqual(z,zeros)) return x; //i,e syndrome is zero.
        
        
        /*
        * initialization Step
        */
        double[] bsc ={1-f,f};
        for(int i = 0; i< bits.size();i++){
            for(int j = 0; j<bits.get(i).size(); j++){
                q.put(bits.get(i).get(j),i, bsc);
            }
        }
        
        for(int iterations = 0; iterations< Max_iterations; iterations++){
            
            horizontalStep();
            verticalStep(); //also updates n_estimate
                  
            boolean[] newSyndrome = computeSyndrome(n_estimate);
           
            //check if syndrome of error pattern estimated is same as syndrome of codeword. 
            //If so add estimated error pattern to received word and return it
            if(isEqual(newSyndrome,z) == true){   
                boolean [] result = new boolean[x.length];
                
                for(int i = 0; i< x.length; i++){
                    result[i] = x[i] ^ n_estimate[i];
                }
                return result;
            }
        }
        
        //Even if doesn't converge in Max_iterations, return something
        boolean [] result = new boolean[x.length];
        for(int i = 0; i< x.length; i++){
                    result[i] = x[i] ^ n_estimate[i];
        }
        return result;
    }
    
    
    /*
    * Horizontal step. Update r. 
    * i.e, send messages to bit nodes from check nodes. Calculate this from messages recieved from bit nodes.
    */   
    private void horizontalStep(){
        for(int m = 0; m< checks.size(); m++){  //Do this for each check node

            int N_m = checks.get(m).size();
            double[][] temp = new double[N_m][2]; //by default this is zeroes.

            //iterate over all possible bit sequences and put probabilities accordingly
            for(int s =0; s< 1<<N_m; s++){

                double product_s = 1;

                //Calculate product so that this need not be computed over and over
                for(int j = 0; j<N_m; j++){
                    int n = checks.get(m).get(j);
                    product_s = product_s* q.get(m,n)[(s>>j)&1];
                }


                //calculate r in to a temperory array
                for(int j = 0; j<N_m; j++){
                    int n = checks.get(m).get(j);
                    int Zm = (z[m] == true)?1:0;
                    if( ((s>>j)&1) == 0){
                        if(Integer.bitCount(s)%2 == Zm){
                            temp[j][(s>>j)&1]+= product_s/q.get(m,n)[0];
                        }
                        else{
                            temp[j][1-((s>>j)&1)]+= product_s/q.get(m,n)[0];  
                        }
                    }
                }
            }

            //Update r
            for(int j = 0; j<N_m; j++){
                int n = checks.get(m).get(j);
                r.put(m,n,temp[j]);
            }

        }
        
    }
    
     /*
     * Vertical step. Update q. i.e, send messages to check nodes from bit nodes
     * Compute these from messages received from checknodes.
     * Also update, n_estimate
     */
    private void verticalStep(){
        for(int n = 0; n<bits.size();n++){
            double product0 = 1;
            double product1 = 1;

            for(int j = 0; j< bits.get(n).size(); j++){
                int m = bits.get(n).get(j);
                product0 = product0* r.get(m,n)[0];
                product1 = product1* r.get(m,n)[1];
            }


            for(int j = 0; j< bits.get(n).size(); j++){
                int m = bits.get(n).get(j);
                double[] temp = {0,0};
                double a = (1-f)*product0/r.get(m,n)[0];
                double b = f*product1/r.get(m,n)[1];
                temp[0] = a/(a+b);
                temp[1] = b/(a+b);

                q.put(m,n, temp);   //update q
            }

            double qn0 = (1-f)*product0;
            double qn1 = f*product1;
            q0[n] = qn0/(qn0+qn1);
            n_estimate[n] = q0[n] < 0.5;    //update n_estimate
        }        
    }
    
    /*
    * Efficient way of computing syndrome based on graph structure
    */
    private boolean[] computeSyndrome(boolean[] x){
        boolean[] Syndrome = new boolean[checks.size()];
        for(int m = 0; m< checks.size(); m++){
            Syndrome[m] = false;
            for(int j = 0; j< checks.get(m).size();j++){
                int n = checks.get(m).get(j);
                Syndrome[m] = Syndrome[m]^x[n];
            }
        }
        return Syndrome;
    }
    
    private boolean isEqual(boolean[] a, boolean[] b){
        if (a.length!=b.length) return false;
        boolean isEqual = true;
        for(int m = 0; m< a.length; m++){
            if(a[m] != b[m]){
                isEqual = false;
                break;
            }
        }
        return isEqual;
    }
    
    /**
     * Illustrates decoding with monitor output
     * @param x
     * @return
     */
    public boolean[] decodeWithDebug(boolean[] x){
        if(x.length != bits.size() ) System.out.println("input must be size"+ bits.size());
        
        //computing syndrome
        boolean isNonzero = false;
        for(int m = 0; m< checks.size(); m++){
            z[m] = false;
            for(int j = 0; j< checks.get(m).size();j++){
                int n = checks.get(m).get(j);
                z[m] = z[m]^x[n];
            }
            isNonzero = isNonzero || z[m];//if one of these is true
        }
        System.out.println("sydrome = "+Arrays.toString(LDPC.toInt(z)));

        
        if(isNonzero == false) return x; //i,e syndrome is zero.
        
        
        //initialization Step
        double[] bsc ={1-f,f};
        //System.out.println(Arrays.toString(bsc));
        for(int i = 0; i< bits.size();i++){
            for(int j = 0; j<bits.get(i).size(); j++){
                q.put(bits.get(i).get(j),i, bsc);
                System.out.println("message from bitnode"+ i +" to checknode"+bits.get(i).get(j)+" - "+Arrays.toString(bsc));
            }
        }
        
        
        for(int iterations = 0; iterations< Max_iterations; iterations++){
            
            /*
            * Horizontal step. Update r from messages to check nodes.
            */
            for(int m = 0; m< checks.size(); m++){

                int N_m = checks.get(m).size();
                double[][] temp = new double[N_m][2]; //by default this is zeroes.

                //iterate over all possible bit sequences and put probabilities accordingly
                for(int s =0; s< 1<<N_m; s++){

                    double product_s = 1;

                    //Calculate product so that this need not be computer over and over
                    for(int j = 0; j<N_m; j++){
                        int n = checks.get(m).get(j);
                        product_s = product_s* q.get(m,n)[(s>>j)&1];
                    }


                    //calculate r in temp
                    for(int j = 0; j<N_m; j++){
                        int n = checks.get(m).get(j);
                        int Zm = (z[m] == true)?1:0;
                        if( ((s>>j)&1) == 0){
                            if(Integer.bitCount(s)%2 == Zm){
                                temp[j][(s>>j)&1]+= product_s/q.get(m,n)[0];
                            }
                            else{
                                temp[j][1-((s>>j)&1)]+= product_s/q.get(m,n)[0];  
                            }
                        }
                    }
                }
                
                //Update r
                for(int j = 0; j<N_m; j++){
                    int n = checks.get(m).get(j);
                    r.put(m,n,temp[j]);
                    System.out.println("message from checknode "+m +" to bit node"+n+" - "+Arrays.toString(temp[j]));
                }
               
            }
            
            /*
            * Vertical step. Update q from messages to bit nodes.
            */
            for(int n = 0; n<bits.size();n++){
                double product0 = 1;
                double product1 = 1;
                
                for(int j = 0; j< bits.get(n).size(); j++){
                    int m = bits.get(n).get(j);
                    product0 = product0* r.get(m,n)[0];
                    product1 = product1* r.get(m,n)[1];
                }
                
                
                for(int j = 0; j< bits.get(n).size(); j++){
                    int m = bits.get(n).get(j);
                    double[] temp = {0,0};
                    double a = (1-f)*product0/r.get(m,n)[0];
                    double b = f*product1/r.get(m,n)[1];
                    temp[0] = a/(a+b);
                    temp[1] = b/(a+b);
                   
                    q.put(m,n, temp);   //update q
                    System.out.println("message from bitnode"+ n +" to checknode "+m+" - "+Arrays.toString(temp));
                }
                
                double qn0 = (1-f)*product0;
                double qn1 = f*product1;
                q0[n] = qn0/(qn0+qn1);
                n_estimate[n] = q0[n] < 0.5;    //update n_estimate
            }
            // Uncomment for illustration
            System.out.println("Current Belief that each bit is 0 - "+Arrays.toString(q0));
            System.out.println(Arrays.toString(n_estimate));
            
            /*
            * checking if our estimate is a codeword.
            * Below is just a better way of multiplying with H using graph formed earlier.
            */
            boolean[] newSyndrome = new boolean[checks.size()];
            boolean isEqual = true;
            for(int m = 0; m< checks.size(); m++){
                newSyndrome[m] = false;
                for(int j = 0; j< checks.get(m).size();j++){
                    int n = checks.get(m).get(j);
                    newSyndrome[m] = newSyndrome[m]^n_estimate[n];
                }
                if(newSyndrome[m] != z[m]){
                    isEqual = false;
                    break;
                }             //no need to compute further if any one bit of syndrome is not same as Z.
            }
            System.out.println("New syndrome with this estimate - " + Arrays.toString(newSyndrome));

            if(isEqual == true){   //i,e syndrome matches
                boolean [] result = new boolean[x.length];
                for(int i = 0; i< x.length; i++){
                    result[i] = x[i] ^ n_estimate[i];
                }
                return result;
            } 
            
            System.out.println("Next iteration - " + iterations);
        }
        
        System.out.println("Maximum iterations reached. Decode fail " + Max_iterations);
        System.out.println("Noise estimate" + Arrays.toString(LDPC.toInt(n_estimate)));
        boolean [] result = new boolean[x.length];
        for(int i = 0; i< x.length; i++){
                    result[i] = x[i] ^ n_estimate[i];
        }
        return result;
    }
     
}
