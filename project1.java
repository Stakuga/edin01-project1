/**
 * Rough implementation of exercise 1 for 
 * EDIN01 Project 1.
 * @author Louis Copland
 * @author Raquel Perez Lopez
 */
import java.lang.Math;
import java.math.BigInteger;
import java.math.BigDecimal;
import java.math.RoundingMode;
import org.mathIT.numbers.BigNumbers;
import org.mathIT.numbers.Factors;
import java.util.*;
import Jama.*;

public class project1 {

    private static ArrayList<Factors> factorlist = new ArrayList<Factors>();
    private static ArrayList<BigInteger> rvalues = new ArrayList<BigInteger>();

    /**
     * Finds the square root of the largest possible 25-digit value, 
     * and determines the number of seconds it would take to test 
     * all values up to that value, which is the worst-case scenario. 
     * Final value is output in terms of hours after conversion.
     */
    public static void ex1() {
        BigInteger biggestvalue = new BigInteger("9999999999999999999999999");
        BigInteger tenmil = new BigInteger("10000000");
        BigInteger valuerooted = biggestvalue.sqrt(); // Possibly round up rather than truncate
        BigInteger numsecs = valuerooted.divide(tenmil);
        BigInteger hourconversion = new BigInteger("3600");
        BigInteger numhours = numsecs.divide(hourconversion);
        System.out.println(numhours);
    }

    public static void ex2() {
        BigInteger biggestvalue = new BigInteger("9999999999999999999999999");
        BigInteger tenmil = new BigInteger("10000000");
        BigInteger valuerooted = biggestvalue.sqrt(); // Possibly round up rather than truncate
        BigDecimal decimalrooted = new BigDecimal(valuerooted);

        // Precompute number of primes using approximation by Prime Number Theorem
        BigDecimal lnx = BigNumbers.ln(decimalrooted);
        RoundingMode rm1;
        rm1 = RoundingMode.valueOf("UP");
        int precision = 20;
        BigDecimal approx = decimalrooted.divide(lnx, precision, rm1);
        BigInteger numofprimes = approx.toBigInteger();
        
        // Run calculation 
        BigInteger numsecs = numofprimes.divide(tenmil);
        BigInteger hourconversion = new BigInteger("3600");
        BigInteger numhours = numsecs.divide(hourconversion);
        System.out.println(numhours);

        // We see that the computational time theoretically would be significantly reduced.
        // But, it would take ages to compute the primes accurately, which could make 
        // the endeavour pointless.
    }

    public static void ex3() {
        BigInteger testint = new BigInteger("592");
        Factors factortree = new Factors(testint);
        //System.out.println(butz);
    }

    public static void biglad() {
        int N = 16637;
        BigInteger bigprime = new BigInteger("29");
        int L = 12;
        int[] primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29};
        HashMap<Integer, Integer> F = new HashMap<Integer, Integer>(10);
        for (int i = 0; i < 10; i++) {
            F.put(primes[i], i);
        }
        
        double[][] M = new double[12][10];


        M[0] = outrow(F, M, N, 2, 3, bigprime);
        M[1] = outrow(F, M, N, 4, 4, bigprime);
        M[2] = outrow(F, M, N, 3, 5, bigprime);
        M[3] = outrow(F, M, N, 4, 5, bigprime);
        M[4] = outrow(F, M, N, 2, 6, bigprime);
        M[5] = outrow(F, M, N, 2, 7, bigprime);
        M[6] = outrow(F, M, N, 6, 10, bigprime);
        M[7] = outrow(F, M, N, 4, 11, bigprime);
        M[8] = outrow(F, M, N, 12, 12, bigprime);
        M[9] = outrow(F, M, N, 4, 13, bigprime);
        M[10] = outrow(F, M, N, 8, 13, bigprime);
        M[11] = outrow(F, M, N, 8, 14, bigprime);

        Matrix em = new Matrix(M);

    }

    public static void test() {
        //int[][] m = {{1,1,0,0},{1,1,0,1},{0,1,1,1},{0,0,1,0},{0,0,0,1}};
        int[][] m = {{1,1,0,1,0,0,1,0,0,0},
        {0,0,0,0,0,1,0,0,0,0},
        {1,0,0,0,1,0,1,0,0,0},
        {0,1,0,1,1,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,1,1},
        {1,0,1,1,0,0,1,0,0,0},
        {0,0,0,0,0,1,1,1,0,0},
        {1,0,0,0,0,0,1,0,0,0},
        {1,0,1,0,1,0,0,0,1,0},
        {1,0,1,0,0,0,0,0,1,0},
        {1,0,0,0,0,0,0,1,0,0},
        {1,1,0,1,0,0,0,1,0,0}};
        int numrows = m.length;
        int numcolumns = m[0].length;
        int[] markedrows = new int[numrows];
        int[] rowsums = new int[numrows];

        for (int i = 0; i < numcolumns; i++) { // for every column index
            for (int j = 0; j < numrows; j++) { // for every row element in said column
                if (m[j][i] == 1) { // if we come across a column element equal to 1
                    markedrows[j] = 1;
                    for (int k = 0; k < numcolumns; k++) { // for each element in the same row
                        if (m[j][k] == 1 && k != i) { // if we have a row element that is not the current column AND is 1
                            // add column i and column k
                            int[] newcolumn = new int[numrows]; // initialise new empty column
                            for (int l = 0; l < numrows; l++) {
                                int number = (m[l][i] + m[l][k]) % 2;
                                newcolumn[l] = number;
                            }
                            for (int l = 0; l < numrows; l++) {
                                m[l][k] = newcolumn[l];
                            }
                        }
                    }
                    break;
                }
            }
        }
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[i].length; j++) {
                System.out.print(m[i][j]);
                System.out.print(" ");
            }
            System.out.println(" ");
        }
        for (int i = 0; i < markedrows.length; i++) {
            System.out.print(markedrows[i]);
            System.out.print(" ");
        }
        System.out.println(" ");

        // make row sums
        for (int i = 0; i < m.length; i++) {
            int sum = 0;
            for (int j = 0; j < m[i].length; j++) {
                sum = sum + m[i][j];
            }
            rowsums[i] = sum;
        }

        for (int i = 0; i < markedrows.length; i++) {
            ArrayList<Integer> rowstomultiply = new ArrayList<Integer>();
            if (markedrows[i] == 0) { // if we find an unmarked row
                rowstomultiply.add(i);
                int[] operatingrow = m[i];
                //System.out.println(i + 1);
                BigInteger righthand = BigInteger.ONE;
                BigInteger lefthand = BigInteger.ONE;
                BigInteger N = new BigInteger("16637");
                
                for (int j = 0; j < operatingrow.length; j++) { // find column with 1
                    if (operatingrow[j] == 1) {
                        for (int k = 0; k < numrows; k++) { // check across rows with same column = to 1
                            if (m[k][j] == 1 && rowsums[k] == 1) { // if its a cool row
                                rowstomultiply.add(k);
                            }
                        }
                    }
                }
                // now, for each marked row we find left and right hand sides
                for (int j = 0; j < rowstomultiply.size(); j++) {
                    int rownumber = rowstomultiply.get(j);
                    BigInteger rightvalue = factorlist.get(rownumber).value();
                    righthand = righthand.multiply(rightvalue);
                    BigInteger leftvalue = rvalues.get(rownumber);
                    lefthand = lefthand.multiply(leftvalue);
                }
                righthand = righthand.sqrt().mod(N);
                lefthand = lefthand.mod(N);

                // now, subtract left from right hand side

                BigInteger difference = righthand.subtract(lefthand);

                // find gcd
                BigInteger gcd = difference.gcd(N);

                if (gcd.compareTo(BigInteger.ONE) != 0) {
                    System.out.println(gcd);
                    BigInteger otherFactor = N.divide(gcd);
                    System.out.println(otherFactor);
                    break;
                }
            }
        }

    }

    private static double[] outrow(HashMap<Integer, Integer> F, double[][] M, int N, int j, int k, BigInteger bigprime) {
        int r = routput(j, k, N);
        double rsquared = Math.pow(r, 2);
        
        int rsquare = (int)rsquared;
        int yee = rsquare % N;
        BigInteger rvalue = new BigInteger(Integer.valueOf(r).toString());
        rvalues.add(rvalue);
        Factors eff = new Factors(yee);
        factorlist.add(eff);

        BigInteger[] values = eff.keySet().toArray(new BigInteger[0]);
        Integer[] exponents = eff.values().toArray(new Integer[0]);

        double[] row = new double[10];

        for(int i = 0; i < exponents.length; i++) {
            if (exponents[i] % 2 != 0) {
                //System.out.println(F.get(values[i].intValue()));
                row[F.get(values[i].intValue())] = 1;
            }
        }

        return row;
    }

    private static int routput(int j, int k, int N) {
        double J = j;
        double K = k;
        double n = N;
        double output = Math.sqrt(K * n) + J;
        int intoutput = (int)output;
        return intoutput;
    }

    private static class Tuple {
        private int j;
        private int k;
        public Tuple(int j, int k) {
            this.j = j;
            this.k = k;
        }
        public int getj() {
            return j;
        }
        public int getk() {
            return k;
        }
    }

    public static void main(String[] args) {
        //ex1();
        //ex2();
        //ex3();
        biglad();
        test();
        return;
    }
}
