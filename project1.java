/**
 * Rough implementation of exercise 3 for 
 * EDIN01 Project 1.
 * @author Louis Copland
 * @author Raquel Perez Lopez
 */
import java.lang.Math;
import java.math.BigInteger;
import org.mathIT.numbers.BigNumbers;
import org.mathIT.numbers.Factors;
import java.util.*;
import Jama.*;

public class project1 {

    private static ArrayList<BigInteger> factorlist = new ArrayList<BigInteger>();
    private static ArrayList<BigInteger> rvalues = new ArrayList<BigInteger>();
    
    //private static BigInteger N = new BigInteger("31741649");
    //private static BigInteger N = new BigInteger("3205837387"); // Factorised with F=2000 and L=F+200
    //private static BigInteger N = new BigInteger("392742364277");
    private static BigInteger N = new BigInteger("145968946107052219367611");

    private static int F = 1000;
    private static int L = F + 50;

    public static void test() {
        int[][] m = produceMatrix();

        int numrows = L;
        int numcolumns = F;
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
                    BigInteger rightvalue = factorlist.get(rownumber);
                    righthand = righthand.multiply(rightvalue);
                    BigInteger leftvalue = rvalues.get(rownumber);
                    lefthand = lefthand.multiply(leftvalue);
                }
                righthand = righthand.sqrt().mod(N); // check this
                lefthand = lefthand.mod(N);

                // now, subtract left from right hand side

                BigInteger difference = righthand.subtract(lefthand);

                // find gcd
                BigInteger gcd = difference.gcd(N);

                if (gcd.compareTo(BigInteger.ONE) != 0 && gcd.compareTo(N) != 0) {
                    System.out.println(gcd);
                    BigInteger otherFactor = N.divide(gcd);
                    System.out.println(otherFactor);
                    break;
                }
            }
        }

    }

    public static int[][] produceMatrix() {
        int[][] matrix = new int[L][F];

        int[] primes = producePrimes();
        BigInteger biggestPrime = new BigInteger(Integer.toString(primes[F - 1]));

        HashMap<Integer, Integer> primeMap = new HashMap<Integer, Integer>(F);
        for (int i = 0; i < F; i++) {
            primeMap.put(primes[i], i);
        }

        int fulfilledRows = 0;



        //test values
        int bad = 0;
        double good = 0.;

        // We want to generate L rows; conversely do something to assign a row L times

        for (int k = 1; k <= Integer.MAX_VALUE; k++) {
            if (fulfilledRows == L) {
                break;
            }
            for (int j = 1; j <= k; j++) {
                if (fulfilledRows == L) {
                    break;
                }


                BigInteger r = routput(j, k, N);
    
                BigInteger rSquaredModN = r.modPow(BigInteger.TWO, N);
    
                BigInteger factorTest = rSquaredModN;
                if (rSquaredModN.compareTo(BigInteger.ONE) != 1) {
                    continue;
                }
                //System.out.println("Testing new r value");
                //System.out.println(rSquaredModN);
                boolean isSmooth = false; // check if number is smooth
    
                int[] exponentArray = new int[primes.length];
                boolean divisible = true;
                while (divisible == true && isSmooth == false) {
                    for (int a = 0; a < primes.length; a++) {
                        //System.out.println("Next factor is " + primes[a]);
                        //System.out.println(factorTest);
                        int primeNumber = primes[a];
                        BigInteger bigIntPrime = new BigInteger(Integer.toString(primeNumber));
                        if (factorTest.mod(bigIntPrime).compareTo(BigInteger.ZERO) == 0) {
                            exponentArray[a]++;
                            //System.out.println("incremenetd exponent array" + a);
                            factorTest = factorTest.divide(bigIntPrime);
                            //System.out.println(factorTest);
                        }
                        if (factorTest.compareTo(BigInteger.ONE) == 0) {
                            isSmooth = true;
                            break;
                        }
                        if (a == (primes.length - 1) && factorTest.mod(bigIntPrime).compareTo(BigInteger.ZERO) != 0) {
                            divisible = false;
                        }
                    }
                }
    
                if (isSmooth == false) {
                    //bad++;
                    //System.out.println("Bad row");
                    //System.out.println(bad);
                    continue;
                }
    
                rvalues.add(r);
                
                /*Factors factors = new Factors(rSquaredModN);
                System.out.println(factors.toString());
                System.out.println(biggestPrime);
                BigInteger[] values = factors.keySet().toArray(new BigInteger[0]);
                Integer[] exponents = factors.values().toArray(new Integer[0]);
    
                if (values[values.length - 1].compareTo(biggestPrime) == 1) {
                    bad = bad + 1.;
                    continue;
                }
                else {
                    good = good + 1.;
                }
    
                factorlist.add(factors);*/
    
                factorlist.add(rSquaredModN);
        
                int[] row = new int[F];
                for (int a = 0; a < row.length; a++) {
                    row[a] = exponentArray[a] % 2;
                }
        
                /*for(int i = 0; i < exponents.length; i++) {
                    if (exponents[i] % 2 != 0) {
                        row[primeMap.get(values[i].intValue())] = 1;
                    }
                }*/
                boolean isSame = true;
                for (int c = 0; c < matrix.length; c++) {
                    int[] matrixRow = matrix[c];
                    int rowSum = 0;
                    int matSum = 0;
                    for (int d = 0; d < row.length; d++) {
                        matSum = matSum + matrixRow[d];
                        rowSum = rowSum + row[d];
                        if (matrixRow[d] != row[d]) {
                            isSame = false;
                        }
                    }
                }
                if (isSame) {
                    continue;
                }
    
                matrix[fulfilledRows] = row;
    
                fulfilledRows = fulfilledRows + 1;
                System.out.println("Row successfully made!");
                System.out.println(fulfilledRows);
            }
        }
        // while we fulfilledRows != L;

        return matrix;
    }

    private static BigInteger routput(int j, int k, BigInteger N) {
        String jString = Integer.toString(j);
        BigInteger jBigInt = new BigInteger(jString);

        String kString = Integer.toString(k);
        BigInteger kBigInt = new BigInteger(kString);

        BigInteger kN = kBigInt.multiply(N);
        BigInteger kNsqrt = kN.sqrt();
        BigInteger output = kNsqrt.add(jBigInt);

        return output;
    }

    /**
     * Produces the first F primes.
     * @return int[] array representing factorbase.
     */
    private static int[] producePrimes() {
        int[] output = new int[F];
        BigInteger prime = BigInteger.ONE;
        for (int i = 0; i < F; i++) {
            prime = prime.nextProbablePrime();
            int nextPrime = prime.intValue();
            output[i] = nextPrime;
        }
        return output;
    }

    public static void main(String[] args) {
        test();

        return;
    }
}
