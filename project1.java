import java.lang.Math;
import java.math.*;
import org.mathIT.numbers.BigNumbers;
import java.util.*;

/**
 * Implementation of exercises 1, 2 and 3 
 * for EDIN01 project 1.
 * Each exercise is run by a corresponding method 
 * that is static to the 'project 1' class.
 * @author Louis Copland
 * @author Raquel Perez Lopez
 */
public class project1 {

    private static ArrayList<BigInteger> factorlist = new ArrayList<BigInteger>();
    private static ArrayList<BigInteger> rvalues = new ArrayList<BigInteger>();
    
    /*
     * N, F and L values are changed or 
     * commented/uncommented as appropriate.
     */

    //private static BigInteger N = new BigInteger("31741649");
    //private static BigInteger N = new BigInteger("3205837387");
    //private static BigInteger N = new BigInteger("392742364277");
    private static BigInteger N = new BigInteger("145968946107052219367611");
    //private static BigInteger N = new BigInteger("92434447339770015548544881401");

    private static int F = 1000;
    private static int L = F + 10;


    /**
     * Finds the square root of the largest possible 25-digit value, 
     * and determines the number of seconds it would take to test 
     * all values up to that value, which is the worst-case scenario. 
     * Final value is output in terms of hours after conversion.
     */
    public static void exercise1() {
        BigInteger ten = new BigInteger("10");
        BigInteger biggestvalue = ten.pow(12);
        BigInteger tenmil = new BigInteger("10000000");
        BigInteger numsecs = biggestvalue.divide(tenmil);
        BigInteger hourconversion = new BigInteger("3600");
        BigInteger numhours = numsecs.divide(hourconversion);
        System.out.println(numhours);
    }

    public static void exercise2() {
        BigInteger biggestvalue = new BigInteger("9999999999999999999999999");
        BigInteger tenmil = new BigInteger("10000000");
        BigInteger valuerooted = biggestvalue.sqrt(); 
        BigDecimal decimalrooted = new BigDecimal(valuerooted);

        // Precompute number of primes using 
        // approximation by Prime Number Theorem
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
    }

    public static void exercise3() {
        gaussian();
    }

    /**
     * Runs Gaussian elimination on the 
     * binary matrix produced by produceMatrix().
     * Solutions to the Gaussian elimination 
     * correspond to relation rows that can be 
     * multiplied together to create perfect squares 
     * of the form x^2 = y^2modN.
     * GCD(y - x, N) is then calculated, and if the 
     * result is unique (i.e. not 1 or N) then a factor 
     * (factor 1) of N is calculated. Factor 2 is also 
     * calculated by dividing N by factor 1, and these 
     * results are printed.
     */
    public static void gaussian() {
        int[][] m = produceMatrix();

        int numrows = L;
        int numcolumns = F;
        int[] markedrows = new int[numrows];
        int[] rowsums = new int[numrows];

        /*
         * This for loop demarcates a Gaussian elimination 
         * method described by the paper 'A Fast Algorithm 
         * for Gaussian Elimination over GF(2)'.
         * https://www.cs.umd.edu/~gasarch/TOPICS/factoring/fastgauss.pdf
         */
        for (int i = 0; i < numcolumns; i++) { 
            for (int j = 0; j < numrows; j++) { 
                if (m[j][i] == 1) { 
                    markedrows[j] = 1;
                    for (int k = 0; k < numcolumns; k++) { 
                        if (m[j][k] == 1 && k != i) { 
                            int[] newcolumn = new int[numrows]; 
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

        /*
         * Determine sums of matrix rows 
         * following initial Gaussian 
         * elimination process.
         * Matrix rows with sum equal to 1 
         * must have only 1 index equal to 1.
         */
        for (int i = 0; i < m.length; i++) {
            int sum = 0;
            for (int j = 0; j < m[i].length; j++) {
                sum = sum + m[i][j];
            }
            rowsums[i] = sum;
        }

        /*
         * The Gaussian elimination process involves 
         * finding 'unmarked' rows, which are then 
         * linearly dependent on some other rows. 
         * These other rows are found, and altogether 
         * these correspond to the relations that can 
         * be multiplied together to produce a relation 
         * of the form x^2 = y^2modN.
         */
        for (int i = 0; i < markedrows.length; i++) {
            ArrayList<Integer> rowstomultiply = new ArrayList<Integer>();
            if (markedrows[i] == 0) {
                rowstomultiply.add(i);
                int[] operatingrow = m[i];
                BigInteger righthand = BigInteger.ONE;
                BigInteger lefthand = BigInteger.ONE;
                
                for (int j = 0; j < operatingrow.length; j++) { 
                    if (operatingrow[j] == 1) {
                        for (int k = 0; k < numrows; k++) { 
                            if (m[k][j] == 1 && rowsums[k] == 1) { 
                                rowstomultiply.add(k);
                            }
                        }
                    }
                }

                /*
                 * Determine x (left-hand value) and y 
                 * (right-hand value) for relation of 
                 * form x^2 = y^2modN.
                 */
                for (int j = 0; j < rowstomultiply.size(); j++) {
                    int rownumber = rowstomultiply.get(j);
                    BigInteger rightvalue = factorlist.get(rownumber);
                    righthand = righthand.multiply(rightvalue);
                    BigInteger leftvalue = rvalues.get(rownumber);
                    lefthand = lefthand.multiply(leftvalue);
                }
                righthand = righthand.sqrt().mod(N); 
                lefthand = lefthand.mod(N);

                // Determine difference between x and y.
                BigInteger difference = righthand.subtract(lefthand);

                // Determine GCD.
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

    /**
     * For some specified factor base F 
     * and L, and some given N, outputs a 
     * binary matrix where each row corresponds 
     * to the factor decomposition of some 
     * value r^2modN that is smooth over F.
     * Each index in each row corresponds to 
     * the exponent value of the factor 
     * decomposition, with odd exponents equal 
     * to 1 and even exponents equal to 0.
     * @return A 2D binary matrix.
     */
    public static int[][] produceMatrix() {
        int[][] matrix = new int[L][F];

        int[] primes = producePrimes();
        BigInteger biggestPrime = new BigInteger(Integer.toString(primes[F - 1]));

        int fulfilledRows = 0;

        /*
         * For varying values of j and k, 
         * numbers are generated and then 'sieved' 
         * over the factor base, i.e. accepted numbers 
         * are those that can be factored over the 
         * factor base.
         * L numbers are produced so that a suitable 
         * number of solutions present following 
         * Gaussian elimination.
         */
        for (int k = 1; k <= Integer.MAX_VALUE; k++) {
            if (fulfilledRows == L) {
                break;
            }
            for (int j = 1; j <= k; j++) {
                if (fulfilledRows == L) {
                    break;
                }

                // Produce values of the form r^2modN.
                BigInteger r = routput(j, k, N);
                BigInteger rSquaredModN = r.modPow(BigInteger.TWO, N);
                BigInteger factorTest = rSquaredModN;

                if (rSquaredModN.compareTo(BigInteger.ONE) != 1) {
                    continue;
                }

                // Check if r^2modN is smooth over factor base.
                boolean isSmooth = false; 
                int[] exponentArray = new int[primes.length];
                boolean divisible = true;
                while (divisible == true && isSmooth == false) {
                    for (int a = 0; a < primes.length; a++) {
                        int primeNumber = primes[a];
                        BigInteger bigIntPrime = new BigInteger(Integer.toString(primeNumber));
                        if (factorTest.mod(bigIntPrime).compareTo(BigInteger.ZERO) == 0) {
                            exponentArray[a]++;
                            factorTest = factorTest.divide(bigIntPrime);
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
                    continue;
                }
    
                rvalues.add(r);
                factorlist.add(rSquaredModN);
        
                // Produce binary row corresponding 
                // to exponents.
                int[] row = new int[F];
                for (int a = 0; a < row.length; a++) {
                    row[a] = exponentArray[a] % 2;
                }
        
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
            }
        }

        return matrix;
    }

    /**
     * Produces r value corresponding to 
     * formula r = (k * N) ^ 0.5 + j.
     * @param j An integer corresponding 
     * to j in the formula.
     * @param k An integer corresponding 
     * to k in the formula.
     * @param j A BigInteger corresponding 
     * to N in the formula.
     * @return BigInteger The formula output.
     */
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

    /**
     * Driver function for program.
     */
    public static void main(String[] args) {
        //exercise1();
        //exercise2();
        exercise3();
        return;
    }
}
