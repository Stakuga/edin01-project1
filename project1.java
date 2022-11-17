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

    private static ArrayList<Factors> factorlist = new ArrayList<Factors>();
    private static ArrayList<BigInteger> rvalues = new ArrayList<BigInteger>();
    
    //private static BigInteger N = new BigInteger("31741649");
    //private static BigInteger N = new BigInteger("3205837387"); // Factorised with F=2000 and L=F+200
    //private static BigInteger N = new BigInteger("392742364277");
    private static BigInteger N = new BigInteger("145968946107052219367611");
    //private static BigInteger N = new BigInteger("87463");

    private static int F = 1000;
    private static int L = 10;

    public static void test() {
        int[][] m = produceMatrix();
        System.out.println("Matrix produced.");

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
                // now, for each unmarked row we find left and right hand sides
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

                if (gcd.compareTo(BigInteger.ONE) != 0 && gcd.compareTo(N) != 0) {
                    System.out.println(gcd);
                    BigInteger otherFactor = N.divide(gcd);
                    System.out.println(otherFactor);
                    System.exit(0);
                    break;
                }
            }
        }
        System.out.println("nup");
    }

    public static int[][] produceMatrix() {
        ArrayList<BigInteger> primes = producePrimes();

        int[][] matrix = new int[primes.size() + L][primes.size()];

        BigInteger biggestPrime = primes.get(primes.size() - 1);

        HashMap<BigInteger, Integer> primeMap = new HashMap<BigInteger, Integer>(primes.size());
        for (int i = 0; i < primes.size(); i++) {
            primeMap.put(primes.get(i), i);
        }

        int fulfilledRows = 0;

        int j = 1;
        int k = 1;

        double good = 0;
        double bad = 0;

        BigInteger[] tSquaredArray = N.sqrtAndRemainder();
        BigInteger a;
        if (tSquaredArray[1].compareTo(BigInteger.ZERO) == 0) {
            a = tSquaredArray[0];
        }
        else {
            a = N.sqrt().add(BigInteger.ONE);
        }

        // We want to generate L rows; conversely do something to assign a row L times

        // while we fulfilledRows != L;
        while (fulfilledRows != (primes.size() + L)) {
            // do some shit to matrix[fulfilledRows]

            BigInteger fA = a.pow(2).subtract(N);

            Factors factors = new Factors(fA);
            Collection<BigInteger> values = factors.keySet();
            BigInteger[] valuesArray = factors.keySet().toArray(new BigInteger[0]);
            BigInteger biggestFactor = valuesArray[valuesArray.length - 1];
            //System.out.println(factors.toString());
            //System.out.println("Biggest factor is " + biggestFactor);
            //System.out.println("Biggest prime is " + biggestPrime);
            if (!(biggestPrime.compareTo(biggestFactor) == 1 | biggestPrime.compareTo(biggestFactor) == 0)) {
                a = a.add(BigInteger.ONE);
                bad++;
                System.out.println(good/(good + bad) * 100 + " good rows is " + fulfilledRows);
                continue;
            }

            Integer[] exponents = factors.values().toArray(new Integer[0]);
            HashSet<BigInteger> factorbase = new HashSet(primes);

            if (factorbase.containsAll(values)) { // then fucken add the cunt
                good++;
                System.out.println(good/(good + bad) * 100 + " good rows is " + fulfilledRows);
                rvalues.add(a);
                factorlist.add(factors);
                int[] row = new int[primes.size()];
                for(int i = 0; i < exponents.length; i++) { // check this
                    if (exponents[i] % 2 != 0) {
                        row[primeMap.get(values.toArray()[i])] = 1;
                    }
                }
                matrix[fulfilledRows] = row;
                a = a.add(BigInteger.ONE);
            }
            else { // it ought to be skipped
                bad++;
                System.out.println(good/(good + bad) * 100 + " good rows is " + fulfilledRows);
                a = a.add(BigInteger.ONE);
                continue;
            }

            fulfilledRows = fulfilledRows + 1;
        }

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
    private static ArrayList<BigInteger> producePrimes() {
        ArrayList<BigInteger> output = new ArrayList<BigInteger>();
        BigInteger prime = BigInteger.ONE;
        for (int i = 0; i < F; i++) {
            prime = prime.nextProbablePrime();
            if (prime.compareTo(BigInteger.TWO) == 0) { // case of prime == 2
                output.add(prime);
                continue;
            }
            BigInteger pMinusOne = prime.subtract(BigInteger.ONE);
            BigInteger pMinusOneOverTwo = pMinusOne.divide(BigInteger.TWO);
            BigInteger fin = N.modPow(pMinusOneOverTwo, prime);
            if (fin.compareTo(BigInteger.ONE) != 0) {
                continue;
            }
            else {
                output.add(prime);
            }
        }

        return output;
    }

    public static void main(String[] args) {
        test();

        return;
    }
}
