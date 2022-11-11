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

public class project1 {

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
        
        int[][] M = new int[12][10];

        int r = routput(3, 5, N);
        double rsquared = Math.pow(r, 2);

        int rsquare = (int)rsquared;
        int yee = rsquare % N;
        
        Factors eff = new Factors(yee);

        if (eff.lastKey().compareTo(bigprime) == 1) {
            System.out.println("Ah fuck");
        }
        else {
            System.out.println(eff.keySet());
        }

        BigInteger[] values = eff.keySet().toArray(new BigInteger[0]);

        for(int i = 0; i < values.length; i++) {
            if (F.containsKey(values[i].intValue())) {
                System.out.println(F.get(values[i].intValue()));
                M[0][F.get(values[i].intValue())] = 1;
            }
        }
        System.out.println(" ");
        for (int i = 0; i < M[0].length; i++) {
            System.out.println(M[0][i]);
        }


        // IN ELSE CLAUSE:
        // produce matrix row all set to 0
        // for each key:
        //      if in F hashmap
        //      get position and set thing in empty matrix row to 1

        int[] fuck = new int[10];


    }

    private static int routput(int j, int k, int N) {
        double J = j;
        double K = k;
        double n = N;
        double output = Math.sqrt(K * n) + J;
        int intoutput = (int)output;
        return intoutput;
    }

    public static void main(String[] args) {
        //ex1();
        //ex2();
        //ex3();
        biglad();
        return;
    }
}
