/**
 * Rough implementation of exercise 1 for 
 * EDIN01 Project 1.
 * @author Louis Copland
 * @author Raquel Perez Lopez
 */
import java.math.BigInteger;
import java.math.BigDecimal;
import java.math.RoundingMode;
import org.mathIT.numbers.BigNumbers;

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

    public static void main(String[] args) {
        ex1();
        ex2();
        return;
    }
}
