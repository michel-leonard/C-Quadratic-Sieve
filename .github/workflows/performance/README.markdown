# Quadratic Sieve Factorization Performance Logs

This repository contains performance logs for the **Quadratic Sieve Factorization** implemented in C. These logs provide insights into the execution of the factorization algorithm with different input sizes and settings.

## Test Parameters

- **Number of test series:** 10
- **Bit range:** 30-bit to 260-bit numbers
- **Samples per series:** 231
- **Total factored samples:** 2310
- **Seed generation:** Based on the millisecond timestamp at test start
- **Verbose mode:** `--verbose 4`

## Purpose of Performance Logs

The logs serve multiple purposes:

1. **Analysis of Factorization Performance**
   - Provides data for assessing execution times and efficiency across different bit sizes.
   
2. **Understanding Verbose Output**
   - Shows detailed breakdowns of factorization steps.
   - Includes maintenance messages and progress markers.
   - Beyond 220-bit numbers, displays `.` progress indicators (one per percentage point).

3. **Debugging and Maintenance**
   - Helps diagnose issues and optimize the factorization process.

## Log Format

Each of the 10 log files corresponds to a test series and follows a structured format, looking like :

```
Quadratic Sieve on 65563611653881364695825618063205274348178793774070863314867727808937.
Suggested multipliers are [1, 73, 57, 7, 13, 65, 17].
Suggested multipliers are [1, 73, 7, 57, 25, 13, 17], so use 1.
N is a 226-bit number, and kN is a 226-bit number using 8 words.
The algorithm use the seed 0 and targets 14336 relations.
The iterative list contains 16, 1280, 2560 and 5120.
The single large-prime variation is being processed under 5488112.
Other params include sieve=95232, error_bits=25, threshold=85 and s=9.
The underlying calculations use 558-bit variables.
Allocated 44 MB of memory with a 1915 KB structure.
The factor base of 14336 suitable primes ends with 334099.
..................................................
..................................................
Linear algebra passed on 1st attempt after 221 iterations.
- New divisor 75317872199314295268989 shown.
- New divisor 5942669066932965286905944776009 shown.
    - Performs GCD within the divisors to get 78901181.
    - This prime factor reduces N to 200-bit.
- New divisor 25543141994241849528169 shown.
        - Notes a perfect power to get 22619059.
        - This prime factor reduces N to 126-bit.
    - Performs GCD within the divisors to get 26792707.
    - This prime factor reduces N to 102-bit.
- New divisor 109057505098129307567131 shown.
    - Performs GCD within the divisors to get 38794853.
    - This prime factor reduces N to 76-bit.
  - Divides N to get 24574439.
  - This prime factor reduces N to 52-bit.
- New divisor 65562439 shown.
  - This prime factor reduces N to 26-bit.
  - And allows us for N to get 42877193.
  - The factorization is complete since it's a prime.
The sieve reported 35259 partials which added 5637 smooth relations.
The algorithm completed with 47296 polynomials and 14347 relations.
Used 23 MB of memory during 134.90 second(s).
Number: 65563611653881364695825618063205274348178793774070863314867727808937
Factors: 22619059^3, 24574439 (prime), 26792707 (prime), 38794853 (prime), 42877193 (prime), 65562439 (prime), 78901181 (prime)
Status: fully factored in 135.30 s
```

