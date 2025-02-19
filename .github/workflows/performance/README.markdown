# Factorization Software Performance Logs

This repository contains performance logs for the **Self-Initializing Quadratic Sieve (SIQS)** implemented in C. These logs provide insights into the execution of the factorization algorithm with different input sizes.

## Test Parameters

- **Number of test series:** 10
- **Bit range:** 30-bit to 260-bit numbers
- **Samples per series:** 231
- **Total factored samples:** 2310
- **Seed generation:** Based on the Unix timestamp
- **Verbose mode:** `--verbose 4`

## Purpose of Performance Logs

The logs serve multiple purposes:

1. **Analysis of Factorization Performance**
   - Provides data for assessing execution times and efficiency across different bit sizes.
   
2. **Understanding Verbose Output**
   - Shows detailed data of factorization steps.
   - Includes maintenance messages and progress markers.
   - Beyond 220-bit numbers, displays `.` progress indicators (one per percentage point).

3. **Debugging and Maintenance**
   - Helps diagnose issues and optimize the factorization process.

## Log Format

Each of the 10 log files corresponds to a test series and follows a structured format, looking like :

```
Quadratic Sieve on 823538134720364217987066388519176747269489.
Suggested multipliers are [1, 41, 29, 11, 65, 14, 89].
Suggested multipliers are [1, 41, 9, 29, 11, 113, 89], so use 1.
N is a 140-bit number, and kN is a 140-bit number using 5 words.
The algorithm use the seed 0 and targets 1702 relations.
The iterative list contains 8, 448, 896 and 1792.
Other params include sieve=31744, error_bits=16, threshold=64 and s=6.
The underlying calculations use 372-bit variables.
Allocated 1 MB of memory with a 260 KB structure.
The factor base of 1792 suitable primes ends with 32189.
Linear algebra passed on 1st attempt after 28 iterations.
- New divisor 96678460721497163 shown.
- New divisor 164394747597064169 shown.
      - Performs GCD within the divisors to get 364097879.
      - This prime factor reduces N to 111-bit.
    - Performs GCD within the divisors to get 451512511.
    - This prime factor reduces N to 83-bit.
    - Prunes the divisors to get 214121333.
    - This prime factor reduces N to 55-bit.
- New divisor 81912161 shown.
  - This prime factor reduces N to 29-bit.
  - And allows us for N to get 285619237.
  - The factorization is complete since it's a prime.
The algorithm completed with 312 polynomials and 1719 relations.
Used 0 MB of memory during 0.20 second(s).
Number: 823538134720364217987066388519176747269489
Factors: 81912161 (prime), 214121333 (prime), 285619237 (prime), 364097879 (prime), 451512511 (prime)
Status: fully factored in 0.21 s
```

