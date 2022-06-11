# C Factorization using Quadratic Sieve

Pure C factorizer using self-initialising  **Quadratic Sieve**.

This ~2500 lines project :

 - is imediately compatible with Microsoft Windows, Linux (no one dependancy)
 - is a sunday handyman **command line** factorizer from 0 to 210 ... 260 bits
 - is built so that you can easily use and test the software
 - use its own "big num" library named **cint** 
 - use **[AVL trees](https://en.wikipedia.org/wiki/AVL_tree)** to organize informations
 - use **[Lanczos Block](https://en.wikipedia.org/wiki/Lanczos_algorithm)**, a pure C iterative matrix eigenvalues finder algorithm
 - use **[Pollard's Rho](https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm)** algorithm to answer under 65 bits

# GNU General Public License

- As an undergraduate student, this project is part of my computer science + maths training.
- This software proposition is from Michel Leonard (student at Université de Franche-Comté, Mon, 11 Jul 2022)
- There is no guarantee of any kind on the software
- C code is shared under the terms of the **GNU General Public License**
- The **main mathematical and logical inspiration source** is located at :
    - [http://web.mit.edu/sage/export/flintqs-0.0.20070817/QS.cpp](http://web.mit.edu/sage/export/flintqs-0.0.20070817/QS.cpp) - **GNU General Public License**

This software implementation would have been impossible without "FLINT: Fast Library for Number Theory" maintained by [William Hart](https://github.com/wbhart).

# Usage
You want to use the software, you may have to compile it :
- Ubuntu provide you a C compiler by the command `sudo apt install build-essential`
- On Windows you can install [MinGW](https://winlibs.com/), it will, like as Ubuntu, provide you a C compiler

With **Terminal** on Ubuntu or **Powershell** on Windows you can go to the downloaded directory ([cd](https://en.wikipedia.org/wiki/Cd_(command)) command) :
- Execute `gcc -Wall -pedantic -O3 main.c -o qs`, GCC will create the right executable for your device

The compilation took a few seconds, you can use the software :
- Execute `./qs 1389767379868103466550952369852270458062404848980444227682523391`

The software will show you its answer. To see the help use the `-h` option.

# RSA factorization
|Bits| Command| Took  |
|--|--|--|
| 130 | `./qs 982374584994591973035454918323883152991`  | 0.1 s
| 150 | `./qs 1179676342138800493769972880071793674490763037`  | 0.2 s
| 170 | `./qs 1107814594796361407351721529681773419585636720020897`  | 1 s
| 190 | `./qs 995209482127497644492758995962031762505687007535997155961`  | 3 s
| 210 | `./qs 1143938601578848425045957187554857460827103223569512762813742971`  | 20 s
| 230 | `./qs 1268631359685752304166485771456206118633849920528997030903788793364933`  | 1 min 40 s
| 250 | `./qs 1401811817899460116600945074728583412740519573015376930481203561750251051823`  | 18 min 20 s

The initial software goal was to **factor 200-bit RSA in 30 seconds**, after which there were fewer situations tested.

Only a few dozen of RSA numbers greater than 230-bit have been tested.

# Fermat numbers factorization
|F| Value | Took  |
|--|--|--|
| 7 | `340282366920938463463374607431768211457`  | 150 ms
| 8 | `115792089237316195423570985008687907853269984665640564039457584007913129639937`  | 28 minutes

Tests like software development were made by laptop Honor MagicBook on Windows (64-bit).\
Largest number factored during development was the 79 digits 8th Fermat Number.


# Mersenne numbers factorization

|Number| Decimal digits | Fully factored in  |
|--|--|--|
| 2^259 - 1 | 78  | 8 min
| 2^360 - 1 | 109  | 3 s
| 2^468 - 1 | 141  | 2 min 30 s

- Mersenne numbers were factored by the trial division algorithm and then completed by the quadratic sieve
- the quadratic sieve for conveinance reject inputs greater than 220-bit, so the option `-limit=230` was used

# Other factorizations

The software is intended as a generalist factorization solution.

## Random

Tens of thousands of random factorization tests have been performed using:
- Single core Debian (Linux) on [OVH-cloud-amd64](https://www.ovhcloud.com/en/vps/)
- PHP7.4 (time measure + output analysis)
- [GMP](https://www.php.net/manual/en/book.gmp.php) (input proposition + answer verification)

I thanks [William Hart](https://github.com/wbhart), a factor has always been shown ... an extract of tested numbers is provided in text files.

## Primes

- quotes around a number means the software know it's not a prime
- the usual answers contain prime numbers

# Testing

Not needed at all for regular use, a basic ~ 100 lines testing feature is available :
- the goal is simply to give you a little background
- test durations are between 30 seconds (130-bit) and 3 minutes (200-bit)
- test start with `./qs test=1` for a crescendo, or `./qs test=160` for 160-bit

# cint

You have access to the source code of **cint**, the lite "big num" library **cint** can :
- take input from an regular integer
- take input from an "arbitrary" long string, in any base from 2 to 62
- perform basic math operations, including nth_root, modular_inverse, is_prime
- output its content as a string in any base from 2 to 62

To perform intermediates computations **cint** does not use global variable, excepting rand it's thread-safe.

# AVL trees

You have access to the source code of **AVL**, a fast self-balancing binary search tree utility  :

- used to store relations and answers
- uses a **height field** and a **parent pointer** for each tree node
- a worst case is "less complex" than `1.44 * log2(number of entries)`
- duplicate keys aren't possible in the tree, and keys remain sorted
- this software implementation use **cint** as tree entries (or keys)
- tested (and fast) with many more keys than quadratic sieve uses

**cint** and **AVL** are initially separate projects from the quadratic sieve.

# The file main.c

This file contains the factor function, the entry point of any factorization :
- just 70 lines, so you can easily check the software structure
- `inner_continuation_condition` function decide if sieving continue or stop
- `outer_continuation_condition` function decide **after data analysis** if algorithm sieve again or returns

# The file fac_utils.c

The file contains utilities that are not specifically intended for a quadratic sieve :
- `c_factor`, the front end factorizer that format the solution
- `log_computation` 
- `multiplication_modulo` 
- `power_modulo` 
- `kronecker` symbol
- `tonelli_shanks`
- `modular_inverse` 
- `mem_aligned` 

# The file fac_quadratic.c

The file contains the quadratic sieve algorithm, and the C99 code it needs.

### preparation_part_2 .. 3

The input **N is** duplicated and called "**kN**" after this function complete :
- multiply **N** by a prime to reach 120 bits
- apply a [KS](https://books.google.fr/books?id=qRt6CwAAQBAJ&lpg=PA328&ots=ryHIirZZQ2&dq=Knuth%20Schroeppel%20analysis&hl=fr&pg=PA328#v=onepage&q=Knuth%20Schroeppel%20analysis&f=false) multiplier to **N**, intended to optimize runtime

|Variable| Information|
|--|--|
| N | Prime factors are removed from **N** and **N** is updated until `N = 1`  |
| kN | Algorithm computes with **kN** which always remains a constant |

### qs_parametrize, preparation_part_4
- define the algorithm parameters, you may try to ameliore
- allocates a block of memory for the quadratic sieve computations
- prepare constants, variables, buffers, data arrays ...

There is a struct inside the **qs_sheet** (or manager) called **mem** :
- it holds the  **base** entry point of the malloced memory
- it holds a **now** void* pointer which is the available memory pointer

**qs_sheet** holds 2 AVL tree manager :
- one to store the relations
- one to store the known divisors of N

At any moment with small precautions you are supposed to be able to store someting in **now**, then to update **now** accordingly to what you stored.  **now** is always supposed to contain only zero until its end, this is the main allocation technique used by this software.  `mem_align` aims to offset any pointer to a nice pointer, wasting a few bits if necessary.

### preparation_part_5 .. 6
Fill the manager's **base** array with prime numbers provided [by](https://stackoverflow.com/a/61895974/18765627) a fast constant expression :
- verify that kN is a square mod prime, ignore it otherwise
- associates the prime with square root of kN mod prime
- associates the prime with its size (log2)
- computes invariants like **D** used to generate the **A** polynomial coefficient

### iteration_analyzer
- normally do nothing (very large majority of tested cases, and all < 220-bit cases)
- may help you if you try to better configure the software (helped me around 230-bit)
- can unblock sieving by randomizing **D**, an important variable

*I have identified that in some rare cases after 220-bit, sieving was blocked (no more relations found) due to a 65-bit result stored into a 64-bit variable. I have identified that because i cleared the 64-th bit of all variables during a test, that's how I was able to generate this case on command.* This function resolved that case and the software "always" completed its execution, a 128-bit support may do it too.

### Polynomial `AX^2 + 2BX + C` coefficients : 

| coefficient | const qualified after  | const qualified until |
|--|--|--|
| A | `iteration_part_2` |`inner_continuation_condition` completion
| B | `iteration_part_5` |for loop final expression
| C | `iteration_part_7` |for loop final expression

Polynomial and related data are computed until `iteration_part_7` complete.
`iteration_part_8` and `iteration_part_9` are used for sieving.

### Search sieve for relations
-  `register_relations` defines the variable **X** , computes the polynomial value in **X** then searches sieve for relations. Data are buffered and sent to `register_relation_kind_1` or `register_relation_kind_2`.

Knowledge is then structured by `register_relations_commons` :
- informations potentially useful are saved into a struct **qs_relation***
- `register_relation_kind_1` immediately build a matrix of **qs_relation***, using AVL tree
- `register_relation_kind_2` combine **qs_relation*** together, using AVL tree

The AVL tree is used as a hashmap.

### inner_continuation_condition
Decides if sieving shoud continue or break.

### lanczos_block
This algorithm is supposed to find matrix eigenvalues.
- the process need memory, **fac_lanczos.c** have its own array builder
- all but answer is zeroed and reusable after computations

### finalization_part_1 .. 2 .. 3

The function exploit `lanczos_block` answer, it use a worker to remove prime factors from N.

In certain cases such as :
```
./qs 51460938795049063955433175628971167803839994111348342302522016010379
```
when N is not fully factored, it tries to complete the factorization with GCDs and perfect power checks.

### outer_continuation_condition
This function is used to decide if the algorithm should return the control, or return to sieving.
- normal case is to return the control, answer was found
- other case is when no answer was found

By default, if no answer was found the algorithm searches 10%, 25% then 50% more relations before giving up.
