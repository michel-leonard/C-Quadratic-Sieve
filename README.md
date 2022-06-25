
# C Factorization using Quadratic Sieve

Pure C factorizer using self-initialising  **Quadratic Sieve**.

This ~2500 lines project :

 - is imediately compatible with Microsoft Windows, Linux (no one dependancy)
 - is a C99 **command line** factorizer from 0 to 300 bits
 - is built so that you can easily use and test the software
 - use its own "big num" library named **cint** 
 - use **[AVL trees](https://en.wikipedia.org/wiki/AVL_tree)** to organize informations
 - use **[Lanczos Block](https://en.wikipedia.org/wiki/Lanczos_algorithm)**, a pure C iterative matrix eigenvalues finder algorithm
 - use **[Pollard's Rho](https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm)** algorithm to answer under 64 bits

# GNU General Public License

- As an undergraduate student, this project is part of my computer science + maths training.
- This software proposition is from Michel Leonard (student at Université de Franche-Comté, Mon, 11 Jul 2022)
- There is no guarantee of any kind on the software
- C code is shared under the terms of the **GNU General Public License**
- The **main mathematical and logical inspiration source** is located at :
    - [http://web.mit.edu/sage/export/flintqs-0.0.20070817/QS.cpp](http://web.mit.edu/sage/export/flintqs-0.0.20070817/QS.cpp) - **GNU General Public License**

This software implementation would have been impossible without "FLINT: Fast Library for Number Theory" maintained by [William Hart](https://github.com/wbhart).

# Usage
If you don't know how to get an executable, try to follow this fast procedure :
- Ubuntu provide you a C compiler by the command `sudo apt install build-essential`
- On Windows you can install [MinGW](https://winlibs.com/), it will, like as Ubuntu, provide you a C compiler

With **Terminal** on Ubuntu or **Powershell** on Windows you can go to the downloaded directory ([cd](https://en.wikipedia.org/wiki/Cd_(command)) command) :
- Execute `gcc -Wall -pedantic -O3 main.c -o qs`, GCC will create the right executable for your device

The compilation took a few seconds, you can use the software :
- Execute `./qs 1389767379868103466550952369852270458062404848980444227682523391`

The software will show you its answer. To see the help use the `-h` option.

# RSA factorization

Factoring RSA numbers from 65 to 130 bits takes about the same time as 130 bits.

|Bits| Command| Took
|--|--|--|
| 130 | `./qs 982374584994591973035454918323883152991`  | 0.1 s
| 150 | `./qs 1179676342138800493769972880071793674490763037`  | 0.2 s
| 170 | `./qs 1107814594796361407351721529681773419585636720020897`  | 1 s
| 190 | `./qs 995209482127497644492758995962031762505687007535997155961`  | 3 s
| 210 | `./qs 1143938601578848425045957187554857460827103223569512762813742971`  | 15 s
| 230 | `./qs 1268631359685752304166485771456206118633849920528997030903788793364933`  | 45 s
| 250 | `./qs 1401811817899460116600945074728583412740519573015376930481203561750251051823`  | 4 min
| 270 | `./qs 1415606447884291776606783262139201189953436249643759632827004228713595295320953939`  | 13 min
| 290 | `./qs 1311212776762431067307390080218400290276233806447008972197996326428961782072889375956641`  | 1 h

The initial software goal was to **factor 200-bit RSA in 30 seconds**, after which there were fewer situations tested.
- Some RSA numbers larger than 250 bits have been tested, they have been factorized.
- During a stress test that lasted 2 hours, the software factored a **300-bit RSA** number.

Aditionally, the software factored the **321-bit** RSA number relating to the "[bank card case](https://www.enseignement.polytechnique.fr/profs/informatique/Eric.Goubault/Cours09/qs.pdf)" during a technical test that laster 12 hours.
# Fermat numbers factorization
|F| Value | Took  |
|--|--|--|
| 7 | `340282366920938463463374607431768211457`  | 150 ms
| 8 | `115792089237316195423570985008687907853269984665640564039457584007913129639937`  | 3 min

These tests like software development were made by laptop Honor MagicBook on Windows (64-bit).\
One of the largest number factored during development was the 79 digits 8th Fermat Number.

# Mersenne numbers factorization

|Number| Decimal digits | Fully factored in  |
|--|--|--|
| 2^259 - 1 | 78  | 1 min 40 s
| 2^360 - 1 | 109  | 3 s
| 2^468 - 1 | 141  | 40 s

- Mersenne numbers were factored by the trial division algorithm and then completed by the quadratic sieve
- the quadratic sieve for conveinance reject inputs greater than 220-bit, so the option `-limit=230` was used

# Other factorizations



The software is intended as a generalist factorization solution.

## Random

Tens of thousands of random factorization tests have been performed using:
- Single core Debian (Linux) on [OVH-cloud-amd64](https://www.ovhcloud.com/en/vps/)
- PHP7.4 (time measure + output analysis)
- [GMP](https://www.php.net/manual/en/book.gmp.php) (input proposition + answer verification)

I again thanks [William Hart](https://github.com/wbhart), a factor has always been shown ... an extract of tested numbers is provided in text files.

## Primes

- the usual answers contain prime numbers
- quotes around a number means the software know it's not a prime

# Testing

Not needed at all for regular use, a basic ~ 100 lines testing feature is available :
- the goal is simply to provide inputs when users test new configurations
- test durations are between 30 seconds (130-bit) and 3 minutes (200-bit)
- the `./qs test=1` test offers a 2 minute crescendo up to 200+ bits
- the `./qs test=160` test offers a minute of 160-bit random odd numbers

The time measurements provided in this page are indicative, not all were taken by the same device.

# Memory

Memory allocations are reasonably sized, so this project passes pointers to **assert**.
- the program will stop if the memory is refused, showing you an error message
- it would be recommended to restart your device if you see this kind of message

Technical : [valgrind](https://valgrind.org/) shows around **10MB** and **80MB** allocated depending on the SIQS input size.

# cint

You have access to the source code of **cint**, the lite "big num" library **cint** is designed to :
- take input from a regular integer
- take input from an "arbitrary" long string in base 2 to 62
- perform basic math operations, including nth_root, modular_inverse, is_prime
- output its content as a string in base from 2 to 62

To perform intermediates computations **cint** does not use global variable, excepting rand it's thread-safe.

# AVL trees

You have access to the source code of **AVL**, a fast self-balancing binary search tree utility  :

- used to store relations and answers
- uses a **height field** and a **parent pointer** for each tree node
- tree can store and retrieve keys in worst case `1.44 * log2 (number of entries)`
- this software implementation use **cint** as tree entries (or keys)
- tested and fast with many more keys than quadratic sieve uses

**cint** and **AVL** are initially separate projects from the quadratic sieve.

# Improvements

Developers can for example replace **cint** with GMP, **AVL** with a hashmap and better reconfigure the software.

# The file main.c

This file has just some lines, so you can easily perform a refactoring.

# The file fac_utils.c

The file contains utilities that are not specifically intended for a quadratic sieve :
- `c_factor`, the front-end factorizer that format the solution
- `pollard_rho` 
- `is_prime_1062961`
- `log_computation` 
- `multiplication_modulo` 
- `power_modulo` 
- `kronecker_symbol`
- `tonelli_shanks`
- `modular_inverse` 
- `mem_aligned` 

# The file fac_quadratic.c

The quadratic sieve file structure is as follows:
- the ~40 lines function that allows to see the algorithm **structure**
- the 2 important loop **conditions**
- the algorithm **parameters**
- the **functions** approximately in the order they are called

### preparation_part_2 .. 3

The input **N is** duplicated and called "**kN**" after this function complete :
- multiply **N** by a prime to reach 120 bits
- apply a [Knuth-Schroeppel](https://books.google.fr/books?id=qRt6CwAAQBAJ&lpg=PA328&ots=ryHIirZZQ2&dq=Knuth%20Schroeppel%20analysis&hl=fr&pg=PA328#v=onepage&q=Knuth%20Schroeppel%20analysis&f=false) multiplier to **N**, intended to optimize runtime

|Variable| Information|
|--|--|
| N | Prime factors are removed from **N** and **N** is updated until `N = 1`  |
| kN | Algorithm computes with **kN** which always remains a constant |

*A **Knuth-Schroeppel** multiplier make the algorithm completes up to 2.5 times faster, 1.45 times on average.*

*The quadratic sieve could factor number fewer than 120-bit without using two multipliers, this possibly interesting challenge is left to any keen developer.*

### qs_parametrize, preparation_part_4

*Good parameters can improve the speed*.

- define the algorithm parameters
- allocates a block of memory for the quadratic sieve computations
- prepare constants, variables, buffers, data arrays ...

There is a struct inside the **qs_sheet** (or manager) called **mem** :
- it holds the  **base** entry point of the allocated memory
- it holds a **now** void* pointer which represent the current available memory

**qs_sheet** holds 3 AVL tree manager :
- one to store the regular relations
- one to store the relations that wait to be paired (partials)
- one to store the known divisors of N

With small precautions you are supposed to be able to store anything in **now**, then to update **now** accordingly to what you stored.  **now** is always supposed to contain only zero until its end. This is the main memory management technique used by this software.  `mem_align` aims to provide aligned pointers, wasting a few bits if necessary.

### preparation_part_5 .. 6
Fill the manager's **base** array with prime numbers provided [by](https://stackoverflow.com/a/61895974/18765627) a constant expression :
- verify that kN is a square mod prime, ignore it otherwise
- associates the prime with square root of kN mod prime
- associates the prime with its size (log2)
- computes invariants like **D** used to generate the **A** polynomial coefficient

### get_started_iteration
- can restore the relations previously saved by Lanczos algorithm
- can be used to perform analysis, save/restore factorization to file or other periodic actions

### Polynomial `AX^2 + 2BX + C` coefficients : 

| coefficient | is a constant after  | is a constant until|
|--|--|--|
| A | `iteration_part_1` |`inner_continuation_condition` completion
| B | `iteration_part_4` |for loop final expression
| C | `iteration_part_6` |for loop final expression

What is a for loop final expression ? `for ([initialization]; [condition]; [final-expression])`.

Polynomial and related data are computed until `iteration_part_6` complete.\
The algorithm prepares data that will help generate polynomial values that will be multiples of the factor base.\
`iteration_part_7` and `iteration_part_8` are used for sieving.

### Search sieve for relations
`register_relations` searches sieve for relations, it reads the sieve. When the function finds interesting to define the variable **X** , calculates the value of the polynomial in **X** then divides the value with the factor base, it thus tries to establish relations.

Buffered knowledge is structured by `register_relation_kind_1` and `register_relation_kind_2` :
- informations potentially useful are saved into a struct **qs_relation***
- `register_relation_kind_1` immediately build a matrix of **qs_relation***, using AVL tree
- `register_relation_kind_2` combine **qs_relation*** together, using AVL tree

The AVL tree is used to identify duplicates and retrieve their data.

### inner_continuation_condition
Decides if sieving shoud continue or break, the relation counter is usually the condition.

### lanczos_block
This algorithm is supposed to find matrix eigenvalues.

- the process need memory, **fac_lanczos.c** have its own array builder
- all memory taken in **mem.now** is zeroed and reusable after calculations
- the process is usually fast with small numbers, it can take many iterations with 300-bit numbers
- a **reduce_matrix** function is used before giving up

### finalization_part_1 .. 2 .. 3

The function exploit `lanczos_block` answer, a worker is used to remove prime factors from N.

In certain cases such as :
```
./qs 51460938795049063955433175628971167803839994111348342302522016010379
```
when N is not fully factored, it completes the factorization with GCDs and perfect power checks.

### outer_continuation_condition
This function is used to decide if the algorithm should return the control, or return to sieving.
- normal case is to return the control, answer was found
- unusual case is when N isn't fully factored (maybe parameters was wrong)

So before giving up, the algorithm searches 10%, 25% then 50% more relations.

# Thank you
There are many people to thank, the following list is not exhaustive :
- Carl Pomerance
- William Hart
- Jason Papadopoulos
- the professors at "University of Franche-Comté" who taught me programming
