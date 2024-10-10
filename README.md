# Factorization software using Quadratic Sieve

This **general-purpose integer factorization software** is released into the public domain. It is provided "as is" without any warranty, express or implied.

## Your Pure C99 Implementation

- simply has no dependencies, not even `math.h`.
- is immediately compatible with Microsoft Windows and Linux.
- use Pollard's Rho algorithm for factoring 64-bit integers.
- comes with a Big Number library for handling large numbers.
- includes a Self-Initializing Quadratic Sieve (SIQS).
- use AVL trees to assist the Quadratic Sieve in data management.
- is capable of batch factor numbers with JSON or CSV output.
- has been tested, produces reproducible results and is thread-safe.

## The Factorization Manager

- reads numbers to be factored from either a file or the command line.
- performs preliminary checks (e.g., trial division, perfect power) before invoking the algorithms.
- verifies the factorization (the product of factors should be equal to the input number).
- returns a `0` exit code when all numbers are fully factored, non-zero otherwise.

## Using the Software

The `factor.exe` file is available for Windows users. Otherwise, to get an executable, follow these steps:

- **Ubuntu:** Install a C compiler with the command: `sudo apt install build-essential`.
- **Windows:** Install [MinGW](https://winlibs.com/), which will provide a **GCC** compiler similar to Ubuntu.

Once you have the compiler, use **Terminal** on Ubuntu or **PowerShell** on Windows to navigate to the downloaded directory using the [`cd`](https://en.wikipedia.org/wiki/Cd_(command)) command:

- Compile the program by executing: `gcc -Wall -pedantic -O3 main.c -o factor`. GCC will create the appropriate executable for your device.

Compilation should take just a few seconds, and then you can use the software:

- Run the executable with a command like: `./factor 2840222005804112554242019945562782034485263657432907`

The software will quickly display the factors:

```
Number: 2840222005804112554242019945562782034485263657432907
Factors: 186565459987 (prime), 123384476380949634131^2
```
## Command Line Options

These options are designed to facilitate your interaction with the software:

| Name | Short | Goal | Example |
|:---:|:---:|:---:|---|
| `--input-file` | `-i` | Factor the entire file. | `-i input.txt` |
| `--output-file` | `-o` | Write results to a file. | `-o output.txt` |
| `--timeout` | `-t` | Interrupt the Quadratic Sieve after a specified duration in seconds. | `-t 1` |
| `--force` | `-f` | Override the default limits: 8191 digits for a number and 220-bit for the Quadratic Sieve.| `-f` |
| `--verbose` | `-v` | Show more or fewer details, such as Quadratic Sieve progress, depending on the verbosity level. The default value is 1. | • `-v 0` <br>• `-v 1` <br>• `-v 2` <br>• `-v 3` <br>• `-v 4` |
| `--output-csv` | `-c` | Format results as CSV. | `-c` |
| `--output-json` | `-j` | Format results as JSON. | `-j` |
| `--output-json-compact` | `-J` | Same as `--output-json`, but disable the pretty print. | `-J` |
| `--output-sql` | `-sql`| Format results as SQL. | `-sql` |
| `--rand-seed` | `-r` | Specify a custom XorShift random seed (64-bit). | `-r 123` |
| `--generate` | `-g` | Designed for testing, this option create a `generated.txt` suitable as input file, containing non-trivial sample numbers. Its parameters, all optional, are `-g <bits-min> <bits-max> <count>`. |• `-g` <br>• `-g 150` <br>• `-g 130 150` <br>• `-g 150 150 1000`|
| `--help` | `-h` | Print the help. |`-h`|

## Quadratic sieve parameters

You can enforce parameters to the Quadratic Sieve :

|Quadratic Sieve option|Typical value range|
|--|--|
|`--qs-multiplier <value>`|from `1` to `150`|
|`--qs-base-size <value>`|from `1000` to `50000`|
|`--qs-large-prime <value>`|from `300000` to `30000000`|
|`--qs-alloc-mb <value>`|from `1` to a few hundred|
|`--qs-sieve <value>`|`31744` by default|
|`--qs-threshold <value>`|from `60` to `110`|
|`--qs-error-bits <value>`|from `15` to `35`|
|`--qs-laziness <value>`|from `80` to `120` (in percentage)|

Navigate through the source code to see their default value and usage.

## Speed of the factorizer

The software is designed to be a **general-purpose factorizer** for integers, also supporting negative numbers. The algorithm can factor large numbers (e.g., a 1500-bit highly composite number or 10,000!), but performance can vary significantly depending on the structure of the number. For every additional 3 decimal digits (or 10 bits) to the input number, the Quadratic Sieve duration roughly doubles.

- **Default target performance**
  - The software can factor **170-bit RSA** integers in a second.
  - The software can factor **220-bit RSA** integers in approximately 30 seconds.

- **Performance with basic tuning (sometimes without tuning)**
  - Successfully factors RSA numbers greater than 250 bits.
  - Factored a **300-bit RSA** number in 2 hours.
  - Factored the **321-bit** RSA number related to the [bank card case](https://www.enseignement.polytechnique.fr/profs/informatique/Eric.Goubault/Cours09/qs.pdf).
  - Factored RSA-100, a 330-bit number.

## Testing the Software

You can use the built-in generator to easily produce various reproducible samples (consistent across all platforms) in a `generated.txt` file. These examples are for Ubuntu:

Example of a testing command that factor 10 numbers of 150-bit based on the seed `1`, taking a few seconds:
```sh
./factor --generate 150 150 10 --rand-seed 1
./factor -i generated.txt -c -v 2 && echo "All the numbers have been fully factored."
```
- You can use `--rand-seed 2` to get another set of sample numbers.
- In Ubuntu, use `--rand-seed $(date +%s%N)` to systematically get another set of sample numbers.
- In Windows the `&&` operator isn't always supported, and the factorizer often has a `.exe` extension.
- When the software returns a `0` exit code, it means that all numbers have been fully factored.

Example of a verbose testing command that factor `1189625923578446633052664796204789022410184682730580575561447`, a 200-bit number:
```sh
./factor --generate 200
./factor -i generated.txt -v 4 && echo "The number have been fully factored."
```

  Example of a testing command that factor sample numbers ranging from 120-bit to 170-bit, taking 20 seconds:
```sh
./factor --generate 120 170
./factor -i generated.txt -v 2 && echo "All the numbers have been fully factored."
```

Example of the testing command that generated the provided `.csv` and `.json` samples of output files. To give you an idea of ​​the environment, it took 11 minutes per file on the machine that tested the software:
```sh
./factor --generate
./factor --input-file generated.txt --output-file example-result.csv --output-csv
./factor --input-file generated.txt --output-file example-result.json --output-json
```
Example of a testing command that shows that the default multiplier is better than 1 (no multiplier) and speeds up execution time by a factor 2:
```sh
./factor --generate 170 170 10
./factor -v 2 -c -i generated.txt
./factor -v 2 -c -i generated.txt --qs-multiplier 1
```

Example of a quick testing command based on the OEIS list of [highly composite numbers](https://oeis.org/A002182/b002182.txt):
```sh
wget -qO- https://oeis.org/A002182/b002182.txt | cut -d' ' -f 2 > oeis-list.txt
./factor -i oeis-list.txt -c && echo "All the numbers have been fully factored."
```

## How Does the Factorization Manager Operate?

- it reads options from the command line.
- it prepares a session capable of handling the largest number.
- it recursively factors each number using the available algorithms.
- it reuses the session to factor all numbers.

Before starting, the input integers are parsed; if an error is detected, the factorization will not begin, and you will be notified via the standard error output. The software attempts to commit the output after each factorization by clearing buffers.

## Memory Management

- The manager allocates a few memory (a session) based on the largest number.
- The Quadratic Sieve manages its memory independently, based on its input.
- Memory used by the Quadratic Sieve is freed as soon as it has factored its input.
- The manager reuse its session, the memory is feeed once all numbers are factored.
- Valgrind was used to check for memory leaks, no leaks were found.

Memory allocations are relatively small (typically a few dozen megabytes) and are passed to `assert`. In case of memory allocation errors (e.g., insufficient available memory), the software will exit, and it is recommended to restart your device.

## Big Numbers

The `cint` library (released "as is", into the public domain, without any warranty, express or implied) is provided for handling large integers. It includes a variety of basic and advanced mathematical functions to support your calculations.

## AVL Binary Search Trees

The `avl` functions (released "as is", into the public domain, without any warranty, express or implied) manage thousands of integers in AVL trees, making Quadratic Sieve relations storage and retrieval efficient and straightforward.

## Block Lanczos

The Block Lanczos algorithm (released "as is", into the public domain, without any warranty, express or implied) solve a large, sparse system of linear equations, and identify dependencies among the rows of the matrix formed by the relations identified by the Quadratic Sieve.

## Files

Your software's source code is organized into the following files:

| Filename                | Purpose                                                                                                               |
|-------------------------|-----------------------------------------------------------------------------------------------------------------------|
| `64-bits-factorization.c` | Contains functions for factoring classical integers, including Pollard's Rho algorithm and a Miller-Rabin primality checker. |
| `avl-trees.c`           | Contains the AVL Trees functions, originally a separate project from the integer factorization project.               |
| `basic-math.c`          | Includes a prime number generator capable of yielding the first 326,983 prime numbers, along with the Kronecker symbol function, Tonelli-Shanks function, and modular inverse function. |
| `big-num.c`             | Contains the `cint` library, also initially a separate project from the integer factorization project.               |
| `block-lanczos.c`       | Provides the linear algebra part of the Quadratic Sieve, including a matrix reducer and a manager for this task.     |
| `headers.h`             | Contains definitions for types, structures, functions, and `DEBUG` macros used throughout the project.               |
| `main.c`                | Contains the program's entry point and transfers control to the factorization manager after about 30 lines.          |
| `manager.c`             | Includes I/O functions, help documentation, shortcuts commonly used by the Quadratic Sieve, and the factorization manager source code. |
| `quadratic-sieve.c`     | Contains the Quadratic Sieve implementation, divided into around thirty functions, with a main function that manages the factorization process in about 30 lines. Search for "debug" in this file to find the debugging elements displayed. |

The software :
- uses a stateful approach, it is thread-safe and no global variables are used.
- uses reproductible **XorShift** random values that depend on the `--rand-seed` (default to `0`).
- produces reproducible results regardless of the order in which the factorizations are batched.

## The file quadratic-sieve.c

The quadratic sieve file structure is as follows:

- the ~40 lines function that allows to see the algorithm **structure**
- the 2 important loop **conditions**
- the algorithm **parameters**
- the **functions** approximately in the order they are called

### preparation_part_2 .. 3

The input **N is** duplicated and called "**kN**" after this function complete :

- multiply **N** by a prime number to reach 120 bits
- apply a multiplier to **N** intended to optimize runtime

The `preparation_part_3_default_proposition` function provide a multiplier that make the algorithm completes **faster** on average. A `preparation_part_3_alternative_proposition` is also provided for the same purpose but is unused.

| Variable | Information |
|----------|-------------|
| N        | Prime factors are removed from **N**, and **N** is updated until `N = 1`. |
| kN       | The algorithm computes with **kN**, which remains constant. |

See of how multipliers impact the execution of `./factor 10146246660331135529849066552891483667878618078493758287`:

| Option used | Makes the algorithm factoring | Comment on the multiplier | Took |
|:---:|:---:|:---:|:---:|
|`--qs-multiplier 7`| `kN` = `7` * `N`|**Selected by the algorithm**|3.1 s|
|`--qs-multiplier 1`| `kN` = `1` * `N`|No multiplier|6.3 s|
|`--qs-multiplier 113`| `kN` = `113` * `N`|Large multipliers are sometimes good|8.2 s|
|`--qs-multiplier 29`| `kN` = `29` * `N`|One of the worst|19.4 s|
|`--qs-multiplier 35`| `kN` = `35` * `N`|One of the worst|26.3 s|

While the proposed method doesn't always select the absolute best multiplier, it finds a balance that improves average runtime by a factor 2. The best multiplier can be up to **10 times faster** than the worst. Can you guess which multiplier achieves this ? For those interested, the **Knuth-Schroeppel analysis** explores why choosing the right multiplier makes such a significant difference ([Learn more](https://www.google.com/search?q=Knuth+Schroeppel+analysis)).

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

With small precautions you are supposed to be able to store anything in **now**, then to update **now** accordingly to what you stored. **now** is always supposed to contain only zero until its end. This is the main memory management technique used by this software. `mem_align` aims to provide aligned pointers, wasting a few bits if necessary.

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
|:---:|:---:|:---:|
| A | `iteration_part_1` |`inner_continuation_condition` completion
| B | `iteration_part_4` |for loop final expression
| C | `iteration_part_6` |for loop final expression

- Polynomial and related data are computed until `iteration_part_6` complete
- The algorithm prepares data that will help generate polynomial values that will be multiples of the factor base
- `iteration_part_7` and `iteration_part_8` are used for sieving

### Search sieve for relations
`register_relations` searches sieve for relations, it reads the sieve. When the function finds interesting to define the variable **X** , calculates the value of the polynomial in **X** then divides the value with the factor base, it thus tries to establish relations.

Buffered knowledge is structured by `register_regular_relation` and `register_partial_relation` :

- informations potentially useful are saved into a struct **qs_relation***
- `register_regular_relation` immediately build a matrix of regular full **qs_relation*** using AVL tree
- `register_partial_relation` combine (single large primes variation) **qs_relation*** together using AVL tree

The AVL trees are used to identify duplicates and retrieve their data.

### inner_continuation_condition
Decides if sieving shoud continue or break, the relation counter is usually the condition.

### lanczos_block

This algorithm find a subset of all exponent vectors such that the sum of their exponent vectors is the zero vector :

- the process need memory, **fac_lanczos.c** have its own array builder
- all memory taken in **mem.now** is zeroed and reusable after calculations
- finding the matrix eigenvalues is usually fast, it can still take 400+ iterations with large N
- the algorithm is probabilistic with high success rate, generally it is lucky enough
- several tries with and without **reduce_matrix** are performed before giving up

The `reduce_matrix` function helps `lanczos_block` and could be called unconditionally, but for now it still first tries not to call it. This was an experience that helped me understand why we should sometimes discard a few relations.

### finalization_part_1 .. 2

The `finalization_part_1` function performs the **square root + GCD** step based on the `lanczos_block` answer, it leads to a factorization of N.

### qs_register_divisor

The `qs_register_divisor` function :
- perform checks (such as primality, perfect power and GCD between known divisors)
- register the appropriate divisors (prime factors of N) to the caller's routine

### outer_continuation_condition
This function is used to decide if the algorithm should return the control, or return to sieving.

- normal case is to return the control, answer was found
- unusual case is when N isn't fully factored (maybe parameters was wrong)

So before giving up, the algorithm searches 10%, 25% then 50% more relations.

## Large Prime Variation

*In short: Given two partial relations with the same large prime, the `register_partial_relation` function merges them to obtain a full relation.*

This step is optional for completing a factorization but can significantly speed up the process. To disable this feature, simply avoid calling the `register_partial_relation` function.

During the execution of `register_relations`, after trial division by the factor base primes (`cint_remove`), if the quotient remains greater than 1 but is less than the large prime bound, the software registers this partial relation with `register_partial_relation`. This function stores the new `qs_relation` in a separate location from full relations.

When the large prime collected matches a prime from another **partial relation** (detected via the AVL tree), the software **combines these two relations** and then calls the `register_regular_relation` function, which registers all regular full relations. At this point, the newly created relation is treated like any other full relation.

|N bits| Relations added immediately | Relations added by single large primes variation | Relations needs
|:---:|:---:|:---:|:---:|
| 150 | 2551 | 3| 2554
| 170 | 3168 | 659 | 3827
| 190 | 3212 | 1786 | 4998
| 210 | 3409 | 3012 | 6421
| 230 | 5659 | 4686 | 10345
| 250 | 10051 | 9793 | 19844
| 270 | 14026 | 15484 | 29510
| 290 | 12820 | 20201 | 33021

## Online Tools

- [Dario Alpern's Integer Factorization Calculator](https://www.alpertron.com.ar/ECM.HTM): This tool processes the factorization directly in your browser.
- [Number Empire Factoring Calculator](https://www.numberempire.com/factoringcalculator.php): This tool performs factorization on the server.


# Thank you
There are many people to thank :

- Carl Pomerance
- Jason Papadopoulos
- William Hart
- my professors at University of Franche-Comté ❤️
- GitHub and SourceForge users reporting issues

## Any issue ?

You can report any issues, and we'll investigate to correct the source files.
