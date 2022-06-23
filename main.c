#include "avl.c"            // the trees.
#include "cint.c"           // the integers.
#include "fac_headers.h"    // factor headers.
#include "fac_utils.c"      // utilities and front-end.
#include "fac_quadratic.c"  // quadratic sieve source.
#include "fac_lanczos.c"    // quadratic sieve Lanczos.
#include "fac_testing.c"    // quadratic sieve tests.

// Why this project use "cint" instead of GMP ?
// - Author search to understand the problematics of 64+ bit integers.
// - Original software goal was to factor 200-bit RSA in 30 seconds.
// - "cint" allow us to see what is sufficient to reach the goal.

static inline void fac_display_verbose(fac_cint **ans);
static inline void fac_display_help(char *name);

int main(int argc, char *argv[]){
	cint N ;
	fac_params config = {0};
	char * n ; // the string to factor in base 10.
	n = fac_fill_params(&config, argc, argv);
	if (config.testing) fac_mini_tests(&config);
	else if (config.help) fac_display_help(argv[0]);
	else if (n) {
		const int bits = 64 + 4 * (int) strlen(n);
		cint_init_by_string(&N, bits, n, 10); // init the number as a cint.
		fac_cint ** answer = c_factor(&N, &config); // execute the routine.
		fac_display_verbose(answer); // print answer.
		free(answer); // release answer memory.
		free(N.mem); // release number memory.
	} else
		fputs("usage : qs [-h] [-s] [number]", stderr);
	return 0 ;
}

static inline void fac_display_verbose(fac_cint ** ans) {
	for(int i = 0; i < 100; ++i)
		putchar(' ');
	putchar('\r');
	char * str = fac_answer_to_string(ans);
	puts(str);
	free(str);
}

static inline void fac_display_help(char *name) {
	char * str = 1 + strrchr(name, '/');
	if (str < name) str = 1 + strrchr(name, '\\');
	if (str < name) str = name;
	puts("=== [ Welcome to the factor function help ] === \n");
	printf(" - use     ./%s 123     to see the factors of 123\n", str);
	printf(" - use     ./%s -test=130     to see a 130-bit factorization test (30 seconds)\n", str);
	printf(" - use     ./%s -limit=150     to define a limit of bit for the quadratic sieve\n", str);
	printf(" - use     ./%s -s [number]    to not see the progress of quadratic sieve\n", str);
	printf(" - numbers around quotes are identified not to be prime\n");
	printf(" - numbers shown have passed some tests like perfect square, perfect cube, primality, trial divisions\n");
	printf(" - best function for > 64-bit numbers is a quadratic sieve largely inspired by a William Hart FLINT implementation\n");
	printf(" - use     gcc -Wall -pedantic -O3 main.c -o qs     to make the software compile with a great optimizer\n\n");
	printf("Bests in math with your factorizer...\n");
	putchar('\n');
}

// Development (Linux) done with gcc version 8.3.0 (Debian 8.3.0-6)
// Development (Microsoft Windows) done with gcc version 11.2.0 (MinGW-W64 x86_64-ucrt-posix-seh)
// In 2022, the common software speed is 1s for 170-bit RSA, 1 min for 230-bit, 10 min for 260-bit.

// ================================================================================================

// Solving of a general modular quadratic equation uses Shanks-Tonelli algorithm for
// all prime divisors of the modulus, lifting the result by Hensel lemma to higher prime
// powers and finally Chinese Remainder Theorem to combine the results for distinct modulus divisors.

// Cut down the computational expenses for Tonelli-Shanks to a single modular exponentiation see :
// R.D. Silverman, The Multiple Polynomial Quadratic Sieve, Math. Comput. 48, (1987), 329-339.

// Hensel's lemma asserts that every factorization of h modulo m into coprime polynomials can be lifted in a
// unique way into a factorization modulo m^k for every k, see : wikipedia.org "Hensel's Lemma"

// See : https://www2.karlin.mff.cuni.cz/~krypto/Implementace_MPQS_SIQS_files/main_file.pdf for interesting precisions.
