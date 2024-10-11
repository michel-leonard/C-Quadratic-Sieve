// This version of the factorization software (all files) is released into the public domain.
// The software is provided "as it" without any warranty, express or implied.

// The goal of the software is to factor integers (from the command line or a file), with no dependency.
// The manager performs trial divisions, checks whether numbers are perfect powers, tests their primality.
// The Quadratic Sieve is used for large numbers (the --force option is available for even larger numbers).
// Output formats JSON, CSV, SQL and Plain are available (in file or in terminal).

#include "avl-trees.c"				// the binary search trees (originally a separate project)
#include "big-num.c"				// the 64+bit integers (originally a separate project)
#include "headers.h"				// the headers
#include "manager.c"				// the factorization manager and i/o utils
#include "basic-math.c"				// the math shortcuts and functions (math.h isn't used)
#include "64-bits-factorization.c"	// the factorization for small numbers (using Pollard's Rho)
#include "quadratic-sieve.c"		// the Quadratic Sieve (essential part of this project)
#include "block-lanczos.c"			// the Lanczos Block algorithm (essential part of this project)

int main(int argc, char *argv[]) {

	// Default state (the integer factorization software don't use global variables).
	state state = {0};

	// Random number generator consistent across platforms.
	state.params.rand.seed = 0x4eee588af4f90e1 ;

	// The software is silent when there is no terminal.
	state.params.verbose = isatty(STDOUT_FILENO) != 0;

	// Read the command line arguments.
	for (int i = 1; i < argc; ++i)
		if (!(i + 3 < argc && read_key_and_3_values(argv + i, &state) && (i += 3)))
			if (!(i + 2 < argc && read_key_and_2_values(argv + i, &state) && (i += 2)))
				if (!(i + 1 < argc && read_key_value(argv + i, &state) && ++i))
					read_flags(argv + i, &state);

	// Ensure the remaining command line arguments are non-zero integers to factorize.
	for (int i = 1; i < argc; ++i)
		if (argv[i] && !validate_string_number(argv[i], &state))
			fprintf(stderr, "%s: Unknown argument '%s'.\n", argv[0], (state.code = 1, argv[i]));

	if (state.params.help)
		print_help(argv[0]); // Option --help or -h was found.
	else if (state.code == 0 && *state.params.generate)
		generate_input_file(&state); // Generate a file containing sample numbers.
	else if (state.code == 0 && prepare_file_descriptors(&state))
		process_multi(argc, argv, &state); // Process the request(s).

	return state.code;
}

// Mon, 09 Sep 2024 06:00:00 GMT - gcc version 8.3.0 on Debian GNU/Linux 10 (buster)
// Compilation is done using "gcc -Wall -pedantic -O3 -std=c99 main.c -o factor"
// Usage: ./factor 51460938795049063955433175628971167803839994111348342302522016010379
