#include <time.h>
#include <sys/time.h>
// Basic ~ 100 lines factorization tester
// use test=1 or test=160 to execute it

static inline void fac_mini_tests(fac_params *m) {
	// Welcome to the testing feature
	// init 5 numbers + a computation sheet
	cint nums[5];
	for (int i = 0; i < 5; ++i)
		cint_init(&nums[i], 2048, 0);
	cint *N = &nums[0], *Q = N + 1;
	cint_sheet *sheet = cint_new_sheet(1 << 10);

	fac_params params = {0}; params.silent = 1 ;
	//---------------

	unsigned sr = add_rand_seed(sheet);
	sr ^= time(0);

	int error_number = 0;

	unsigned nth, bits = m->testing >= 3 && m->testing <= 220 ? m->testing : 150;
	int seconds = bits < 140 ? 30 : bits < 160 ? 60 : bits < 180 ? 120 : 180;
	const char *seconds_str = bits < 140 ? "30 seconds" : bits < 160 ? "minute" : bits < 180 ? "2 minutes" : "3 minutes";
	if (m->testing < 2) printf("-- crescendo : your %s factorization test -- \n\n", seconds_str);
	else printf("-- %3d-bit : your %s factorization test -- \n\n", bits, seconds_str);

	double timeout = 1e6 * seconds;
	double chronometer = 0;
	for ( nth = 1; chronometer < timeout && !error_number && nth <= 1000; ++nth) {
		if (m->testing < 2) bits = nth > 2 ? nth : 3 ;

		int trial_max = (bits > 50) * 100000;
		{
			retry :
			cint_random_bits(N, bits), *N->mem |= trial_max != 0;
			if (cint_is_prime(sheet, N, 2))
				goto retry; // a prime number isn't submitted

			for (int n = 3; n < trial_max; n += 2)
				if (is_prime_4669921(n))
					if (cint_reinit(Q, n), cint_remove(sheet, N, Q))
						goto retry; // a 'trial' divisible number isn't submitted
		}

		char *str = cint_to_string(N, 10);
		printf("%2d. %s = ", nth, str);
		free(str);
		fflush(stdout);

		cint_dup(Q, N);

		struct timeval tv;
		gettimeofday(&tv, 0), chronometer -= tv.tv_sec * 1e6 + tv.tv_usec;
		fac_cint **factors = c_factor(N, &params);
		assert(factors); // answer is not null
		gettimeofday(&tv, 0), chronometer += tv.tv_sec * 1e6 + tv.tv_usec;

		sr += (unsigned) timeout;
		srand(sr);

		for (int i = 0; factors[i]; ++i) {
			str = cint_to_string(&factors[i]->cint, 10);
			// a power of the factor must be removable from number.
			const int powers = (int) cint_remove(sheet, N, &factors[i]->cint);
			switch (powers) {
				case 0 : error_number |= 2; printf(" [%s] ", str); break;
				case 1 : printf(factors[i + 1] ? "%s * " : "%s", str); break;
				default: printf(factors[i + 1] ? "(%s ^ %d) * " : "(%s ^ %d)", str, powers); break;
			}
			free(str);
			fflush(stdout);
			const int n_bits = (int) cint_count_bits(&factors[i]->cint);
			if (n_bits != factors[i]->bits) error_number |= 4;

			if (n_bits > 1) {
				// the factor must be 'prime' with 3 Miller-Rabin iterations.
				const int prime = (int) cint_is_prime(sheet, &factors[i]->cint, 3);
				if (prime == 0) error_number |= 8;
			}

		}

		// number is equal to 1 after all factors removed.
		if (cint_count_bits(N) != 1) error_number |= 16;

		putchar('\n');
		free(factors);
	}

	if (error_number == 0 && m->testing > 1) puts("Prime numbers wasn't submitted, the software has checked the answers.");
	if (error_number & 2) puts("number wasn't a multiple of the factor");
	if (error_number & 4) puts("bit count of the factor was wrong");
	if (error_number & 8) puts("factor of number wasn't 'prime'");
	if (error_number & 16) puts("number wasn't correctly factored");
	if (error_number) printf("qs rand seed was %u\n", params.qs_rand_seed);

	chronometer /= 1e6 ;
	printf("Thank you, technical chronometer displays    %6.3f s   .\n", chronometer);
	if (m->testing >= 2) printf("On average a %3d-bit factorization take      %6.3f s   .\n",  bits, chronometer / nth);

	// Clear the numbers + the computation sheet.
	cint_clear_sheet(sheet);
	for (int i = 0; i < 5; ++i)
		free(nums[i].mem);
}
