//  GNU General Public License

//      As an undergraduate student, this project is part of my computer science + maths training.
//          - This software proposition is from Michel Leonard (student at Université de Franche-Comté, Mon, 11 Jul 2022)
//          - There is no guarantee of any kind on the software
//          - C code is shared under the terms of the GNU General Public License
//          - The main mathematical and logical inspiration source is located at :
//              http://web.mit.edu/sage/export/flintqs-0.0.20070817/QS.cpp - GNU General Public License

//  This software implementation would have been impossible without "FLINT: Fast Library for Number Theory" maintained by William Hart.

// Quadratic sieve main routine
int quadratic_sieve(fac_caller *caller) {
	// The routine is invoked by a caller, and :
	// - must factor the caller's "number"
	// - must register answers to the caller's routines
	// - can use caller's resources (vars + configurations)
	// - resources are enough to copy 400-bit into vars

	int limit = (int) caller->params->qs_limit;
	if (limit == 0) limit = 220; // the default limit is such that the QS answer in seconds.
	if (caller->number->bits < 64 || caller->number->bits > limit) return 0;

	qs_sheet qs = {0};
	preparation_part_1(&qs, caller);
	preparation_part_2(&qs);
	preparation_part_3(&qs);
	qs_parametrize(&qs);
	preparation_part_4(&qs);
	preparation_part_5(&qs);
	preparation_part_6(&qs);
	do {
		do {
			get_started_iteration(&qs);
			iteration_part_1(&qs, &qs.poly.D, &qs.poly.A);
			iteration_part_2(&qs, &qs.poly.A, &qs.poly.B);
			iteration_part_3(&qs, &qs.poly.A, &qs.poly.B);
			for (qs_sm i = 0, addi, *corr; i < qs.poly.gray_max && qs.n_bits != 1; ++i, ++qs.poly.curves) {
				addi = iteration_part_4(&qs, i, &corr, &qs.poly.B);
				iteration_part_5(&qs, &qs.constants.kN, &qs.poly.B);
				iteration_part_6(&qs, &qs.constants.kN, &qs.poly.A, &qs.poly.B, &qs.poly.C);
				iteration_part_7(&qs, addi, corr);
				iteration_part_8(&qs, addi, corr);
				register_relations(&qs, &qs.poly.A, &qs.poly.B, &qs.poly.C);
			}
		} while (inner_continuation_condition(&qs));
		// Lanczos may "do more" than read-only to the relations, in this case a snapshot will be available.
		finalization_part_1(&qs, lanczos_block(&qs));
		finalization_part_2(&qs);
	} while (outer_continuation_condition(&qs));
	const int res = finalization_part_3(&qs);
	free(qs.mem.base);
	return res;
}

// Quadratic sieve main condition 1 : user configuration
// continue to search relations or break ? res is the answer.
int inner_continuation_condition(qs_sheet *qs) {
	int res = 1;
	res &= qs->n_bits != 1 ; // the bit count of N may have changed.
	res &= qs->relations.length.now < qs->relations.length.needs; // the condition.
	if (qs->caller->params->silent == 0) {
		const double rel_begin = (double) qs->relations.length.now, rel_end = (double) qs->relations.length.needs ;
		fac_display_progress("Quadratic sieve", 100. * rel_begin / rel_end); // progress isn't linear
	}
	return res;
}

// Quadratic sieve main condition 2 : user configuration
// return to sieving or stop the algorithm ? res is the answer.
int outer_continuation_condition(qs_sheet *qs) {
	int res = qs->sieve_again_perms-- > 0; // avoid infinite loop.
	res &= qs->divisors.total_primes < qs->sieve_again_perms; // search prime factors.
	res &= qs->n_bits != 1 ; // the bit count of N isn't 1.
	if (res) {
		qs_sm new_needs = qs->relations.length.needs;
		new_needs += new_needs >> (1 + qs->sieve_again_perms);
		if (qs->caller->params->silent == 0)
			printf("quadratic sieve targets %u more relations\n", new_needs - qs->relations.length.needs);
		qs->relations.length.needs = new_needs ;
	}
	return res;
}

// Quadratic sieve parameters : user configuration
qs_sm linear_param_resolution(const double v[][2], const qs_sm point) {
	qs_sm res, i, j ;
	if (v[1][0] == 0)
		res = (qs_sm) v[0][1];
	else {
		for (j = 1; v[j + 1][0] && point > v[j][0]; ++j);
		i = j - 1;
		if (v[i][0] > point) res = (qs_sm) v[i][1];
		else if (v[j][0] < point) res = (qs_sm) v[j][1];
		else {
			const double a = (v[j][1] - v[i][1]) / (v[j][0] - v[i][0]);
			const double b = v[i][1] - a * v[i][0];
			res = (qs_sm) (a * point + b);
		}
	}
	return res + (res > 512) * (16 - res % 16) ;
}

void qs_parametrize(qs_sheet *qs) {
	const qs_sm bits = (qs_sm) cint_count_bits(qs->caller->vars); // kN
	qs->kn_bits = bits; // input was adjusted, there is at least 115-bit.

	// params as { bits, value } take the extremal value if bits exceed.
	static const double param_base_size [][2]= { {110, 800}, {130, 1500}, {210, 6e3}, {240, 12e3}, {260, 24e3}, {320, 36e3}, {0} };
	qs->base.length = linear_param_resolution(param_base_size, bits);

	static const double param_laziness [][2]= {{110, 90}, {180, 100}, {220, 100}, {250, 110}, {0} };
	// collecting more/fewer relations than recommended (used to verify "sieve again" feature).
	qs->relations.length.needs = qs->base.length * linear_param_resolution(param_laziness, bits) / 100 ;

	static const double param_m_value [][2]= { {110, 1}, {190, 4}, {0} };
	qs->m.length = (qs->m.length_half = 31744 * linear_param_resolution(param_m_value, bits)) << 1;

	qs->m.cache_size = 126976 ; // algorithm reaches "M length" by steps of "cache size".

	static const double param_error [][2]= { {110, 13}, {300, 33}, {0} };
	qs->error_bits = linear_param_resolution(param_error, bits);

	static const double param_threshold [][2]= { {110, 63}, {220, 78}, {300, 99}, {0} };
	qs->threshold.value = linear_param_resolution(param_threshold, bits);

	static const double param_alloc [][2]= { {1e3, 2}, {3e3, 8}, {5e3, 20}, {15e3, 80}, {36e3, 256}, {0} };
	qs->mem.bytes_allocated = linear_param_resolution(param_alloc, qs->base.length) << 20 ; // { base size, overall needs of megabytes }

	qs->sieve_again_perms = 3; // Sieve again up to 3 times before giving up

	// Iterative list
	qs->iterative_list[0] = 1; // one
	static const double param_first_prime [][2]= { {170, 8}, {210, 12}, {300, 30}, {0} };
	qs->iterative_list[1] = linear_param_resolution(param_first_prime, bits); // first
	const qs_sm large_base = qs->base.length > 5000 ? 5000 : qs->base.length;
	qs->iterative_list[2] = large_base >> 2; // medium
	qs->iterative_list[3] = large_base >> 1; // mid
	qs->iterative_list[4] = large_base; // sec
	qs->iterative_list[5] = qs->base.length; // factor base size

	const double last_prime_in_base = qs->base.length * 2.5 * log_computation(qs->base.length);
	qs->relations.large_prime = (qs_md) (last_prime_in_base * last_prime_in_base * 1.25);

	// Other parameters (they have never been changed)
	qs->s.values.double_value = (qs->s.values.defined = (qs->s.values.subtract_one = bits / 28) + 1) << 1;
	qs->poly.gray_max = 1 << (qs->s.values.defined - 3); // computing the roots of f(x) once for all these polynomials.

	// The algorithm itself completes its configuration during the last preparation part.
}

// Quadratic sieve source : algorithm
void preparation_part_1(qs_sheet *qs, fac_caller *caller) {
	// initializing with the caller's resources
	qs->caller = caller;
	qs->calc = caller->calc;
	qs->n_bits = caller->number->bits;
}

void preparation_part_2(qs_sheet *qs) {
	cint * N = &qs->caller->number->cint, * kN = qs->caller->vars, *ADJUSTOR = kN + 1 ;
	// Input is "transparently" adjusted by a prime number to measure at least 115-bit.
	static const int prime_generator[] = {
			9, 7, 5, 3, 17, 27, 3, 1, 29, 3, 21, 7, 17, 15,
			9, 43, 35, 15, 29, 3, 11, 3, 11, 15, 17, 25, 53,
			31, 9, 7, 23, 15, 27, 15, 29, 7, 59, 15, 5, 21,
			69, 55, 21, 21, 5, 159, 3, 81, 9, 69, 131, 33, 15 };
	const qs_sm bits = qs->n_bits;
	if (bits < 115) {
		qs->adjustor = bits < 115 ? ((qs_md)1 << (124 - bits)) + prime_generator[115 - bits] : 1 ;
		simple_int_to_cint(ADJUSTOR, qs->adjustor);
		cint_mul(N, ADJUSTOR, kN);
	} else
		qs->adjustor = 1, cint_dup(kN, N);
}

void preparation_part_3(qs_sheet *qs) {
	// the use of a multiplier can improve the run time of the algorithm.
	qs_sm mul = (qs_sm) qs->caller->params->qs_multiplier ;
	if (mul == 0)
		mul = preparation_part_3_michel(qs);
	qs->multiplier = mul ;
	cint *kN = qs->caller->vars, *MUL = kN + 1, *N = kN + 2 ;
	if (qs->multiplier > 1) {
		simple_int_to_cint(MUL, qs->multiplier);
		cint_dup(N, kN);
		cint_mul(MUL, N, kN);
	}
}

static inline qs_sm preparation_part_3_michel(qs_sheet *qs) {
	// The function propose a multiplier for N, that should (on average) allow accumulating relations faster.
	qs_sm mul[] = {1, 2, 3, 4, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43}, n_mul = (qs_sm)(sizeof(mul) / sizeof(qs_sm));
	cint *N = qs->caller->vars, *PRIME = N + 1, *Q = N + 2, *R = N + 3;
	qs_sm i, j;
	double score[n_mul];
	const qs_sm n_mod_8 = (qs_sm) (*N->mem % 8);
	for (i = j = 0; i < n_mul; ++i) {
		const qs_sm x = n_mod_8 * mul[i] % 8;
		if (x == 1 || x == 7) {
			mul[j] = mul[i];
			score[j] = 1.38629436 - log_computation((double) mul[i]) / 2.0;
			++j;
		}
	}
	if (j)
		for (n_mul = j, i = 3; i < 500 ; i += 2)
			if (is_prime_4669921(i)) {
				const double intake = log_computation((double) i) / (i - 1);
				simple_int_to_cint(PRIME, i);
				cint_div(qs->calc, N, PRIME, Q, R);
				const int kronecker = kronecker_symbol((qs_sm) simple_cint_to_int(R), i);
				for (j = 0; j < n_mul; ++j)
					score[j] += intake + intake * kronecker * kronecker_symbol(mul[j], i);
			}
	for (i = 1, j = 0; i < n_mul; ++i)
		if (score[i] > 1.075 * score[0])
			score[0] = score[j = i];
	return mul[j];
}

void preparation_part_4(qs_sheet *qs) {
	void *mem;
	mem = qs->mem.base = calloc(1, qs->mem.bytes_allocated);
	assert(mem);

	if (qs->caller->params->qs_rand_seed) srand(qs->rand_seed = qs->caller->params->qs_rand_seed);
	else qs->caller->params->qs_rand_seed = qs->rand_seed = add_rand_seed(&mem);

	// kN was computed into the caller's courtesy memory, now the QS has parametrized and "allocated"
	const size_t kn_size = qs->caller->vars[0].end - qs->caller->vars[0].mem + 1 ;
	// the quadratic sieve variables can store at most kN ^ 2 in terms of bits
	const size_t vars_size = kn_size << 1 ;
	const size_t buffers_size = qs->base.length + (qs->iterative_list[1] << 1);
	// more relation pointers than "guessed" are available (sieve again feature).
	const size_t relations_size = (qs->base.length < qs->relations.length.needs ? qs->relations.length.needs : qs->base.length) * 9 / 4 ;

	{
		// list of the numbers used by the algorithm
		cint * const n[] = {
				&qs->vars.N,
				// polynomial
				&qs->poly.A,
				&qs->poly.B,
				&qs->poly.C,
				&qs->poly.D,
				// temporary
				&qs->vars.TEMP[0], &qs->vars.TEMP[1], &qs->vars.TEMP[2], &qs->vars.TEMP[3], &qs->vars.TEMP[4],
				// "my" are not used by the algorithm, they serve to facilitate development checks.
				&qs->vars.MY[0], &qs->vars.MY[1], &qs->vars.MY[2], &qs->vars.MY[3], &qs->vars.MY[4],
				// relations finder
				&qs->vars.X,
				&qs->vars.KEY,
				&qs->vars.VALUE,
				&qs->vars.CYCLE,
				// a factor of N
				&qs->vars.FACTOR,
				// constants
				&qs->constants.kN,
				&qs->constants.ONE,
				&qs->constants.M_HALF,
				&qs->constants.LARGE_PRIME,
				&qs->constants.MULTIPLIER,
				0,
		};
		for (int i = 0; n[i]; ++i) {
			n[i]->mem = n[i]->end = mem_aligned(mem) ;
			mem = n[i]->mem + (n[i]->size = vars_size);
		}
	}

	cint_dup(&qs->vars.N, &qs->caller->number->cint);
	cint_dup(&qs->constants.kN, qs->caller->vars);

	simple_int_to_cint(&qs->constants.ONE, 1);
	simple_int_to_cint(&qs->constants.M_HALF, qs->m.length_half);
	simple_int_to_cint(&qs->constants.LARGE_PRIME, qs->relations.large_prime);
	simple_int_to_cint(&qs->constants.MULTIPLIER, qs->multiplier);

	// Allocates "s" rows
	qs->s.data = mem_aligned(mem);
	mem = mem_aligned(qs->s.data + qs->s.values.defined);
	for (qs_sm i = 0; i < qs->s.values.defined; ++i) {
		simple_inline_cint(&qs->s.data[i].B_terms, kn_size, &mem); // also "s" more cint
		qs->s.data[i].A_inv_double_value_B_terms = mem;
		mem = mem_aligned(qs->s.data[i].A_inv_double_value_B_terms + qs->base.length);
	}
	qs->s.A_indexes = mem_aligned(mem); // the indexes of the prime numbers that compose A

	// Allocates "base length" rows
	qs->base.data = mem_aligned(qs->s.A_indexes + qs->s.values.double_value);
	qs->m.positions[0] = mem_aligned(qs->base.data + qs->base.length);
	qs->m.positions[1] = mem_aligned(qs->m.positions[0] + qs->base.length);
	qs->m.sieve = mem_aligned(qs->m.positions[1] + qs->base.length);
	qs->m.sieve[qs->m.length] = 0xFF ; // the end of the sieve evaluates to "true" under any appropriated mask.
	qs->m.flags = mem_aligned(qs->m.sieve + qs->m.length + sizeof(uint64_t));
	// buffer[0] is zeroed after use, buffer[1] isn't supposed zeroed.
	qs->buffer[0] = mem_aligned(qs->m.flags + qs->base.length);
	qs->buffer[1] = mem_aligned(qs->buffer[0] + buffers_size);

	// Other allocations
	qs->relations.length.reserved = (qs_sm) relations_size ;
	// Lanczos Block has a part of memory, it takes a "lite" snapshot before throwing relations.
	qs->lanczos.snapshot = mem_aligned(qs->buffer[1] + buffers_size) ;
	qs->relations.data = mem_aligned(qs->lanczos.snapshot + relations_size);
	qs->divisors.data = mem_aligned(qs->relations.data + relations_size);
	qs->mem.now = mem_aligned(qs->divisors.data + 512);

	const qs_sm n_trees = (qs_sm) (sizeof(qs->uniqueness) / sizeof(struct avl_manager));
	for (qs_sm i = 0; i < n_trees; ++i) {
		// the trees are used to identify duplicates (relations, partials, factors of N)
		qs->uniqueness[i].inserter_argument = &qs->mem.now;
		qs->uniqueness[i].inserter = &avl_cint_inserter;
		qs->uniqueness[i].comparator = (int (*)(const void *, const void *)) &h_cint_compare;
		// they use default sign-less comparator.
	}
}

void preparation_part_5(qs_sheet *qs) {
	static const double inv_ln_2 = 1.44269504088896340736;
	cint *A = qs->vars.TEMP, *B = A + 1, *C = A + 2;
	qs_sm i = 0, prime;
	// The factor base contain the multiplier different from 2.
	if (qs->multiplier != 2)
		qs->base.data[i].size = (qs_sm) (.35 + inv_ln_2 * log_computation(qs->base.data[i].num = qs->multiplier)), ++i;

	// Then the number 2.
	qs->base.data[i].num = 2, qs->base.data[i].size = 1;
	qs->base.data[i].kN_sqrt_mod_prime = *qs->constants.kN.mem % 8 == 1 || *qs->constants.kN.mem % 8 == 7, ++i;

	// Then the prime numbers for which kN is a quadratic residue modulo.
	for (prime = 3; i < qs->base.length; prime += 2)
		if (is_prime_4669921(prime)) {
			simple_int_to_cint(A, prime);
			cint_div(qs->calc, &qs->constants.kN, A, B, C);
			const qs_sm kn_mod_prime = (qs_sm) simple_cint_to_int(C);
			qs->base.data[i].kN_sqrt_mod_prime = tonelli_shanks(kn_mod_prime, prime);
			if (qs->base.data[i].kN_sqrt_mod_prime) {
				qs->base.data[i].num = prime;
				qs->base.data[i].size = (qs_sm) (.35 + inv_ln_2 * log_computation(prime)), ++i;
			} // 2.5 * (base size) * ln(base size) is close to the largest prime number in factor base.
		}
}

void preparation_part_6(qs_sheet *qs) {
	// completes the configuration by the algorithm itself.
	// computes D, a template for the A polynomial coefficient.
	const qs_sm s = qs->s.values.defined ;
	qs_sm i, min;
	qs->poly.span_half = (qs->poly.span = qs->base.length / s / s / 2) >> 1;
	cint *kN = qs->vars.TEMP, *TMP = kN + 1, *R = kN + 2;
	cint_dup(kN, &qs->constants.kN);
	cint_left_shifti(kN, 1);
	cint_sqrt(qs->calc, kN, TMP, R);
	cint_div(qs->calc, TMP, &qs->constants.M_HALF, &qs->poly.D, R);
	qs->poly.d_bits = (qs_sm) cint_count_bits(&qs->poly.D);
	cint_nth_root(qs->calc, &qs->poly.D, s, R); assert(s >= 3); // use the s-th root of D.
	const qs_sm root = (qs_sm) simple_cint_to_int(R) ;
	for (i = 1; assert(i < qs->base.length), qs->base.data[i].num <= root; ++i);
	assert(i >= qs->poly.span);
	for (min = i - qs->poly.span_half, i *= i; i / min < qs->poly.span + min; --min);
	qs->poly.min = min ;
	//
	i = qs->poly.span_half + min ; i = 1 + i * i / min ;
	if (i < qs->poly.span + min) i = qs->poly.span + min ;
	if (qs->iterative_list[3] < i) qs->iterative_list[3] = i ;
	assert(qs->iterative_list[3] <= qs->iterative_list[4]);
}

void get_started_iteration(qs_sheet *qs) {
	if (qs->lanczos.snapshot[0].relation) {
		// the operation is fast, it shouldn't happen in average case.
		// it restores the relations reduced by Lanczos algorithm.
		qs_sm i ;
		for(i = 0; qs->lanczos.snapshot[i].relation; ++i) {
			qs->relations.data[i] = qs->lanczos.snapshot[i].relation;
			qs->relations.data[i]->Y.length = qs->lanczos.snapshot[i].y_length;
			qs->lanczos.snapshot[i].relation = 0 ;
		}
		qs->relations.length.now = i ;
	}
	// D is randomized if algorithm remarks that no relation accumulates (the software tester didn't remove it)
	if (qs->relations.length.prev == qs->relations.length.now && qs->poly.curves)
		cint_random_bits(&qs->poly.D, qs->poly.d_bits);
	qs->relations.length.prev = qs->relations.length.now;
}

void iteration_part_1(qs_sheet * qs, const cint * D, cint * A) {
	// A is a "random" product of "s" distinct prime numbers from the factor base.
	cint * X = qs->vars.TEMP, * Y = X + 1, *TMP ;
	qs_sm a, b = qs->poly.min + qs->poly.span_half, i = 0, j;
	if (qs->s.values.defined & 1) TMP = A, A = X, X = TMP ;
	// swap pointers so the last multiplication completes inside the A variable.
	simple_int_to_cint(A, 1);
	for (a = 0, b *= b; a < qs->s.values.subtract_one; ++a) {
		if (a & 1) i = b / (i + qs->poly.min) - (qs_sm) rand_upto(10);
		else i = qs->poly.span_half + (qs_sm) rand_upto(qs->poly.span_half) + qs->poly.min;
		for (j = 0; j < a; j = i == qs->s.data[j].prime_index ? ++i, 0 : j + 1);
		qs->s.data[a].prime_index = i; // the selected divisor of A wasn't already present in the product.
		simple_int_to_cint(Y, qs->base.data[i].num);
		cint_mul(A, Y, X), TMP = A, A = X, X = TMP ;
	}
	// it uses a prime number from the factor base that complete A, which must be close to D.
	cint_div(qs->calc, D, A, X, Y);
	const qs_sm d_over_a = (qs_sm) simple_cint_to_int(X);
	for (i = qs->base.data[0].num != 2 ; qs->base.data[i].num <= d_over_a; ++i);
	for (j = 0; j < qs->s.values.subtract_one; j = i == qs->s.data[j].prime_index ? ++i, 0 : j + 1);
	qs->s.data[qs->s.values.subtract_one].prime_index = i ;
	simple_int_to_cint(Y, qs->base.data[i].num);
	cint_mul(A, Y, X); // generated A values should always be distinct, and "A" will no longer change.
	assert(X == &qs->poly.A);
}

void iteration_part_2(qs_sheet * qs, const cint * A, cint * B) {
	cint *X = qs->vars.TEMP, *PRIME = X + 1, *Y = X + 2, *R = X + 3;
	qs_sm i, *pen = qs->s.A_indexes;
	cint_erase(B);
	for (i = 0; i < qs->s.values.defined; ++i) {
		const qs_sm prime = qs->base.data[qs->s.data[i].prime_index].num;
		// write [index of prime number, power] of the A factors into buffer.
		*pen++ = qs->s.data[i].prime_index, *pen++ = 1;
		qs->s.data[i].prime_squared = (qs_md)prime * (qs_md)prime ;
		simple_int_to_cint(PRIME, prime);
		cint_div(qs->calc, A, PRIME, X, R), assert(R->mem == R->end); // div exact.
		cint_div(qs->calc, X, PRIME, Y, R);
		qs->s.data[i].A_over_prime_mod_prime = (qs_sm) simple_cint_to_int(R);
		qs_md x = modular_inverse(qs->s.data[i].A_over_prime_mod_prime, prime);
		x = x * qs->base.data[qs->s.data[i].prime_index].kN_sqrt_mod_prime % prime;
		simple_int_to_cint(X, x > prime >> 1 ? prime - x : x);
		cint_mul(A, X, Y);
		cint_div(qs->calc, Y, PRIME, &qs->s.data[i].B_terms, R), assert(R->mem == R->end); // div exact.
		cint_addi(B, &qs->s.data[i].B_terms);
	}
}

void iteration_part_3(qs_sheet * qs, const cint * A, const cint * B) {
	cint *Q = qs->vars.TEMP, *R = Q + 1, *PRIME = Q + 2;
	qs_md i, j, x, y;
	for (i = 0; i < qs->base.length; ++i) {
		// prepare the "roots" and "A_inv_double_value_B_terms". The algorithm will be able
		// to fill 2 ** (s - 3) sieves by using these values and adding "prime sizes".
		const qs_sm prime = qs->base.data[i].num;
		simple_int_to_cint(PRIME, prime);
		cint_div(qs->calc, A, PRIME, Q, R);
		const qs_sm a_mod_prime = (qs_sm) simple_cint_to_int(R) ;
		cint_div(qs->calc, B, PRIME, Q, R) ;
		const qs_sm b_mod_prime = (qs_sm) simple_cint_to_int(R) ;
		const qs_sm a_inv_double_value = modular_inverse(a_mod_prime, prime) << 1 ;
		// Arithmetic shifts "<<" and ">>" performs multiplication or division by powers of two.
		x = y = prime;
		x += qs->base.data[i].kN_sqrt_mod_prime;
		y -= qs->base.data[i].kN_sqrt_mod_prime;
		x -= b_mod_prime;
		x *= a_inv_double_value >> 1;
		y *= a_inv_double_value ;
		x += qs->m.length_half ;
		x %= prime ;
		y += x ;
		y %= prime ;
		qs->base.data[i].root[0] = (qs_sm) x ;
		qs->base.data[i].root[1] = (qs_sm) y ;
		for (j = 0; j < qs->s.values.defined; ++j) {
			cint_div(qs->calc, &qs->s.data[j].B_terms, PRIME, Q, R);
			const qs_md b_term = simple_cint_to_int(R);
			qs->s.data[j].A_inv_double_value_B_terms[i] = (qs_sm)(a_inv_double_value * b_term % prime);
		}
	}
	// The next function operates over "B_terms" multiplied by 2.
	for (i = 0; i < qs->s.values.defined; cint_left_shifti(&qs->s.data[i++].B_terms, 1));
}

qs_sm iteration_part_4(const qs_sheet * qs, const qs_sm nth_curve, qs_sm ** corr, cint *B) {
	qs_sm i, gray_act; // the Gray code in "nth_curve" indicates which "B_term" to consider.
	for (i = 0; nth_curve >> i & 1; ++i);
	if (gray_act = (nth_curve >> i & 2) != 0, gray_act)
		cint_addi(B, &qs->s.data[i].B_terms) ;
	else // and which action to perform.
		cint_subi(B, &qs->s.data[i].B_terms) ;
	*corr = qs->s.data[i].A_inv_double_value_B_terms;
	return gray_act; // B values generated here should always be distinct.
}

void iteration_part_5(qs_sheet *  qs, const cint * kN, const cint * B) {
	cint *P = qs->vars.TEMP, *Q = P + 1, *R_kN = P + 2, *R_B = P + 3, *TMP = P + 4;
	for (qs_sm a = 0; a < qs->s.values.defined; ++a) {
		const qs_sm i = qs->s.data[a].prime_index;
		const qs_tmp prime = qs->base.data[i].num ;
		simple_int_to_cint(P, qs->s.data[a].prime_squared);
		cint_div(qs->calc, B, P, Q, R_B);
		cint_div(qs->calc, kN, P, Q, R_kN);
		if (B->nat < 0) cint_addi(R_B, P); // if B is negative.
		const qs_tmp rem_b = (qs_tmp) simple_cint_to_int(R_B);
		const qs_tmp rem_kn = (qs_tmp) simple_cint_to_int(R_kN);
		qs_tmp s ; // the remainders are modulo the prime number squared.
		if (rem_b >> 31){
			// the multiplication overflows.
			cint_mul(R_B, R_B, TMP);
			cint_subi(TMP, R_kN);
			simple_int_to_cint(P, (qs_md) prime);
			cint_div(qs->calc, TMP, P, Q, R_B);
			s = (qs_tmp) simple_cint_to_int(Q);
			if (Q->nat < 0) s = -s ;
		} else {
			// the system multiplies.
			s = rem_b * rem_b - rem_kn;
			s /= prime ;
		}
		//
		qs_tmp bezout = (rem_b % prime) * (qs_tmp) qs->s.data[a].A_over_prime_mod_prime ;
		bezout = (qs_tmp) modular_inverse((qs_sm) (bezout % prime), (qs_sm) prime);
		//
		s = (qs_tmp) qs->m.length_half - s * bezout ;
		s %= prime ;
		s += (s < 0) * prime ;
		qs->base.data[i].root[0] = (qs_sm) s;
		qs->base.data[i].root[1] = (qs_sm) -1;
	}
}

void iteration_part_6(qs_sheet *qs, const cint *kN, const cint *A, const cint *B, cint *C) {
	cint *TMP = qs->vars.TEMP, *R = TMP + 1;
	cint_mul(B, B, TMP); // (B * B) % A = kN % A
	cint_subi(TMP, kN); // C = (B * B - kN) / A
	cint_div(qs->calc, TMP, A, C, R), assert(R->mem == R->end); // div exact.
}

void iteration_part_7(qs_sheet * qs, const qs_sm gray_addi, const qs_sm * restrict corr) {
	// Sieve for larger prime numbers.
	memset(qs->m.sieve, 0, qs->m.length * sizeof(*qs->m.sieve));
	memset(qs->m.flags, 0, qs->base.length * sizeof(*qs->m.flags));
	uint8_t * restrict end = qs->m.sieve + qs->m.length, *p_0, *p_1;
	for(qs_sm i = qs->iterative_list[3], j = qs->iterative_list[4]; i < j; ++i) {
			const qs_sm prime = qs->base.data[i].num, size = qs->base.data[i].size, co = gray_addi ? prime - corr[i] : corr[i];
			qs->base.data[i].root[0] += co; if (qs->base.data[i].root[0] >= prime) qs->base.data[i].root[0] -= prime;
			qs->base.data[i].root[1] += co; if (qs->base.data[i].root[1] >= prime) qs->base.data[i].root[1] -= prime;
			p_0 = qs->m.sieve + qs->base.data[i].root[0];
			p_1 = qs->m.sieve + qs->base.data[i].root[1];
			for (; end > p_0 && end > p_1;)
				*p_0 += size, p_0 += prime, *p_1 += size, p_1 += prime;
			*p_0 += (end > p_0) * size, *p_1 += (end > p_1) * size;
		}
	for(qs_sm i = qs->iterative_list[4], j = qs->iterative_list[5]; i < j; ++i){
		const qs_sm prime = qs->base.data[i].num, size = qs->base.data[i].size, co = gray_addi ? prime - corr[i] : corr[i] ;
		qs->base.data[i].root[0] += co; if (qs->base.data[i].root[0] >= prime) qs->base.data[i].root[0] -= prime;
		qs->base.data[i].root[1] += co; if (qs->base.data[i].root[1] >= prime) qs->base.data[i].root[1] -= prime;
		for(p_0 = qs->m.sieve + qs->base.data[i].root[0]; end > p_0; qs->m.flags[i] |= 1 << ((p_0 - qs->m.sieve) & 7), *p_0 += size, p_0 += prime);
		for(p_1 = qs->m.sieve + qs->base.data[i].root[1]; end > p_1; qs->m.flags[i] |= 1 << ((p_1 - qs->m.sieve) & 7), *p_1 += size, p_1 += prime);
	}
}

void iteration_part_8(qs_sheet * qs, const qs_sm gray_addi, const qs_sm *  corr) {
	// Sieving means taking an interval [−M/2, +M/2] and determining for
	// which X in [−M/2, +M/2] a given prime number divides AX^2 + 2BX + C.
	uint8_t * chunk_begin = qs->m.sieve, *chunk_end = chunk_begin;
	uint8_t * sieve_end = chunk_begin + qs->m.length ;
	qs_sm *buffer = qs->buffer[0], walk_idx, * walk = buffer;
	// Since the previous function, the check is performed for the prime numbers of the factor base.
	for(qs_sm i = 0; i < qs->iterative_list[3]; ++i)
		if (qs->base.data[i].root[1] != (qs_sm) -1) {
			*walk++ = i ; // the current prime number isn't a factor of A.
			const qs_sm prime = qs->base.data[i].num, co = gray_addi ? prime - corr[i] : corr[i] ;
			qs->base.data[i].root[0] += co; if (qs->base.data[i].root[0] >= prime) qs->base.data[i].root[0] -= prime;
			qs->base.data[i].root[1] += co; if (qs->base.data[i].root[1] >= prime) qs->base.data[i].root[1] -= prime;
			qs->m.positions[0][i] = chunk_begin + qs->base.data[i].root[0];
			qs->m.positions[1][i] = chunk_begin + qs->base.data[i].root[1];
		}
	for(walk_idx = 0; buffer[walk_idx] < qs->iterative_list[1]; ++walk_idx);
	do{ // iterates step by step until the entire sieve is filled.
		walk = buffer + walk_idx ;
		chunk_end = chunk_end + qs->m.cache_size < sieve_end ? chunk_end + qs->m.cache_size : sieve_end;
		do{
			const qs_sm size = qs->base.data[*walk].size, prime = qs->base.data[*walk].num, times = 4 >> (*walk > qs->iterative_list[2]) ;
			uint8_t ** const p_0 = qs->m.positions[0] + *walk, ** const p_1 = qs->m.positions[1] + *walk;
			const ptrdiff_t diff = *p_1 - *p_0 ;
			for(const uint8_t * const bound = chunk_end - prime * times; bound > *p_0;)
				for(qs_sm i = 0; i < times; ++i)
					**p_0 += size, *(*p_0 + diff) += size, *p_0 += prime;
			for(; *p_0 < chunk_end && *p_0 + diff < chunk_end;)
				**p_0 += size, *(*p_0 + diff) += size, *p_0 += prime;
			*p_1 = *p_0 + diff ;
			if (*p_0 < chunk_end) **p_0 += size, *p_0 += prime;
			if (*p_1 < chunk_end) **p_1 += size, *p_1 += prime;
		} while(*++walk);
	} while(chunk_begin = chunk_end, chunk_begin < sieve_end);
	memset(qs->buffer[0], 0, (walk - qs->buffer[0]) * sizeof(*walk));
}

int qs_register_factor(qs_sheet * qs){
	// the current function have "permission" to update N, not kN.
	// it returns -1 if the factorization is completed, when N = 1.
	cint * F = &qs->vars.FACTOR ;
	int i, res = h_cint_compare(F, &qs->constants.ONE) > 0 && h_cint_compare(F, &qs->vars.N) < 0 ;
	if (res) {
		F->nat = 1 ; // absolute value of the factor.
		const struct avl_node *node = avl_at(&qs->uniqueness[2], F);
		if (qs->uniqueness[2].affected)
			for (i = 0; i < 2 && qs->n_bits != 1; ++i) {
				const int is_prime = cint_is_prime(qs->calc, F, -1);
				if (is_prime) {
					const int power = (int) cint_remove(qs->calc, &qs->vars.N, F);
					assert(power); // 200-bit RSA take about 10,000,000+ "duplications".
					fac_push(qs->caller, F, 1, power, 0);
					++qs->divisors.total_primes;
					// if "n_bits" is greater than one the functions must continue to search factors.
					qs->n_bits = (qs_sm) cint_count_bits(&qs->vars.N);
					if (qs->n_bits == 1) res = -1;
					else cint_dup(F, &qs->vars.N);
				} else if (i++ == 0)
					// the given number isn't prime, it's only registered with the known divisors of N.
					qs->divisors.data[qs->divisors.length++] = node->key;
			}
	}
	return res ;
}

void register_relations(qs_sheet * qs, const cint * A, const cint * B, const cint * C) {
	cint *  TMP = qs->vars.TEMP, * K = &qs->vars.KEY, * V = &qs->vars.VALUE ;
	qs_sm m_idx, idx, mod;
	// iterates the values of X in [-M/2, +M/2].
	for (m_idx = 0; m_idx < qs->m.length; ++m_idx)
		// the trick that permit to step faster using typecast and 0xC0C0C0C0C0C0C0C0 isn't implemented.
		if (qs->m.sieve[m_idx] >= qs->threshold.value) {
			// over the threshold, compute f(X) and check candidate for smoothness.
			simple_int_to_cint(&qs->vars.X, m_idx);
			cint_subi(&qs->vars.X, &qs->constants.M_HALF); // X = "current index" - M/2
			cint_mul(A, &qs->vars.X, TMP); // TMP = AX
			cint_addi(TMP, B); // TMP = AX + B
			cint_dup(qs->vars.MY, TMP);
			cint_dup(K, TMP); // Key = copy of AX + B
			cint_addi(TMP, B); // TMP = AX + 2B
			cint_mul(TMP, &qs->vars.X, V); // V = AX^2 + 2BX
			cint_addi(V, C); // Value = f(X) = AX^2 + 2BX + C
			// it should hold that kN = (Key * Key) - (A * Value)
			V->nat = 1 ; // absolute value
			qs_sm target_bits = (qs_sm) cint_count_bits(V) - qs->error_bits;
			qs_sm removed_bits = 0, * restrict pen = qs->buffer[1];
			//  writes the pairs [index of the prime number, power found in V].
			if (qs->base.data[0].num != 1) {
				simple_int_to_cint(TMP, qs->base.data[0].num);
				*pen++ = 0; // remove powers of the multiplier.
				*pen = (qs_sm) cint_remove(qs->calc, V, TMP);
				if (*pen) removed_bits += *pen++ * qs->base.data[0].size; else --pen;
			}
			for (idx = 1; idx < qs->iterative_list[1]; ++idx)
				if (qs->base.data[idx].root[1] == (qs_sm) -1 || (mod = m_idx % qs->base.data[idx].num, mod == qs->base.data[idx].root[0] || mod == qs->base.data[idx].root[1])) {
					simple_int_to_cint(TMP, qs->base.data[idx].num);
					// for a given prime number of the factor base, "remove" returns
					// the numbers of powers that was present in V, and V is updated.
					*pen++ = idx;
					*pen = (qs_sm) cint_remove(qs->calc, V, TMP);
					if (*pen) removed_bits += *pen++ * qs->base.data[idx].size; else --pen;
				}
			if (removed_bits + qs->m.sieve[m_idx] >= target_bits) {
				// there is a chance to register a new relation.
				for (removed_bits = 0, target_bits = qs->m.sieve[m_idx]; idx < qs->iterative_list[4] && removed_bits < target_bits; ++idx)
					if (qs->base.data[idx].root[1] == (qs_sm) -1 || (mod = m_idx % qs->base.data[idx].num, mod == qs->base.data[idx].root[0] || mod == qs->base.data[idx].root[1])) {
						simple_int_to_cint(TMP, qs->base.data[idx].num);
						*pen++ = idx;
						*pen = (qs_sm) cint_remove(qs->calc, V, TMP);
						if (*pen) removed_bits += *pen++ * qs->base.data[idx].size; else --pen;
					}
				for (const uint8_t mask = 1 << (m_idx & 7); idx < qs->iterative_list[5] && removed_bits < target_bits; ++idx)
					if (qs->m.flags[idx] & mask)
						if (mod = m_idx % qs->base.data[idx].num, mod == qs->base.data[idx].root[0] || mod == qs->base.data[idx].root[1]) {
							simple_int_to_cint(TMP, qs->base.data[idx].num);
							*pen++ = idx;
							*pen = (qs_sm) cint_remove(qs->calc, V, TMP);
							if (*pen) removed_bits += *pen++ * qs->base.data[idx].size; else --pen;
						}
				const qs_sm * restrict const prime_indexes_and_powers[4] = {
						qs->s.A_indexes, // really factoring A * f(X), so commit outstanding A factors.
						qs->s.A_indexes + qs->s.values.double_value,
						qs->buffer[1],
						pen,
				};
				if (h_cint_compare(V, &qs->constants.ONE) == 0)
					register_relation_kind_1(qs, K, prime_indexes_and_powers);
					// #1 warwick.ac.uk, not a relation, perhaps a partial ?
				else if (h_cint_compare(V, &qs->constants.LARGE_PRIME) < 0)
					register_relation_kind_2(qs, K, V, prime_indexes_and_powers);
			}
		}
}

void register_relation_kind_1(qs_sheet * qs, const cint * KEY, const qs_sm * const restrict args[4]) {
	struct avl_node *node = avl_at(&qs->uniqueness[0], KEY);
	if (node->value)
		return; // duplicates at this stage are ignored.
	struct qs_relation * rel = qs->mem.now;
	qs_sm i, j ;
	// a relation must be swappable for Lanczos Block reducing.
	qs->mem.now = rel + 1 ;
	rel->X = node->key; // constant X is stored by the node key.
	rel->Y.data = qs->mem.now; // Y data has a bounded length.
	const size_t y_length = (args[1] - args[0] + args[3] - args[2]) >> 1 ;
	rel->axis.Z.data = rel->Y.data + y_length; // writes Z ahead.
	for (i = 0; i < 4;) {
		// processes the given column arrays.
		const qs_sm * restrict idx = args[i++], * restrict const end_index = args[i++];
		for (; idx < end_index; idx += 2) {
			const qs_sm power = *(idx + 1) ;
			if (power & 1) {
				// remove from Y the indexes of the prime numbers that are already listed (factors of A).
				for (j = 0; j < rel->Y.length && rel->Y.data[j] != *idx; ++j);
				if (j == rel->Y.length) // add, the index wasn't present.
					rel->Y.data[rel->Y.length++] = *idx;
				else // or remove.
					memmove(rel->Y.data + j, rel->Y.data + j + 1, (--rel->Y.length - j) * sizeof(*rel->Y.data));
			}
			for (j = 0; j < power; ++j)
				rel->axis.Z.data[rel->axis.Z.length++] = *idx;
		}
	}
	qs->mem.now = rel->axis.Z.data + rel->axis.Z.length; // Z length is known.
	int verified = 0 ;
	if (rel->Y.length > qs->s.values.defined) {
		// it often passes.
		cint *A = qs->vars.TEMP, *B = A + 1;
		simple_int_to_cint(A, 1);
		for (j = 0; j < rel->axis.Z.length; ++j) {
			simple_int_to_cint(B, qs->base.data[rel->axis.Z.data[j]].num);
			cint_mul_modi(qs->calc, A, B, &qs->constants.kN);
		}
		cint_mul_mod(qs->calc, rel->X, rel->X, &qs->constants.kN, B);
		verified = !cint_compare(A, B) || (cint_addi(A, B), !cint_compare(A, &qs->constants.kN));
	}
	if (verified){
		node->value = qs->relations.data[qs->relations.length.now] = rel;
		qs->mem.now = rel->axis.Z.data + rel->axis.Z.length;
		rel->id = ++qs->relations.length.now; // Keep the relation
	} else {
		char * open = (char*) rel, * close = qs->mem.now ;
		qs->mem.now = memset(open, 0, close - open); // Throw
	}
}

void register_relation_kind_2(qs_sheet * qs, const cint * KEY, const cint * VALUE, const qs_sm * const restrict args[4]) {
	// this function can register many relations during long factorizations.
	if (qs->kn_bits < 150)
		return; // but does not have time to "start" with a small N.
	// searches 2 different KEY sharing the same VALUE.
	struct avl_node *node = avl_at(&qs->uniqueness[1], VALUE);
	struct qs_relation *old, *new;
	cint * BEZOUT = 0;
	old = node->value;
	if (old) {
		if (old->id) return; // the "partials" have no id.
		if (old->X == 0) return; // the value is already marked as "ignored".
		if (old->axis.next) return; // accepting all "next" without caring would reduce the "chance".
		for (;old && h_cint_compare(KEY, old->X); old = old->axis.next);
		if (old) return; // same KEY already registered.
		old = node->value;
		if (old->axis.next == 0) {
			cint *A = qs->vars.TEMP, *B = A + 1;
			if (qs->multiplier != 1)
				if (cint_gcd(qs->calc, VALUE, &qs->constants.MULTIPLIER, A), *A->mem != 1){
					old->X = 0; // shouldn't be related so close to the multiplier.
					return;
				}
			// the given duplicate isn't 3rd, 4th .. it's the 2nd ... so compute BEZOUT
			cint_modular_inverse(qs->calc, VALUE, &qs->constants.kN, A);
			if (A->mem == A->end) {
				old->X = 0; // no solution to the linear congruence.
				cint_gcd(qs->calc, VALUE, &qs->constants.kN, &qs->vars.FACTOR);
				cint_div(qs->calc, &qs->vars.N, &qs->vars.FACTOR, A, B);
				if (B->mem == B->end) // found a small factor of N ?
					qs_register_factor(qs);
				return; // nothing.
			} else
				BEZOUT = A;
		}
	}

	new = mem_aligned(qs->mem.now);
	qs->mem.now = new + 1;
	new->X = qs->mem.now;

	if (BEZOUT) {
		// BEZOUT is stored directly after X, the new "partial" become the root of the linked list.
		qs->mem.now = new->X + 2;
		simple_dup_cint(new->X, KEY, &qs->mem.now);
		simple_dup_cint(new->X + 1, BEZOUT, &qs->mem.now);
		new->axis.next = old, node->value = new = old = new;
	} else {
		qs->mem.now = new->X + 1;
		simple_dup_cint(new->X, KEY, &qs->mem.now);
		if (old) {
			for (; old->axis.next; old = old->axis.next);
			old->axis.next = new, old = node->value;
		} else node->value = new;
	}

	// data buffered isn't persistent, it may later be needed, so it's copied.
	qs_sm * data = new->Y.data = mem_aligned(qs->mem.now);
	new->Y.length = (qs_sm) (args[1] - args[0]);
	memcpy(data, args[0], new->Y.length * sizeof(*data));
	memcpy(data + new->Y.length, args[2], (args[3] - args[2]) * sizeof(*data));
	new->Y.length += (qs_sm) (args[3] - args[2]);
	qs->mem.now = new->Y.data + new->Y.length;

	if (old) {
		BEZOUT = old->X + 1 ; // the modular inverse was stored here.
		cint_mul_mod(qs->calc, new->X, BEZOUT, &qs->constants.kN, &qs->vars.CYCLE);
		do {
			if (old != new) {
				// combines, it registers a regular relation using the 2 buffers.
				cint_mul_mod(qs->calc, &qs->vars.CYCLE, old->X, &qs->constants.kN, &qs->vars.KEY);
				qs_sm * restrict begin = qs->buffer[0], * restrict pen = begin;
				data = memset(qs->buffer[1], 0, qs->base.length * sizeof(*data));
				for (qs_sm i = 0; i < new->Y.length; i += 2)
					data[new->Y.data[i]] += new->Y.data[i + 1];
				for (qs_sm i = 0; i < old->Y.length; i += 2)
					data[old->Y.data[i]] += old->Y.data[i + 1];
				for (qs_sm i = 0; i < qs->base.length; ++i)
					if (data[i]) // writes [index of the prime number, power]
						*pen++ = i, *pen++ = data[i];
				args = (const qs_sm * restrict const[4]){ begin, pen, 0, 0, };
				register_relation_kind_1(qs, &qs->vars.KEY, args);
				memset(begin, 0, (char *) pen - (char *) begin); // zeroed.
			}
		} while ((old = old->axis.next));
		// the linked list can handle 2+ entries, but their more complex combinations isn't implemented.
	}
}

void finalization_part_1(qs_sheet * qs, const uint64_t * restrict const lanczos_answer) {
	const uint64_t mask = *lanczos_answer, * restrict null_rows = lanczos_answer + 1;
	// Lanczos "linear algebra" answer is simply "mask followed by null_rows", with read-only.
	if (mask == 0 || null_rows == 0)
		return;
	cint *A = qs->vars.TEMP, *B = A + 1, *C = A + 2, *D = A + 3;
	qs_sm * restrict power_of_primes;
	for(qs_sm row = 0; row < 64 && qs->n_bits != 1; ++row)
		if (mask >> row & 1){
			simple_int_to_cint(A, 1), simple_int_to_cint(B, 1), simple_int_to_cint(C, 1);
			power_of_primes = memset(qs->buffer[1], 0, qs->base.length * sizeof(*power_of_primes));
			for (qs_sm i = 0; i < qs->relations.length.now; ++i)
				if (null_rows[i] >> row & 1) {
					// algorithm must retrieve both the X and Z values, which
					// are related to the Y values provided to Lanczos Block.
					const struct qs_relation * restrict const rel = qs->relations.data[i];
					cint_mul_modi(qs->calc, A, rel->X, &qs->vars.N);
					for (qs_sm j = 0; j < rel->axis.Z.length; ++j)
						++power_of_primes[rel->axis.Z.data[j]];
				}
			for (qs_sm i = 0; i < qs->base.length; ++i)
				if (power_of_primes[i]){
					// powers are even ... square root ...
					simple_int_to_cint(B, qs->base.data[i].num);
					simple_int_to_cint(D, power_of_primes[i] >> 1);
					cint_pow_modi(qs->calc, B, D, &qs->vars.N);
					cint_mul_modi(qs->calc, C, B, &qs->vars.N);
				}
			h_cint_subi(C, A);
			if (C->mem != C->end) {
				cint_gcd(qs->calc, &qs->vars.N, C, &qs->vars.FACTOR);
				// The RSA-100 number has been factored by the software in 2022.
				// A. K. Lenstra announced this factorization completed in 1991.
				qs_register_factor(qs);
			}
		}
}

void finalization_part_2(qs_sheet * qs) {
	if (qs->n_bits == 1)
		return;

	// Algorithm checked Lanczos answer but N is still greater than one.
	// Perform basic checks until no new divisor can be discovered.

	cint * F = &qs->vars.FACTOR, **di = qs->divisors.data, *Q = qs->vars.TEMP, *R = Q + 1 ;
	qs_sm i, j, k, count = 0 ;

	for(i = qs->divisors.processing_index, j = qs->divisors.length; i < j; ++i)
		for (k = 1 + i; k < j; ++k)
			if (cint_gcd(qs->calc, di[i], di[k], F), qs_register_factor(qs) == -1)
				return; // register "gcd of new divisors with old divisors"
	do {
		for (; i = j, j = qs->divisors.length, i != j;)
			for (; i < j; ++i)
				for (k = 0; k < i; ++k)
					if (cint_gcd(qs->calc, di[i], di[k], F), qs_register_factor(qs) == -1)
						return; // register "gcd of new divisors with old divisors"

		for (i = qs->divisors.processing_index; i < j; ++i)
			if (fac_any_root_check(qs->caller, di[i], &qs->vars.FACTOR, R))
				if (qs_register_factor(qs) == -1)
					return; // register "perfect power of new divisors"

		if (count != qs->divisors.total_primes)
			if (fac_any_root_check(qs->caller, &qs->vars.N, &qs->vars.FACTOR, R))
				if (qs_register_factor(qs) == -1)
					return; // register "prefect root of N"

		count = qs->divisors.total_primes ;
		qs->divisors.processing_index = j ;
	} while(j != qs->divisors.length);

}

int finalization_part_3(qs_sheet * qs) {
	// Usually do nothing, (N = 1)
	// Otherwise push something non-trivial to the caller's routine
	int res = qs->n_bits == 1 ;
	if (res == 0){
		cint * F = &qs->vars.FACTOR;
		for(qs_sm i = 0; i < qs->divisors.length; ++i) {
			if (qs->caller->params->silent == 0) {
				char *str = cint_to_string(qs->divisors.data[i], 10);
				printf("- quadratic sieve can't silently ignore %s\n", str);
				free(str);
			}
			if (h_cint_compare(&qs->vars.N, qs->divisors.data[i]) > 0) {
				const int power = (int) cint_remove(qs->calc, &qs->vars.N, qs->divisors.data[i]);
				if (power) {
					assert(power == 1);
					fac_push(qs->caller, F, 0, 1, 1);
				}
			}
		}
		res = h_cint_compare(&qs->vars.N, &qs->caller->number->cint) ;
		if (res) // res is true if QS was able to decompose N.
			if (h_cint_compare(&qs->vars.N, &qs->constants.ONE))
				fac_push(qs->caller, &qs->vars.N, -1, 1, 1);
	}
	return res;
}

// "Pure C99 factorizer using self-initialising Quadratic Sieve."
// compiling with "gcc -Wall -pedantic -O3 main.c -o qs" can speed up by a factor 2 or 3.
