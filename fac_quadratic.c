//  GNU General Public License

//      As an undergraduate student, this project is part of my computer science + maths training.
//          - This software proposition is from Michel Leonard (student at Université de Franche-Comté, Mon, 11 Jul 2022)
//          - There is no guarantee of any kind on the software
//          - C code is shared under the terms of the GNU General Public License
//          - The main mathematical and logical inspiration source is located at :
//              http://web.mit.edu/sage/export/flintqs-0.0.20070817/QS.cpp - GNU General Public License

//  This software implementation would have been impossible without "FLINT: Fast Library for Number Theory" maintained by William Hart.

// Quadratic sieve main routine
static int quadratic_sieve(fac_caller *caller) {
	// The routine is invoked by a caller, and :
	// - must factor the caller's "number"
	// - must register answers to the caller's routines
	// - can use caller's resources (variables + configurations)
	// - resources are enough to copy 400-bit into variables

	int limit = (int) caller->params->qs_limit;
	if (limit == 0) limit = 220; // the default qs_limit is such that the QS answer in seconds.
	if (caller->number->bits < 64 || caller->number->bits > limit) return 0;

	qs_sheet qs = {0};
	preparation_part_1(&qs, caller);
	preparation_part_2(&qs);
	preparation_part_3(&qs);
	// adjustor and multiplier are known, parametrize.
	qs_parametrize(&qs);
	preparation_part_4(&qs);
	preparation_part_5(&qs);
	qs.the.min = preparation_part_6(&qs, &qs.variables.D);
	do {
		do {
			get_started_iteration(&qs);
			iteration_part_1(&qs, &qs.variables.D, &qs.variables.A);
			iteration_part_2(&qs, &qs.variables.A, &qs.variables.B);
			iteration_part_3(&qs, &qs.variables.A, &qs.variables.B);
			for (qs_sm i = 0, addi, *corr; i < qs.curves.max && qs.n_bits != 1; ++i, ++qs.curves.now) {
				addi = iteration_part_4(&qs, i, &corr, &qs.variables.B);
				iteration_part_5(&qs, &qs.constants.kN, &qs.variables.B);
				iteration_part_6(&qs, &qs.constants.kN, &qs.variables.A, &qs.variables.B, &qs.variables.C);
				iteration_part_7(&qs, addi, corr);
				iteration_part_8(&qs, addi, corr);
				register_relations(&qs, &qs.variables.A, &qs.variables.B, &qs.variables.C);
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
static inline int inner_continuation_condition(qs_sheet *qs) {
	int res = 1;
	res &= qs->n_bits != 1 ; // the bit count of N is updated when algorithm registered a factor.
	res &= qs->relations.length.now < qs->relations.length.needs; // the condition.
	if (qs->caller->params->silent == 0) {
		const double rel_begin = (double) qs->relations.length.now, rel_end = (double) qs->relations.length.needs ;
		fac_display_progress("Quadratic sieve", 100. * rel_begin / rel_end); // progress isn't linear
	}
	return res;
}

// Quadratic sieve main condition 2 : user configuration
// return to sieving or stop the algorithm ? res is the answer.
static inline int outer_continuation_condition(qs_sheet *qs) {
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
static inline qs_sm linear_param_resolution(const double v[][2], const qs_sm point) {
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

static inline void qs_parametrize(qs_sheet *qs) {
	const qs_sm bits = (qs_sm) cint_count_bits(qs->caller->vars); // kN
	qs->kn_bits = bits; // input was adjusted so there is at least 115-bit.

	// params as { bits, value } take the extremal value if bits exceed.
	static const double param_base_size [][2]= { {110, 800}, {130, 1500}, {210, 4500}, {240, 9000}, {250, 15000}, {290, 25000}, {0} };
	qs->base.length = linear_param_resolution(param_base_size, bits);

	static const double param_laziness [][2]= {{110, 90}, {190, 100}, {220, 100}, {250, 120}, {0} };
	// collecting more/fewer relations than recommended (used to verify "sieve again" feature).
	qs->relations.length.needs = qs->base.length * linear_param_resolution(param_laziness, bits) / 100 ;

	static const double param_m_value [][2]= { {110, 1}, {190, 4}, {0} };
	qs->m.double_value = (qs->m.value = 31744 * linear_param_resolution(param_m_value, bits)) << 1;

	static const double param_error [][2]= { {110, 13}, {300, 33}, {0} };
	qs->error_bits = linear_param_resolution(param_error, bits);

	static const double param_threshold [][2]= { {110, 63}, {220, 78}, {300, 102}, {0} };
	qs->threshold.value = linear_param_resolution(param_threshold, bits);

	static const double param_alloc [][2]= { {1e3, 2}, {3e3, 8}, {5e3, 20}, {15e3, 80}, {25e3, 140}, {0} };
	qs->mem.bytes_allocated = linear_param_resolution(param_alloc, qs->base.length) << 20 ; // { Base size, Megabytes }

	qs->sieve_again_perms = 3; // Sieve again up to 3 times before giving up

	// Other parameters
	qs->s.values.double_value = (qs->s.values.defined = (qs->s.values.subtract_one = bits / 28) + 1) << 1;
	qs->curves.max = (1 << qs->s.values.subtract_one) - 2;

	{
		// Iterative list (not always fully ordered)
		qs->iterative_list[0] = 1; // one
		static const double param_first_prime [][2]= { {170, 8}, {210, 12}, {300, 30}, {0} };
		qs->iterative_list[1] = linear_param_resolution(param_first_prime, bits); // first
		const qs_sm large_base = qs->base.length > 5000 ? 5000 : qs->base.length;
		qs->iterative_list[2] = large_base >> 2; // medium
		qs->iterative_list[3] = large_base >> 1; // mid
		qs->iterative_list[4] = large_base; // sec
		qs->iterative_list[5] = qs->base.length; // factor base size
		static const double param_large_prime [][2]= { {110, 25e4}, {170, 1e6}, {240, 3e7}, {280, 2e8}, {0} };
		qs->iterative_list[6] =  linear_param_resolution(param_large_prime, bits); // large prime
	}

	qs->the.span_half = (qs->the.span = qs->base.length / qs->s.values.defined / qs->s.values.defined / 2) >> 1;
}

// Quadratic sieve source : algorithm
static inline void preparation_part_1(qs_sheet *qs, fac_caller *caller) {
	// initializing with the caller's resources
	qs->caller = caller;
	qs->calc = caller->calc;
	qs->n_bits = caller->number->bits;
}

static inline void preparation_part_2(qs_sheet *qs) {
	cint * N = &qs->caller->number->cint, * kN = qs->caller->vars, *ADJUSTOR = kN + 1 ;
	// Input is "transparently" adjusted by a prime to measure at least 115-bit.
	static const int prime_generator[] = {
			9, 7, 5, 3, 17, 27, 3, 1, 29, 3, 21, 7, 17, 15,
			9, 43, 35, 15, 29, 3, 11, 3, 11, 15, 17, 25, 53,
			31, 9, 7, 23, 15, 27, 15, 29, 7, 59, 15, 5, 21,
			69, 55, 21, 21, 5, 159, 3, 81, 9, 69, 131, 33, 15 };
	const qs_sm bits = qs->n_bits;
	if (bits < 115) {
		const qs_md adjustor = bits < 115 ? ((qs_md)1 << (124 - bits)) + prime_generator[115 - bits] : 1;
		simple_int_to_cint(ADJUSTOR, adjustor);
		cint_mul(N, ADJUSTOR, kN);
	} else
		cint_dup(kN, N);
}

static inline void preparation_part_3(qs_sheet *qs) {
	// the following part can speed up the factorization by a factor 2.5.
	// so it is possible to compare propositions.
	qs_sm mul = (qs_sm) qs->caller->params->qs_multiplier ;
	if (mul == 0) // like everyone, listen propositions.
		mul = errno ? preparation_part_3_proposition(qs) : preparation_part_3_original(qs);
	qs->knuth_schroeppel = mul ;
	// relations accumulate 1.46 times faster (on average)
	if (qs->knuth_schroeppel > 1){
		cint *kN = qs->caller->vars, *A = kN + 1, *B = kN + 2 ;
		simple_int_to_cint(A, qs->knuth_schroeppel);
		cint_dup(B, kN);
		cint_mul(A, B, kN);
	}
}

static inline qs_sm preparation_part_3_proposition(qs_sheet *qs) {
	// The function propose a Knuth-Schroeppel multiplier for N, intended to optimize QS runtime.
	cint * kN = qs->caller->vars, *A = kN + 1, *B = kN + 2, *C = kN + 3;
	static const int mul[] = {1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19, 21, 22, 23, 26, 29, 30, 31, 33, 34, 35, 37, 38, 39, 41, 42, 43, 46, 47, 51, 53, 55, 57, 58, 59, 61, 62, 65, 66, 67, 69, 70, 71, 73}, n_mul = sizeof(mul) / sizeof(*mul);
	int a, b, c;
	double score[n_mul], intake;
	b = (int) (*kN->mem & 7);
	for (a = 0; a < n_mul; ++a) {
		const int mod = b * mul[a] & 7;
		score[a] = -.5 * log_computation(mul[a]) + (mod & 1) * 0.3465736 * (mod == 1 ? 4 : mod == 5 ? 2 : 1);
	}
	// scan primes contributions until they are greater than a convenient qs_limit.
	b = 2048 ;
	for (a = 3; a < b; a += 2)
		if (is_prime_4669921(a)) {
			simple_int_to_cint(A, a);
			intake = log_computation(a) / (a - 1);
			cint_div(qs->calc, kN, A, B, C);
			const qs_sm d = (qs_sm) simple_cint_to_int(C); // kN mod prime
			for (c = 0; c < n_mul; ++c) {
				const qs_sm e = mul[c] * d % a ;
				score[c] += e ? power_modulo(e, (a - 1) >> 1, a) == 1 ? 2 * intake : 0 : intake;
			}
		}
	for (b = 1, c = 0; b < n_mul; ++b)
		if (score[b] > score[0])
			score[0] = score[c = b];
	return mul[c];
}

static inline qs_sm preparation_part_3_original(qs_sheet *qs) {
	// The function propose a Knuth-Schroeppel multiplier for N, intended to optimize QS runtime.
	static const qs_sm mul[] = {1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43}, n_mul = (qs_sm)(sizeof(mul) / sizeof(qs_sm));
	cint *kN = qs->caller->vars, *A = kN + 1, *B = kN + 2, *C = kN + 3;
	qs_sm a, b, c;
	double score[n_mul];
	const qs_md n_mod_8 = *kN->mem % 8;
	for (b = 0; b < n_mul; ++b) {
		qs_md calc = n_mod_8 * mul[b] % 8;
		if (calc == 1) score[b] = 1.38629436;
		else if (calc == 5) score[b] = 0.69314718;
		else score[b] = 0.34657359 ;
		score[b] -= log_computation((double) mul[b]) / 2.0;
	}
	// the time spent here can vary with the importance of the factor base.
	static const double max_prime_bound [][2]= {{120, 2e4}, {215, 11e4}, {240, 2e5}, {260, 35e5}, {0} };
	// another proposal would be b = 10,000.
	b = linear_param_resolution(max_prime_bound, qs->n_bits);
	for (a = 3; a < b ; a += 2)
		if (is_prime_4669921(a)) {
			const double intake = log_computation((double) a) / a;
			simple_int_to_cint(A, a);
			cint_div(qs->calc, kN, A, B, C);
			const int kronecker = kronecker_symbol((qs_sm) simple_cint_to_int(C), a);
			for (c = 0; c < n_mul; ++c)
				score[c] += intake + intake * kronecker * kronecker_symbol(mul[c], a);
		}
	for (b = 1, c = 0; b < n_mul; ++b)
		if (score[b] > score[0])
			score[0] = score[c = b];
	return mul[c];
}

static inline void preparation_part_4(qs_sheet *qs) {
	void *mem;
	mem = qs->mem.base = calloc(1, qs->mem.bytes_allocated);
	assert(mem);
	qs->caller->params->qs_rand_seed = 1 ;
	if (qs->caller->params->qs_rand_seed) srand(qs->rand_seed = qs->caller->params->qs_rand_seed);
	else qs->caller->params->qs_rand_seed = qs->rand_seed = add_rand_seed(&mem);

	// kN was computed into the caller's courtesy memory, now the QS has parametrized and "allocated"
	const size_t kn_size = qs->caller->vars[0].end - qs->caller->vars[0].mem + 1 ;
	// the quadratic sieve variables are able to temporarily hold at most kN ^ 2
	const size_t vars_size = kn_size << 1 ;
	const size_t buffers_size = qs->base.length + (qs->iterative_list[1] << 1);

	{
		// list of numbers used by the algorithm
		cint * const n[] = {
				&qs->variables.N,
				// polynomial
				&qs->variables.A,
				&qs->variables.B,
				&qs->variables.C,
				&qs->variables.D,
				// temporary
				&qs->variables.TEMP[0],
				&qs->variables.TEMP[1],
				&qs->variables.TEMP[2],
				&qs->variables.TEMP[3],
				&qs->variables.TEMP[4],
				// applicative
				&qs->variables.X,
				&qs->variables.KEY,
				&qs->variables.VALUE,
				&qs->variables.CYCLE,
				// a factor of N
				&qs->variables.FACTOR,
				// constants
				&qs->constants.kN,
				&qs->constants.ONE,
				&qs->constants.M,
				&qs->constants.SMALL_PRIME,
				&qs->constants.LARGE_PRIME,
				0,
		};
		for (int i = 0; n[i]; ++i) {
			n[i]->mem = n[i]->end = mem_aligned(mem) ;
			mem = n[i]->mem + (n[i]->size = vars_size);
		}
	}

	cint_dup(&qs->variables.N, &qs->caller->number->cint);
	cint_dup(&qs->constants.kN, qs->caller->vars);

	simple_int_to_cint(&qs->constants.ONE, 1);
	simple_int_to_cint(&qs->constants.SMALL_PRIME, qs->iterative_list[2]);
	simple_int_to_cint(&qs->constants.LARGE_PRIME, qs->iterative_list[6]);
	simple_int_to_cint(&qs->constants.M, qs->m.value);

	// Allocates "base length" rows
	qs->base.data = mem;

	// Allocates "s" rows
	qs->s.data = mem_aligned(qs->base.data + qs->base.length);
	mem = mem_aligned(qs->s.data + qs->s.values.defined);
	for (qs_sm i = 0; i < qs->s.values.defined; ++i) {
		simple_inline_cint(&qs->s.data[i].B_terms, kn_size, &mem); // also "s" more cint
		qs->s.data[i].A_inv_2B = mem;
		mem = mem_aligned(qs->s.data[i].A_inv_2B + qs->base.length);
	}

	// Other allocations
	qs->relations.length.reserved = qs->base.length < qs->relations.length.needs ? qs->relations.length.needs : qs->base.length;
	qs->relations.length.reserved <<= 1 ; // 2 * more relations than first guessed are available, hard qs_limit.
	// Lanczos Block has its memory to take a "lite" snapshot before removing relations.
	qs->lanczos.snapshot = mem_aligned(mem) ;
	qs->relations.data = mem_aligned(qs->lanczos.snapshot + qs->relations.length.reserved);
	qs->others.a_divisors = mem_aligned(qs->relations.data + qs->relations.length.reserved);
	qs->others.buffer[0] = mem_aligned(qs->others.a_divisors + qs->s.values.double_value);
	qs->others.buffer[1] = mem_aligned(qs->others.buffer[0] + buffers_size);
	// - the buffer[0] is zeroed after use.
	// - the buffer[1] is left as is after use.
	qs->others.pos[0] = mem_aligned(qs->others.buffer[1] + buffers_size);
	qs->others.pos[1] = mem_aligned(qs->others.pos[0] + qs->base.length);
	qs->others.flags = mem_aligned(qs->others.pos[1] + qs->base.length);
	qs->others.sieve = mem_aligned(qs->others.flags + qs->base.length);
	qs->divisors.data = mem_aligned(qs->others.sieve + qs->m.double_value + 4);
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

static inline void preparation_part_5(qs_sheet *qs) {
	static const double inv_ln_2 = 1.44269504088896340736;
	cint *A = qs->variables.TEMP, *B = A + 1, *C = A + 2;
	qs_sm i = 0, prime;
	// The factor base contains [ multiplier, 2, primes... ]
	if (qs->knuth_schroeppel != 2)
		qs->base.data[i].size = (qs_sm) (.35 + inv_ln_2 * log_computation(qs->base.data[i].num = qs->knuth_schroeppel)), ++i;

	qs->base.data[i].num = 2, qs->base.data[i].size = 1;
	qs->base.data[i].sqrt = *qs->constants.kN.mem % 8 == 1 || *qs->constants.kN.mem % 8 == 7;
	++i;

	for (prime = 3; i < qs->base.length; prime += 2)
		if (is_prime_4669921(prime)) {
			simple_int_to_cint(A, prime);
			cint_div(qs->calc, &qs->constants.kN, A, B, C);
			qs->base.data[i].sqrt = tonelli_shanks((qs_sm) simple_cint_to_int(C), prime);
			if (qs->base.data[i].sqrt) {
				qs->base.data[i].num = prime; // Solution to the congruence exists.
				qs->base.data[i].size = (qs_sm) (.35 + inv_ln_2 * log_computation(prime));
				++i;
			}
		}
}

static inline qs_sm preparation_part_6(qs_sheet *qs, cint *D) {
	cint *A = qs->variables.TEMP, *B = A + 1, *C = A + 2, *E = A + 3;
	qs_sm a, b, min;
	cint_dup(A, &qs->constants.kN);
	cint_left_shifti(A, 1);
	cint_sqrt(qs->calc, A, B, C);
	cint_div(qs->calc, B, &qs->constants.M, D, E);
	qs->d_bits = (qs_sm) cint_count_bits(D);
	cint_nth_root(qs->calc, D, qs->s.values.defined, C);
	assert(cint_count_bits(C) < 16);
	for (a = (qs_sm) simple_cint_to_int(C), min = 1; assert(min < qs->base.length), a >= qs->base.data[min].num; ++min);
	for (b = min * min, min -= qs->the.span_half, assert(b > min); b / min < qs->the.span + min; --min);
	return min;
}

static inline void get_started_iteration(qs_sheet *qs) {
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
	if (qs->relations.length.prev == qs->relations.length.now && qs->curves.now)
		cint_random_bits(&qs->variables.D, qs->d_bits);
	qs->relations.length.prev = qs->relations.length.now;
}

static inline void iteration_part_1(qs_sheet * qs, const cint * D, cint * A) {
	// D is given as a constant, the function must compute A to be close to D.
	// A is a product of "s" random distinct prime numbers from the factor base.
	cint * X = qs->variables.TEMP, * Y = X + 1, *TMP ;
	qs_sm a, b = qs->the.min + qs->the.span_half, idx = 0;
	if (qs->s.values.defined & 1) TMP = A, A = X, X = TMP ;
	// swap pointers so the final multiplication ends up inside A.
	simple_int_to_cint(A, 1);
	for (a = 0, b *= b; a < qs->s.values.subtract_one; ++a) {
		if (a & 1)
			idx = b / (idx + qs->the.min) - (qs_sm) rand_upto(10);
		else idx = qs->the.span_half + (qs_sm) rand_upto(qs->the.span_half) + qs->the.min;
		for (qs_sm i = -1; i != a;)
			for (++idx, i = 0; i < a && qs->s.data[i].a_divisor_idx != idx; ++i);
		qs->s.data[a].a_divisor_idx = idx;
		simple_int_to_cint(Y, qs->base.data[idx].num);
		cint_mul(A, Y, X), TMP = A, A = X, X = TMP ;
	}
	// searches the prime number greater than D/A that isn't in the product.
	cint_div(qs->calc, D, A, X, Y);
	const qs_sm d_over_a = (qs_sm) simple_cint_to_int(X);
	for (idx = qs->base.data[0].num != 2 ; qs->base.data[idx].num <= d_over_a; ++idx);
	for (qs_sm i = 0; i < qs->s.values.subtract_one; ++i)
		if (idx == qs->s.data[i].a_divisor_idx)
			++idx, i = -1 ; // can't use this prime number twice.
	qs->s.data[qs->s.values.subtract_one].a_divisor_idx	= idx ;
	simple_int_to_cint(Y, qs->base.data[idx].num);
	cint_mul(A, Y, X); // A coefficient can be const qualified now.
	TMP = X, X = A, A = TMP ; assert(A == &qs->variables.A);
}

static inline void iteration_part_2(qs_sheet * qs, const cint * A, cint * B) {
	cint *C = qs->variables.TEMP, *D = C + 1, *E = C + 2, *F = C + 3;
	qs_sm a, b, c, *data = qs->others.a_divisors;
	cint_erase(B);
	for (a = 0; a < qs->s.values.defined; ++a) {
		// write [power, index] of the A divisors into buffer.
		*data++ = 1, *data++ = qs->s.data[a].a_divisor_idx;
		b = qs->base.data[qs->s.data[a].a_divisor_idx].num;
		simple_int_to_cint(D, b);
		cint_div(qs->calc, A, D, E, F), assert(F->mem == F->end); // div exact.
		cint_div(qs->calc, E, D, C, F);
		qs->s.data[a].a_mod_prime = (qs_sm) simple_cint_to_int(F);
		c = (qs_sm) modular_inverse(qs->s.data[a].a_mod_prime, b);
		c = (qs_sm) multiplication_modulo(c, qs->base.data[qs->s.data[a].a_divisor_idx].sqrt, b);
		simple_int_to_cint(F, c > b >> 1 ? b > c ? b - c : c - b : c);
		cint_mul(A, F, E);
		cint_div(qs->calc, E, D, &qs->s.data[a].B_terms, F);
		cint_addi(B, &qs->s.data[a].B_terms);
	}
}

static inline void iteration_part_3(qs_sheet * qs, const cint * A, const cint * B) {
	cint *C = qs->variables.TEMP, *D = C + 1, *E = C + 2, *F = C + 3;
	qs_tmp i, j, x, y;
	for (i = 0; i < qs->base.length; ++i) {
		const qs_sm prime = qs->base.data[i].num;
		simple_int_to_cint(C, prime);
		cint_div(qs->calc, A, C, D, E);
		const qs_sm inv_a_prime_mul_2 = (qs_sm) modular_inverse((qs_sm) simple_cint_to_int(E), prime) << 1 ;
		cint_div(qs->calc, B, C, D, E) ;
		x = y = prime;
		x += qs->base.data[i].sqrt;
		y -= qs->base.data[i].sqrt;
		x -= (qs_tmp) simple_cint_to_int(E);
		x *= inv_a_prime_mul_2 >> 1;
		y *= inv_a_prime_mul_2 ;
		x += qs->m.value ;
		x %= prime ;
		y += x ;
		y %= prime ;
		qs->base.data[i].sol[0] = (qs_sm) x ;
		qs->base.data[i].sol[1] = (qs_sm) y ;
		// it was using intermediates signed numbers.
		for (j = 0; j < qs->s.values.defined; ++j) {
			cint_div(qs->calc, &qs->s.data[j].B_terms, C, E, F);
			qs->s.data[j].A_inv_2B[i] = (qs_sm) multiplication_modulo(simple_cint_to_int(F), inv_a_prime_mul_2, prime);
		}
	}
	for (i = 0; i < qs->s.values.defined; cint_left_shifti(&qs->s.data[i++].B_terms, 1));
	// "B terms" were multiplied by 2 for the next function, they are const now.
}

static inline qs_sm iteration_part_4(const qs_sheet * qs, const qs_sm curves, qs_sm ** corr, cint *B) {
	qs_sm i, gray_addi;
	for (i = 0; curves >> i & 1; ++i); // which "B_term" is considered ? the Gray code in "curves" indicates.
	if (gray_addi = (curves >> i & 2) != 0, gray_addi) // which action to perform ? the Gray code in "curves" indicates.
		cint_addi(B, &qs->s.data[i].B_terms) ;
	else
		cint_subi(B, &qs->s.data[i].B_terms) ;
	*corr = qs->s.data[i].A_inv_2B;
	return gray_addi; // All the B generated are distinct, a testing AVL tree must not recognize any duplicate.
}

static inline void iteration_part_5(qs_sheet *  qs, const cint * KN, const cint * B) {
	cint *A = qs->variables.TEMP, *C = A + 1, *D = A + 2, *E = A + 3;
	qs_sm a, b;
	qs_md c;
	qs_tmp s; // temporarily signed.
	for (a = 0; a < qs->s.values.defined; ++a) {
		const qs_sm id = qs->s.data[a].a_divisor_idx, p = qs->base.data[id].num;
		simple_int_to_cint(A, (qs_md) p * (qs_md) p);
		cint_div(qs->calc, KN, A, C, D);
		cint_div(qs->calc, B, A, C, E);
		c = simple_cint_to_int(E);
		b = (qs_sm) multiplication_modulo(c, qs->s.data[a].a_mod_prime, p);
		b = (qs_sm) modular_inverse(E->nat < 0 ? p - b : b, p) ;
		s = (qs_tmp) simple_cint_to_int(D);
		if ((s | c) >> 31){
			// "cint" is slower but the value may be greater than 64-bit
			cint_mul(E, E, C);
			cint_subi(C, D);
			simple_int_to_cint(A, (qs_md) p);
			cint_div(qs->calc, C, A, D, E);
			s = (qs_tmp) simple_cint_to_int(D) ;
			c = D->nat < 0;
		} else {
			// the value can often fit in 64-bit
			s = (qs_tmp) c * (qs_tmp) c - s;
			s /= p;
			c = s < 0;
		}
		s = (qs_tmp) multiplication_modulo(c ? (qs_md) -(s % p) : (qs_md) s, b, p);
		s = (qs_tmp) (c ? s + qs->m.value : -s + qs->m.value) % p;
		s += p * c;
		qs->base.data[id].sol[0] = (qs_sm) s;
		qs->base.data[id].sol[1] = (qs_sm) -1;
	}
}

static inline void iteration_part_6(qs_sheet *qs, const cint *N, const cint *A, const cint *B, cint *C) {
	cint *D = qs->variables.TEMP, *E = D + 1;
	cint_mul(B, B, D); // (B * B) % A = N % A
	cint_subi(D, N);
	cint_div(qs->calc, D, A, C, E);
	assert(E->mem == E->end); // div exact.
}

static inline void iteration_part_7(qs_sheet * qs, const qs_sm gray_addi, const qs_sm * restrict corr) {
	memset(qs->others.sieve, 0, qs->m.double_value);
	memset(qs->others.flags, 0, qs->base.length);
	uint8_t * restrict end = qs->others.sieve + qs->m.double_value, *p_0, *p_1;
	*end = 255;
	for(qs_sm i = qs->iterative_list[3], j = qs->iterative_list[4]; i < j; ++i){
		const qs_sm prime = qs->base.data[i].num, size = qs->base.data[i].size, co = gray_addi ? prime - corr[i] : corr[i] ;
		for(qs->base.data[i].sol[0] += co; qs->base.data[i].sol[0] >= prime; qs->base.data[i].sol[0] -= prime);
		for(qs->base.data[i].sol[1] += co; qs->base.data[i].sol[1] >= prime; qs->base.data[i].sol[1] -= prime);
		p_0 = qs->others.sieve + qs->base.data[i].sol[0] ;
		p_1 = qs->others.sieve + qs->base.data[i].sol[1] ;
		for(; end > p_0 && end > p_1;)
			*p_0 += size, p_0 += prime, *p_1 += size, p_1 += prime ;
		*p_0 += (end > p_0) * size, *p_1 += (end > p_1) * size;
	}
	for(qs_sm i = qs->iterative_list[4], j = qs->iterative_list[5]; i < j; ++i){
		const qs_sm prime = qs->base.data[i].num, size = qs->base.data[i].size, co = gray_addi ? prime - corr[i] : corr[i] ;
		qs->base.data[i].sol[0] += co; for(; qs->base.data[i].sol[0] >= prime; qs->base.data[i].sol[0] -= prime);
		qs->base.data[i].sol[1] += co; for(; qs->base.data[i].sol[1] >= prime; qs->base.data[i].sol[1] -= prime);
		for(p_0 = qs->others.sieve + qs->base.data[i].sol[0]; end > p_0; qs->others.flags[i] |= 1 << ((p_0 - qs->others.sieve) & 7), *p_0 += size, p_0 += prime);
		for(p_1 = qs->others.sieve + qs->base.data[i].sol[1]; end > p_1; qs->others.flags[i] |= 1 << ((p_1 - qs->others.sieve) & 7), *p_1 += size, p_1 += prime);
	}
}

static inline void iteration_part_8(qs_sheet * qs, const qs_sm gray_addi, const qs_sm *  corr) {
	// this operation takes time.
	uint8_t * chunk_open = qs->others.sieve, * chunk_close = chunk_open;
	uint8_t * sieve_close = chunk_open + qs->m.double_value ;
	qs_sm *buffer = qs->others.buffer[0], walk_idx, * walk = buffer;
	// initialize
	for(qs_sm i = 0; i < qs->iterative_list[3]; ++i)
		if (qs->base.data[i].sol[1] != (qs_sm) -1) {
			*walk++ = i ;
			const qs_sm prime = qs->base.data[i].num, co = gray_addi ? prime - corr[i] : corr[i] ;
			qs->base.data[i].sol[0] += co; for(; qs->base.data[i].sol[0] >= prime; qs->base.data[i].sol[0] -= prime);
			qs->base.data[i].sol[1] += co; for(; qs->base.data[i].sol[1] >= prime; qs->base.data[i].sol[1] -= prime);
			qs->others.pos[0][i] = chunk_open + qs->base.data[i].sol[0];
			qs->others.pos[1][i] = chunk_open + qs->base.data[i].sol[1];
		}
	for(walk_idx = 0; buffer[walk_idx] < qs->iterative_list[1]; ++walk_idx);
	// iterates until the next chunk "begin" where the sieve "ends".
	do{
		walk = buffer + walk_idx ;
		chunk_close = sieve_close;
		// initially, the process was advancing by chunks using the following :
		// chunk_close = chunk_close + CACHE_SIZE < sieve_close ? chunk_close + CACHE_SIZE : sieve_close;
		do{
			const qs_sm size = qs->base.data[*walk].size, prime = qs->base.data[*walk].num, times = 4 >> (*walk > qs->iterative_list[2]) ;
			uint8_t ** const p_0 = qs->others.pos[0] + *walk, ** const p_1 = qs->others.pos[1] + *walk;
			const ptrdiff_t diff = *p_1 - *p_0 ;
			for(const uint8_t * const bound = chunk_close - prime * times; bound > *p_0;)
				for(qs_sm i = 0; i < times; ++i)
					**p_0 += size, *(*p_0 + diff) += size, *p_0 += prime;
			for(; *p_0 < chunk_close && *p_0 + diff < chunk_close;)
				**p_0 += size, *(*p_0 + diff) += size, *p_0 += prime;
			*p_1 = *p_0 + diff ;
			if (*p_0 < chunk_close) **p_0 += size, *p_0 += prime;
			if (*p_1 < chunk_close) **p_1 += size, *p_1 += prime;
		} while(*++walk);
	} while(chunk_open = chunk_close, chunk_open < sieve_close);
	memset(qs->others.buffer[0], 0, (walk - qs->others.buffer[0]) * sizeof(walk));
}

static int qs_register_factor(qs_sheet * qs){
	// the function have "permission" to update N.
	// returns -1 if the algorithm should stop, accept any divisor of N.
	cint * F = &qs->variables.FACTOR ;
	int i, res = h_cint_compare(F, &qs->constants.ONE) > 0 && h_cint_compare(F, &qs->variables.N) < 0 ;
	if (res) {
		F->nat = 1 ; // absolute value of the factor.
		const struct avl_node *node = avl_at(&qs->uniqueness[2], F);
		if (qs->uniqueness[2].affected) {
			for (i = 0; i < 2 && qs->n_bits != 1; ++i) {
				const int prime = cint_is_prime(qs->calc, F, -1);
				if (prime) {
					const int power = (int) cint_remove(qs->calc, &qs->variables.N, F);
					assert(power);
					// 200-bit RSA take about 10,000,000+ "duplications", so it will perform the last.
					fac_push(qs->caller, F, 1, power, 0);
					++qs->divisors.total_primes;
					// functions can check if "n_bits" is greater than one
					// if they believe the information has not yet been given above.
					qs->n_bits = (qs_sm) cint_count_bits(&qs->variables.N);
					if (qs->n_bits == 1)
						res = -1;
					else cint_dup(F, &qs->variables.N);
				} else if (i++ == 0)
						qs->divisors.data[qs->divisors.length++] = node->key;
			}
		}
	}
	return res ;
}

static inline void register_relation_kind_2(qs_sheet * qs, const cint * KEY, const cint * VALUE, const qs_sm * args[4]) {
	if (qs->kn_bits < 150)
		return; // this kind of relation is ignored for small inputs.

	// the function searches 2 different KEY sharing the same VALUE.
	struct avl_node *node = avl_at(&qs->uniqueness[1], VALUE);
	struct qs_relation *old, *new;
	cint *BEZOUT = 0;
	qs_sm a, b, *data, *open, *close;
	old = node->value;
	if (old) {
		if (old->id) return; // normally there is no collision.
		if (old->X == 0) return; // modular inverse (p, number) already failed.
		if (old->axis.next) return; // accepting all "next" without caring would reduce the chance to find a nontrivial Lanczos solution.
		for (; old; old = old->axis.next)
			if (h_cint_compare(KEY, old->X) == 0) return; // same KEY already registered.
		old = node->value;
		if (old->axis.next == 0) {
			// this duplicate isn't 3rd, 4th .. it's the 2nd ... so compute BEZOUT
			cint *A = qs->variables.TEMP, *B = A + 1;
			cint_modular_inverse(qs->calc, VALUE, &qs->constants.kN, A);
			if (A->mem == A->end) {
				// no solution to the linear congruence.
				old->X = 0;
				cint_div(qs->calc, &qs->constants.kN, VALUE, A, B);
				a = (qs_sm) simple_cint_to_int(B);
				b = (qs_sm) simple_cint_to_int(VALUE);
				for (; (a %= b) && (b %= a););
				simple_int_to_cint(A, a | b);
				cint_div(qs->calc, &qs->variables.N, A, &qs->variables.FACTOR, B);
				if (B->mem == B->end) // found a factor of N by computing this GCD ("sheer luck")
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
		// BEZOUT is stored directly after X, the newcomer become the root of the linked list.
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

	// data buffered isn't persistent, so it's copied.
	data = new->Y.data = mem_aligned(qs->mem.now);
	new->Y.length = (qs_sm) (args[1] - args[0]);
	memcpy(data, args[0], new->Y.length * sizeof(*data));
	memcpy(data + new->Y.length, args[2], (args[3] - args[2]) * sizeof(*data));
	new->Y.length += (qs_sm) (args[3] - args[2]);
	qs->mem.now = new->Y.data + new->Y.length;

	if (old) {
		int n_new_relations = 0 ;
		BEZOUT = old->X + 1 ; // the modular inverse was stored here.
		cint_mul_mod(qs->calc, new->X, BEZOUT, &qs->constants.kN, &qs->variables.CYCLE);
		do {
			if (old != new) {
				++n_new_relations ;
				// combines, it registers a regular relation using the 2 buffers.
				cint_mul_mod(qs->calc, &qs->variables.CYCLE, old->X, &qs->constants.kN, &qs->variables.KEY);
				open = close = qs->others.buffer[0];
				data = memset(qs->others.buffer[1], 0, qs->base.length * sizeof(*data));
				for (a = 0; a < new->Y.length; a += 2)
					data[new->Y.data[a + 1]] += new->Y.data[a];
				for (a = 0; a < old->Y.length; a += 2)
					data[old->Y.data[a + 1]] += old->Y.data[a];
				for (a = 0; a < qs->base.length; ++a)
					if (data[a])
						*close++ = data[a], *close++ = a;
				args = (const qs_sm *[4]){ open, close, 0, 0, };
				register_relation_kind_1(qs, &qs->variables.KEY, args);
				memset(open, 0, (char *) close - (char *) open); // zeroed.
			}
		} while ((old = old->axis.next));
		// the linked list can handle many entries, but their combinations isn't implemented.
		assert(n_new_relations == 1);
	}
}

static inline void register_relation_kind_1(qs_sheet * qs, const cint * KEY, const qs_sm * args[4]) {
	struct avl_node *node = avl_at(&qs->uniqueness[0], KEY);
	if (node->value)
		return; // duplicates at this stage are ignored.
	struct qs_relation * rel = qs->mem.now;
	qs_sm i, j ;
	// create a new relation (they must be swappable for Lanczos Block reducing).
	qs->mem.now = rel + 1 ;
	rel->X = node->key; // constant X is stored by the node key.
	rel->Y.data = qs->mem.now; // data Y has a known bounded length.
	rel->axis.Z.data = rel->Y.data + (args[1] - args[0]) + (args[3] - args[2]); // writes Z ahead.
	cint_dup(qs->variables.TEMP, &qs->variables.A);
	for (i = 0; i < 4;) {
		// processes the given column arrays.
		const qs_sm * restrict data = args[i++], * restrict data_end = args[i++];
		for (; data < data_end; data += 2) {
			if (*data & 1) {
				// remove index of the prime numbers that are divisors of both "A" and "AX^2 + 2BX + C".
				for (j = 0; j < rel->Y.length && rel->Y.data[j] != *(data + 1); ++j);
				if (j == rel->Y.length)
					// add, the index wasn't present
					rel->Y.data[rel->Y.length++] = *(data + 1);
				else
					// or remove
					memmove(rel->Y.data + j, rel->Y.data + j + 1, (--rel->Y.length - j) * sizeof(*rel->Y.data));
			}
			for (j = 0; j < *data; ++j)
				rel->axis.Z.data[rel->axis.Z.length++] = *(data + 1);
		}
	}
	// Z length is known.
	qs->mem.now = rel->axis.Z.data + rel->axis.Z.length;
	int verified = 0 ;
	if (verified == 0){
		// it often passes but should be verified.
		cint *A = qs->variables.TEMP, *B = A + 1;
		simple_int_to_cint(A, 1);
		for (i = 0; i < rel->axis.Z.length; ++i) {
			simple_int_to_cint(B, qs->base.data[rel->axis.Z.data[i]].num);
			cint_mul_modi(qs->calc, A, B, &qs->constants.kN);
		}
		cint_mul_mod(qs->calc, rel->X, rel->X, &qs->constants.kN, B);
		verified = !cint_compare(A, B) || (cint_addi(A, B), !cint_compare(A, &qs->constants.kN));
	}
	if (verified){
		// Keep this relation
		node->value = qs->relations.data[qs->relations.length.now] = rel;
		qs->mem.now = rel->axis.Z.data + rel->axis.Z.length;
		rel->id = ++qs->relations.length.now;
	} else {
		// Forget all
		char * open = (char*) rel, * close = qs->mem.now ;
		qs->mem.now = memset(open, 0, close - open);
	}
}

static inline void register_relations(qs_sheet * qs, const cint * A, const cint * B, const cint * C) {
	cint *  TMP = qs->variables.TEMP, * K = &qs->variables.KEY, * V = &qs->variables.VALUE ;
	qs_sm a, b, bits, extra, mod, v_1, v_2, verification, *data;
	for (a = 0; a < qs->m.double_value; ++a)
		if (qs->others.sieve[a] >= qs->threshold.value) {
			simple_int_to_cint(&qs->variables.X, a);
			cint_subi(&qs->variables.X, &qs->constants.M); // X = (a - M)
			cint_mul(A, &qs->variables.X, TMP); // TMP = AX
			cint_addi(TMP, B); // TMP = AX + B
			cint_dup(K, TMP); // K = copy of AX + B
			cint_addi(TMP, B); // TMP = AX + 2B
			cint_mul(TMP, &qs->variables.X, V); // V = AX^2 + 2BX
			cint_addi(V, C); // V = AX^2 + 2BX + C
			V->nat = 1 ;
			bits = (qs_sm) cint_count_bits(V) - qs->error_bits;
			extra = 0, verification = 1;
			data = qs->others.buffer[1]; // prepare data into buffer for the next function.
			for (b = 0; b < 2; ++b)
				if (qs->base.data[b].num != 1) {
					simple_int_to_cint(TMP, qs->base.data[b].num);
					*data = (qs_sm) cint_remove(qs->calc, V, TMP);
					if (*data) {
						extra += b ? *data : qs->base.data[b].size;
						*++data = b;
						++data;
					}
				}
			for (b = 2; b < qs->iterative_list[1] && verification; ++b) {
				v_1 = (v_2 = 0, qs->base.data[b].sol[1] == (qs_sm) -1) || (v_2 = 1, mod = a % qs->base.data[b].num, mod == qs->base.data[b].sol[0] || mod == qs->base.data[b].sol[1]);
				if (v_1) {
					simple_int_to_cint(TMP, qs->base.data[b].num);
					// remove is the action of updating a number so that it isn't a multiple of another number anymore.
					*data = (qs_sm) cint_remove(qs->calc, V, TMP);
					verification &= v_2 <= *data;
					if (*data) {
						*++data = b, ++data;
						extra += qs->base.data[b].size;
					}
				}
			}
			if (qs->others.sieve[a] += extra, qs->others.sieve[a] >= bits) {
				for (; b < qs->iterative_list[4] && extra < qs->others.sieve[a] && verification; ++b) {
					v_1 = (v_2 = 0, qs->base.data[b].sol[1] == (qs_sm) -1) || (v_2 = 1, mod = a % qs->base.data[b].num, mod == qs->base.data[b].sol[0] || mod == qs->base.data[b].sol[1]);
					if (v_1) {
						simple_int_to_cint(TMP, qs->base.data[b].num);
						*data = (qs_sm) cint_remove(qs->calc, V, TMP);
						verification &= v_2 <= *data;
						if (*data || v_2) {
							*++data = b, ++data;
							extra += qs->base.data[b].size;
						}
					}
				}
				const uint8_t mask = 1 << (a & 7);
				for (b = qs->iterative_list[4]; b < qs->base.length && extra < qs->others.sieve[a] && verification; ++b)
					if (qs->others.flags[b] & mask)
						if (mod = a % qs->base.data[b].num, mod == qs->base.data[b].sol[0] || mod == qs->base.data[b].sol[1]) {
							simple_int_to_cint(TMP, qs->base.data[b].num);
							*data = (qs_sm) cint_remove(qs->calc, V, TMP);
							verification &= *data != 0;
							++data, *data++ = b;
							extra += qs->base.data[b].size;
						}
				if (verification) {
					const qs_sm * args[4] = {
							qs->others.a_divisors,
							qs->others.a_divisors + qs->s.values.double_value,
							qs->others.buffer[1],
							data,
					};
					if (h_cint_compare(V, &qs->constants.SMALL_PRIME) < 0) {
						register_relation_kind_1(qs, K, args);
					} else if (h_cint_compare(V, &qs->constants.LARGE_PRIME) < 0)
						register_relation_kind_2(qs, K, V, args);
				}
				// If verification is false : "An abnormally unusable data row has been ignored".
			}
		}
}

static inline void finalization_part_1(qs_sheet * qs, const uint64_t * lanczos_answer) {
	const uint64_t mask = *lanczos_answer, * null_rows = lanczos_answer + 1;
	// Lanczos "linear algebra" answer is simply "mask followed by null_rows", with read-only.
	if (mask == 0 || null_rows == 0)
		return;
	cint *A = qs->variables.TEMP, *B = A + 1, *C = A + 2, *D = A + 3;
	qs_md a, b, c;
	qs_sm *power_of_primes;
	for(c = 0; c < 64; ++c)
		if (mask >> c & 1){
			simple_int_to_cint(A, 1), simple_int_to_cint(B, 1), simple_int_to_cint(C, 1);
			power_of_primes = memset(qs->others.buffer[1], 0, qs->base.length * sizeof(*power_of_primes));
			for (a = 0; a < qs->relations.length.now; ++a)
				if (null_rows[a] >> c & 1) {
					const struct qs_relation * restrict const rel = qs->relations.data[a];
					cint_mul_modi(qs->calc, A, rel->X, &qs->variables.N);
					for (b = 0; b < rel->axis.Z.length; ++b)
						++power_of_primes[rel->axis.Z.data[b]];
				}
			for (b = 0; b < qs->base.length; ++b)
				if (power_of_primes[b]){
					simple_int_to_cint(B, qs->base.data[b].num);
					simple_int_to_cint(D, power_of_primes[b] >> 1); // power is even, take square root
					cint_pow_modi(qs->calc, B, D, &qs->variables.N);
					cint_mul_modi(qs->calc, C, B, &qs->variables.N);
				}
			h_cint_subi(C, A);
			if (C->mem != C->end) {
				C->nat = 1; // absolute value
				cint_gcd(qs->calc, &qs->variables.N, C, &qs->variables.FACTOR); // gcd
				if(qs_register_factor(qs) == -1)
					break;
			}
		}
}

static inline void finalization_part_2(qs_sheet * qs) {
	if (qs->n_bits == 1)
		return;

	// Algorithm has finalized but N is still greater than one.
	// Perform checks until no new divisor can be discovered.

	cint * F = &qs->variables.FACTOR, **di = qs->divisors.data, *Q = qs->variables.TEMP, *R = Q + 1 ;
	qs_sm i, j, k, count = 0 ;

	for(i = qs->divisors.processing_index, k = qs->divisors.length; i < k; ++i)
		for (j = 1 + i; j < k; ++j)
			if (cint_gcd(qs->calc, di[i], di[j], F), qs_register_factor(qs) == -1)
				return; // register "gcd of new divisors with old divisors"
	do {
		for (; i = k, k = qs->divisors.length, i != k;)
			for (; i < k; ++i)
				for (j = 0; j < i; ++j)
					if (cint_gcd(qs->calc, di[i], di[j], F), qs_register_factor(qs) == -1)
						return; // register "gcd of new divisors with old divisors"

		for (i = qs->divisors.processing_index; i < k; ++i)
			if (fac_any_root_check(qs->caller, di[i], &qs->variables.FACTOR, R))
				if (qs_register_factor(qs) == -1)
					return; // register "perfect power of new divisors"

		if (count != qs->divisors.total_primes)
			if (fac_any_root_check(qs->caller, &qs->variables.N, &qs->variables.FACTOR, R))
				if (qs_register_factor(qs) == -1)
					return; // register "prefect root of N"

		count = qs->divisors.total_primes ;
		qs->divisors.processing_index = k ;
	} while(k != qs->divisors.length);

}

static inline int finalization_part_3(qs_sheet * qs) {
	// Usually do nothing, (N = 1)
	// Otherwise push something non-trivial to the caller's routine
	int res = qs->n_bits == 1 ;
	if (res == 0){
		cint * F = &qs->variables.FACTOR;
		for(qs_sm i = 0; i < qs->divisors.length; ++i) {
			if (qs->caller->params->silent == 0) {
				char *str = cint_to_string(qs->divisors.data[i], 10);
				printf("- quadratic sieve can't silently ignore %s\n", str);
				free(str);
			}
			if (h_cint_compare(&qs->variables.N, qs->divisors.data[i]) > 0) {
				const int power = (int) cint_remove(qs->calc, &qs->variables.N, qs->divisors.data[i]);
				if (power) {
					assert(power == 1);
					fac_push(qs->caller, F, 0, 1, 1);
				}
			}
		}
		res = h_cint_compare(&qs->variables.N, &qs->caller->number->cint) ;
		if (res) // res is true if QS was able to decompose N.
			if (h_cint_compare(&qs->variables.N, &qs->constants.ONE))
				fac_push(qs->caller, &qs->variables.N, -1, 1, 1);
	}
	return res;
}

// "Pure C99 factorizer using self-initialising Quadratic Sieve."
// compiling with "gcc -Wall -pedantic -O3 main.c -o qs" can speed up by a factor 2 or 3.
