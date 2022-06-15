//  GNU General Public License

//      As an undergraduate student, this project is part of my computer science + maths training.
//          - This software proposition is from Michel Leonard (student at Université de Franche-Comté, Mon, 11 Jul 2022)
//          - There is no guarantee of any kind on the software
//          - C code is shared under the terms of the GNU General Public License
//          - The main mathematical and logical inspiration source is located at :
//              http://web.mit.edu/sage/export/flintqs-0.0.20070817/QS.cpp - GNU General Public License

//  This software implementation would have been impossible without "FLINT: Fast Library for Number Theory" maintained by William Hart.

// Quadratic sieve manager
static int quadratic_sieve(fac_caller *caller) {
	// The routine is invoked by a caller, and :
	// - must factor the caller's "number"
	// - must register answers to the caller's routines
	// - can use caller's resources (vars + configurations)
	// - resources are enough to copy 300-bit into variables
	int limit = caller->params->limit;
	if (limit == 0) limit = 220;
	if (caller->number->bits < 65 || caller->number->bits > limit) return 0;

	qs_sheet qs = {0};
	preparation_part_1(&qs, caller);
	preparation_part_2(&qs, qs.vars.N, qs.constants.A, qs.constants.kN);
	preparation_part_3(&qs, qs.constants.kN, qs.constants.M);
	// adjustor and qs_multiplier of number are defined, parametrize.
	qs_parametrize(&qs);
	preparation_part_4(&qs);
	preparation_part_5(&qs);
	qs.info.the_min = preparation_part_6(&qs, qs.vars.D);
	do {
		do {
			iteration_analyzer(&qs);
			iteration_part_1(&qs, qs.vars.A);
			cint_dup(qs.vars.A, iteration_part_2(&qs, qs.vars.A, qs.vars.D));
			iteration_part_3(&qs, qs.vars.A, qs.vars.B);
			iteration_part_4(&qs, qs.vars.A, qs.vars.B);
			for (qs_sm i = 1, add, *corr; i < qs.info.poly_max; ++i, ++qs.analyzer.curves) {
				add = iteration_part_5(&qs, i, &corr, qs.vars.B);
				iteration_part_6(&qs, qs.constants.kN, qs.vars.B);
				iteration_part_7(&qs, qs.constants.kN, qs.vars.A, qs.vars.B, qs.vars.C);
				iteration_part_8(&qs, add, corr);
				iteration_part_9(&qs, add, corr);
				register_relations(&qs, qs.vars.A, qs.vars.B, qs.vars.C);
			}
		} while (inner_continuation_condition(&qs));
		// Lanczos may read-only (small N) or write to the matrix (larger N).
		finalization_part_1(&qs, lanczos_block(&qs));
		finalization_part_2(&qs);
	} while (outer_continuation_condition(&qs));
	const int res = finalization_part_3(&qs);
	free(qs.mem.base);
	return res;
}

// Quadratic sieve main condition 1
static inline int inner_continuation_condition(qs_sheet *qs) {
	// Used to decide if the inner loop (sieving loop) should continue to search relations or break.
	int res = 1;
	res &= qs->relations.length.now < qs->relations.length.expected; // the condition.
	if (qs->caller->params->silent == 0) {
		const int rel_begin = (int) qs->relations.length.now, rel_end = (int) qs->relations.length.expected;
		fac_display_progress("Quadratic sieve", rel_begin * 100 / rel_end);
	}
	return res;
}

// Quadratic sieve main condition 2
static inline int outer_continuation_condition(qs_sheet *qs) {
	// The algorithm has normally finished its work, the prime factors if they exist have been identified.
	// Used to decide whether the algorithm returns to the inner loop (sieving) or returns control.
	int res = qs->analyzer.retry_perms-- > 0; // avoid infinite loop.
	res &= qs->divisors.total_primes < qs->analyzer.retry_perms; // search at least (2.. 1.. 0) prime factors.
	if (res) {
		// puts("quadratic sieve need more relations");
		// the new parameter is to collect a little more relations.
		qs->relations.length.expected = qs->relations.length.now + (qs->relations.length.now >> (1 + qs->analyzer.retry_perms));
	}
	return res;
}

// Quadratic sieve parameters
static inline qs_sm linear_param_resolution(const double v[][2], const qs_sm bits) {
	qs_sm res, i, j = 1 ;
	for(; v[j + 1][0] && bits > v[j][0]; ++j);
	i = j - 1 ;
	if (v[i][0] > bits) res = (qs_sm) v[i][1] ;
	else if (v[j][0] < bits) res = (qs_sm) v[j][1] ;
	else {
		const double a = (v[j][1] - v[i][1]) / (v[j][0] - v[i][0]);
		const double b = v[i][1] - a * v[i][0];
		res = (qs_sm) (a * bits + b);
	}
	return res + (res > 512) * (16 - res % 16) ;
}

static inline void qs_parametrize(qs_sheet *qs) {
	const qs_sm bits = (qs_sm) cint_count_bits(qs->constants.kN);
	qs->info.kn_bits = bits; // input was adjusted so there is at least 115-bit.

	// params are { bits, value }
	// static const double param_base_size   [][2]= { {110, 800}, {130, 1500}, {200, 3200}, {260, 15000}, {290, 30000}, {0} };
	static const double param_base_size   [][2]= { {130, 847}, {150, 1131}, {170, 1970}, {190, 3336}, {210, 5520}, {230, 8947}, {250, 14244}, {270, 22320}, {290, 34479}, {0} };
	// static const double param_base_size   [][2]= { {120, 800}, {140, 1200}, {160, 1800}, {180, 2500}, {200, 4000}, {220, 6000}, {240, 10000}, {260, 25000}, {290, 35000}, {0} };
	qs->base.length = qs->relations.length.expected = linear_param_resolution(param_base_size, bits); // given by (int) exp(.32 * sqrt(log($n) * log(log($n))))

	static const double param_m_value     [][2]= { {110, 64e3}, {200, 256e3}, {0} };
	qs->info.m.value = linear_param_resolution(param_m_value, bits);

	static const double param_m_alloc     [][2]= { {110, 1e7}, {270, 1e8}, {290, 2e8}, {0} };
	qs->info.total_bytes_allocated = linear_param_resolution(param_m_alloc, bits);

	static const double param_error      [][2]= { {110, 14}, {290, 32}, {0} };
	qs->info.error_bits = linear_param_resolution(param_error, bits);

	//static const double param_threshold   [][2]= { {110, 63}, {220, 77}, {300, 102}, {0} };
	static const double param_threshold   [][2]= { {120, 64}, {140, 66}, {160, 70}, {180, 73}, {200, 80}, {220, 87}, {240, 93}, {260, 100}, {290, 114}, {0} };
	qs->info.threshold = linear_param_resolution(param_threshold, bits);


	// Other parameters
	qs->analyzer.retry_perms = 3; // Sieve again 3 times before giving up.
	qs->s.values.double_value = (qs->s.values.defined = (qs->s.values.subtract_one = bits / 28) + 1) << 1;
	qs->info.poly_max = (1 << qs->s.values.subtract_one) - 1;
	qs->info.cache_block_size = 32000;

	// Computations
	qs->info.p_list[0] = 1; // one

	static const double param_first_prime [][2]= { {110, 8}, {210, 11}, {260, 25}, {290, 28}, {0} };
	qs->info.p_list[1] = linear_param_resolution(param_first_prime, bits); // first

	qs->info.p_list[2] = 768; // medium
	qs->info.p_list[3] = qs->base.length < 2048 ? qs->base.length : 2048; // mid
	qs->info.p_list[4] = qs->base.length < 5120 ? qs->base.length : 5120; // sec
	qs->info.p_list[5] = qs->base.length; // factor base size

	static const double param_large_prime [][2]= { {110, 25e4}, {170, 1e6}, {210, 1e7}, {240, 3e7}, {290, 3e8}, {0} };
	qs->info.p_list[6] =  linear_param_resolution(param_large_prime, bits); // large

	qs->info.the_span = qs->base.length / qs->s.values.defined / qs->s.values.defined / 2;
	assert(qs->info.cache_block_size <= qs->info.m.value);
	{
		qs->info.m.double_value = qs->info.m.value << 1;
		const qs_sm q = qs->info.m.double_value / qs->info.cache_block_size;
		const qs_sm r = qs->info.m.double_value % qs->info.cache_block_size;
		qs->info.m.q = q;
		qs->info.m.r = r;
		qs->info.m.n_reps = 1 + (q > 0) * ((q > 1) * (q - 2) + 1 + (r != 0));
		qs->info.m.divided = qs->info.m.double_value / sizeof(uint64_t);
	}
}

// Quadratic sieve source
static inline void preparation_part_1(qs_sheet *qs, fac_caller *caller) {
	// Algorithm is initialized with the caller's resources
	qs->caller = caller;
	qs->calc = caller->calc;
	qs->constants.A = caller->vars; // adjustor
	qs->constants.M = caller->vars + 1; // qs_multiplier
	qs->constants.kN = caller->vars + 2;
	qs->vars.N = caller->vars + 3;
	qs->vars.temp = caller->vars + 4;
	cint_dup(qs->vars.N, &caller->number->cint);
	qs->info.n_bits = caller->number->bits;
}

static inline void preparation_part_2(qs_sheet *qs, const cint * N, cint * ADJUSTOR, cint * kN) {
	// Input is transparently adjusted by a prime to measure at least 115-bit.
	static const int prime_generator[] = {
			9, 7, 5, 3, 17, 27, 3, 1, 29, 3, 21, 7, 17, 15,
			9, 43, 35, 15, 29, 3, 11, 3, 11, 15, 17, 25, 53,
			31, 9, 7, 23, 15, 27, 15, 29, 7, 59, 15, 5, 21,
			69, 55, 21, 21, 5, 159, 3, 81, 9, 69, 131, 33,};
	const qs_sm bits = qs->info.n_bits;
	if (bits < 115) {
		const qs_md adjustor = bits < 115 ? (1LLU << (124 - bits)) + prime_generator[115 - bits] : 1;
		simple_int_to_cint(ADJUSTOR, adjustor);
		cint_mul(N, ADJUSTOR, kN);
	} else {
		cint_reinit(ADJUSTOR, 1);
		cint_dup(kN, N);
	}
}

static inline int preparation_part_3(qs_sheet *qs, cint * kN, cint * MULTIPLIER) {
	cint *A = qs->vars.temp, *B = A + 1, *C = A + 2, *D = A + 3;
	static const int mul[] = {1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,}, n_mul = sizeof(mul) / sizeof(*mul);
	int a, b, c, d, res;
	if (qs->caller->params->qs_multiplier)
		res = qs->caller->params->qs_multiplier ;
	else {
		double factors[n_mul], logarithm;
		b = (int) (*kN->mem % 8);
		for (a = 0; a < n_mul; ++a) {
			int numerator = b * mul[a] % 8;
			numerator = 1 << (numerator == 1 ? 4 : numerator == 5 ? 2 : 1);
			factors[a] = log_computation((double) numerator / (double) mul[a]) / 2;
		}
		b = 10 + (int) (cint_count_bits(kN) << 1); // Function of the bit count of number.
		for (a = 3; b; a += 2)
			if (is_prime_4669921(a)) {
				logarithm = log_computation(a) / a;
				simple_int_to_cint(B, a);
				cint_div(qs->calc, kN, B, C, D);
				d = kronecker_symbol(simple_cint_to_int(D), a);
				for (c = 0; c < n_mul; ++c)
					factors[c] += logarithm + logarithm * d * kronecker_symbol(mul[c], a);
				--b;
			}
		for (a = 1, c = 0; a < n_mul; ++a)
			if (factors[a] > factors[0])
				factors[0] = factors[c = a];
		res = mul[c];
	}
	simple_int_to_cint(MULTIPLIER, res);
	if (res > 1)
		cint_dup(B, kN), cint_mul(B, MULTIPLIER, kN);
	return res;
}

static inline void preparation_part_4(qs_sheet *qs) {
	void *mem;
	mem = qs->mem.base = calloc(1, qs->info.total_bytes_allocated);
	assert(mem);
	qs->info.s_rand = mix_rand_seed(&mem);

	// cint are 16 integers ... contiguously allocated, not larger than kN, not resized during execution
	cint *num = qs->vars.temp = mem_aligned(mem);
	mem = num + 16;
	for (int i = 0; i < 16; ++i)
		simple_inline_cint(&num[i], qs->constants.kN->size, &mem);

	// 5 are left to compute anything

	// 5 are named vars
	qs->vars.A = num + 5; qs->vars.B = num + 6;
	qs->vars.C = num + 7, qs->vars.D = num + 8;
	qs->vars.FACTOR = num + 9 ;

	// the caller's memory was used to initialize the algorithm
	// now N and kN are copied to the algorithm memory block
	cint_dup(num + 10, qs->vars.N), qs->vars.N = num + 10;
	cint_dup(num + 11, qs->constants.kN), qs->constants.kN = num + 11;

	// 4 are named constants
	qs->constants.ONE = num + 12, qs->constants.LOWER = num + 13 ;
	qs->constants.M_2 = num + 14, qs->constants.UPPER = num + 15 ;

	simple_int_to_cint(qs->constants.ONE, 1);
	simple_int_to_cint(qs->constants.LOWER, 1000);
	simple_int_to_cint(qs->constants.UPPER, qs->info.p_list[6]);
	simple_int_to_cint(qs->constants.M_2, qs->info.m.value);

	// Allocates "base length" rows
	qs->base.data = mem;
	mem = qs->base.data + qs->base.length;

	// Allocates "s" rows
	qs->s.data = mem;
	mem = qs->s.data + qs->s.values.defined;
	for (size_t i = 0; i < qs->s.values.defined; ++i) {
		simple_inline_cint(&qs->s.data[i].B_terms, qs->constants.kN->size, &mem); // also "s" more cint
		qs->s.data[i].A_inv_2B = mem;
		mem = qs->s.data[i].A_inv_2B + qs->base.length;
	}

	// Other allocations
	const size_t relations_reserved = qs->relations.length.expected << 2 ;
	qs->relations.data = mem_aligned(mem); // 4 * more relations than first guessed are available in array, hard limit.
	qs->others.sm_buffer = mem_aligned(qs->relations.data + relations_reserved); // Small buffer, sized for A-invariants
	const size_t medium_buffer_size = qs->base.length + (qs->info.p_list[1] << 1);
	qs->others.md_uncleared_buffer = mem_aligned(qs->others.sm_buffer + qs->s.values.double_value); // Medium buffer, not cleared after usage
	qs->others.md_cleared_buffer = mem_aligned(qs->others.md_uncleared_buffer + medium_buffer_size); // Medium buffer, cleared after usage
	qs->others.offsets[0] = mem_aligned(qs->others.md_cleared_buffer + medium_buffer_size);
	qs->others.offsets[1] = mem_aligned(qs->others.offsets[0] + qs->base.length);
	qs->others.flags = mem_aligned(qs->others.offsets[1] + qs->base.length);
	qs->others.sieve = mem_aligned(qs->others.flags + qs->base.length);
	qs->divisors.data = mem_aligned(qs->others.sieve + qs->info.m.double_value + 4);
	qs->mem.now = mem_aligned(qs->divisors.data + 512);

	struct avl_manager *trees[2] = {&qs->relations.tree, &qs->divisors.tree,};
	for (int i = 0; i < 2; ++i) {
		trees[i]->inserter_argument = &qs->mem.now;
		trees[i]->inserter = &avl_cint_inserter;
		trees[i]->comparator = (int (*)(const void *, const void *)) &h_cint_compare; // use default sign-less comparator.
	}

}

static inline void preparation_part_5(qs_sheet *qs) {
	static const double inv_ln_2 = 1.44269504088896340736;
	cint *A = &qs->vars.temp[0], *B = A + 1, *C = A + 2;
	qs_md a = 0, multiplier = simple_cint_to_int(qs->constants.M);

	if (multiplier != 2)
		qs->base.data[a].size = (qs_sm) (.35 + inv_ln_2 * log_computation(qs->base.data[a].num = multiplier)), ++a;

	qs->base.data[a].num = 2, qs->base.data[a].size = 1;
	qs->base.data[a].sqrt = *qs->constants.kN->mem % 8 == 1 || *qs->constants.kN->mem % 8 == 7;
	++a;

	for (qs_sm i = 3; a < qs->base.length; i += 2)
		if (is_prime_4669921(i)) {
			cint_div(qs->calc, qs->constants.kN, cint_immediate(qs->calc, i), B, C);
			qs->base.data[a].sqrt = tonelli_shanks(simple_cint_to_int(C), i);
			if (qs->base.data[a].sqrt) {
				qs->base.data[a].num = i; // Solution to the congruence exists.
				qs->base.data[a].size = (qs_sm) (.35 + inv_ln_2 * log_computation(i));
				++a;
			}
		}
}

static inline qs_sm preparation_part_6(qs_sheet *qs, cint *res) {
	cint *A = &qs->vars.temp[0], *B = A + 1, *C = A + 2, *D = A + 3;
	qs_sm a, b, min;
	cint_dup(A, qs->constants.kN);
	cint_left_shifti(A, 1);
	cint_sqrt(qs->calc, A, B, C);
	cint_div(qs->calc, B, qs->constants.M_2, res, D);
	cint_nth_root(qs->calc, res, qs->s.values.defined, C);
	assert(cint_count_bits(C) < 16);
	for (a = (qs_sm) simple_cint_to_int(C), min = 1; assert(min < qs->base.length), a >= qs->base.data[min].num; ++min);
	for (b = min * min, min -= qs->info.the_span >> 1, assert(b > min); b / min < qs->info.the_span + min; --min);
	return min;
}

static inline void iteration_analyzer(qs_sheet *qs) {
	qs->analyzer.blank += qs->relations.length.now == qs->analyzer.progress;
	if (qs->analyzer.blank == qs->s.values.defined) {
		qs->analyzer.blank = 0;
		// multiple iterations was performed without new relations, so unblock sieving.
		// a rare case I have encountered after 230-bit, it may be solved by 128-bit support...
		cint_random_bits(qs->vars.D, cint_count_bits(qs->vars.D));
	}
	qs->analyzer.progress = qs->relations.length.now;
}

static inline void iteration_part_1(qs_sheet *qs, cint *A) {
	cint *B = &qs->vars.temp[0], *C = B + 1, *D;
	qs_sm a, b = 0, c, d = qs->info.the_span >> 1, f = qs->info.the_min + d;
	if (qs->s.values.defined & 1) D = B, B = A, A = D;
	cint_reinit(B, 1);
	for (f *= f, a = 0; a < qs->s.values.defined; D = B, B = A, A = D, ++a) {
		if (a & 1)
			b = f / (b + qs->info.the_min) - (qs_sm) rand_upto(10) - qs->info.the_min;
		else b = d + (qs_sm) rand_upto(d);
		for (c = -1; c != a;)
			for (++b, c = 0; c < a && qs->s.data[c].a_ind != b; ++c);
		qs->s.data[a].a_ind = b;
		simple_int_to_cint(C, qs->base.data[b + qs->info.the_min].num);
		cint_mul(B, C, A);
	}
	for (a = 0; a < qs->s.values.defined; qs->s.data[a++].a_ind += qs->info.the_min);
}

static inline cint *iteration_part_2(qs_sheet *qs, const cint *A, const cint *D) {
	cint *B = &qs->vars.temp[0], *C = B + 1;
	qs_sm a, b, c = qs->s.values.subtract_one;
	cint_div(qs->calc, D, A, B, C);
	b = simple_cint_to_int(B);
	for (a = 1; b >= qs->base.data[a].num; ++a);
	for (b = -1; b != c; ++a)
		for (b = 0; b < c && qs->s.data[b].a_ind != a; ++b);
	qs->s.data[c].a_ind = --a;
	simple_int_to_cint(B, qs->base.data[a].num);
	cint_mul(A, B, C);
	return C;
}

static inline void iteration_part_3(qs_sheet *qs, const cint *A, cint *B) {
	cint *C = &qs->vars.temp[0], *D = C + 1, *E = C + 2, *F = C + 3;
	qs_sm a, b, c, *buffer = qs->others.sm_buffer;
	cint_reinit(B, 0);
	for (a = 0; a < qs->s.values.defined; ++a) {
		*buffer++ = 1, *buffer++ = qs->s.data[a].a_ind; // write later required A-invariants into sm_buffer.
		b = qs->base.data[qs->s.data[a].a_ind].num;
		simple_int_to_cint(D, b);
		cint_div(qs->calc, A, D, E, F);
		cint_div(qs->calc, E, D, C, F);
		qs->s.data[a].a_mod_p = simple_cint_to_int(F);
		c = modular_inverse(qs->s.data[a].a_mod_p, b);
		c = multiplication_modulo(c, qs->base.data[qs->s.data[a].a_ind].sqrt, b);
		simple_int_to_cint(F, c > b >> 1 ? b > c ? b - c : c - b : c);
		cint_mul(A, F, E);
		cint_div(qs->calc, E, D, &qs->s.data[a].B_terms, F);
		cint_addi(B, &qs->s.data[a].B_terms);
	}
}

static inline void iteration_part_4(qs_sheet *qs, const cint *A, const cint *B) {
	cint *C = &qs->vars.temp[0], *D = C + 1, *E = C + 2, *F = C + 3;
	qs_sm a, b;
	qs_md_tmp_si s; // temporarily signed.
	for (a = 0; a < qs->base.length; ++a) {
		const qs_sm c = qs->base.data[a].num;
		simple_int_to_cint(C, c);
		cint_div(qs->calc, A, C, D, E);
		qs->base.data[a].A_inv = modular_inverse(simple_cint_to_int(E), c);
		for (b = 0; b < qs->s.values.defined; ++b) {
			cint_div(qs->calc, &qs->s.data[b].B_terms, C, E, F);
			qs->s.data[b].A_inv_2B[a] = multiplication_modulo(simple_cint_to_int(F), qs->base.data[a].A_inv << 1, c);
		}
		cint_div(qs->calc, B, C, D, E);
		s = c;
		s += qs->base.data[a].sqrt;
		s -= (qs_md_tmp_si) simple_cint_to_int(E);
		s *= qs->base.data[a].A_inv;
		s += (qs_md_tmp_si) qs->info.m.value;
		qs->base.data[a].sol[0] = (qs_sm) (s % c);
		assert(qs->base.data[a].sol[0] < c);
		s = qs->base.data[a].sqrt;
		s -= c;
		s *= qs->base.data[a].A_inv << 1;
		s = -s % c;
		qs->base.data[a].sol[1] = (qs_sm) (s + qs->base.data[a].sol[0]);
		assert(qs->base.data[a].sol[1] < c + c);
	}
}

static inline qs_sm iteration_part_5(const qs_sheet *qs, const qs_sm curves, qs_sm **corr, cint *B) {
	qs_sm a, b, act;
	for (a = 0; a < qs->s.values.defined && !(curves >> a & 1); ++a);
	void (*const fn)(cint *, const cint *) = (act = curves >> a & 2) ? &cint_addi : &cint_subi;
	for (b = 0; b < 2; ++b)
		fn(B, &qs->s.data[a].B_terms);
	*corr = qs->s.data[a].A_inv_2B;
	return act;
}

static inline void iteration_part_6(qs_sheet *qs, const cint *KN, const cint *B) {
	cint *A = &qs->vars.temp[0], *C = A + 1, *D = A + 2, *E = A + 3;
	qs_sm a, b;
	qs_md c;
	qs_md_tmp_si s; // temporarily signed.
	for (a = 0; a < qs->s.values.defined; ++a) {
		const qs_sm id = qs->s.data[a].a_ind, p = qs->base.data[id].num;
		simple_int_to_cint(A, (qs_md) p * (qs_md) p);
		cint_div(qs->calc, KN, A, C, D);
		cint_div(qs->calc, B, A, C, E);
		s = (qs_md_tmp_si) simple_cint_to_int(D);
		c = simple_cint_to_int(E);
		b = multiplication_modulo(c, qs->s.data[a].a_mod_p, p);
		b = modular_inverse(E->nat > 0 ? b : p - b, p) % p;
		s = (qs_md_tmp_si) (c * c - s) / p;
		c = s > 0;
		s = (qs_md_tmp_si) multiplication_modulo(c ? s : -s, b, p);
		s = (qs_md_tmp_si) (c ? -s + qs->info.m.value : s + qs->info.m.value);
		s %= p;
		s += p * (s < 0);
		qs->base.data[id].sol[0] = (qs_sm) s;
		qs->base.data[id].sol[1] = -1U;
	}
}

static inline void iteration_part_7(qs_sheet *qs, const cint *N, const cint *A, const cint *B, cint *C) {
	cint *D = &qs->vars.temp[0], *E = D + 1;
	cint_mul(B, B, D);
	cint_subi(D, N);
	cint_div(qs->calc, D, A, C, E);
	assert(E->mem == E->end); // div exact.
}

static inline void iteration_part_8(qs_sheet *qs, const qs_sm add, const qs_sm *corr) {
	qs_sm a, b, c, d;
	memset(qs->others.sieve, 0, qs->info.m.double_value);
	memset(qs->others.flags, 0, qs->base.length);
	uint8_t *pos[2], *end = qs->others.sieve + qs->info.m.double_value;
	*end = 255;
	for (a = 3, b = 1 + a; a < 5; ++a, ++b)
		for (c = qs->info.p_list[a]; c < qs->info.p_list[b]; ++c) {
			for (d = 0; d < 2; pos[d] = qs->others.sieve + qs->base.data[c].sol[d], ++d) {
				qs->base.data[c].sol[d] += add ? -corr[c] + qs->base.data[c].num : corr[c];
				for (; qs->base.data[c].sol[d] >= qs->base.data[c].num; qs->base.data[c].sol[d] -= qs->base.data[c].num);
			}
			for (; end > pos[0] && end > pos[1];)
				for (d = 0; d < 2; ++d)
					*pos[d] += qs->base.data[c].size, pos[d] += qs->base.data[c].num;
			if (a == 2) {
				for (d = 0; d < 2; ++d)
					if (end > pos[d])
						*pos[d] += qs->base.data[c].size;
			} else
				for (d = 0; d < 2; ++d)
					for (; end > pos[d];)
						qs->others.flags[c] |= 1 << ((pos[d] - qs->others.sieve) & 7), *pos[d] += qs->base.data[c].size, pos[d] += qs->base.data[c].num;
		}
}

static inline void iteration_part_9(qs_sheet *qs, const qs_sm add, const qs_sm *corr) {
	qs_sm a, b, c, d, e;
	ptrdiff_t diff;
	uint8_t *pos[2], *curr, *end, *tmp;
	for (a = 0; a < qs->info.m.n_reps; ++a) {
		const int cond_1 = a == 0;
		const int cond_2 = qs->info.m.n_reps == 1 || a != qs->info.m.n_reps - 1;
		curr = qs->others.sieve + qs->info.cache_block_size * a;
		end = curr + (cond_2 ? qs->info.cache_block_size : qs->info.m.r);
		for (b = 0, c = 1 + b; b < 3; ++b, ++c) {
			const qs_sm times = 1 << (3 - b);
			if (b || cond_1) {
				for (d = qs->info.p_list[b]; d < qs->info.p_list[c]; ++d)
					if (b == 2 || qs->base.data[d].sol[1] != -1U) {
						if (cond_1) {
							for (e = 0; e < 2; ++e) {
								qs->base.data[d].sol[e] += add ? -corr[d] + qs->base.data[d].num : corr[d];
								for (; qs->base.data[d].sol[e] >= qs->base.data[d].num; qs->base.data[d].sol[e] -= qs->base.data[d].num);
							}
						}
						if (b) {
							for (e = 0; e < 2; ++e)
								pos[e] = cond_1 ? curr + qs->base.data[d].sol[e] : qs->others.offsets[e][d];
							diff = pos[1] - pos[0];
							tmp = end - qs->base.data[d].num * times;
							for (; tmp > pos[0];)
								for (e = 0; e < times; ++e)
									*pos[0] += qs->base.data[d].size, *(pos[0] + diff) += qs->base.data[d].size, pos[0] += qs->base.data[d].num;
							if (b == 1) {
								for (; tmp = pos[0] + diff, end > pos[0] && end > tmp;)
									*pos[0] += qs->base.data[d].size, *tmp += qs->base.data[d].size, pos[0] += qs->base.data[d].num;
								pos[1] = pos[0] + diff;
							} else
								for (pos[1] = pos[0] + diff; end > pos[0] && end > pos[1];)
									for (e = 0; e < times; ++e)
										*pos[e] += qs->base.data[d].size, pos[e] += qs->base.data[d].num;
							for (e = 0; e < 2; ++e)
								if (end > pos[e])
									*pos[e] += qs->base.data[d].size, pos[e] += qs->base.data[d].num;
							if (cond_2)
								for (e = 0; e < 2; ++e)
									qs->others.offsets[e][d] = pos[e];
						}
					}
			}
		}
	}
}

static inline void register_relations(qs_sheet *qs, const cint *A, const cint *B, const cint *C) {
	cint *D = &qs->vars.temp[0], *E = D + 1, *KEY = D + 2, *VALUE = D + 3;
	qs_sm a, b, bits, extra, mod, v_1, v_2, verification, *data;
	for (a = 0; a < qs->info.m.double_value; ++a)
		if(qs->others.sieve[a] >= qs->info.threshold){
			simple_int_to_cint(D, a);
			cint_subi(D, qs->constants.M_2); // X
			cint_mul(A, D, E); // AX
			cint_addi(E, B); // AX + B
			cint_dup(KEY, E); // copy of AX + B
			cint_addi(E, B); // AX + 2B
			cint_mul(E, D, VALUE); // AX^2 + 2BX
			cint_addi(VALUE, C); // AX^2 + 2BX + C
			KEY->nat = VALUE->nat = 1 ; // abs values
			bits = (qs_sm) cint_count_bits(VALUE) - qs->info.error_bits;
			extra = 0, verification = 1;
			data = qs->others.md_uncleared_buffer; // buffered data may be used by next function.
			for (b = 0; b < 2; ++b) {
				if (qs->base.data[b].num != 1) {
					simple_int_to_cint(D, qs->base.data[b].num);
					*data = cint_remove(qs->calc, VALUE, D);
					if (*data) {
						extra += b ? *data : qs->base.data[b].size;
						*++data = b;
						++data;
					}
				}
			}
			for (b = 2; b < qs->info.p_list[1] && verification; ++b) {
				v_1 = (v_2 = 0, qs->base.data[b].sol[1] == -1U) || (v_2 = 1, mod = a % qs->base.data[b].num, mod == qs->base.data[b].sol[0] || mod == qs->base.data[b].sol[1]);
				if (v_1) {
					simple_int_to_cint(D, qs->base.data[b].num);
					*data = cint_remove(qs->calc, VALUE, D);
					verification &= v_2 <= *data;
					if (*data) {
						*++data = b, ++data;
						extra += qs->base.data[b].size;
					}
				}
			}
			if (qs->others.sieve[a] += extra, qs->others.sieve[a] >= bits) {
				for (; b < qs->info.p_list[4] && extra < qs->others.sieve[a] && verification; ++b) {
					v_1 = (v_2 = 0, qs->base.data[b].sol[1] == -1U) || (v_2 = 1, mod = a % qs->base.data[b].num, mod == qs->base.data[b].sol[0] || mod == qs->base.data[b].sol[1]);
					if (v_1) {
						simple_int_to_cint(D, qs->base.data[b].num);
						*data = cint_remove(qs->calc, VALUE, D);
						verification &= v_2 <= *data;
						if (*data || v_2) {
							*++data = b, ++data;
							extra += qs->base.data[b].size;
						}
					}
				}
				const uint8_t mask = 1 << (a & 7);
				for (b = qs->info.p_list[4]; b < qs->base.length && extra < qs->others.sieve[a] && verification; ++b)
					if (qs->others.flags[b] & mask) {
						if (mod = a % qs->base.data[b].num, mod == qs->base.data[b].sol[0] || mod == qs->base.data[b].sol[1]) {
							simple_int_to_cint(D, qs->base.data[b].num);
							*data = cint_remove(qs->calc, VALUE, D);
							verification &= *data != 0;
							++data, *data++ = b;
							extra += qs->base.data[b].size;
						}
					}
				if (verification) {
					if (h_cint_compare(VALUE, qs->constants.ONE) == 0) {
						qs_sm *open_1 = qs->others.sm_buffer;
						qs_sm *close_1 = open_1 + qs->s.values.double_value;
						qs_sm *open_2 = qs->others.md_uncleared_buffer;
						register_relation_kind_1(qs, KEY, open_1, close_1, open_2, data);
					} else if (h_cint_compare(VALUE, qs->constants.UPPER) < 0)
						register_relation_kind_2(qs, data, KEY, VALUE);
				}
				// If verification is false : "An abnormally unusable data row has been ignored during relation registering".
			}
		}
}

static int qs_register_factor(qs_sheet *qs) {
	// returns -1 if the algorithm should stop, accept any divisor of N.
	cint * F = qs->vars.FACTOR ;
	int i, res = h_cint_compare(F, qs->constants.ONE) && h_cint_compare(F, qs->vars.N) ;
	if (res) {
		const struct avl_node *node = avl_at(&qs->divisors.tree, F);
		if (qs->divisors.tree.affected) {
			fac_cint *ans = &qs->caller->factor;
			for (i = 0; i < 2; ++i) {
				ans->prime = cint_is_prime(qs->calc, F, -1);
				if (ans->prime) {
					cint_dup(&ans->cint, F);
					// 200-bit RSA take about 10,000,000+ "duplications", so perform the last.
					ans->power = qs->caller->number->power * (int) cint_remove(qs->calc, qs->vars.N, F);
					assert(ans->power);
					fac_push(qs->caller, ans, 1);
					++qs->divisors.total_primes;
					if (i || h_cint_compare(qs->vars.N, qs->constants.ONE) == 0)
						i = 1, res = -1;
					else cint_dup(F, qs->vars.N);
				} else {
					if (i == 0) qs->divisors.data[qs->divisors.length++] = node->key;
					break;
				}
			}
		}
	}
	return res ;
}

static inline void process_column_array(struct qs_relation * restrict rel, const qs_sm * restrict ptr){
	qs_sm a, b;
	if (*ptr & 1) {
		for (a = 0; a < rel->Y.length; ++a)
			if (rel->Y.data[a] == *(ptr + 1)) {
				for (b = a; b + 1 < rel->Y.length; ++b)
					rel->Y.data[b] = rel->Y.data[b + 1];
				--rel->Y.length;
				a = -2U;
			}
		if (a != -1U)
			rel->Y.data[rel->Y.length++] = *(ptr + 1);
	}
	for (a = 0; a < *ptr; ++a)
		rel->axis.Z.data[rel->axis.Z.length++] = *(ptr + 1);
}

static inline void register_relation_kind_1(qs_sheet *qs, const cint *KEY, qs_sm *p_1, const qs_sm *p_2, qs_sm *p_3, const qs_sm *p_4) {
	struct avl_node *node = avl_at(&qs->relations.tree, KEY);
	if (node->value)
		return; // duplicates at this stage are ignored.
	char *open = qs->mem.now = mem_aligned(qs->mem.now), *close;
	assert(open + (1 << 21) < (char*)qs->mem.base + qs->info.total_bytes_allocated);
	// between "open" and "close" data will be stored (commit) or be zeroed (rollback).
	struct qs_relation *rel = qs->mem.now;
	qs->mem.now = rel + 1 ;
	// create a new relation.
	rel->X = node->key; // constant X is const-stored by the node key.
	rel->Y.data = qs->mem.now; // data Y has a known length which only decrease.
	rel->axis.Z.data = rel->Y.data + (p_2 - p_1) + (p_4 - p_3); // writes Z ahead.
	for (; p_1 < p_2; p_1 += 2) process_column_array(rel, p_1);
	for (; p_3 < p_4; p_3 += 2) process_column_array(rel, p_3);
	qs->mem.now = rel->axis.Z.data + rel->axis.Z.length;
	cint *A = &qs->vars.temp[0], *B = A + 1, *C = A + 2, *D = A + 3;
	cint_reinit(A, 1);
	for (qs_sm a = 0; a < rel->axis.Z.length; ++a) {
		simple_int_to_cint(B, qs->base.data[rel->axis.Z.data[a]].num);
		cint_mul_modi(qs->calc, A, B, qs->constants.kN);
	}
	cint_mul(rel->X, rel->X, C);
	cint_div(qs->calc, C, qs->constants.kN, D, B);
	if (cint_compare(A, B) && (cint_addi(A, B), cint_compare(A, qs->constants.kN))) {
		close = qs->mem.now;
		qs->mem.now = memset(open, 0, close - open);
	} else {
		// relation saved, commit, memory updated.
		node->value = qs->relations.data[qs->relations.length.now] = rel;
		qs->mem.now = rel->axis.Z.data + rel->axis.Z.length;
		rel->id = ++qs->relations.length.now;
	}
}

static inline void register_relation_kind_2(qs_sheet *qs, const qs_sm *data_end, const cint *KEY, const cint *VALUE) {
	if (qs->info.kn_bits < 160)
		return; // not faster, not slower during tests.

	// the function searches 2 different KEY sharing the same VALUE.
	// keys and values are stored in the same tree, there is no collisions.
	struct avl_node *node = avl_at(&qs->relations.tree, VALUE);
	struct qs_relation *old, *new;
	cint *BEZOUT = 0;
	qs_sm a = 0, b, *data, *open, *close;
	old = node->value;
	if (old) {
		if (old->id)
			return; // normally no collision.
		if (old->X == 0)
			return; // modular inverse (p, number) failed.
		for (; old; old = old->axis.next)
			if (h_cint_compare(KEY, old->X) == 0)
				return; // same KEY already registered.
		old = node->value;
		if (old->axis.next == 0) {
			// this duplicate isn't 3rd, 4th .. it's the 2nd ... so compute BEZOUT
			cint *A = &qs->vars.temp[0], *B = A + 1;
			cint_modular_inverse(qs->calc, VALUE, qs->constants.kN, A);
			if (A->mem == A->end) {
				// no solution to the linear congruence.
				old->X = 0;
				cint_div(qs->calc, qs->constants.kN, VALUE, A, B);
				a = simple_cint_to_int(B);
				b = simple_cint_to_int(VALUE);
				for (; (a %= b) && (b %= a););
				simple_int_to_cint(A, a |= b);
				cint_div(qs->calc, qs->vars.N, A, qs->vars.FACTOR, B);
				if (B->mem == B->end)
					qs_register_factor(qs);
				// nothing.
				return;
			} else BEZOUT = A;
		}
	}

	new = mem_aligned(qs->mem.now);
	qs->mem.now = new + 1;
	new->X = qs->mem.now;

	if (BEZOUT) {
		// BEZOUT is stored directly after X, newcomer become the root of the linked list.
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
	a = qs->s.values.double_value;
	memcpy(data, qs->others.sm_buffer, a * sizeof(*data));
	data += a, a = data_end - qs->others.md_uncleared_buffer;
	memcpy(data, qs->others.md_uncleared_buffer, a * sizeof(*data));
	new->Y.length = data + a - new->Y.data;
	qs->mem.now = new->Y.data + new->Y.length;

	if (old) {
		cint *A = &qs->vars.temp[0], *B = A + 1, *C = A + 2, *NEW_KEY = A + 3;
		BEZOUT = old->X + 1; // BEZOUT was stored here.
		cint_mul(BEZOUT, new->X, B);
		for (; old->axis.next; old = old->axis.next)
			if (old != new) {
				// combines.
				cint_mul(B, old->X, A);
				cint_div(qs->calc, A, qs->constants.kN, C, NEW_KEY);
				// it's using 2 buffers.
				open = close = qs->others.md_cleared_buffer;
				data = memset(qs->others.md_uncleared_buffer, 0, qs->base.length * sizeof(*data));
				for (a = 0; a < new->Y.length; a += 2)
					data[new->Y.data[a + 1]] += new->Y.data[a];
				for (a = 0; a < old->Y.length; a += 2)
					data[old->Y.data[a + 1]] += old->Y.data[a];
				for (a = 0; a < qs->base.length; ++a)
					if (data[a])
						*close++ = data[a], *close++ = a;
				register_relation_kind_1(qs, NEW_KEY, open, close, 0, 0);
				memset(open, 0, (char *) close - (char *) open); // zeroed.
			}
	}
}

static inline void finalization_part_1(qs_sheet *qs, const uint64_t *lanczos_answer) {
	const uint64_t mask = *lanczos_answer, *null_rows = lanczos_answer + 1;
	// lanczos answer isn't a struct, it's "mask followed by null_rows".
	if (mask == 0 || null_rows == 0)
		return;

	cint *A = &qs->vars.temp[0], *B = A + 1, *C = A + 2, *D = A + 3;
	qs_md a, b, c;
	qs_sm *power_of_primes;
	//for (a = b = 0; a < 64; ++a) b += (lanczos_mask & 1LLU << a) != 0; assert(b);
	for (c = 0; c < 64; ++c) {
		for (; !(mask & 1LLU << c); ++c);
		cint_reinit(A, 1), cint_reinit(B, 1), cint_reinit(C, 1);
		power_of_primes = memset(qs->others.md_uncleared_buffer, 0, qs->base.length * sizeof(*power_of_primes));
		for (a = b = 0; a < qs->relations.length.now; ++a) {
			if (null_rows[a] & 1LLU << c) {
				const struct qs_relation * restrict const rel = qs->relations.data[a];
				cint_mul_modi(qs->calc, A, rel->X, qs->vars.N);
				for (b = 0; b < rel->axis.Z.length; ++b)
					++power_of_primes[rel->axis.Z.data[b]];
			}
		}
		for (b = 0; b < qs->base.length; ++b)
			if (power_of_primes[b]){
				simple_int_to_cint(B, qs->base.data[b].num);
				if (power_of_primes[b] > 1) {
					simple_int_to_cint(D, power_of_primes[b] >> 1);
					cint_pow_modi(qs->calc, B, D, qs->vars.N);
				}
				cint_mul_modi(qs->calc, C, B, qs->vars.N);
			}
		h_cint_subi(C, A);
		if (C->mem != C->end) {
			C->nat = 1; // ABS(C)
			cint_gcd(qs->calc, qs->vars.N, C, qs->vars.FACTOR);
			if (qs_register_factor(qs) == -1)
				break;
		}
	}
}

static inline void finalization_part_2(qs_sheet *qs) {
	if (cint_count_bits(qs->vars.N) == 1)
		return;

	cint * F = qs->vars.FACTOR, **di = qs->divisors.data, *Q = qs->vars.temp, *R = Q + 1 ;
	qs_sm i, j, k, count = 0 ;

	for(i = qs->divisors.processing_index, k = qs->divisors.length; i < k; ++i)
		for (j = 1 + i; j < k; ++j)
			if (cint_gcd(qs->calc, di[i], di[j], F), qs_register_factor(qs) == -1)
				return; // gcd of new divisors with old divisors
	do {
		for (; i = k, k = qs->divisors.length, i != k;)
			for (; i < k; ++i)
				for (j = 0; j < i; ++j)
					if (cint_gcd(qs->calc, di[i], di[j], F), qs_register_factor(qs) == -1)
						return; // gcd of new divisors with old divisors

		for (i = qs->divisors.processing_index; i < k; ++i)
			if (fac_any_root_check(qs->caller, di[i], qs->vars.FACTOR, R))
				if (qs_register_factor(qs) == -1)
					return; // perfect power of new divisors

		if (count != qs->divisors.total_primes)
			if (fac_any_root_check(qs->caller, qs->vars.N, qs->vars.FACTOR, R))
				if (qs_register_factor(qs) == -1)
					return; // if N changed, check prefect power of N

		count = qs->divisors.total_primes ;
		qs->divisors.processing_index = k ;
	} while(k != qs->divisors.length); // until no new divisor

}

static inline int finalization_part_3(qs_sheet *qs) {
	// Usually do nothing, N equals 1
	// Otherwise push something non-trivial to the caller's routine
	int res = h_cint_compare(qs->vars.N, qs->constants.ONE) == 0 ;
	if (res == 0){
		fac_cint *ans = &qs->caller->factor;
		ans->prime = 0 ;
		ans->power = qs->caller->number->power ;
		for(qs_sm i = 0; i < qs->divisors.length; ++i) {
			if (qs->caller->params->silent == 0) {
				char *str = cint_to_string(qs->divisors.data[i], 10);
				printf("- quadratic sieve can't silently ignore %s\n", str);
				free(str);
			}
			if (h_cint_compare(qs->vars.N, qs->divisors.data[i]) > 0) {
				const int power = (int) cint_remove(qs->calc, qs->vars.N, qs->divisors.data[i]);
				if (power) {
					assert(power == 1);
					cint_dup(&ans->cint, qs->divisors.data[i]);
					fac_push(qs->caller, ans, 0);
				}
			}
		}
		res = h_cint_compare(qs->vars.N, &qs->caller->number->cint) ;
		if (res) // res is true if QS was able to decompose N.
			if (h_cint_compare(qs->vars.N, qs->constants.ONE)){
				cint_dup(&ans->cint, qs->vars.N);
				fac_push(qs->caller, ans, 0);
			}
	}
	return res;
}
