// Pure C99 Self-Initializing Quadratic Sieve (SIQS) released "as it" into the public domain, without any warranty, express or implied.

int factorization_quadratic_sieve(state * state, int bits) {
	// The state contain the request (a number N) and the Quadratic Sieve :
	// - must register (ideally prime) factors of N (using the state)
	// - must update N (present in the state) accordingly to the registered factors
	// - can use resources present in the state (temporary variables + parameters)
	// - return a non-zero value if a factor (not necessarily prime) was registered

	if (bits < 65 || (220 < bits && !state->params.force))
		return 0; // The execution time is doubled for each 10-bit increase.

	// On the testing machine 260 bits take 12 minutes, 270-bit take 24 minutes.

	qs_sheet qs = {0};
	preparation_part_1(&qs, state, bits);
	preparation_part_2(&qs);
	preparation_part_3(&qs);
	qs_parametrize(&qs);
	preparation_part_4(&qs);
	preparation_part_5(&qs);
	preparation_part_6(&qs);
	do {
		do {
			// Keep randomly trying various polynomials.
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
		// Analyzes all observations made by the algorithm.
		finalization_part_1(&qs, lanczos_block(&qs));
	} while (outer_continuation_condition(&qs));
	const int res = finalization_part_2(&qs);
	free(qs.mem.base);
	return res;
}

// Quadratic sieve main condition 1. Continue to search relations or break ?
int inner_continuation_condition(qs_sheet *qs) {
	int answer = 1 ;
	if (qs->time.tick % 16 == 1 && (qs->time.end || qs->state->params.verbose))
		qs->time.now = get_time_ms();
	answer &= qs->n_bits != 1 ; // the bit count of N may have changed.
	answer &= (qs->relations.length.peak = qs->relations.length.now) < qs->relations.length.needs; // the condition.
	answer &= !qs->time.end || qs->time.tick % 16 == 1 || qs->time.now < qs->time.end ;
	if (qs->state->params.verbose) {
		const double rel_begin = (double) qs->relations.length.now, rel_end = (double) qs->relations.length.needs ;
		/*if (qs->time.tick % 16 == 1)
			time_remaining(qs);*/
		display_progress("Quadratic Sieve", 100. * rel_begin / rel_end); // progress isn't linear.
	}
	return answer;
}

// Quadratic sieve main condition 2. Return to sieving or stop the algorithm ?
int outer_continuation_condition(qs_sheet *qs) {
	int answer = 1 ;
	answer &= qs->sieve_again_perms-- > 0; // avoid infinite loop.
	answer &= qs->divisors.total_primes < qs->sieve_again_perms; // search prime factors.
	answer &= qs->n_bits != 1 ; // the bit count of N isn't 1.
	answer &= !qs->time.end || get_time_ms() < qs->time.end ;
	if (answer) {
		qs_sm new_needs = qs->relations.length.needs;
		new_needs += new_needs >> (1 + qs->sieve_again_perms);
		DEBUG_CRITICAL("Quadratic Sieve targets %u more relations.\n", new_needs - qs->relations.length.needs);
		qs->relations.length.needs = new_needs ;
	}
	return answer;
}

void qs_parametrize(qs_sheet *qs) {

	const qs_sm bits = qs->kn_bits; // N adjusted has at least 115-bit.
	qs->kn_bits = (qs_sm) cint_count_bits(qs->state->session.tmp); // kN may be slight larger.

	DEBUG_NORMAL("N is a %u-bit number, and kN is a %u-bit number using %u words.\n", (qs_sm) cint_count_bits(&qs->state->session.num), qs->kn_bits, (unsigned)(qs->state->session.tmp->end - qs->state->session.tmp->mem));

	DEBUG_NORMAL("Default parameters are to factor a %u-bit number.\n", bits);

	qs_md tmp ;
	// params as { bits, value } take the extremal value if bits exceed.
	static const double param_base_size [][2]= { {135, 1300}, {165, 4200}, {200, 10000}, {260, 20000}, {330, 55000}, {0} };
	qs->base.length = (tmp = qs->state->params.qs_base_size) ? tmp : linear_param_resolution(param_base_size, bits);
	DEBUG_NORMAL("The factor base will include %u prime numbers.\n", qs->base.length);

	static const double param_laziness [][2]= {{150, 95}, {250, 100}, {0} };
	// collecting more/fewer relations than recommended (in percentage).
	qs->relations.length.needs = qs->base.length * ((tmp = qs->state->params.qs_laziness) ? tmp : linear_param_resolution(param_laziness, bits)) / 100;
	DEBUG_NORMAL("The algorithm use the seed %" PRIu64 " and targets %u relations.\n", qs->state->params.rand.custom, qs->relations.length.needs);

	static const double param_m_value [][2]= { {120, 1}, {330, 6}, {0} };
	qs->m.length = (qs->m.length_half = (qs->state->params.qs_sieve ? qs->state->params.qs_sieve : 31744) * linear_param_resolution(param_m_value, bits)) << 1;

	qs->m.cache_size = 95232 ; // algorithm reaches "M length" by steps of "cache size".

	static const double param_error [][2]= { {120, 15}, {330, 35}, {0} };
	// tolerance for small errors in the smoothness test during relation collection.
	qs->error_bits = (tmp = qs->state->params.qs_error_bits) ? tmp : linear_param_resolution(param_error, bits);

	static const double param_threshold [][2]= { {120, 60}, {330, 110}, {0} };
	qs->threshold.value = (tmp = qs->state->params.qs_threshold) ? tmp : linear_param_resolution(param_threshold, bits);

	static const double param_alloc [][2]= { {130, 992}, {140, 1280}, {150, 2176}, {160, 3584}, {170, 7168}, {180, 12288}, {190, 14336}, {190, 14336}, {200, 24576}, {210, 30720}, {220, 40960}, {230, 49152}, {240, 57344}, {250, 67584}, {260, 81920}, {270, 98304}, {280, 114688}, {290, 122880}, {300, 139264}, {310, 163840}, {320, 196608}, {330, 229376}, {0} };
	qs->mem.bytes_allocated = (tmp = qs->state->params.qs_alloc_mb) ? tmp << 20 : linear_param_resolution(param_alloc, qs->kn_bits) << 10;

	qs->sieve_again_perms = 3; // Sieve again up to 3 times before giving up

	// Iterative list
	qs->iterative_list[0] = 1; // one
	static const double param_first_prime [][2]= { {170, 8}, {210, 12}, {320, 40}, {0} };
	qs->iterative_list[1] = linear_param_resolution(param_first_prime, bits); // first
	const qs_sm large_base = 5120 < qs->base.length ? 5120 : qs->base.length;
	qs->iterative_list[2] = large_base >> 2; // medium
	qs->iterative_list[3] = large_base >> 1; // mid
	qs->iterative_list[4] = large_base; // sec

	DEBUG_NORMAL("Iterative list contain %u, %u, %u and %u.\n", qs->iterative_list[1], qs->iterative_list[2], qs->iterative_list[3], qs->iterative_list[4]);

	const qs_md last_prime_in_base = (qs_md) (qs->base.length * 2.5 * log_computation(qs->base.length));
	qs->relations.too_large_prime = (tmp = qs->state->params.qs_large_prime) ? tmp : last_prime_in_base << 4;

	DEBUG_NORMAL("Primes greater than or equal to %" PRIu64 " are ignored.\n", qs->relations.too_large_prime);

	qs->s.values.double_value = (qs->s.values.defined = (qs->s.values.subtract_one = bits / 28) + 1) << 1;
	qs->poly.gray_max = 1 << (qs->s.values.defined - 3); // computing the roots of f(x) once for all these polynomials.

	DEBUG_NORMAL("Other params include sieve=%u, error_bits=%u, threshold=%u and s=%u.\n", qs->m.length_half, qs->error_bits, qs->threshold.value, qs->s.values.defined);

	// The algorithm itself completes its configuration during the last preparation part.
	assert(qs->s.values.defined >= 3);

}

// Quadratic sieve parameters (utility, dev configuration)
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

// Quadratic sieve source (algorithm)
void preparation_part_1(qs_sheet *qs, state *state, int bits) {
	// initializing (until kN is computed) with the caller's resources.
	qs->state = state;
	DEBUG_NORMAL("\nQuadratic Sieve for %s.\n", simple_cint_string(state, &state->session.num));
	qs->sheet = state->session.sheet;
	qs->seed = state->params.rand.seed;
	qs->n_bits = qs->kn_bits = bits;
	if (qs->time.start = get_time_ms(), qs->state->params.timeout)
		qs->time.end = qs->time.start + 1000 * qs->state->params.timeout;
}

void preparation_part_2(qs_sheet *qs) {
	// Not N, but kN is transparently adjusted by a prime number to measure at least 115-bit.
	cint * N = &qs->state->session.num, * kN = qs->state->session.tmp, *ADJUSTOR = kN + 1 ;
	static const int prime_generator[] = {
			9, 7, 5, 3, 17, 27, 3, 1, 29, 3, 21, 7, 17, 15,
			9, 43, 35, 15, 29, 3, 11, 3, 11, 15, 17, 25, 53,
			31, 9, 7, 23, 15, 27, 15, 29, 7, 59, 15, 5, 21,
			69, 55, 21, 21, 5, 159, 3, 81, 9, 69, 131, 33, 15 };
	const qs_sm bits = (qs_sm) qs->n_bits;
	if (bits < 115) {
		qs->adjustor = (1LLU << (124 - bits)) + prime_generator[115 - bits] ;
		simple_int_to_cint(ADJUSTOR, qs->adjustor);
		cint_mul(N, ADJUSTOR, kN);
		qs->kn_bits = (qs_sm) cint_count_bits(kN);
		DEBUG_NORMAL("Input (%u bits) is multiplied by %" PRIu64 " to reach %u bits.\n", bits, qs->adjustor, qs->kn_bits);
	} else
		qs->adjustor = 1, cint_dup(kN, N);
}

void preparation_part_3(qs_sheet *qs) {
	// Frequently select a small multiplier that will outperform a factorization done without a multiplier.
	// After it, the algorithm will factor kN instead of N, where k is a small constant named "multiplier".
	const qs_sm mul = (qs_sm) qs->state->params.qs_multiplier ;
	if (mul){
		DEBUG_NORMAL("The multiplier is %u.\n", mul);
		qs->multiplier = mul ;
	} else {
		const size_t total_best = 7;
		qs_sm best[total_best];
		for (int i = qs->state->params.verbose < 2; i < 2; ++i) {
			if (i)
				preparation_part_3_default_proposition(qs, best, total_best);
			else
				preparation_part_3_alternative_proposition(qs, best, total_best);
			DEBUG_NORMAL("%s", "Suggested multipliers are [");
			for (size_t j = 0; j < total_best - 1; ++j)
				DEBUG_NORMAL("%u, ", best[j]);
			DEBUG_NORMAL("%u]%s", best[total_best - 1], i ? "" : ".\n");
		}
		qs->multiplier = *best;
		DEBUG_NORMAL(", so use %u.\n", *best);
	}
	if (qs->multiplier > 1) {
		cint *kN = qs->state->session.tmp, *MUL = kN + 1, *N = kN + 2 ;
		simple_int_to_cint(MUL, qs->multiplier);
		cint_dup(N, kN);
		cint_mul(MUL, N, kN);
	}
}

static inline void preparation_part_3_default_proposition(qs_sheet *qs, qs_sm *caller_res, const size_t caller_res_len) {
	// Choose a multiplier that make the input more favorable for smoothness
	// over the future factor base, and lead to faster relation gathering.
	struct {
		qs_sm mul;
		double score;
	} res[128] ; // 127 is the greatest multiplier.

	cint *N = qs->state->session.tmp, *PRIME = N + 1, *Q = N + 2, *R = N + 3;
	const double log_2 = 0.6931471805599453;
	const size_t len = sizeof(res) / sizeof(*res) - 1;
	for (qs_sm i = 0; i < len ; ++i) {
		res[i].mul = i + 1;
		res[i].score = -0.5 * log_computation(res[i].mul);
		switch (*N->mem * res[i].mul % 8) {
			// Special case against 2, the first prime number.
			case 3 : case 7 : res[i].score += 0.5 * log_2; break;
			case 5 : res[i].score += 1.0 * log_2; break;
			case 1 : res[i].score += 3.0 * log_2; break;
		}
	}

	for (qs_sm prime = 3; prime < 312; prime += 2)
		if (is_prime_4669913(prime)) {
			// Normal case against the odd primes.
			simple_int_to_cint(PRIME, prime);
			cint_div(qs->sheet, N, PRIME, Q, R);
			const qs_sm n_mod_prime = (qs_sm) simple_cint_to_int(R);
			const double intake = 2.0 / (prime - 1) * log_computation(prime);
			const int kronecker = kronecker_symbol(n_mod_prime, prime);
			for (qs_sm i = 0; i < len; ++i)
				if (kronecker * kronecker_symbol(res[i].mul, prime) == 1)
					res[i].score += intake;
		}

	// Sort the results.
	for (int i = 0; i < len; ++i)
		for (int j = 1 + i; j < len; ++j)
			if (res[i].score < res[j].score)
				res[len] = res[i], res[i] = res[j], res[j] = res[len];

	for(int i = 0; i < caller_res_len; ++i)
		caller_res[i] = res[i].mul ;
}

static inline void preparation_part_3_alternative_proposition(qs_sheet *qs, qs_sm *caller_res, const size_t caller_res_len) {
	// Choose a multiplier that make the input more favorable for smoothness
	// over the future factor base, and lead to faster relation gathering.
	struct {
		qs_sm mul;
		double score;
	} res[] = {{1, 0}, {2, 0}, {3, 0}, {5, 0}, {6, 0}, {7, 0}, {10, 0}, {11, 0}, {13, 0}, {14, 0}, {15, 0}, {17, 0}, {19, 0}, {21, 0}, {22, 0}, {23, 0}, {26, 0}, {29, 0}, {30, 0}, {31, 0}, {33, 0}, {34, 0}, {35, 0}, {37, 0}, {38, 0}, {39, 0}, {41, 0}, {42, 0}, {43, 0}, {46, 0}, {47, 0}, {51, 0}, {53, 0}, {55, 0}, {57, 0}, {58, 0}, {59, 0}, {61, 0}, {62, 0}, {65, 0}, {66, 0}, {67, 0}, {69, 0}, {70, 0}, {71, 0}, {73, 0}, {79, 0}, {83, 0}, {89, 0}, {97, 0}, {101, 0}, {103, 0}, {107, 0}, {109, 0}, {0, 0}};

	cint *N = qs->state->session.tmp, *TMP = N + 1, *Q = N + 2, *R = N + 3;
	const double log_2 = 0.6931471805599453;
	const size_t len = sizeof(res) / sizeof(*res) - 1;
	for (qs_sm i = 0; i < len; ++i) {
		res[i].score = -0.5 * log_computation(res[i].mul);
		switch (*N->mem * res[i].mul % 8) {
			// Special case against 2, the first prime number.
			case 3 : case 7 : res[i].score += 0.5 * log_2; break;
			case 5 : res[i].score += 1.0 * log_2; break;
			case 1 : res[i].score += 3.0 * log_2; break;
		}
	}

	for (qs_sm prime = 3; prime < 504; prime += 2)
		if (is_prime_4669913(prime)) {
			// Normal case against the odd primes.
			simple_int_to_cint(TMP, prime);
			cint_div(qs->sheet, N, TMP, Q, R);
			const qs_sm n_mod_prime = (qs_sm) simple_cint_to_int(R);
			const double intake = log_computation(prime) / (prime - 1);
			for (qs_sm i = 0; i < len; ++i) {
				const qs_sm kn_mod_prime = n_mod_prime * res[i].mul % prime;
				if (kn_mod_prime == 0)
					res[i].score += intake;
				else if (kronecker_symbol(kn_mod_prime, prime) == 1)
					res[i].score += 2.0 * intake;
			}
		}

	// Sort the results.
	for (int i = 0; i < len; ++i)
		for (int j = 1 + i; j < len; ++j)
			if (res[i].score < res[j].score)
				res[len] = res[i], res[i] = res[j], res[j] = res[len];

	for(int i = 0; i < caller_res_len; ++i)
		caller_res[i] = res[i].mul ;
}

void preparation_part_4(qs_sheet *qs) {
	void *mem;
	mem = qs->mem.base = calloc(1, qs->mem.bytes_allocated);
	assert(mem);

	// kN was computed into the caller's memory, now the QS has parametrized and allocated
	const size_t kn_size = qs->state->session.tmp[0].end - qs->state->session.tmp[0].mem + 1 ;
	// the quadratic sieve variables can store at most kN ^ 2 in terms of bits
	const size_t vars_size = kn_size << 1 ;

	DEBUG_NORMAL("Quadratic Sieve use %u-bit variables for computations.\n", (unsigned)(vars_size * cint_exponent));

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
				// "MY" is used to facilitate development controls
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
				&qs->constants.TOO_LARGE_PRIME,
				&qs->constants.MULTIPLIER,
				0,
		};
		for (int i = 0; n[i]; ++i) {
			n[i]->mem = n[i]->end = mem_aligned(mem) ;
			mem = n[i]->mem + (n[i]->size = vars_size);
		}
	}

	cint_dup(&qs->vars.N, &qs->state->session.num);
	cint_dup(&qs->constants.kN, qs->state->session.tmp);

	simple_int_to_cint(&qs->constants.ONE, 1);
	simple_int_to_cint(&qs->constants.M_HALF, qs->m.length_half);
	simple_int_to_cint(&qs->constants.TOO_LARGE_PRIME, qs->relations.too_large_prime);
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
	qs->m.sieve[qs->m.length] = 0xFF ; // the end of the sieve evaluates to "true" under any "truthy" mask.
	qs->m.flags = mem_aligned(qs->m.sieve + qs->m.length + sizeof(uint64_t));
	// Usage: buffer[0] is zeroed after use, buffer[1] isn't supposed zeroed after use.
	qs->buffer[0] = mem_aligned(qs->m.flags + qs->base.length);
	qs->buffer[1] = mem_aligned(qs->buffer[0] + buffers_size);

	// Other allocations
	qs->relations.length.capacity = (qs_sm) relations_size ;
	// Lanczos Block has a part of memory, it takes a "lite" snapshot before throwing relations.
	qs->lanczos.snapshot = mem_aligned(qs->buffer[1] + buffers_size) ;
	qs->relations.data = mem_aligned(qs->lanczos.snapshot + relations_size);
	qs->divisors.data = mem_aligned(qs->relations.data + relations_size);
	qs->mem.now = mem_aligned(qs->divisors.data + 127);

	const qs_sm n_trees = (qs_sm) (sizeof(qs->uniqueness) / sizeof(struct avl_manager));
	for (qs_sm i = 0; i < n_trees; ++i) {
		// the trees are used to identify duplicates (relations, partials, factors of N)
		qs->uniqueness[i].inserter_argument = &qs->mem.now;
		qs->uniqueness[i].inserter = &avl_cint_inserter;
		qs->uniqueness[i].comparator = (int (*)(const void *, const void *)) &h_cint_compare;
		// they use default sign-less comparator.
	}
	DEBUG_NORMAL("Allocated %u MB of memory with a %u KB structure.\n", qs->mem.bytes_allocated >> 20, (unsigned)((char*)qs->mem.now - (char*)qs->mem.base) >> 10);
}

void preparation_part_5(qs_sheet *qs) {
	// Prepare the factor base (a set of small prime numbers used to find smooth numbers).
	static const double inv_ln_2 = 1.4426950408889634;
	cint *A = qs->vars.TEMP, *B = A + 1, *C = A + 2;
	qs_sm i = 0, prime;

	// the factor base contain the multiplier if different from 2.
	if (qs->multiplier != 2)
		qs->base.data[i].size = (qs_sm) (.35 + inv_ln_2 * log_computation(qs->base.data[i].num = qs->multiplier)), ++i;

	// then it contains the number 2.
	qs->base.data[i].num = 2, qs->base.data[i].size = 1;
	qs->base.data[i].sqrt_kN_mod_prime = *qs->constants.kN.mem % 8 == 1 || *qs->constants.kN.mem % 8 == 7, ++i;

	// then the prime numbers for which kN is a quadratic residue modulo.
	for (prime = 3; i < qs->base.length; prime += 2)
		if (is_prime_4669913(prime)) {
			simple_int_to_cint(A, prime);
			cint_div(qs->sheet, &qs->constants.kN, A, B, C);
			const qs_sm kn_mod_prime = (qs_sm) simple_cint_to_int(C);
			qs->base.data[i].sqrt_kN_mod_prime = tonelli_shanks(kn_mod_prime, prime);
			// Update: if a prime is a prime factor of the multiplier, it's included too.
			if (qs->base.data[i].sqrt_kN_mod_prime) {
				qs->base.data[i].num = prime;
				qs->base.data[i].size = (qs_sm) (.35 + inv_ln_2 * log_computation(prime)), ++i;
			}
		}
	// 2.5 * (base size) * ln(base size) is close to the largest prime number in factor base.
	qs->base.largest = qs->base.data[i - 1].num ;
	DEBUG_NORMAL("The largest prime in the factor base is %u.\n", qs->base.largest);
}

void preparation_part_6(qs_sheet *qs) {
	// completes the configuration by the algorithm itself.
	// computes D : a template for the A polynomial coefficient.
	qs_sm i, min, span;
	const qs_sm s = qs->s.values.defined ;
	qs->poly.span.half = (span = qs->base.length / (s * (s + s))) >> 1;
	cint *kN = qs->vars.TEMP, *TMP = kN + 1, *R = kN + 2;
	cint_dup(kN, &qs->constants.kN);
	cint_left_shifti(kN, 1);
	cint_sqrt(qs->sheet, kN, TMP, R);
	cint_div(qs->sheet, TMP, &qs->constants.M_HALF, &qs->poly.D, R);
	qs->poly.d_bits = (qs_sm) cint_count_bits(&qs->poly.D);
	cint_nth_root(qs->sheet, &qs->poly.D, s, R); // use the s-th root of D.
	const qs_sm root = (qs_sm) simple_cint_to_int(R) ;
	for (i = 1; assert(i < qs->base.length), qs->base.data[i].num <= root; ++i);
	assert(i >= span);
	for (min = i - qs->poly.span.half, i *= i; i / min < span + min; --min);
	qs->poly.span.x_1 = min ;
	qs->poly.span.x_2 = min + qs->poly.span.half ;
	qs->poly.span.x_3 = qs->poly.span.x_2 * qs->poly.span.x_2 ;
	assert(qs->poly.span.x_1 < qs->base.length);
	assert(qs->poly.span.x_2 < qs->base.length);
}

void get_started_iteration(qs_sheet *qs) {
	if (qs->lanczos.snapshot[0].relation) {
		// the operation is fast, it shouldn't happen in average case.
		// it restores the relations reduced by the linear algebra step that failed.
		qs_sm i ;
		for(i = 0; qs->lanczos.snapshot[i].relation; ++i) {
			qs->relations.data[i] = qs->lanczos.snapshot[i].relation;
			qs->relations.data[i]->Y.length = qs->lanczos.snapshot[i].y_length;
			qs->lanczos.snapshot[i].relation = 0 ;
		}
		assert(qs->relations.length.prev == i) ; // When the matrix is reduced, the original (prev) length is stored.
		DEBUG_NORMAL("Restore the relations previously reduced for the linear algebra to a size of %u.\n", i);
		qs->relations.length.now = i ;
	}
	//  Increase the tick value in each iteration of the algorithm.
	if (++qs->time.tick % 32 == 0) {
		if (qs->relations.length.prev == qs->relations.length.now) {
			// The algorithm notices that no relations accumulates, and reacts to unblock the situation.
			// Common issues include insufficient smooth relations, incorrect bounds for the factor base.
			DEBUG_NORMAL("Relation counter remained %u, no relation accumulates, D is randomized.\n", qs->relations.length.now);
			cint_random_bits(&qs->poly.D, qs->poly.d_bits, &qs->seed);
			*qs->poly.D.mem |= 1; // Shouldn't happen, D becomes a randomized odd number.
		}
		qs->relations.length.prev = qs->relations.length.now;
	}
}

void iteration_part_1(qs_sheet * qs, const cint * D, cint * A) {
	qs_sm n_tries = 0 ; // several attempts may rarely be necessary.
	retry:;
	// A is a "random" product of "s" distinct prime numbers from the factor base.
	cint * X = qs->vars.TEMP, * Y = X + 1, *TMP, *_A = A ;
	qs_sm a, b, i = 0, j;
	if (qs->s.values.defined & 1) TMP = A, A = X, X = TMP ;
	// swap pointers so the last multiplication completes inside the A variable.
	simple_int_to_cint(A, 1);
	for (a = 0, b = qs->poly.span.x_3; a < qs->s.values.subtract_one; ++a) {
		if (a & 1) i = b / (i + qs->poly.span.x_1) - (qs_sm) xor_rand(&qs->seed, 0, 9);
		else i = qs->poly.span.x_2 + (qs_sm) xor_rand(&qs->seed, 0, qs->poly.span.half);
		for (j = 0; j < a; j = i == qs->s.data[j].prime_index ? ++i, 0 : j + 1);
		qs->s.data[a].prime_index = i; // the selected divisor of A wasn't already present in the product.
		simple_int_to_cint(Y, qs->base.data[i].num);
		cint_mul(A, Y, X), TMP = A, A = X, X = TMP;
	}
	// a prime number from the factor base completes A, which must be close to D.
	cint_div(qs->sheet, D, A, X, Y);
	const qs_sm d_over_a = (qs_sm) simple_cint_to_int(X);
	for (i = qs->base.data[0].num != 2 ; qs->base.data[i].num <= d_over_a; ++i);
	for (j = 0; j < qs->s.values.subtract_one; j = i == qs->s.data[j].prime_index ? ++i, 0 : j + 1);
	if (qs->base.length <= i) {
		DEBUG_CRITICAL("Quadratic Sieve inconsistency, %u isn't less than the factor base size.\n", i);
		assert(++n_tries <= 9);
		A = _A ; // shouldn't happen.
		goto retry;
	}
	qs->s.data[qs->s.values.subtract_one].prime_index = i ;
	simple_int_to_cint(Y, qs->base.data[i].num);
	cint_mul(A, Y, X); // generated A values should always be distinct, "A" no longer change.
	assert(X == &qs->poly.A);
}

void iteration_part_2(qs_sheet * qs, const cint * A, cint * B) {
	cint *X = qs->vars.TEMP, *PRIME = X + 1, *Y = X + 2, *R = X + 3;
	qs_sm i, *pen = qs->s.A_indexes;
	cint_erase(B);
	for (i = 0; i < qs->s.values.defined; ++i) {
		const qs_sm idx = qs->s.data[i].prime_index, prime = qs->base.data[idx].num;
		if (idx >= qs->iterative_list[3])
			qs->iterative_list[3] = 8 + idx - idx % 8 ;
		// write [index of prime number, power] of the A factors into buffer.
		*pen++ = qs->s.data[i].prime_index, *pen++ = 1;
		qs->s.data[i].prime_squared = (qs_md)prime * (qs_md)prime ;
		simple_int_to_cint(PRIME, prime);
		cint_div(qs->sheet, A, PRIME, X, R), assert(R->mem == R->end); // div exact.
		cint_div(qs->sheet, X, PRIME, Y, R);
		qs->s.data[i].A_over_prime_mod_prime = (qs_sm) simple_cint_to_int(R);
		qs_md x = modular_inverse(qs->s.data[i].A_over_prime_mod_prime, prime);
		x = x * qs->base.data[qs->s.data[i].prime_index].sqrt_kN_mod_prime % prime;
		simple_int_to_cint(X, x > prime >> 1 ? prime - x : x);
		cint_mul(A, X, Y);
		cint_div(qs->sheet, Y, PRIME, &qs->s.data[i].B_terms, R), assert(R->mem == R->end); // div exact.
		cint_addi(B, &qs->s.data[i].B_terms);
	}
}

void iteration_part_3(qs_sheet * qs, const cint * A, const cint * B) {
	cint *Q = qs->vars.TEMP, *R = Q + 1, *PRIME = Q + 2;
	qs_md i, j, x, y;
	for (i = 0; i < qs->base.length; ++i) {
		// prepare the "roots" and "A_inv_double_value_B_terms". The algorithm will
		// fill 2 ** (s - 3) sieves by using these values and adding "prime sizes".
		const qs_sm prime = qs->base.data[i].num;
		simple_int_to_cint(PRIME, prime);
		cint_div(qs->sheet, A, PRIME, Q, R);
		const qs_sm a_mod_prime = (qs_sm) simple_cint_to_int(R) ;
		cint_div(qs->sheet, B, PRIME, Q, R) ;
		const qs_sm b_mod_prime = (qs_sm) simple_cint_to_int(R) ;
		const qs_sm a_inv_double_value = modular_inverse(a_mod_prime, prime) << 1 ;
		// Arithmetic shifts "<<" and ">>" performs multiplication or division by powers of two.
		x = y = prime;
		x += qs->base.data[i].sqrt_kN_mod_prime;
		y -= qs->base.data[i].sqrt_kN_mod_prime;
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
			// compute the roots update value for all "s".
			cint_div(qs->sheet, &qs->s.data[j].B_terms, PRIME, Q, R);
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
		cint_div(qs->sheet, B, P, Q, R_B);
		cint_div(qs->sheet, kN, P, Q, R_kN);
		if (B->nat < 0) cint_addi(R_B, P); // if B is negative.
		const qs_tmp rem_b = (qs_tmp) simple_cint_to_int(R_B);
		const qs_tmp rem_kn = (qs_tmp) simple_cint_to_int(R_kN);
		qs_tmp s ; // both remainders are modulo the prime number squared.
		if (rem_b < 0xb504f334) {
			// the multiplication is straightforward.
			s = rem_b * rem_b - rem_kn;
			s /= prime ;
		} else {
			// the common multiplication would overflow.
			cint_mul(R_B, R_B, TMP);
			cint_subi(TMP, R_kN);
			simple_int_to_cint(P, (qs_md) prime);
			cint_div(qs->sheet, TMP, P, Q, R_B);
			s = (qs_tmp) simple_cint_to_int(Q);
			if (Q->nat < 0) s = -s ;
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
	cint_div(qs->sheet, TMP, A, C, R), assert(R->mem == R->end); // div exact.
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
	for(qs_sm i = qs->iterative_list[4], j = qs->base.length; i < j; ++i){
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
			const qs_tmp diff = *p_1 - *p_0 ;
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

int qs_register_divisor(qs_sheet *qs) {
	// Register a new divisor of N, and update N when a prime factor is discovered.
	// Return -1 when the factorization is completed.
	int res = 0, pow;
	cint *F = &qs->vars.FACTOR, *Q = qs->state->session.tmp + 3, *R = Q + 1;
	if (h_cint_compare(F, &qs->constants.ONE) > 0 && h_cint_compare(F, &qs->vars.N) < 0) {
		F->nat = 1; // Absolute value ensure this number is greater than 1 and lower than N.
		struct avl_node * node ;
		next:
		if (node = avl_at(&qs->uniqueness[2], F), qs->uniqueness[2].affected) {
			// This new divisor of N wasn't present in the AVL tree.
			if (++res == 1)
				DEBUG_MATH("- New divisor %s shown.\n", simple_cint_string(qs->state, F));
			else
				DEBUG_MATH("  - So let's verify %s.\n", simple_cint_string(qs->state, F));
			if (cint_is_prime(qs->sheet, F, -1)) {
				pow = (int) cint_remove(qs->sheet, &qs->vars.N, F);
				assert(pow);
				manager_add_factor(qs->state, F, pow, 1);
				++qs->divisors.total_primes;
				qs->n_bits = (qs_sm) cint_count_bits(&qs->vars.N);
				if (qs->n_bits == 1) {
					DEBUG_MATH("%s", "  - The factorization is completed.\n");
					res = -1;
				} else {
					cint_dup(F, &qs->vars.N);
					// 200-bit RSA take about 10,000,000+ "duplications".
					DEBUG_MATH("  - This prime number reduces N to %u-bit.\n", qs->n_bits);
					goto next;
				}
			} else if ((pow = factorization_any_root_checker(qs->state, F, Q, R))) {
				cint_dup(F, Q);
				DEBUG_MATH("  - Express it using a power of %d.\n", pow);
				goto next;
			} else if (cint_div(qs->sheet, &qs->vars.N, F, Q, R), cint_is_prime(qs->sheet, Q, -1)) {
				cint_dup(F, Q);
				qs->divisors.data[qs->divisors.length++] = node->key;
				DEBUG_MATH("%s", "  - Just divide N by this number.\n");
				goto next;
			}  else {
				for (qs_sm i = 0; i < qs->divisors.length; ++i) {
					cint_gcd(qs->sheet, F, qs->divisors.data[i], Q);
					if (cint_is_prime(qs->sheet, Q, -1)) {
						cint_dup(F, Q);
						DEBUG_MATH("  - Perform GCD with the known divisor%s.\n", 1 < qs->divisors.length ? "s" : "");
						goto next;
					}
				}
				// Simply register a known divisor of N.
				qs->divisors.data[qs->divisors.length++] = node->key;
				DEBUG_MATH("%s", "  - Store this divisor for further use.\n");
			}
		}
	}
	return res;
}

void register_relations(qs_sheet * qs, const cint * A, const cint * B, const cint * C) {
	cint *  TMP = qs->vars.TEMP, * K = &qs->vars.KEY, * V = &qs->vars.VALUE ;
	qs_sm m_idx, idx, mod;
	// iterates the values of X in [-M/2, +M/2].
	for (m_idx = 0; m_idx < qs->m.length; ++m_idx)
		if (qs->m.sieve[m_idx] >= qs->threshold.value) {
			// over the threshold, compute f(X) and check candidate for smoothness.
			simple_int_to_cint(&qs->vars.X, m_idx);
			cint_subi(&qs->vars.X, &qs->constants.M_HALF); // X = "current index" - M/2
			cint_mul(A, &qs->vars.X, TMP); // TMP = AX
			cint_addi(TMP, B); // TMP = AX + B
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
				*pen = (qs_sm) cint_remove(qs->sheet, V, TMP);
				if (*pen) removed_bits += *pen++ * qs->base.data[0].size; else --pen;
			}
			for (idx = 1; idx < qs->iterative_list[1]; ++idx)
				if (qs->base.data[idx].root[1] == (qs_sm) -1 || (mod = m_idx % qs->base.data[idx].num, mod == qs->base.data[idx].root[0] || mod == qs->base.data[idx].root[1])) {
					simple_int_to_cint(TMP, qs->base.data[idx].num);
					// for a given prime number of the factor base, "remove" returns
					// the numbers of powers that was present in V, and V is updated.
					*pen++ = idx;
					*pen = (qs_sm) cint_remove(qs->sheet, V, TMP);
					if (*pen) removed_bits += *pen++ * qs->base.data[idx].size; else --pen;
				}
			if (removed_bits + qs->m.sieve[m_idx] >= target_bits) {
				// there is a chance to register a new relation.
				for (removed_bits = 0, target_bits = qs->m.sieve[m_idx]; idx < qs->iterative_list[4] && removed_bits < target_bits; ++idx)
					if (qs->base.data[idx].root[1] == (qs_sm) -1 || (mod = m_idx % qs->base.data[idx].num, mod == qs->base.data[idx].root[0] || mod == qs->base.data[idx].root[1])) {
						simple_int_to_cint(TMP, qs->base.data[idx].num);
						*pen++ = idx;
						*pen = (qs_sm) cint_remove(qs->sheet, V, TMP);
						if (*pen) removed_bits += *pen++ * qs->base.data[idx].size; else --pen;
					}
				for (const uint8_t mask = 1 << (m_idx & 7); idx < qs->base.length && removed_bits < target_bits; ++idx)
					if (qs->m.flags[idx] & mask)
						if (mod = m_idx % qs->base.data[idx].num, mod == qs->base.data[idx].root[0] || mod == qs->base.data[idx].root[1]) {
							simple_int_to_cint(TMP, qs->base.data[idx].num);
							*pen++ = idx;
							*pen = (qs_sm) cint_remove(qs->sheet, V, TMP);
							if (*pen) removed_bits += *pen++ * qs->base.data[idx].size; else --pen;
						}
				const qs_sm * restrict const prime_indexes_and_powers[4] = {
						qs->s.A_indexes, // really factoring A * f(X), commit outstanding A factors.
						qs->s.A_indexes + qs->s.values.double_value,
						qs->buffer[1],
						pen,
				};
				if (h_cint_compare(V, &qs->constants.ONE) == 0)
					register_regular_relation(qs, K, prime_indexes_and_powers);
				else if (154 < qs->kn_bits && h_cint_compare(V, &qs->constants.TOO_LARGE_PRIME) < 0)
					//  Store it until another partial share the same large prime.
					register_partial_relation(qs, K, V, prime_indexes_and_powers);
			}
		}
}

void register_regular_relation(qs_sheet * qs, const cint * KEY, const qs_sm * const restrict args[4]) {
	struct avl_node *node = avl_at(&qs->uniqueness[0], KEY);
	if (node->value)
		return; // duplicates at this stage are ignored.
	struct qs_relation * rel = qs->mem.now;
	qs_sm i, j ;
	qs->mem.now = rel + 1 ; // a relation must be swappable for Lanczos Block reducing.
	rel->X = node->key; // constant X is stored by the node key.
	rel->Y.data = qs->mem.now; // Y data length only decreases.
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
			cint_mul_modi(qs->sheet, A, B, &qs->constants.kN);
		}
		cint_mul_mod(qs->sheet, rel->X, rel->X, &qs->constants.kN, B);
		verified = !cint_compare(A, B) || (cint_addi(A, B), !cint_compare(A, &qs->constants.kN));
	}
	if (verified){
		node->value = qs->relations.data[qs->relations.length.now] = rel;
		qs->mem.now = rel->axis.Z.data + rel->axis.Z.length;
		rel->id = ++qs->relations.length.now; // Keep the relation
	} else {
		DEBUG_NORMAL("Discard the relation at index %u.\n", qs->relations.length.now);
		char * open = (char*) rel, * close = qs->mem.now ;
		qs->mem.now = memset(open, 0, close - open); // Throw
	}
}

void register_partial_relation(qs_sheet * qs, const cint * KEY, const cint * VALUE, const qs_sm * const restrict args[4]) {
	// searches 2 different KEY sharing the same VALUE.
	struct avl_node *node = avl_at(&qs->uniqueness[1], VALUE);
	struct qs_relation *old, *new;
	cint * BEZOUT = 0;
	old = node->value;
	if (old) {
		if (old->X == 0) return; // the value is already marked as "ignored".
		if (old->axis.next) return; // accepting all "next" without caring reduce the chance.
		for (; old && h_cint_compare(KEY, old->X); old = old->axis.next);
		if (old) return; // same KEY already registered.
		old = node->value;
		if (old->axis.next == 0) {
			// there is an "old" using the same VALUE, and it has no "next" yet.
			cint *A = qs->vars.TEMP, *B = A + 1;
			if (qs->multiplier != 1)
				if (cint_gcd(qs->sheet, VALUE, &qs->constants.MULTIPLIER, A), cint_compare_char(A, 1)){
					old->X = 0; // VALUE shouldn't be related so close to the multiplier.
					return;
				}
			// so compute BEZOUT.
			cint_modular_inverse(qs->sheet, VALUE, &qs->constants.kN, A);
			if (A->mem == A->end) {
				old->X = 0; // no solution to the linear congruence.
				cint_gcd(qs->sheet, VALUE, &qs->constants.kN, &qs->vars.FACTOR);
				cint_div(qs->sheet, &qs->vars.N, &qs->vars.FACTOR, A, B);
				if (B->mem == B->end) // found a small factor of N ?
					qs_register_divisor(qs);
				return; // nothing.
			} else
				BEZOUT = A;
		}
	}

	new = mem_aligned(qs->mem.now);
	qs->mem.now = new + 1;
	new->X = qs->mem.now;

	if (BEZOUT) {
		// BEZOUT is stored directly after the new X, like in an array.
		qs->mem.now = new->X + 2;
		simple_dup_cint(new->X, KEY, &qs->mem.now);
		simple_dup_cint(new->X + 1, BEZOUT, &qs->mem.now);
		// The 2nd newcomer become the root of the linked list.
		new->axis.next = old, node->value = new = old = new;
	} else {
		// All but the 2nd have no special treatment.
		qs->mem.now = new->X + 1; // they come at the end of the linked list.
		simple_dup_cint(new->X, KEY, &qs->mem.now);
		if (old) {
			for (; old->axis.next; old = old->axis.next);
			old->axis.next = new, old = node->value;
		} else node->value = new;
	}

	// data buffered isn't persistent, it may be needed, so it's copied.
	qs_sm * data = new->Y.data = mem_aligned(qs->mem.now);
	new->Y.length = (qs_sm) (args[1] - args[0]);
	memcpy(data, args[0], new->Y.length * sizeof(*data));
	memcpy(data + new->Y.length, args[2], (args[3] - args[2]) * sizeof(*data));
	new->Y.length += (qs_sm) (args[3] - args[2]);
	qs->mem.now = new->Y.data + new->Y.length;

	if (old) {
		BEZOUT = old->X + 1 ; // the modular inverse was stored here.
		cint_mul_mod(qs->sheet, new->X, BEZOUT, &qs->constants.kN, &qs->vars.CYCLE);
		do {
			if (old != new) {
				// combines, it registers a regular relation using the 2 buffers.
				cint_mul_mod(qs->sheet, &qs->vars.CYCLE, old->X, &qs->constants.kN, &qs->vars.KEY);
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
				register_regular_relation(qs, &qs->vars.KEY, args);
				++qs->relations.length.by_partial;
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
	cint *X = qs->vars.TEMP, *Y = X + 1, *TMP = X + 2, *POW = X + 3;
	qs_sm * restrict power_of_primes;
	for(qs_sm row = 0; row < 64 && qs->n_bits != 1; ++row)
		if (mask >> row & 1){
			// use the Fermat's (1607 - 1665) method to compute a factorization of N.
			simple_int_to_cint(X, 1), simple_int_to_cint(TMP, 1), simple_int_to_cint(Y, 1);
			power_of_primes = memset(qs->buffer[1], 0, qs->base.length * sizeof(*power_of_primes));
			for (qs_sm i = 0; i < qs->relations.length.now; ++i)
				if (null_rows[i] >> row & 1) {
					const struct qs_relation * restrict const rel = qs->relations.data[i];
					// The algorithm must retrieve the "X" and "Z" relation fields
					// related to the "Y" field initially submitted to Lanczos Block.
					cint_mul_modi(qs->sheet, X, rel->X, &qs->vars.N);
					for (qs_sm j = 0; j < rel->axis.Z.length; ++j)
						++power_of_primes[rel->axis.Z.data[j]];
				}
			for (qs_sm i = 0; i < qs->base.length; ++i)
				if (power_of_primes[i]){
					// powers are even ... square root ...
					simple_int_to_cint(TMP, qs->base.data[i].num);
					simple_int_to_cint(POW, power_of_primes[i] >> 1);
					cint_pow_modi(qs->sheet, TMP, POW, &qs->vars.N);
					cint_mul_modi(qs->sheet, Y, TMP, &qs->vars.N);
				}
			h_cint_subi(Y, X);
			if (Y->mem != Y->end) {
				cint_gcd(qs->sheet, &qs->vars.N, Y, &qs->vars.FACTOR);
				// 100 digits RSA number has been factored by the software in 2022.
				if (qs_register_divisor(qs) == -1)
					break;
			}
		}
}

int finalization_part_2(qs_sheet *qs) {
	if (qs->n_bits != 1 && qs->divisors.length) {
		cint *F = &qs->vars.FACTOR;
		// N isn't fully factored, register divisor(s) of N from the smallest to the caller's routine.
		qsort(qs->divisors.data, qs->divisors.length, sizeof(*qs->divisors.data), (int (*)(const void *, const void *)) h_cint_compare);
		for (qs_sm i = 0; i < qs->divisors.length; ++i)
			if (h_cint_compare(&qs->vars.N, qs->divisors.data[i]) > 0) {
				if (qs->state->params.verbose)
					DEBUG_CRITICAL("Quadratic Sieve identified the divisor %s.\n", simple_cint_string(qs->state, qs->divisors.data[i]));
				const int power = (int) cint_remove(qs->sheet, &qs->vars.N, qs->divisors.data[i]);
				if (power) {
					manager_add_factor(qs->state, F, power, -1);
					qs->n_bits = (qs_sm) cint_count_bits(&qs->vars.N);
				}
			}
		DEBUG_NORMAL("Quadratic Sieve reduced N to %s.\n", simple_cint_string(qs->state, &qs->vars.N));
	}

	DEBUG_NORMAL("The sieve found %u smooth relations using %u polynomials.\n", qs->relations.length.peak, qs->poly.curves);
	if (qs->uniqueness[1].count)
		DEBUG_NORMAL("%u partials were used to add %u smooth relations.\n", (unsigned)qs->uniqueness[1].count, qs->relations.length.by_partial);
	else
		DEBUG_NORMAL("%s", "No partial relation were used.\n");
	DEBUG_NORMAL("Used %u MB of memory during %.02f second(s).\n", (unsigned)((char*)qs->mem.now - (char*) qs->mem.base) >> 20, 0.001 * (get_time_ms() - qs->time.start));

	const int res = cint_compare(&qs->state->session.num, &qs->vars.N) != 0;
	if (res) // Update the input number present in the state (we worked on a copy)
		cint_dup(&qs->state->session.num, &qs->vars.N);
	return res;
}

// compiling with "gcc -Wall -pedantic -O3 -std=c99 main.c -o factor" can speed up by a factor 2 or 3.
