// Front-End Factor manager
static inline fac_cint **c_factor(const cint *N, fac_params *config) {
	void *mem;
	fac_caller m = {0};
	const int input_bits = (int) cint_count_bits(N);
	int alloc_bits = (1 + (input_bits >> 9)) << 13;
	// initially allocates 8 Kb for each 512-bit chunk.
	mem = m.mem.base = calloc(1, alloc_bits);
	assert(mem);
	add_rand_seed(&mem);
	if (config) m.params = config;
	else m.params = mem, mem = m.params + 1;
	m.calc = cint_new_sheet((1 + (input_bits >> 10)) << 10);
	assert(m.calc);

	// prepare 5 variables.
	const size_t vars_size = 500 * (1 + input_bits / 500) / cint_exponent;
	for (int i = 0; i < 5; ++i)
		simple_inline_cint(&m.vars[i], vars_size, &mem);

	simple_inline_cint(&m.trial.cint, vars_size, &mem);

	// prepare a working array.
	const int max_factors = input_bits / 10 + 32;
	m.questions.data = mem, mem = m.questions.data + max_factors ;
	m.answers.data = mem, mem = m.answers.data + max_factors ;

	m.number = &m.questions.data[m.questions.index++] ;
	simple_dup_cint(&m.number->cint, N, &mem);
	m.number->bits = input_bits ;
	m.number->power = 1, m.number->prime = -1 ;

	m.mem.now = mem ;
	// iterates the array until it's empty, begin with the input N.
	// functions must not push their input to the stack, they return 0 instead.
	do {
		m.number = &m.questions.data[--m.questions.index];
		int res =   fac_special_cases(&m)
		            || fac_trial_division(&m, 1)
		            || fac_perfect_checker(&m)
		            || fac_primality_checker(&m)
		            || fac_pollard_rho_64_bits(&m)
		            || quadratic_sieve(&m)
		            || fac_trial_division(&m, 2);
		if (res == 0)
			fac_push(&m, &m.number->cint, 0, 1, 0);
	} while (m.questions.index);



	// answer goes into an appropriately sized memory allocation.
	size_t bytes = 0 ;
	for(unsigned i = 0; i < m.answers.index; ++i)
		bytes += m.answers.data[i].cint.end - m.answers.data[i].cint.mem + 1 ;
	bytes *= sizeof(h_cint_t);
	bytes += (sizeof(fac_cint) + sizeof(fac_cint*)) * (m.answers.index + 1) ;
	fac_cint  ** res = mem = calloc(1, bytes);
	assert(mem);

	qsort(m.answers.data, m.answers.index, sizeof(fac_cint), &fac_sort_result);

	mem = res + m.answers.index + 1 ;
	for(unsigned i = 0; i < m.answers.index; ++i) {
		fac_cint * factor = &m.answers.data[i] ;
		res[i] = mem, mem = res[i] + 1 ;
		res[i]->power = factor->power ;
		res[i]->prime = factor->prime ;
		res[i]->bits = (int) cint_count_bits(&factor->cint) ;
		simple_inline_cint(&res[i]->cint, factor->cint.size, &mem);
		cint_dup(&res[i]->cint, &factor->cint);
	}
	free(m.mem.base);
	cint_clear_sheet(m.calc);
	return res;
}

static inline int fac_special_cases(fac_caller *m) {
	int res = m->number->bits < 3 ;
	if (res && m->answers.index == 0) {
		const int prime = m->number->bits > 1;
		fac_push(m, &m->number->cint, prime, 1, 0);
	}
	return res ;
}

static inline int fac_trial_division(fac_caller *m, const int level) {
	cint * F = m->vars ;
	int res = (*m->number->cint.mem & 1) == 0 ; // remove power of 2.
	if (m->trial.done_up_to == 0){
		if (res) {
			simple_int_to_cint(F, 2);
			const int power = (int) cint_count_zeros(&m->number->cint);
			fac_push(m, F, 1, power, 0);
			cint_right_shifti(&m->number->cint, power);
		}
		m->trial.done_up_to = 1 ;
	}

	unsigned bound;
	if (m->number->bits <= 64) bound = 1024;
	else if (bound = 4669921, level == 1)
		for (int i = 0; i < 250; bound >>= (m->number->bits < i) * 1, i += 30);

	for (; (m->trial.done_up_to += 2) < bound;) {
		if (is_prime_4669921(m->trial.done_up_to)) {
			simple_int_to_cint(F, m->trial.done_up_to);
			const unsigned power = cint_remove(m->calc, &m->number->cint, F);
			if (power){
				fac_push(m, F, 1, (int) power, 0);
				++res;
			}
		}
	}
	simple_int_to_cint(&m->trial.cint, m->trial.done_up_to);
	if (res) fac_push(m, &m->number->cint, -1, 1, 1);
	return res ;
}

static inline int fac_any_root_check(fac_caller * m, const cint *N, cint *Q, cint *R){
	// Can normally say if a number is a perfect power.
	// It takes in account the trial divisions initially done.
	int res = 0 ;
	const int max_root = 30 ;
	for(int nth = 2; nth < max_root ; ++nth)
		if (is_prime_4669921(nth)) {
			cint_nth_root_remainder(m->calc, N, nth, Q, R);
			if (R->mem == R->end){
				res = nth ;
				break;
			}
			if (h_cint_compare(Q, &m->trial.cint) <= 0)
				break;
		}
	return res ;
}

static inline int fac_perfect_checker(fac_caller *m) {
	assert(m->number->bits > 2);
	cint *Q = m->vars, *R = Q + 1;
	int power = fac_any_root_check(m, &m->number->cint, Q, R);
	if (power)
		fac_push(m, Q, -1, power, 1);
	return power;
}

static inline int fac_primality_checker(fac_caller *m) {
	m->number->prime = cint_is_prime(m->calc, &m->number->cint, m->number->bits > 2048 ? 1 : -1);
	if (m->number->prime)
		fac_push(m, &m->number->cint, 1, 1, 0);
	return m->number->prime;
}

static inline int fac_pollard_rho_64_bits(fac_caller *m) {
	int res = m->number->bits > 0 && m->number->bits < 65;
	if (res) {
		// Perform a Pollard's Rho test, this function can't complete with a prime number.
		qs_md n[2] = {simple_cint_to_int(&m->number->cint), 1,}; // number and its factor.
		for (size_t limit = 7; n[1] == 1 || n[0] == n[1]; ++limit) {
			if (m->params->silent == 0)
				fac_display_progress("Pollard Rho", 100. * (double) limit / 21);
			size_t a = -1, b = 2, c = limit;
			qs_md d, e = rand_64(), f = 1;
			for (n[1] = 1, d = e %= n[0]; n[1] == 1 && --c; e = d, b <<= 1, a = -1) {
				for (; n[1] |= f, n[1] == 1 && ++a != b;) {
					d = multiplication_modulo(d, d, n[0]);
					for (++d, d *= d != n[0], f = n[0], n[1] = d > e ? d - e : e - d; (n[1] %= f) && (f %= n[1]););
				}
			}
		}
		n[0] /= n[1];
		cint * F = m->vars ;
		for (int i = 0; i < 2; ++i) {
			simple_int_to_cint(F, n[i]);
			fac_push(m, F, -1, 1, 1);
		}
	}
	return res;
}

// functions submit factors of N, they don't push N itself with "forward".
static void fac_push(fac_caller *m, const cint * num, const int prime, const int power, const int forward) {
	// the product of "stack last" and "stack next" must remain N.
	fac_cint * row ;
	if (forward){
		row = &m->questions.data[m->questions.index++];
		const size_t needed_size = num->end - num->mem + 1;
		if (row->cint.size < needed_size)
			simple_inline_cint(&row->cint, needed_size, &m->mem.now);
		row->bits = (int) cint_count_bits(num);
	} else {
		row = &m->answers.data[m->answers.index++] ;
		simple_inline_cint(&row->cint, num->end - num->mem + 1, &m->mem.now);
	}
	cint_dup(&row->cint, num);
	row->prime = prime;
	row->power = power * m->number->power;
	assert(row->power);
}

// Math
static inline int is_prime_4669921(const qs_sm n) {
// Used to iterate over primes, there are 326,983 prime numbers under 4,669,921, very fast.
	return ((n > 1) & ((n < 6) * 42 + 0x208A2882) >> n % 30 && (n < 49 || (n % 7 && n % 11 && n % 13 && n % 17 && n % 19 && n % 23 && n % 29 && (n < 961 || (n % 31 && n % 37 && n % 41 && n % 43 && n % 47 && n % 53 && n % 59 && n % 61 && n % 67 && (n < 5041 || (n % 71 && n % 73 && n % 79 && n % 83 && n % 89 && n % 97 && n % 101 && n % 103 && n % 107 && (n < 11881 || (n % 109 && n % 113 && n % 127 && n % 131 && n % 137 && n % 139 && n % 149 && n % 151 && n % 157 && (n < 26569 || (n % 163 && n % 167 && n % 173 && n % 179 && n % 181 && n % 191 && n % 193 && n % 197 && n % 199 && (n < 44521 || (n % 211 && n % 223 && n % 227 && n % 229 && n % 233 && n % 239 && n % 241 && n % 251 && n % 257 && (n < 69169 || (n % 263 && n % 269 && n % 271 && n % 277 && n % 281 && n % 283 && n % 293 && n % 307 && n % 311 && (n < 97969 || (n % 313 && n % 317 && n % 331 && n % 337 && n % 347 && n % 349 && n % 353 && n % 359 && n % 367 && (n < 139129 || (n % 373 && n % 379 && n % 383 && n % 389 && n % 397 && n % 401 && n % 409 && n % 419 && n % 421 && (n < 185761 || (n % 431 && n % 433 && n % 439 && n % 443 && n % 449 && n % 457 && n % 461 && n % 463 && n % 467 && (n < 229441 || (n % 479 && n % 487 && n % 491 && n % 499 && n % 503 && n % 509 && n % 521 && n % 523 && n % 541 && (n < 299209 || (n % 547 && n % 557 && n % 563 && n % 569 && n % 571 && n % 577 && n % 587 && n % 593 && n % 599 && (n < 361201 || (n % 601 && n % 607 && n % 613 && n % 617 && n % 619 && n % 631 && n % 641 && n % 643 && n % 647 && (n < 426409 || (n % 653 && n % 659 && n % 661 && n % 673 && n % 677 && n % 683 && n % 691 && n % 701 && n % 709 && (n < 516961 || (n % 719 && n % 727 && n % 733 && n % 739 && n % 743 && n % 751 && n % 757 && n % 761 && n % 769 && (n < 597529 || (n % 773 && n % 787 && n % 797 && n % 809 && n % 811 && n % 821 && n % 823 && n % 827 && n % 829 && (n < 703921 || (n % 839 && n % 853 && n % 857 && n % 859 && n % 863 && n % 877 && n % 881 && n % 883 && n % 887 && (n < 822649 || (n % 907 && n % 911 && n % 919 && n % 929 && n % 937 && n % 941 && n % 947 && n % 953 && n % 967 && (n < 942841 || (n % 971 && n % 977 && n % 983 && n % 991 && n % 997 && n % 1009 && n % 1013 && n % 1019 && n % 1021 && (n < 1062961 || (n % 1031 && n % 1033 && n % 1039 && n % 1049 && n % 1051 && n % 1061 && n % 1063 && n % 1069 && n % 1087 && (n < 1190281 || (n % 1091 && n % 1093 && n % 1097 && n % 1103 && n % 1109 && n % 1117 && n % 1123 && n % 1129 && n % 1151 && (n < 1329409 || (n % 1153 && n % 1163 && n % 1171 && n % 1181 && n % 1187 && n % 1193 && n % 1201 && n % 1213 && n % 1217 && (n < 1495729 || (n % 1223 && n % 1229 && n % 1231 && n % 1237 && n % 1249 && n % 1259 && n % 1277 && n % 1279 && n % 1283 && (n < 1661521 || (n % 1289 && n % 1291 && n % 1297 && n % 1301 && n % 1303 && n % 1307 && n % 1319 && n % 1321 && n % 1327 && (n < 1852321 || (n % 1361 && n % 1367 && n % 1373 && n % 1381 && n % 1399 && n % 1409 && n % 1423 && n % 1427 && n % 1429 && (n < 2053489 || (n % 1433 && n % 1439 && n % 1447 && n % 1451 && n % 1453 && n % 1459 && n % 1471 && n % 1481 && n % 1483 && (n < 2211169 || (n % 1487 && n % 1489 && n % 1493 && n % 1499 && n % 1511 && n % 1523 && n % 1531 && n % 1543 && n % 1549 && (n < 2411809 || (n % 1553 && n % 1559 && n % 1567 && n % 1571 && n % 1579 && n % 1583 && n % 1597 && n % 1601 && n % 1607 && (n < 2588881 || (n % 1609 && n % 1613 && n % 1619 && n % 1621 && n % 1627 && n % 1637 && n % 1657 && n % 1663 && n % 1667 && (n < 2785561 || (n % 1669 && n % 1693 && n % 1697 && n % 1699 && n % 1709 && n % 1721 && n % 1723 && n % 1733 && n % 1741 && (n < 3052009 || (n % 1747 && n % 1753 && n % 1759 && n % 1777 && n % 1783 && n % 1787 && n % 1789 && n % 1801 && n % 1811 && (n < 3323329 || (n % 1823 && n % 1831 && n % 1847 && n % 1861 && n % 1867 && n % 1871 && n % 1873 && n % 1877 && n % 1879 && (n < 3568321 || (n % 1889 && n % 1901 && n % 1907 && n % 1913 && n % 1931 && n % 1933 && n % 1949 && n % 1951 && n % 1973 && (n < 3916441 || (n % 1979 && n % 1987 && n % 1993 && n % 1997 && n % 1999 && n % 2003 && n % 2011 && n % 2017 && n % 2027 && (n < 4116841 || (n % 2029 && n % 2039 && n % 2053 && n % 2063 && n % 2069 && n % 2081 && n % 2083 && n % 2087 && n % 2089 && (n < 4405801 || (n % 2099 && n % 2111 && n % 2113 && n % 2129 && n % 2131 && n % 2137 && n % 2141 && n % 2143 && n % 2153)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))));
}

static double log_computation(const double n) {
	// Basic logarithm computation, believing you can't include <math.h>.
	static const double euler = 2.7182818284590452354;
	unsigned a = 0, d;
	double b, c, e, f;
	if (n > 0) {
		for (c = n < 1 ? 1 / n : n; (c /= euler) > 1; ++a);
		c = 1 / (c * euler - 1), c = c + c + 1, f = c * c, b = 0;
		for (d = 1, c /= 2; e = b, b += 1 / (d * c), b - e > 0.00001;) {
			d += 2, c *= f;
		}
	} else { b = (n == 0) / 0.; }
	return n < 1 ? -(a + b) : a + b;
}

static inline qs_md multiplication_modulo(qs_md a, qs_md b, const qs_md mod) {
#ifdef __SIZEOF_INT128__
	return (__uint128_t) a * (__uint128_t) b % (__uint128_t) mod ;
#else
	// Return (a * b) % mod, avoiding overflow errors while doing modular multiplication.
	qs_md res = 0, tmp;
	for (b %= mod; a; a & 1 ? b >= mod - res ? res -= mod : 0, res += b : 0, a >>= 1, (tmp = b) >= mod - b ? tmp -= mod : 0, b += tmp);
	return res % mod;
#endif
}

static inline qs_md power_modulo(qs_md n, qs_md exp, const qs_md mod) {
	// Return (n ^ exp) % mod
	qs_md res = 1;
	for (n %= mod; exp; exp & 1 ? res = multiplication_modulo(res, n, mod) : 0, n = multiplication_modulo(n, n, mod), exp >>= 1);
	return res;
}

static inline int kronecker_symbol(qs_md a, qs_md b) {
	static const int s[8] = {0, 1, 0, -1, 0, -1, 0, 1};
	qs_md c;
	int res = (int) (a | b);
	if (a && b)
		if (res & 1) {
			for (c = 0; !(b & 1); ++c, b >>= 1);
			// When b is odd Jacobi and Kronecker symbols are identical, in factorization algorithms b is often the prime number.
			// When b is an odd prime number, Jacobi symbol is equal to the Legendre symbol.
			for (res = c & 1 ? s[a & 7] : 1; a; c & 1 ? res *= s[b & 7] : 0, a & b & 2 ? res = -res : 0, c = b % a, b = a, a = c)
				for (c = 0; !(a & 1); ++c, a >>= 1);
			res = b == 1 ? res : 0;
		} else res = 0;
	else res = res == 1;
	return res;
}

static qs_md tonelli_shanks(const qs_md n, const unsigned mod) {
	// return root such that (root * root) % mod congruent to n % mod.
	// return 0 if no solution to the congruence exists.
	// mod is assumed odd prime, if mod = 2 then res is (n & 7 == 1 || n & 7 == 7).
	const qs_sm a = n % mod;
	qs_md res = kronecker_symbol(a, mod) == 1, b, c, d, e, f, g, h;
	if (res)
		switch (mod & 7) {
			case 3 : case 7 :
				res = power_modulo(a, (mod + 1) >> 2, mod);
				break;
			case 5 :
				res = power_modulo(a, (mod + 3) >> 3, mod);
				if (multiplication_modulo(res, res, mod) != a) {
					b = power_modulo(2, (mod - 1) >> 2, mod);
					res = multiplication_modulo(res, b, mod);
				}
				break;
			default :
				if (a == 1)
					res = 1;
				else {
					for (c = mod - 1, d = 2; d < mod && power_modulo(d, c >> 1, mod) != c; ++d);
					for (e = 0; !(c & 1); ++e, c >>= 1);
					f = power_modulo(a, c, mod);
					b = power_modulo(d, c, mod);
					for (h = 0, g = 0; h < e; h++) {
						d = power_modulo(b, g, mod);
						d = multiplication_modulo(d, f, mod);
						d = power_modulo(d, 1 << (e - h - 1), mod);
						if (d == mod - 1)
							g += 1 << h;
					}
					f = power_modulo(a, (c + 1) >> 1, mod);
					b = power_modulo(b, g >> 1, mod);
					res = multiplication_modulo(f, b, mod);
				}
		}
	return res;
}

static qs_md modular_inverse(qs_md ra, qs_md rb) {
	// Return a modular multiplicative inverse of n with respect to the modulus.
	// Return 0 if the linear congruence has no solutions.
	// The answer is also called "u1" in the context of extended Euclidean algorithm.
	qs_md rc, sa = 1, sb = 0, sc, i = 0;
	if (rb > 1)
		do {
			rc = ra % rb;
			sc = sa - (ra / rb) * sb;
			sa = sb, sb = sc;
			ra = rb, rb = rc;
		} while (++i, rc);
	sa *= (i *= ra == 1) != 0;
	sa += (i & 1) * sb;
	return sa;
}

static inline qs_md rand_64() {
	// Last bit of rand is very random.
	qs_md res = 0;
	for (qs_sm i = 65; --i; res <<= 1, res |= rand() & 1);
	return res;
}

static inline qs_md rand_upto(const qs_md limit) {
	return rand_64() % limit;
}

static inline unsigned add_rand_seed(void *addr) {
	// Take addresses of memory addresses as argument.
	// Addresses are used together to seed the rand, not too often.
	static unsigned seed;
	if (addr) {
		seed ^= *(unsigned *) addr + (unsigned) (uintptr_t) &errno;
		seed = power_modulo(seed + 1, seed - 3, -5);
		srand(seed);
	}
	return seed;
}

// Cint shortcuts
static inline void simple_inline_cint(cint *N, const size_t size, void **mem) {
	// Fixed size cint is inlined , mem is updated accordingly.
	N->mem = N->end = (h_cint_t *) *mem;
	*mem = N->mem + (N->size = size);
}

static inline void simple_dup_cint(cint *A, const cint *B, void **mem) {
	// Duplicates cint using the given memory, which is updated accordingly.
	// It uses the minimal size, the duplicate is not resizable.
	A->mem = A->end = (h_cint_t *) *mem;
	cint_dup(A, B);
	A->size = A->end - A->mem + 1;
	*mem = A->mem + A->size;
}

static inline void simple_int_to_cint(cint *num, qs_md value) {
	// Pass the given 64-bit number into the given cint.
	for (cint_erase(num); value; *num->end++ = (h_cint_t) (value & cint_mask), value >>= cint_exponent);
}

static inline qs_md simple_cint_to_int(const cint *num) {
	// Return the value of a cint as a 64-bit integer.
	qs_md res = 0;
	for (h_cint_t *ptr = num->end; ptr > num->mem; res = (qs_md) (res * cint_base + *--ptr));
	return res;
}

// Avl
static inline struct avl_node *avl_cint_inserter(void *args, const void *key_to_save) {
	// it expects as result a new node containing the given constant key.
	void * mem = *(void**) args ;
	struct avl_node *res = mem;
	res->key = (cint *) (res + 1);
	mem = (cint *) res->key + 1;
	simple_dup_cint(res->key, key_to_save, &mem);
	assert(res->value == 0);
	*(void**)args = mem ;
	return res;
}

// System
static inline void *mem_aligned(void *ptr) {
	// Default alignment of the return value is 64.
	char *res __attribute__((aligned(64)));
	res = (char *) ptr + (64 - (uintptr_t) ptr) % 64;
	return res;
}

// Misc.
static inline int fac_apply_custom_param(const char *a, const char *b, int length, unsigned *val) {
	int res = memcmp(a, b, length) == 0;
	if (res) {
		for (; *b && !(*b >= '1' && *b <= '9'); ++b);
		for (*val = 0; *b && !(*val >> 26); ++b)
			*val = *val * 10 + *b - '0';
		if (*val == 0) *val = 1 ;
	}
	return res;
}

static inline char *fac_fill_params(fac_params *params, int argc, char **args) {
	char *n = 0;
	for (int i = 1; i < argc; ++i) {
		char *s = args[i];
		for (; *s && !(*s >= '1' && *s <= '9') && !(*s >= 'a' && *s <= 'z'); ++s);
		if (*s >= 'a' && *s <= 'z') {
			int a =
					fac_apply_custom_param("qs_limit=", s, 1, &params->qs_limit) // add command line parameters...
					|| fac_apply_custom_param("testing=", s, 1, &params->testing)
					|| fac_apply_custom_param("silent=", s, 1, &params->silent)
					|| fac_apply_custom_param("multiplier=", s, 1, &params->qs_multiplier)
					|| fac_apply_custom_param("rand=", s, 1, &params->qs_rand_seed)
					|| fac_apply_custom_param("help=", s, 1, &params->help);
			assert(a >= 0);
		} else if (n == 0) {
			for (n = s; *n >= '0' && *n <= '9'; ++n);
			n = *n == 0 ? s : 0;
		}
	}
	return n;
}

static char *fac_answer_to_string(fac_cint **ans) {
	// Should return a string to represent the given answer.
	size_t bytes  = 0;
	int i, j ;
	for(i = 0; ans[i]; ++i){
		bytes += (size_t)(1. + 0.30102999566 * ans[i]->bits) ;
		if (ans[i]->power > 1){
			for(j = 2; ans[i]->power >> j; ++j);
			bytes += (size_t) (6. + 0.30102999566 * j);
		}
		bytes += (ans[i]->prime < 1) << 2 ;
	}
	bytes += 3 * i - 2 ;
	char * str, *res = str = malloc(bytes);
	assert(res);
	for(*str = 0, i = 0; ans[i]; ++i){
		char * s = cint_to_string(&ans[i]->cint, 10);
		if (ans[i]->power > 1) strcat(str++, "(");
		if (ans[i]->prime < 1) strcat(str, " \""), str += 2;
		str += sprintf(str, "%s", s) ;
		if (ans[i]->prime < 1) strcat(str, "\" "), str += 2;
		if (ans[i]->power > 1)
			str += sprintf(str, " ^ %d", ans[i]->power);
		if (ans[i]->power > 1) strcat(str++, ")");
		if (ans[i + 1]) str += sprintf(str, " * ");
		free(s);
	}
	assert(str <= res + bytes);
	return res;
}

static inline void fac_display_progress(const char *name, double percentage) {
	// There are functions that print their progress.
	printf("%s at %.02f %%...", name, percentage);
	putchar('\r');
	fflush(stdout);
}

static inline int fac_sort_result(const void * lhs, const void * rhs) {
	const fac_cint * L = lhs, *R = rhs;
	return h_cint_compare(&L->cint, &R->cint);
}
