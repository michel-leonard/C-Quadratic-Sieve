#ifndef CINT_MASTER
#define CINT_MASTER

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>

// memory is supposed provided by the system, allocations are passed to "assert".
// cint use "computation sheets" instead of global vars.

// the functions name that terminates by "i" means immediate, in place.
// the functions name that begin by "h_" means intended for internal usage.

typedef int64_t h_cint_t;

static const h_cint_t cint_exponent = 4 * sizeof(h_cint_t) - 1;
static const h_cint_t cint_base = (h_cint_t)1 << cint_exponent;
static const h_cint_t cint_mask = cint_base - 1;
static const char *cint_alpha = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

typedef struct {
	h_cint_t *mem; // reserved data pointer.
	h_cint_t *end;
	h_cint_t nat; // -1 = negative, +1 = positive, (zero is a positive)
	size_t size;
} cint;

typedef struct {
	size_t immediate_state ;
	cint temp[10];
} cint_sheet;

static cint_sheet * cint_new_sheet(const size_t bits) {
	// a computation sheet is required by function needing temporary vars.
	cint_sheet * sheet = calloc(1, sizeof(cint_sheet));
	assert(sheet);
	const size_t num_size = 2 + bits / cint_exponent;
	if (sheet) for (size_t i = 0; i < 10; ++i) {
		sheet->temp[i].nat = 1 ;
		sheet->temp[i].mem = sheet->temp[i].end = calloc(num_size, sizeof(h_cint_t));
		assert(sheet->temp[i].mem);
		sheet->temp[i].size = num_size;
	}
	return sheet ;
}

__attribute__((unused)) static inline void cint_negate(cint * N){
	N->nat *= 1 - ((N->mem != N->end) << 1) ;
}

void cint_clear_sheet(cint_sheet *sheet) {
	for (size_t i = 0; i < 10; ++i)
		free(sheet->temp[i].mem);
	free(sheet);
}

static size_t cint_count_bits(const cint * num) {
	size_t res = 0;
	if (num->end != num->mem) {
		for (; *(num->end - 1) >> ++res;);
		res += (num->end - num->mem - 1) * cint_exponent;
	}
	return res;
}

static size_t cint_count_zeros(const cint * num){
	// it returns the total of "right shifts" it takes to turn "num" odd.
	size_t res = 0, i;
	h_cint_t * ptr ;
	for (ptr = num->mem; ptr < num->end && !*ptr; ++ptr, res += cint_exponent);
	for (i = 0; !(*ptr >> i & 1); ++i);
	return res + i ;
}

__attribute__((unused)) static inline int cint_compare_char(const cint * N, const h_cint_t val){
	const h_cint_t res = *N->mem + *(N->mem + 1) - val ;
	return (res > 0) - (res < 0);
}

static inline int h_cint_compare(const cint * lhs, const cint * rhs) {
	h_cint_t res = (h_cint_t) ((lhs->end - lhs->mem) - (rhs->end - rhs->mem));
	if (res == 0 && rhs->end != rhs->mem)
		for (const h_cint_t *l = lhs->end, *r = rhs->end; !(res = *--l - *--r) && l != lhs->mem;);
	return (res > 0) - (res < 0);
}

static inline int cint_compare(const cint * lhs, const cint * rhs) {
	// compare the sign first, then the data
	int res = (int)(lhs->nat - rhs->nat);
	if (res == 0) res = (int) lhs->nat * h_cint_compare(lhs, rhs);
	return res;
}

static void cint_init(cint * num, size_t bits, long long int value) {
	num->size = bits / cint_exponent;
	num->size += 8 - num->size % 4 ;
	num->end = num->mem = calloc(num->size, sizeof(*num->mem));
	assert(num->mem);
	if((num->nat = 1 - ((value < 0) << 1)) < 0) value = -value;
	for (; value; *num->end = (h_cint_t)(value % cint_base), value /= cint_base, ++num->end);
}

static inline void cint_erase(cint * num) {
	num->nat = 1, num->end = memset(num->mem, 0, (num->end - num->mem) * sizeof(h_cint_t));
}

static void cint_reinit(cint * num, long long int value) {
	// it's like an initialization, but without memory allocation
	num->end = memset(num->mem, 0, (num->end - num->mem) * sizeof(h_cint_t));
	if ((num->nat = 1 - ((value < 0) << 1)) < 0) value = -value;
	for (; value; *num->end = (h_cint_t)(value % cint_base), value /= cint_base, ++num->end);
}

__attribute__((unused)) static inline cint * cint_immediate(cint_sheet * sheet, const long long int value){
	cint * res = &sheet->temp[8 + (sheet->immediate_state++ & 1)];
	cint_reinit(res, value);
	return res ;
}

static void cint_reinit_by_string(cint *num, const char *str, const int base) {
	cint_erase(num);
	for (; *str && memchr(cint_alpha, *str, base) == 0; num->nat *= 1 - ((*str++ == '-') << 1));
	for (h_cint_t *p; *str; *num->mem += (h_cint_t) ((char *) memchr(cint_alpha, *str++, base) - cint_alpha), num->end += *num->end != 0)
		for (p = num->end; --p >= num->mem; *(p + 1) += (*p *= base) >> cint_exponent, *p &= cint_mask);
	for (h_cint_t *p = num->mem; p < num->end; *(p + 1) += *p >> cint_exponent, *p++ &= cint_mask);
	num->end += *num->end != 0, num->mem != num->end || (num->nat = 1);
}

static char *cint_to_string(const cint *num, const int base_out) {
	// Very interesting function, only problem is that it can be slow.
	h_cint_t a, b, *c = num->end;
	size_t d, e = 1;
	char *s = malloc(2);
	assert(s);
	strcpy(s, "0");
	for (; --c >= num->mem;) {
		for (a = *c, d = e; d;) {
			b = (h_cint_t) ((char *) memchr(cint_alpha, s[--d], base_out) - cint_alpha), b = b * cint_base + a;
			s[d] = cint_alpha[b % base_out];
			a = b / base_out;
		}
		for (; a; s = realloc(s, ++e + 1), assert(s), memmove(s + 1, s, e), *s = cint_alpha[a % base_out], a /= base_out);
	}
	if (num->nat < 0)
		s = realloc(s, e + 2), assert(s), memmove(s + 1, s, e + 1), *s = '-';
	return s;
}

__attribute__((unused)) static inline void cint_init_by_string(cint *num, const size_t bits, const char *str, const int base) {
	cint_init(num, bits, 0), cint_reinit_by_string(num, str, base);
}

static void cint_reinit_by_double(cint *num, const double value) {
	// sometimes tested, it worked.
	cint_erase(num);
	uint64_t memory;
	memcpy(&memory, &value, sizeof(value));
	uint64_t ex = (memory << 1 >> 53) - 1023, m_1 = 1ULL << 52;
	if (ex < 1024) {
		h_cint_t m_2 = 1 << ex % cint_exponent;
		num->nat *= (value > 0) - (value < 0);
		num->end = 1 + num->mem + ex / cint_exponent;
		h_cint_t *n = num->end;
		for (*(n - 1) |= m_2; --n >= num->mem; m_2 = cint_base)
			for (; m_2 >>= 1;)
				if (m_1 >>= 1)
					*n |= m_2 * ((memory & m_1) != 0);
				else return;
	}
}

__attribute__((unused)) static double cint_to_double(const cint *num) {
	// sometimes tested, it worked.
	uint64_t memory = cint_count_bits(num) + 1022, m_write = 1ULL << 52, m_read = 1 << memory % cint_exponent;
	double res = 0;
	memory <<= 52;
	for (h_cint_t *n = num->end; --n >= num->mem; m_read = 1LL << cint_exponent)
		for (; m_read >>= 1;)
			if (m_write >>= 1)
				memory |= m_write * ((*n & m_read) != 0);
			else
				n = num->mem, m_read = 0;
	memcpy(&res, &memory, sizeof(memory));
	return (double) num->nat * res;
}

__attribute__((unused)) static inline void cint_init_by_double(cint *num, const size_t size, const double value) { cint_init(num, size, 0), cint_reinit_by_double(num, value); }

static void cint_dup(cint *to, const cint *from) {
	// duplicate number (no verification about overlapping or available memory, caller must check)
	const size_t b = from->end - from->mem, a = to->end - to->mem;
	memcpy(to->mem, from->mem, b * sizeof(*from->mem));
	to->end = to->mem + b;
	to->nat = from->nat;
	if (b < a) memset(to->end, 0, (a - b) * sizeof(*to->mem));
}

static void cint_rescale(cint *num, const size_t bits) {
	// rarely tested, it should allow to resize a number transparently.
	size_t curr_size = num->end - num->mem ;
	size_t new_size = 1 + bits / cint_exponent ;
	new_size = new_size + 8 - new_size % 8 ;
	if (curr_size < new_size)
		cint_erase(num), curr_size = 0;
	num->mem = realloc(num->mem, new_size * sizeof(h_cint_t));
	assert(num->mem);
	if (num->size < new_size)
		memset(num->mem + num->size, 0, (new_size - num->size) * sizeof(h_cint_t));
	num->end = num->mem + curr_size ;
	num->size = new_size ;
}

static inline cint * h_cint_tmp(cint_sheet * sheet, const int id, const cint * least){
	// request at least the double of "least" to allow performing multiplication then modulo...
	size_t required = (1 + least->end - least->mem) << 1 ;
	if (sheet->temp[id].size <= required) {
		required = (1 + ((required * cint_exponent) >> 10)) << 10 ;
		cint_rescale(&sheet->temp[id], required);
	}
	return &sheet->temp[id] ;
}

static void h_cint_addi(cint *lhs, const cint *rhs) {
	// perform an addition (without caring of the sign)
	h_cint_t *l = lhs->mem;
	for (const h_cint_t *r = rhs->mem; r < rhs->end;)
		*l += *r++, *(l + 1) += *l >> cint_exponent, *l++ &= cint_mask;
	for (; *l & cint_base; *(l + 1) += *l >> cint_exponent, *l++ &= cint_mask);
	if (rhs->end - rhs->mem > lhs->end - lhs->mem)
		lhs->end = lhs->mem + (rhs->end - rhs->mem);
	lhs->end += *lhs->end != 0;
}

static void h_cint_subi(cint *lhs, const cint *rhs) {
	// perform a subtraction (without caring about the sign, it performs high subtract low)
	h_cint_t a = 0, cmp, *l, *r, *e, *o;
	if (lhs->mem == lhs->end)
		cint_dup(lhs, rhs);
	else if (rhs->mem != rhs->end) {
		cmp = h_cint_compare(lhs, rhs);
		if (cmp) {
			if (cmp < 0) l = lhs->mem, r = rhs->mem, e = rhs->end, lhs->nat = -lhs->nat;
			else l = rhs->mem, r = lhs->mem, e = lhs->end;
			for (o = lhs->mem; r < e; *o = *r++ - *l++ - a, a = (*o & cint_base) != 0, *o++ &= cint_mask);
			for (*o &= cint_mask, o += a; --o > lhs->mem && !*o;);
			lhs->end = 1 + o;
		} else cint_erase(lhs);
	}
}

// regular functions, they care of the input sign
static inline void cint_addi(cint *lhs, const cint *rhs) { lhs->nat == rhs->nat ? h_cint_addi(lhs, rhs) : h_cint_subi(lhs, rhs); }

static inline void cint_subi(cint *lhs, const cint *rhs) { lhs->nat == rhs->nat ? lhs->nat = -lhs->nat, h_cint_subi(lhs, rhs), lhs->mem == lhs->end || (lhs->nat = -lhs->nat), (void) 0 : h_cint_addi(lhs, rhs); }

static void cint_left_shifti(cint *num, const size_t bits) {
	// execute a left shift immediately over the input, for any amount of bits (no verification about available memory)
	if (num->end != num->mem) {
		const size_t a = bits / cint_exponent, b = bits % cint_exponent, c = cint_exponent - b;
		if (a) {
			memmove(num->mem + a, num->mem, (num->end - num->mem + 1) * sizeof(h_cint_t));
			memset(num->mem, 0, a * sizeof(h_cint_t));
			num->end += a;
		}
		if (b) for (h_cint_t *l = num->end, *e = num->mem + a; --l >= e; *(l + 1) |= *l >> c, *l = *l << b & cint_mask);
		num->end += *(num->end) != 0;
	}
}

static void cint_right_shifti(cint *num, const size_t bits) {
	size_t a = bits / cint_exponent, b = bits % cint_exponent, c = cint_exponent - b;
	if (num->end - a > num->mem) {
		if (a) {
			if(num->mem + a > num->end) a = num->end - num->mem;
			memmove(num->mem, num->mem + a, (num->end - num->mem) * sizeof(h_cint_t));
			memset(num->end -= a, 0, a * sizeof(h_cint_t));
		}
		if (b) for (h_cint_t *l = num->mem; l < num->end; *l = (*l >> b | *(l + 1) << c) & cint_mask, ++l);
		if(num->end != num->mem ) num->end -= *(num->end - 1) == 0, num->end == num->mem && (num->nat = 1);
	} else cint_erase(num);
}

static void cint_mul(const cint *lhs, const cint *rhs, cint *res) {
	// the multiplication (no Karatsuba Algorithm, it's the "slow" multiplication)
	h_cint_t *l, *r, *o, *p;
	cint_erase(res);
	if (lhs->mem != lhs->end && rhs->mem != rhs->end) {
		res->nat = lhs->nat * rhs->nat, res->end += (lhs->end - lhs->mem) + (rhs->end - rhs->mem) - 1;
		for (l = lhs->mem, p = res->mem; l < lhs->end; ++l)
			for (r = rhs->mem, o = p++; r < rhs->end; *(o + 1) += (*o += *l * *r++) >> cint_exponent, *o++ &= cint_mask);
		res->end += *res->end != 0;
	}
}

static void cint_powi(cint_sheet *sheet, cint *n, const cint *exp) {
	// read the exponent bit by bit to perform the "fast" exponentiation in place.
	if (n->mem != n->end) {
		size_t bits = cint_count_bits(exp);
		switch (bits) {
			case 0 : cint_reinit(n, n->mem != n->end); break;
			case 1 : break;
			default:;
				cint *a = h_cint_tmp(sheet, 0, n);
				cint *b = h_cint_tmp(sheet, 1, n), *res = n, *tmp;
				cint_erase(a), *a->end++ = 1;
				h_cint_t mask = 1;
				for (const h_cint_t *ptr = exp->mem;;) {
					if (*ptr & mask)
						cint_mul(a, n, b), tmp = a, a = b, b = tmp;
					if (--bits) {
						cint_mul(n, n, b), tmp = n, n = b, b = tmp;
						mask <<= 1;
						if (mask == cint_base) mask = 1, ++ptr;
					} else break;
				}
				if(res != a) cint_dup(res, a);
		}
	}
}

static inline void cint_pow(cint_sheet *sheet, const cint *n, const cint *exp, cint * res) {
	cint_dup(res, n);
	cint_powi(sheet, res, exp);
}

__attribute__((unused)) static void cint_binary_div(const cint *lhs, const cint *rhs, cint *q, cint *r) {
	// the original division algorithm, it doesn't take any temporary variable.
	cint_erase(r);
	if (rhs->end == rhs->mem)
		for (q->nat = lhs->nat * rhs->nat, q->end = q->mem; q->end < q->mem + q->size; *q->end++ = cint_mask); // DBZ
	else {
		h_cint_t a = h_cint_compare(lhs, rhs);
		if (a) {
			cint_erase(q);
			if (a > 0) {
				h_cint_t *l = lhs->end, *k, *qq = q->mem + (lhs->end - lhs->mem);
				for (; --qq, --l >= lhs->mem;)
					for (a = cint_base; a >>= 1;) {
						for (k = r->end; --k >= r->mem; *(k + 1) |= (*k <<= 1) >> cint_exponent, *k &= cint_mask);
						*r->mem += (a & *l) != 0, r->end += *r->end != 0;
						h_cint_compare(r, rhs) >= 0 ? h_cint_subi(r, rhs), *qq |= a : 0;
					}
				q->end += (lhs->end - lhs->mem) - (rhs->end - rhs->mem), q->end += *q->end != 0;
				q->nat = rhs->nat * lhs->nat, (r->end == r->mem) || (r->nat = lhs->nat); // lhs = q * rhs + r
			} else cint_dup(r, lhs);
		} else cint_reinit(q, rhs->nat * lhs->nat);
	}
}

static void h_cint_div_approx(const cint *lhs, const cint *rhs, cint *res) {
	// the division approximation algorithm (answer isn't always exact)
	h_cint_t x, bits = h_cint_compare(lhs, rhs), *o = rhs->end, *p;
	if (bits == 0)
		cint_erase(res), *res->end++ = 1, res->nat = lhs->nat * rhs->nat;
	else if (bits < 0)
		cint_erase(res);
	else {
		cint_dup(res, lhs);
		res->nat *= rhs->nat;
		x = *--o, --o < rhs->mem || (x = x << cint_exponent | *o);
		for (bits = 0; cint_mask < x; x >>= 1, ++bits);
		cint_right_shifti(res, (rhs->end - rhs->mem - 1) * cint_exponent + (bits > 0) * (bits - cint_exponent));
		p = res->end - 3 > res->mem ? res->end - 3 : res->mem;
		for (o = res->end; --o > p; *(o - 1) += (*o % x) << cint_exponent, *o /= x);
		*o /= x;
		res->end -= *(res->end - 1) == 0;
	}
}

static void cint_div(cint_sheet * sheet, const cint *lhs, const cint *rhs, cint *q, cint *r) {
	// The combined division algorithm, it uses the approximation algorithm, "fast" with small inputs.
	assert(rhs->mem != rhs->end);
	cint_erase(q);
	const int cmp = h_cint_compare(lhs, rhs);
	if (cmp < 0)
		cint_dup(r, lhs);
	else if(cmp){
		if (lhs->end < lhs->mem + 3  && rhs->end < rhs->mem + 3){
			// System native division.
			cint_erase(r);
			const h_cint_t a = *lhs->mem | *(lhs->mem + 1) << cint_exponent, b = *rhs->mem | *(rhs->mem + 1) << cint_exponent;
			*q->mem = a / b, *r->mem = a % b;
			if (*q->mem & ~cint_mask) { *++q->end = *q->mem >> cint_exponent, *q->mem &= cint_mask; } q->end += *q->end != 0;
			if (*r->mem & ~cint_mask) { *++r->end = *r->mem >> cint_exponent, *r->mem &= cint_mask; } r->end += *r->end != 0;
		}
		else if(rhs->end == rhs->mem + 1){
			// Special cased "divide by a single word".
			h_cint_t i ;
			cint_erase(r);
			q->end = q->mem + (i = lhs->end - lhs->mem - 1);
			if (lhs->mem[i] < *rhs->mem)
				*r->mem = lhs->mem[i--];
			for(;i >= 0;){
				const h_cint_t tmp = (*r->mem << cint_exponent) | lhs->mem[i];
				q->mem[i--] = tmp / *rhs->mem; *r->mem = tmp % *rhs->mem;
			}
			q->end += *q->end != 0;
			r->end += *r->end != 0;
		} else {
			// Regular division for larger numbers.
			cint *a = h_cint_tmp(sheet, 0, lhs), *b = h_cint_tmp(sheet, 1, lhs);
			cint_dup(r, lhs);
			for (; h_cint_div_approx(r, rhs, b), b->mem != b->end;)
				cint_addi(q, b), cint_mul(b, rhs, a), h_cint_subi(r, a);
			if (r->end != r->mem && r->nat != lhs->nat) // lhs = q * rhs + r
				cint_reinit(b, q->nat), h_cint_subi(q, b), h_cint_subi(r, rhs);
		}
	} else cint_erase(r), *q->end++ = 1 ;
	if (lhs->nat != rhs->nat) // Signs
		q->nat = q->mem == q->end ? 1 : -1, r->nat = r->mem == r->end ? 1 : lhs->nat;
}

__attribute__((unused)) static inline void cint_mul_mod(cint_sheet *sheet, const cint *lhs, const cint *rhs, const cint *mod, cint *res) {
	cint *a = h_cint_tmp(sheet, 2, res), *b = h_cint_tmp(sheet, 3, res);
	cint_mul(lhs, rhs, a);
	cint_div(sheet, a, mod, b, res);
}

static inline void cint_mul_modi(cint_sheet * sheet, cint * lhs, const cint * rhs, const cint * mod){
	cint *a = h_cint_tmp(sheet, 2, lhs), *b = h_cint_tmp(sheet, 3, lhs);
	cint_mul(lhs, rhs, a);
	cint_div(sheet, a, mod, b, lhs);
}

static inline void cint_pow_modi(cint_sheet *sheet, cint *n, const cint *exp, const cint *mod) {
	// same as "power" algorithm, difference is that it take the modulo as soon as possible.
	if (n->mem != n->end) {
		size_t bits = cint_count_bits(exp);
		switch (bits) {
			case 0 : cint_reinit(n, n->mem != n->end); break;
			case 1 : break;
			default:;
				cint *a = h_cint_tmp(sheet, 2, n);
				cint *b = h_cint_tmp(sheet, 3, n);
				cint *c = h_cint_tmp(sheet, 4, n);
				cint_erase(a), *a->end++ = 1;
				h_cint_t mask = 1;
				for (const h_cint_t *ptr = exp->mem;;) {
					if (*ptr & mask)
						cint_mul(a, n, b), cint_div(sheet, b, mod, c, a);
					if (--bits) {
						cint_mul(n, n, b), cint_div(sheet, b, mod, c, n);
						mask <<= 1;
						if (mask == cint_base) mask = 1, ++ptr;
					} else break;
				}
				cint_dup(n, a);
		}
	}
}

__attribute__((unused)) static void cint_pow_mod(cint_sheet *sheet, const cint *n, const cint *exp, const cint *mod, cint *res) {
	cint_dup(res, n);
	cint_pow_modi(sheet, res, exp, mod);
}

static void cint_gcd(cint_sheet * sheet, const cint * lhs, const cint * rhs, cint * gcd){
	// the basic GCD algorithm, by frontal divisions.
	if (rhs->mem == rhs->end)
		cint_dup(gcd, lhs), gcd->nat = 1;
	else {
		cint *A = h_cint_tmp(sheet, 2, lhs),
				*B = h_cint_tmp(sheet, 3, lhs),
				*C = h_cint_tmp(sheet, 4, lhs),
				*TMP, *RES = gcd;
		cint_dup(gcd, lhs);
		cint_dup(A, rhs);
		gcd->nat = A->nat = 1 ;
		for (; A->mem != A->end;) {
			cint_div(sheet, gcd, A, B, C);
			TMP = gcd, gcd = A, A = C, C = TMP;
		}
		gcd->nat = 1 ;
		if (RES != gcd) cint_dup(RES, gcd);
	}

}

__attribute__((unused)) static void cint_binary_gcd(cint_sheet * sheet, const cint * lhs, const cint * rhs, cint * gcd){
	// a binary GCD algorithm.
	if (lhs->mem == lhs->end) cint_dup(gcd, rhs);
	else if(rhs->mem == rhs->end) cint_dup(gcd, lhs);
	else {
		cint *tmp = h_cint_tmp(sheet, 0, lhs),
		*swap, *res = gcd ;
		cint_dup(gcd, lhs), gcd->nat = 1;
		cint_dup(tmp, rhs), tmp->nat = 1;
		const size_t a = cint_count_zeros(lhs), b = cint_count_zeros(rhs);
		for (size_t c = a > b ? b : a;; cint_right_shifti(tmp, cint_count_zeros(tmp))) {
			if (h_cint_compare(gcd, tmp) > 0)
				swap = gcd, gcd = tmp, tmp = swap;
			h_cint_subi(tmp, gcd);
			if (tmp->mem == tmp->end) {
				if (res != gcd)
					cint_dup(res, gcd);
				cint_left_shifti(res, c);
				break;
			}
		}
	}
}

static unsigned cint_remove(cint_sheet * sheet, cint *N, const cint *F) {
	// remove all occurrences of the factor from the input, and return the count.
	size_t res = 0;
	if (N->end != N->mem && F->end != F->mem)
		switch ((*N->mem == 1 && N->end == N->mem + 1) | (*F->mem == 1 && F->end == F->mem + 1) << 1) {
		case 1 : break; // it asks remove other than [-1, 1] but N is [-1, 1].
		case 2 : // it asks remove [-1, 1], so remove one occurrence if N != 0.
		case 3 : res = N->mem != N->end; if (res) N->nat *= F->nat; break;
		default:;
			cint *A = h_cint_tmp(sheet, 2, N), *B = h_cint_tmp(sheet, 3, N);
				// divides N by the factor until there is a remainder
			for (cint *tmp; cint_div(sheet, N, F, A, B), B->mem == B->end; tmp = N, N = A, A = tmp, ++res);
			if (res & 1) cint_dup(A, N);
	}
	return res ;
}

static void cint_sqrt(cint_sheet * sheet, const cint *num, cint *res, cint *rem) {
	// original square root algorithm.
	cint_erase(res), cint_dup(rem, num); // answer ** 2 + rem = num
	if (num->nat > 0 && num->end != num->mem) {
		cint *a = h_cint_tmp(sheet, 0, num), *b = h_cint_tmp(sheet, 1, num);
		cint_erase(a), *a->end++ = 1;
		cint_left_shifti(a, cint_count_bits(num) & ~1);
		for (; a->mem != a->end;) {
			cint_dup(b, res);
			h_cint_addi(b, a);
			cint_right_shifti(res, 1);
			if (h_cint_compare(rem, b) >= 0)
				h_cint_subi(rem, b), h_cint_addi(res, a);
			cint_right_shifti(a, 2);
		}
	}
}

static void cint_cbrt(cint_sheet * sheet, const cint *num, cint *res, cint *rem) {
	// original cube root algorithm.
	cint_erase(res), cint_dup(rem, num); // answer ** 3 + rem = num
	if (num->mem != num->end) {
		cint *a = h_cint_tmp(sheet, 0, num), *b = h_cint_tmp(sheet, 1, num);
		for (size_t c = cint_count_bits(num) / 3 * 3; c < -1U; c -= 3) {
			cint_left_shifti(res, 1);
			cint_dup(a, res);
			cint_left_shifti(a, 1);
			h_cint_addi(a, res);
			cint_mul(a, res, b);
			++*b->mem; // "b" isn't "normalized", it should work for an addition.
			h_cint_addi(b, a);
			cint_dup(a, rem);
			cint_right_shifti(a, c);
			if (h_cint_compare(a, b) >= 0)
				cint_left_shifti(b, c), h_cint_subi(rem, b), cint_erase(b), *b->end++ = 1, h_cint_addi(res, b);
		}
		res->nat = num->nat;
	}
}

static void cint_nth_root(cint_sheet * sheet, const cint *num, const unsigned nth, cint *res) {
	// original nth-root algorithm, it does not try to decompose "nth" into prime factors.
	switch(nth){
		case 0 : cint_reinit(res, num->end == num->mem + 1 && *num->mem == 1) ; break;
		case 1 : cint_dup(res, num); break;
		case 2 : cint_sqrt(sheet, num, res, h_cint_tmp(sheet, 2, num)); break;
		case 3 : cint_cbrt(sheet, num, res, h_cint_tmp(sheet, 2, num)); break;
		default:
			if (num->end > num->mem + 1 || *num->mem > 1) {
				cint *a = h_cint_tmp(sheet, 2, num),
				*b = h_cint_tmp(sheet, 3, num),
				*c = h_cint_tmp(sheet, 4, num),
				*d = h_cint_tmp(sheet, 5, num),
				*e = h_cint_tmp(sheet, 6, num), *r = res, *tmp;
				cint_erase(a), *a->end++ = 1, cint_erase(d), *d->end++ = 1;
				cint_left_shifti(a, (cint_count_bits(num) + nth - 1) / nth);
				h_cint_addi(r, d), cint_reinit(d, nth - 1), cint_reinit(e, nth);
				do {
					tmp = a, a = r, r = tmp, cint_dup(a, num);
					for (unsigned count = nth; --count && (cint_div(sheet, a, r, b, c), tmp = a, a = b, b = tmp, a->mem != a->end););
					cint_mul(r, d, b);
					h_cint_addi(b, a);
					cint_div(sheet, b, e, a, c);
				} while (h_cint_compare(a, r) < 0);
				r == res ? (void) 0 : cint_dup(res, tmp == a ? a : r);
				res->nat = nth & 1 ? num->nat : 1;
			} else cint_dup(res, num);
	}
}

static void cint_nth_root_remainder(cint_sheet * sheet, const cint *num, const unsigned nth, cint *res, cint * rem){
	// nth-root algorithm don't provide the remainder, so here it computes the remainder.
	if (nth == 2) cint_sqrt(sheet, num, res, rem);
	else if(nth == 3) cint_cbrt(sheet, num, res, rem);
	else {
		cint_nth_root(sheet, num, nth, res);
		cint *a = h_cint_tmp(sheet, 2, num);
		cint_reinit(a, nth);
		cint_pow(sheet, res, a, rem);
		cint_subi(rem, num);
	}
}

static void cint_random_bits(cint *num, size_t bits) {
	// provide a random number with exactly the number of bits asked.
	// Normally no one bit more, no one less.
	int i = 0;
	cint_erase(num);
	for (; bits; ++num->end)
		for (i = 0; bits && i < cint_exponent; ++i, --bits)
			*num->end = *num->end << 1 | (rand() & 1);
	if (i) *(num->end - 1) |= 1 << (i - 1);
}

static void cint_modular_inverse(cint_sheet * sheet, const cint * lhs, const cint * rhs, cint * res){
	// original modular inverse algorithm, answer is also called "u1" in extended Euclidean algorithm context.
	if (*rhs->mem > 1 || rhs->end > rhs->mem + 1){
		cint *a = h_cint_tmp(sheet, 2, rhs),
		*b = h_cint_tmp(sheet, 3, rhs),
		*c = h_cint_tmp(sheet, 4, rhs),
		*d = h_cint_tmp(sheet, 5, rhs),
		*e = h_cint_tmp(sheet, 6, rhs),
		*f = h_cint_tmp(sheet, 7, rhs), *tmp, *out = res;
		cint_dup(a, lhs), cint_dup(b, rhs), cint_erase(res), *res->end++ = 1, cint_erase(e);
		a->nat = b->nat = 1 ;
		int i = 0 ;
		do{
			cint_div(sheet, a, b, c, d);
			cint_mul(c, e, f);
			cint_dup(c, res);
			cint_subi(c, f);
			tmp = a, a = b, b = d, d = tmp ;
			tmp = res, res = e, e = c, c = tmp;
		} while(++i, (d->mem == d->end) == (b->mem == b->end));
		if (a->end == a->mem + 1 && *a->mem == 1){
			if (i & 1) cint_addi(res, e);
		} else cint_erase(res);
		if (out != res) cint_dup(out, res);
	} else cint_erase(res);
}

int cint_is_prime(cint_sheet *sheet, const cint *N, int iterations) {
	// is N is considered as a prime number ? the function returns 0 or 1.
	// if the number of Miller-Rabin iteration is negative, the function decides for the caller.
	int res;
	if (*N->mem < 961 && N->mem + 1 >= N->end) {
		int n = (int) *N->mem; // Small numbers for which Miller-Rabin is not used.
		res = (n > 1) & ((n < 6) * 42 + 0x208A2882) >> n % 30 && (n < 49 || (n % 7 && n % 11 && n % 13 && n % 17 && n % 19 && n % 23 && n % 29));
	} else if(res = (*N->mem & 1) != 0, res && iterations) {
		cint *A = h_cint_tmp(sheet, 5, N),
				*B = h_cint_tmp(sheet, 6, N),
				*C = h_cint_tmp(sheet, 7, N);
		size_t a, b, bits = cint_count_bits(N), rand_mod = bits - 3;
		if (iterations < 0) // decides for the caller ... 16 ... 8 ... 4 ... 2 ...
			iterations = 2 << ((bits < 128) + (bits < 256) + (bits < 1024));
		cint_dup(A, N);
		cint_erase(B), *B->end++ = 1;
		cint_subi(A, B);
		cint_dup(C, A); // C = (N - 1)
		a = cint_count_zeros(C);
		cint_right_shifti(C, a); // divides C by 2 until C is odd
		for (bits = 2; iterations-- && res;) {
			cint_random_bits(B, bits); // take any appropriate number
			bits = 3 + *B->mem % rand_mod ;
			cint_pow_modi(sheet, B, C, N); // raise to the power C mod N
			if (*B->mem != 1 || B->end != B->mem + 1) {
				for (b = a; b-- && (res = h_cint_compare(A, B));)
					cint_mul_modi(sheet, B, B, N);
				res = !res;
			} // only a prime number can hold (res = 1) forever
		}
	}
	return res;
}

#endif
