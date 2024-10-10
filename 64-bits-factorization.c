// 64-bit integer factorization algorithm released "as it" into the public domain, without any warranty, express or implied.

// The "worker" algorithm uses the Pollard Rho method when trial division isn't enough
// to fully factor the number, and Miller-Rabin identifies that the number is not prime.

u64 mul_mod(u64 a, u64 b, const u64 mod) {
	u64 res = 0, c; // return (a * b) % mod, avoiding overflow errors while doing modular multiplication.
	for (b %= mod; a; a & 1 ? b >= mod - res ? res -= mod : 0, res += b : 0, a >>= 1, (c = b) >= mod - b ? c -= mod : 0, b += c);
	return res % mod;
}

u64 pow_mod(u64 n, u64 exp, const u64 mod) {
	u64 res = 1; // return (n ^ exp) % mod.
	for (n %= mod; exp; exp & 1 ? res = mul_mod(res, n, mod) : 0, n = mul_mod(n, n, mod), exp >>= 1);
	return res;
}

int is_prime_64_bits(u64 n) {
	// Perform a Miller-Rabin test, it should be a deterministic version.
	static const u64 bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
	static const int n_bases = (int) sizeof(*bases);
	for (int i = 0; i < n_bases; ++i)
		if (n % bases[i] == 0)
			return n == bases[i];
	if (n < bases[n_bases - 1] * bases[n_bases - 1])
		return 1 < n;
	// Depending on the size of the number, we don't need to test all the bases.
	int lim = n < 2152302898747 ? n < 25326001 ? n < 2047 ? 1 : n < 1373653 ? 2 : 3 : n < 3215031751 ? 4 : 5 : n < 341550071728321 ? n < 3474749660383 ? 6 : 7 : n < 3825123056546413051 ? 9 : 12, res = 1, a = 0;
	u64 b, c;
	for (b = n - 1; ~b & 1; b >>= 1, ++a);
	for (int i = 0; i < lim && res; ++i)
		if (c = pow_mod(bases[i], b, n), c != 1) {
			for (int d = a; d-- && (res = c + 1 != n);)
				c = mul_mod(c, c, n);
			res = !res;
		}
	return res;
}

u64 pollard_rho(const u64 n, uint64_t *seed) {
	// Factorize N using the given seed to explore different sequences.
	u64 gcd = 1, a, b, c, i = 0, j = 1, x = 1, y = xor_rand(seed, 1, n - 1);
	for (; gcd == 1; ++i) {
		if (i == j) {
			if (j >> 18)
				break; // Timeout.
			j <<= 1;
			x = y;
		}
		a = y, b = y;
		for (y = 0; a; a & 1 ? b >= n - y ? y -= n : 0, y += b : 0, a >>= 1, (c = b) >= n - b ? c -= n : 0, b += c);
		y = (1 + y) % n;
		for (a = y > x ? y - x : x - y, b = n; (a %= b) && (b %= a););
		gcd = a | b;
	}
	return gcd;
}

u64 nth_root(const u64 n, const u64 nth) {
	// Use an iterative approach for finding the nth-roots.
	u64 a = n, b, c, r = nth ? n + (n > 1) : n == 1 ;
	for (; a < r; b = a + (nth - 1) * r, a = b / nth)
		for (r = a, a = n, c = nth - 1; c && (a /= r); --c);
	return r;
}

u64 square_extraction(u64 *n, int *pow) {
	u64 root = 1;
	if (3 < *n)
		while (root = nth_root(*n, 2), *n == root * root)
			*n = root, *pow <<= 1;
	return 65522U * 65522U < *n ? 65522 : root + 1;
}

void fac_64_worker(state *state, u64 n, fac64_row *rows) {
	uint64_t seed = state->params.rand.custom ^ state->params.rand.seed;
	if (3 < n) {
		int pow = 1;
		if (~n & 1) {
			// Powers of two are removed.
			*rows = (fac64_row) {2, 0};
			do {
				++(*rows).power;
				n >>= 1;
			} while (~n & 1);
			++rows;
		}
		if (8 < n) {
			// The number is odd.
			u64 limit = square_extraction(&n, &pow);
			// Ensure the number has no 16-bit factor by trial division.
			for (u64 prime = 3; prime < limit; prime += 2)
				if (n % prime == 0) {
					int p = 0;
					do ++p, n /= prime;
					while (n % prime == 0);
					*rows++ = (fac64_row) {prime, p * pow};
					limit = square_extraction(&n, &pow);
				}
			while (n >> 32 && !is_prime_64_bits(n)) {
				u64 x, y;
				if (x = nth_root(n, 3), x * x * x == n) // Maybe a cube ?
					*rows++ = (fac64_row) {(n = 1, x), pow * 3};
				else {
					// The number has 2 or 3 prime factors greater than 65536.
					while (x = pollard_rho(n, &seed), x == 1 || x == n);
					n /= x;
					if (x >> 32) {
						y = nth_root(x, 2);
						if (y * y == x)
							// Pollard's Rho produced a square.
							*rows++ = (fac64_row) {y, pow << 1};
						else if (is_prime_64_bits(x))
							// Pollard's Rho produced a prime.
							*rows++ = (fac64_row) {x, pow};
						else {
							// Pollard's Rho produced a composite number.
							while (y = pollard_rho(x, &seed), y == 1 || y == x);
							*rows++ = (fac64_row) {x / y, pow};
							*rows++ = (fac64_row) {y, pow};
						}
					} else if (n % x)
						// Pollard's Rho produced a prime.
						*rows++ = (fac64_row) {x, pow};
					else
						// Pollard's Rho produced a prime that divides the number twice.
						n /= x, *rows++ = (fac64_row) {x, pow + 1};
				}
			}
		}
		if (n != 1)
			*rows++ = (fac64_row) {n, pow};
	} else if (n)
		*rows++ = (fac64_row) {n, 1};
	*rows = (fac64_row) {0};
}
