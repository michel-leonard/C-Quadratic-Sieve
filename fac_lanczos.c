//  GNU General Public License

//      As an undergraduate student, this project is part of my computer science + maths training.
//          - This software proposition is from Michel Leonard (student at Université de Franche-Comté, Mon, 11 Jul 2022)
//          - There is of course no guarantee of any kind on the software
//          - C code is shared under the terms of the GNU General Public License
//          - The main mathematical and logical inspiration source is located at :
//              http://web.mit.edu/sage/export/flintqs-0.0.20070817/QS.cpp - GNU General Public License

//  This software implementation would have been impossible without "FLINT: Fast Library for Number Theory" maintained by William Hart.

// results of operations are last AND non-const arguments.

static inline void lanczos_mul_MxN_Nx64(const qs_sheet *qs, const uint64_t *X, uint64_t *Y) {
	assert(Y != X);
	memset(Y, 0, qs->relations.length.expected * sizeof(uint64_t));
	for (qs_sm a = 0, b ; a < qs->relations.length.now; ++a) {
		struct qs_relation *const rel = qs->relations.data + a;
		for (b = 0; b < rel->Y.length; ++b)
			Y[rel->Y.data[b]] ^= X[a];
	}
}

static inline void lanczos_mul_trans_MxN_Nx64(const qs_sheet *qs, const uint64_t *X, uint64_t *Y) {
	assert(Y != X);
	for (qs_sm a = 0, b; a < qs->relations.length.now; ++a) {
		struct qs_relation *const rel = qs->relations.data + a;
		for (Y[a] = 0, b = 0; b < rel->Y.length; ++b)
			Y[a] ^= X[rel->Y.data[b]];
	}
}

static void lanczos_mul_64xN_Nx64(const qs_sheet *qs, const uint64_t *X, const uint64_t *Y, uint64_t *Z, uint64_t *T) {
	assert(X != Z && Y != Z);
	qs_sm a, b, c, d;
	memset(Z, 0, 256 * 8 * sizeof(*Z));
	memset(T, 0, 64 * sizeof(*T));
	for(a = 0; a < qs->relations.length.now; ++a) {
		const uint64_t tmp = X[a]; // read while writing ?!
		for (b = 0, c = 0; c < 64; c += 8, b += 256)
			Z[b + (tmp >> c & 0xff)] ^= Y[a];
	}
	for (a = 0; a < 8; ++a, ++T) {
		uint64_t tmp[8] = {0};
		for (b = 0; b < 256; ++b)
			if (b >> a & 1)
				for(c = d = 0; c < 8; ++c, d += 256)
					tmp[c] ^= Z[b + d];
		for (b = 0, c = 0; b < 8; ++b, c += 8)
			T[c] = tmp[b];
	}

}

static uint64_t lanczos_find_non_singular_sub(const uint64_t *t, const uint64_t * last_s, uint64_t *s, uint64_t last_dim, uint64_t *w) {
	uint64_t i, j, dim, cols[64];
	uint64_t M[64][2], mask, *row_i, *row_j, m_0, m_1;
	for (i = 0; i < 64; ++i) {
		M[i][0] = t[i], M[i][1] = 1LLU << i;
	}
	mask = 0;
	for (i = 0; i < last_dim; ++i)
		mask |= 1LLU << (cols[63 - i] = last_s[i]);
	for (i = j = 0; i < 64; ++i)
		if (!(mask & (1LLU << i)))
			cols[j++] = i;
	for (i = dim = 0; i < 64; ++i) {
		mask = 1LLU << (cols[i]);
		row_i = M[cols[i]];
		for (j = i; j < 64; ++j) {
			row_j = M[cols[j]];
			if (row_j[0] & mask) {
				m_0 = row_j[0];
				m_1 = row_j[1];
				row_j[0] = row_i[0];
				row_j[1] = row_i[1];
				row_i[0] = m_0;
				row_i[1] = m_1;
				break; // i = j = 64 ;
			}
		}
		if (j < 64) {
			for (j = 0; j < 64; ++j) {
				row_j = M[cols[j]];
				if (row_i != row_j && (row_j[0] & mask)) {
					row_j[0] ^= row_i[0], row_j[1] ^= row_i[1];
				}
			}
			s[dim++] = cols[i];
			continue;
		}
		for (j = i; j < 64; ++j) {
			row_j = M[cols[j]];
			if (row_j[1] & mask) {
				m_0 = row_j[0];
				m_1 = row_j[1];
				row_j[0] = row_i[0];
				row_j[1] = row_i[1];
				row_i[0] = m_0;
				row_i[1] = m_1;
				break; // i = j = 64 ;
			}
		}
		if (j == 64) {
			// submatrix is not invertible
			return 0;
		}

		for (j = 0; j < 64; ++j) {
			row_j = M[cols[j]];
			if (row_i != row_j && (row_j[1] & mask)) {
				row_j[0] ^= row_i[0], row_j[1] ^= row_i[1];
			}
		}

		row_i[0] = row_i[1] = 0;
	}
	for (i = 0; i < 64; ++i)
		w[i] = M[i][1];
	mask = 0;
	for (i = 0; i < dim; ++i)
		mask |= 1LLU << s[i];
	for (i = 0; i < last_dim; ++i)
		mask |= 1LLU << last_s[i];
	if (mask != -1LLU) {
		// not all columns used
		return 0;
	}
	return dim;
}

static void lanczos_mul_Nx64_64x64_acc(qs_sheet *qs, const uint64_t *X, const uint64_t *Y, uint64_t *Z, uint64_t *T) {
	qs_sm a, b, c, d ;
	for (a = 0; a < 8; Y += 8, Z += 256, ++a)
		for (b = 0; b < 256; ++b)
			for (c = Z[b] = 0, d = b; d; d >>= 1, ++c)
				if (d & 1)
					Z[b] ^= Y[c];
	for (a = 0, Z -= 2048; a < qs->relations.length.now; ++a)
		for(b = c = 0; b < 64; b += 8, c += 256) {
			const uint64_t w = X[a];
			T[a] ^= Z[c + (w >> b & 0xff)];
		}
}

static void lanczos_mul_64x64_64x64(const uint64_t *X, const uint64_t *Y, uint64_t *Z) {
	uint64_t a, b, c, d, tmp[64];
	for (a = 0; a < 64; tmp[a++] = c) {
		for (b = 0, c = 0, d = X[a]; d; d >>= 1, ++b)
			if (d & 1)
				c ^= Y[b];
	}
	memcpy(Z, tmp, sizeof(tmp));
}

static void lanczos_transpose_vector(qs_sheet *qs, const uint64_t *X, uint64_t **Y) {
	uint64_t a, b, c, d, * Z ; // Z is zeroed during iteration, nobody above would notice.
	Z = memcpy(qs->mem.now, X, qs->relations.length.now * sizeof(*X));
	for (a = 0; a < qs->relations.length.now; ++a)
		for (b = 0, c = a >> 6, d = 1LLU << (a % 64); Z[a]; Z[a] >>= 1, ++b)
			if (Z[a] & 1)
				Y[b][c] |= d;
}

static void lanczos_combine_cols(qs_sheet *qs, uint64_t *x, uint64_t *v, uint64_t *ax, uint64_t *av) {
	int64_t i, j, k, bit_pos, col, col_words, num_deps ;
	uint64_t mask, *matrix[128], *amatrix[128], *tmp;
	char * ptr_1 = qs->mem.now, *ptr_2 ;
	num_deps = 64 << (v && av);
	col_words = (qs->relations.length.now + 63) / 64;
	for (i = 0; i < num_deps; ++i) {
		matrix[i] = qs->mem.now;
		amatrix[i] =  matrix[i] + col_words;
		qs->mem.now = amatrix[i] + col_words;
	}
	ptr_2 = qs->mem.now ;
	lanczos_transpose_vector(qs, x, matrix);
	lanczos_transpose_vector(qs, ax, amatrix);
	if (num_deps == 128) {
		lanczos_transpose_vector(qs, v, matrix + 64);
		lanczos_transpose_vector(qs, av, amatrix + 64);
	}
	for (i = bit_pos = 0; i < num_deps && bit_pos < qs->relations.length.now; ++bit_pos) {
		mask = 1LLU << (bit_pos % 64);
		col = bit_pos / 64;
		for (j = i; j < num_deps; ++j)
			if (amatrix[j][col] & mask) {
				tmp = matrix[i];
				matrix[i] = matrix[j];
				matrix[j] = tmp;
				tmp = amatrix[i];
				amatrix[i] = amatrix[j];
				amatrix[j] = tmp;
				break;
			}
		if (j == num_deps)
			continue;
		for (++j; j < num_deps; ++j)
			if (amatrix[j][col] & mask)
				for (k = 0; k < col_words; ++k) {
					amatrix[j][k] ^= amatrix[i][k];
					matrix[j][k] ^= matrix[i][k];
				}
		++i;
	}

	for (j = 0; j < qs->relations.length.now; ++j) {
		uint64_t word = 0;
		col = j / 64;
		mask = 1LLU << (j % 64);
		for (k = i; k < 64; ++k) {
			if (matrix[k][col] & mask) {
				word |= 1LLU << k;
			}
		}
		x[j] = word;
	}
	qs->mem.now = memset(ptr_1, 0, ptr_2 - ptr_1);
}

static inline void lanczos_build_array(qs_sheet *qs, uint64_t *** target, const size_t n_rows, const size_t n_cols){
	*target = mem_aligned(qs->mem.now);
	qs->mem.now = mem_aligned(*target + n_rows) ;
	for(size_t i = 0; i < n_rows; ++i) {
		(*target)[i] = qs->mem.now, qs->mem.now = mem_aligned((*target)[i] + n_cols);
	}
}

static inline uint64_t *lanczos_block_worker(qs_sheet *qs) {
	const uint64_t n_cols = qs->relations.length.now, v_size = qs->relations.length.expected;
	uint64_t **md, **xl, **sm, *tmp, i, iter, dim_0, dim_1, mask_0, mask_1 ;
	char *ptr_1, *ptr_2;
	lanczos_build_array(qs, &md, 6, v_size);
	lanczos_build_array(qs, &sm, 13, 64);
	lanczos_build_array(qs, &xl, 2, 1 << 17);
	for (i = 0; i < 64; ++i)
		sm[12][i] = i;
	dim_0 = 0;
	dim_1 = 64;
	mask_1 = (uint64_t) -1;
	iter = 0;
	for (i = 0; i < qs->relations.length.now; ++i)
		md[1][i] = (uint64_t) rand_64();
	memcpy(md[0], md[1], v_size * sizeof(uint64_t));
	lanczos_mul_MxN_Nx64(qs, md[1], xl[1]);
	lanczos_mul_trans_MxN_Nx64(qs, xl[1], md[1]);
	memcpy(xl[0], md[1], v_size * sizeof(uint64_t));
	for (; ++iter;) {
		lanczos_mul_MxN_Nx64(qs, md[1], xl[1]);
		lanczos_mul_trans_MxN_Nx64(qs, xl[1], md[4]);
		lanczos_mul_64xN_Nx64(qs, md[1], md[4], xl[1], sm[3]);
		lanczos_mul_64xN_Nx64(qs, md[4], md[4], xl[1], sm[5]);
		for (i = 0; i < 64 && !(sm[3][i]); ++i);
		if (i == 64) {
			break;
		}
		dim_0 = lanczos_find_non_singular_sub(sm[3], sm[12], sm[11], dim_1, sm[0]);

		if (dim_0 == 0) {
			break;
		}
		mask_0 = 0;
		for (i = 0; i < dim_0; ++i)
			mask_0 |= 1LLU << sm[11][i];
		for (i = 0; i < 64; ++i)
			sm[7][i] = (sm[5][i] & mask_0) ^ sm[3][i];
		lanczos_mul_64x64_64x64(sm[0], sm[7], sm[7]);
		for (i = 0; i < 64; ++i)
			sm[7][i] ^= 1LLU << i;
		lanczos_mul_64x64_64x64(sm[1], sm[3], sm[8]);
		for (i = 0; i < 64; ++i)
			sm[8][i] &= mask_0;
		lanczos_mul_64x64_64x64(sm[4], sm[1], sm[9]);
		for (i = 0; i < 64; ++i)
			sm[9][i] ^= 1LLU << i;
		lanczos_mul_64x64_64x64(sm[2], sm[9], sm[9]);
		for (i = 0; i < 64; ++i)
			sm[10][i] = ((sm[6][i] & mask_1) ^ sm[4][i]) & mask_0;
		lanczos_mul_64x64_64x64(sm[9], sm[10], sm[9]);
		for (i = 0; i < qs->relations.length.now; ++i)
			md[4][i] &= mask_0;
		lanczos_mul_Nx64_64x64_acc(qs, md[1], sm[7], xl[1], md[4]);
		lanczos_mul_Nx64_64x64_acc(qs, md[3], sm[8], xl[1], md[4]);
		lanczos_mul_Nx64_64x64_acc(qs, md[2], sm[9], xl[1], md[4]);
		lanczos_mul_64xN_Nx64(qs, md[1], xl[0], xl[1], sm[7]);
		lanczos_mul_64x64_64x64(sm[0], sm[7], sm[7]);
		lanczos_mul_Nx64_64x64_acc(qs, md[1], sm[7], xl[1], md[0]);
		tmp = md[2], md[2] = md[3], md[3] = md[1], md[1] = md[4], md[4] = tmp;
		tmp = sm[2], sm[2] = sm[1], sm[1] = sm[0], sm[0] = tmp ;
		tmp = sm[4], sm[4] = sm[3], sm[3] = tmp;
		tmp = sm[6], sm[6] = sm[5], sm[5] = tmp;
		memcpy(sm[12], sm[11], 64 * sizeof(int64_t));
		mask_1 = mask_0;
		dim_1 = dim_0;
		assert(iter < 365);
	}

	if(dim_0) {
		lanczos_mul_MxN_Nx64(qs, md[0], md[3]);
		lanczos_mul_MxN_Nx64(qs, md[1], md[2]);
		lanczos_combine_cols(qs, md[0], md[1], md[3], md[2]);
		lanczos_mul_MxN_Nx64(qs, md[0], md[1]);
		for (i = 0; i < n_cols; ++i)
			if (md[1][i]) // shouldn't happen
				md[0] = 0, i = n_cols;
	} else
		md[0] = 0; // linear algebra failed

	ptr_1 = (char*) md[1];
	ptr_2 = qs->mem.now ; // try to clear all but the answer.
	qs->mem.now = memset(ptr_1, 0, ptr_2 - ptr_1);
	return *md;
}

static inline uint64_t *lanczos_block(qs_sheet *qs) {
	uint64_t *res;
	// algorithm do not write to the relations, it only reads, save original lengths.
	const qs_sm initial_expected = qs->relations.length.expected, initial_now = qs->relations.length.now;
	for (qs_sm i = 0; i < 3; ++i) {
		// perform the setup.
		qs->relations.length.expected = qs->relations.length.now > qs->base.length ? qs->relations.length.now : qs->base.length;
		res = lanczos_block_worker(qs);
		if (res)
			break;
		// some relations are thrown if it fails, excepting the first time.
		qs->relations.length.expected -= i << 2, qs->relations.length.now -= i << 2;
	}
	// restore original lengths.
	qs->relations.length.expected = initial_expected, qs->relations.length.now = initial_now;
	return res;
}
