#ifndef QS_HEADERS
#define QS_HEADERS

#include <errno.h>
#include <stddef.h>

// Quadratic sieve integers, used by all factor functions needing large enough unsigned integers.

typedef uint32_t qs_sm; // small size (32-bit),   the normal native integer used in this implementation.
typedef uint64_t qs_md; // medium size (64-bit),  the large native integer used in this implementation.
typedef int64_t qs_md_tmp_si; // medium size,  signed for intermediates computations.

typedef struct {
	cint cint ;
	int power ;
	int prime ;
	int bits ;
} fac_cint;

typedef struct{
	int limit ;
	int testing ;
	int silent ;
	int help ;
	int qs_multiplier ;
	int has_threads;
} fac_params;

typedef struct {

	struct {
		void * base ;
		void * now ;
		void * end ;
	} mem;

	fac_params * params ;

	struct{
		cint * cint;
		int done_up_to ;
	} trial;

	fac_cint * number ;
	fac_cint factor ;

	cint vars[10];
	cint_sheet * calc ;

	struct {
		fac_cint * data ;
		int index ;
	} questions;

	struct {
		fac_cint * data ;
		int index ;
	} answers;

} fac_caller;

// Front-End Factor manager
static inline fac_cint **c_factor(const cint *, fac_params *);
static inline int fac_special_cases(fac_caller *);
static inline int fac_trial_division(fac_caller *, int);
static inline int fac_perfect_checker(fac_caller *);
static inline int fac_primality_checker(fac_caller *);
static inline int fac_pollard_rho_64_bits(fac_caller *);
static inline void fac_push(fac_caller *, const fac_cint *, int);

// Math
static inline int is_prime_4669921(qs_sm n);
static double log_computation(double);
static inline qs_md multiplication_modulo(qs_md, qs_md, qs_md);
static inline qs_md power_modulo(qs_md, qs_md, qs_md);
static qs_md tonelli_shanks(qs_md, unsigned);
static qs_md modular_inverse(qs_md, qs_md);
static inline qs_md rand_64();
static inline qs_md rand_upto(qs_md );
static inline unsigned mix_rand_seed(void *);

// Cint shortcuts
static inline void simple_inline_cint(cint *N, size_t, void **);
static inline void simple_dup_cint(cint *, const cint *, void **);
static inline void simple_int_to_cint(cint *, qs_md);
static inline qs_md simple_cint_to_int(const cint *);

// Avl;
static inline struct avl_node *avl_cint_inserter(void *, const void *);

// System
static inline void *mem_aligned(void *);

// Misc.
static inline int fac_apply_custom_param(const char *, const char *, int, int *);
static inline char *fac_fill_params(fac_params *, int, char **);
static char *fac_answer_to_string(fac_cint **);
static inline void fac_display_progress(const char *, double);
static inline int fac_sort_result(const void * , const void *);


// Quadratic sieve structures

struct qs_relation {
	qs_sm id ; // only definitive relations have a non-zero id.
	cint *X;
	struct {
		qs_sm *data;
		qs_sm length ;
		qs_sm snapshot ;
	} Y;
	union {
		struct {
			qs_sm *data;
			qs_sm length;
		} Z;
		struct qs_relation * next ;
	} axis ;
	// axis :
	// -  "Z" field is used by standard definitive relations.
	// - "next" is used by data that wait to be paired, it uses a linked list instead of a "Z" field.
};

typedef struct {

	qs_sm n_threads ;

	fac_caller * caller ;

	struct {
		cint N ;
		cint A ;
		cint B ;
		cint C ;
		cint D ;
		cint FACTOR ;
		cint X ;
		cint KEY ;
		cint VALUE ;
		cint TEMP[5];
	} vars;

	struct{
		cint kN ;
		cint ONE ;
		cint UPPER ;
		cint M_2 ;
	} constants;

	qs_sm knuth_schroppel ;

	struct {
		qs_sm bytes_allocated;
		void * base ;
		void * now ;
	} mem;

	cint_sheet * calc ;

	// vars and parameters to the algorithm

		qs_sm n_bits;
		qs_sm kn_bits;
		qs_sm d_bits;

		qs_sm cache_block_size;
		struct {
			qs_sm value;
			qs_sm double_value;
			qs_sm divided;
			qs_sm q;
			qs_sm r;
			qs_sm n_reps;
		} m;
		struct{
			qs_sm the_span ;
			qs_sm the_span_half ;
			qs_sm the_min ;
		} mini;
		qs_sm p_list[10];
		qs_sm error_bits;
		qs_sm threshold;
		qs_sm poly_max;
		qs_sm s_rand;

		qs_sm sieve_again_perms;
		qs_sm curves;

	// useful data sharing same length
	struct {
		struct {
			qs_sm num;
			qs_sm size;
			qs_sm A_inv;
			qs_sm sqrt;
			qs_sm tmp;
			qs_sm sol[2];
		} *data;
		size_t length;
	} base;

	// useful data sharing same length
	struct {
		struct {
			qs_sm defined;
			qs_sm subtract_one;
			qs_sm double_value;
		} values;
		struct {
			cint B_terms;
			qs_sm *A_inv_2B;
			qs_sm a_mod_p;
			qs_sm a_ind;
		} *data;
	} s;

	// useful data (special)
	struct {
		qs_sm *sm_buffer; // proportional to the value of "s" (small)
		qs_sm *md_uncleared_buffer; // proportional to "length of factor base" (medium or large)
		qs_sm *md_cleared_buffer; // proportional to "length of factor base" (medium or large)
		uint8_t *sieve;
		uint8_t **offsets[2];
		uint8_t *flags;
	} others;

	// 3 unicity trees : [ relations, cycle finder, divisors of N, ]
	struct avl_manager unicity[3] ;

	// data analysis made by algorithm after sieving
	struct {
		qs_sm reduced_by_lanczos ;
		struct qs_relation **data;
		struct {
			qs_sm now ;
			qs_sm needs ;
			qs_sm allocated ;
		} length;
	} relations;

	struct {
		// divisors of N are kept (symbolically) in a flat array
		qs_sm processing_index ;
		qs_sm total_primes ;
		qs_sm length ;
		cint ** data ;
	} divisors;

} qs_sheet;

// Quadratic sieve functions
static inline qs_sm linear_param_resolution(const double [][2], qs_sm);
static inline void qs_parametrize(qs_sheet *);
static int quadratic_sieve(fac_caller *);
static inline void preparation_part_1(qs_sheet *, fac_caller *);
static inline void preparation_part_2(qs_sheet *);
static inline void preparation_part_3(qs_sheet *);
static inline void preparation_part_4(qs_sheet *);
static inline void preparation_part_5(qs_sheet *);
static inline qs_sm preparation_part_6(qs_sheet *, cint *);
static inline void get_started_iteration(qs_sheet *qs);
static inline void iteration_part_1(qs_sheet *, cint *);
static inline cint *iteration_part_2(qs_sheet *, const cint *, const cint *);
static inline void iteration_part_3(qs_sheet *, const cint *, cint *);
static inline void iteration_part_4(qs_sheet *, const cint *, const cint *);
static inline qs_sm iteration_part_5(const qs_sheet *, qs_sm, qs_sm **, cint *);
static inline void iteration_part_6(qs_sheet *, const cint *, const cint *);
static inline void iteration_part_7(qs_sheet *, const cint *, const cint *, const cint *, cint *);
static inline void iteration_part_8(qs_sheet *, qs_sm , const qs_sm *);
static inline void iteration_part_9(qs_sheet *, qs_sm , const qs_sm *);
static inline void register_relations(qs_sheet *, const cint *, const cint *, const cint *);
static inline int qs_register_factor(qs_sheet *);
static inline void process_column_array(struct qs_relation *, const qs_sm *);
static inline void register_relation_kind_1(qs_sheet *, const cint *, qs_sm *, const qs_sm *, qs_sm *, const qs_sm *);
static inline void register_relation_kind_2(qs_sheet *, const qs_sm *, const cint *, const cint *);
static inline void finalization_part_1(qs_sheet *qs, const uint64_t *);
static inline void finalization_part_2(qs_sheet *qs);
static inline int inner_continuation_condition(qs_sheet *);
static inline int outer_continuation_condition(qs_sheet *);
static inline int finalization_part_3(qs_sheet *qs);

// Quadratic sieve Lanczos part
static inline void lanczos_mul_MxN_Nx64(const qs_sheet *, const uint64_t *, qs_sm, uint64_t *);
static inline void lanczos_mul_trans_MxN_Nx64(const qs_sheet *, const uint64_t *, uint64_t *);
static void lanczos_mul_64xN_Nx64(const qs_sheet *, const uint64_t *, const uint64_t *, uint64_t *, uint64_t *);
static uint64_t lanczos_find_non_singular_sub(const uint64_t *, const uint64_t *, uint64_t *, uint64_t, uint64_t *);
static void lanczos_mul_Nx64_64x64_acc(qs_sheet *, const uint64_t *, const uint64_t *, uint64_t *, uint64_t *);
static void lanczos_mul_64x64_64x64(const uint64_t *, const uint64_t *, uint64_t *);
static void lanczos_transpose_vector(qs_sheet *, const uint64_t *, uint64_t **);
static void lanczos_combine_cols(qs_sheet *, uint64_t *, uint64_t *, uint64_t *, uint64_t *);
static inline void lanczos_build_array(qs_sheet *, uint64_t ***, size_t, size_t);
static inline uint64_t *lanczos_block_worker(qs_sheet *);
static inline uint64_t * lanczos_block(qs_sheet *);

#endif //QS_HEADERS
