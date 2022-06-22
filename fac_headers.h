#ifndef QS_HEADERS
#define QS_HEADERS

#include <errno.h>
#include <stddef.h>

// Quadratic sieve integers, used by all functions needing "large enough" unsigned integers.

typedef uint32_t qs_sm; // small size (32-bit),   the normal native integer used in this implementation.
typedef uint64_t qs_md; // medium size (64-bit),  the large native integer used in this implementation.
typedef int64_t qs_md_tmp_si; // medium size,  signed for intermediates computations.

// The factorization manager calls the quadratic sieve with an appropriated input
// Quadratic sieve can assume N is 64+ bits, isn't trial divisible, isn't a perfect power

typedef struct {
	cint cint ;
	int power ;
	int prime ;
	int bits ;
} fac_cint;

typedef struct{
	unsigned limit ;
	unsigned testing ;
	unsigned silent ;
	unsigned help ;
	unsigned qs_multiplier ;
	unsigned qs_rand_seed ;
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

// Quadratic sieve structures

struct qs_relation {
	qs_sm id ; // definitive relations have a non-zero id.
	cint *X;
	struct {
		qs_sm *data;
		qs_sm length ;
	} Y;
	union {
		struct {
			qs_sm *data;
			qs_sm length;
		} Z;
		struct qs_relation * next ;
	} axis ;
	// axis :
	// -  "Z" field is used by definitive relations.
	// - "next" is used by data that wait to be paired, it uses a linked list instead of a "Z" field.
};

typedef struct {

	// reference to the caller
	fac_caller *caller;

	// computation sheet
	cint_sheet *calc;

	// numbers updated
	struct {
		cint N;
		cint A;
		cint B;
		cint C;
		cint D;
		cint FACTOR;
		cint X;
		cint KEY;
		cint VALUE;
		cint CYCLE;
		cint TEMP[5];
	} variables;

	// constants
	struct {
		cint kN;
		cint ONE;
		cint SMALL_PRIME;
		cint LARGE_PRIME;
		cint M;
	} constants;

	// system
	struct {
		qs_sm bytes_allocated;
		void *base;
		void *now;
	} mem;

	// parameters and miscellaneous variables
	qs_sm knuth_schroeppel;
	qs_sm n_bits;
	qs_sm kn_bits;
	qs_sm d_bits;
	struct {
		qs_sm value;
		qs_sm double_value;
		qs_sm divided;
		qs_sm q;
		qs_sm r;
		qs_sm n_reps;
	} m;
	struct {
		qs_sm span;
		qs_sm span_half;
		qs_sm min;
	} the;
	qs_sm iterative_list[10];
	qs_sm error_bits;
	struct{
		qs_sm value ;
	}threshold;
	qs_sm poly_max;
	qs_sm rand_seed;
	qs_sm curves;
	qs_sm sieve_again_perms;

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
		qs_sm *a_invariants; // proportional to the value of "s" (small)
		qs_sm *buffer[2]; // proportional to "length of factor base" (medium or large)
		uint8_t *sieve;
		uint8_t **pos[2];
		uint8_t *flags;
	} others;

	// uniqueness trees : [ relations, cycle finder, divisors of N, ]
	struct avl_manager uniqueness[3];

	// data collection made by algorithm
	struct {
		struct qs_relation **data;
		struct {
			qs_sm now;
			qs_sm prev;
			qs_sm needs;
			qs_sm max_known;
			qs_sm reserved;
		} length;
	} relations;

	// pointers to the divisors of N are kept in a flat array
	struct {
		qs_sm processing_index;
		qs_sm total_primes;
		qs_sm length;
		cint **data;
	} divisors;

	// lanczos has its own struct
	struct {
		qs_sm safe_length ;
		struct {
			struct qs_relation *relation;
			qs_sm y_length;
		} * snapshot ;
	} lanczos;

} qs_sheet;

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
static inline unsigned add_rand_seed(void *);

// Cint shortcuts
static inline void simple_inline_cint(cint *, size_t, void **);
static inline void simple_dup_cint(cint *, const cint *, void **);
static inline void simple_int_to_cint(cint *, qs_md);
static inline qs_md simple_cint_to_int(const cint *);

// Avl;
static inline struct avl_node *avl_cint_inserter(void *, const void *);

// System
static inline void *mem_aligned(void *);

// Misc.
static inline int fac_apply_custom_param(const char *, const char *, int, unsigned *);
static inline char *fac_fill_params(fac_params *, int, char **);
static char *fac_answer_to_string(fac_cint **);
static inline void fac_display_progress(const char *, double);
static inline int fac_sort_result(const void * , const void *);

// Quadratic sieve functions
static inline qs_sm linear_param_resolution(const double [][2], qs_sm);
static inline void qs_parametrize(qs_sheet *);
static int quadratic_sieve(fac_caller *);
static inline void preparation_part_1(qs_sheet *, fac_caller *);
static inline void preparation_part_2(qs_sheet *);
//
static inline void preparation_part_3(qs_sheet *);
static inline qs_sm preparation_part_3_original(qs_sheet *);
static inline qs_sm preparation_part_3_proposition(qs_sheet *);
//
static inline void preparation_part_4(qs_sheet *);
static inline void preparation_part_5(qs_sheet *);
static inline qs_sm preparation_part_6(qs_sheet *, cint *);
static inline void get_started_iteration(qs_sheet *);
static inline void iteration_part_1(qs_sheet *, cint *);
static inline void iteration_part_2(qs_sheet *, const cint *, cint *);
static inline void iteration_part_3(qs_sheet *, const cint *, cint *);
static inline void iteration_part_4(qs_sheet *, const cint *, const cint *);
static inline qs_sm iteration_part_5(const qs_sheet *, qs_sm, qs_sm **, cint *);
static inline void iteration_part_6(qs_sheet *, const cint *, const cint *);
static inline void iteration_part_7(qs_sheet *, const cint *, const cint *, const cint *, cint *);
static inline void iteration_part_8(qs_sheet *, qs_sm , const qs_sm *);
static inline void iteration_part_9(qs_sheet *, qs_sm , const qs_sm *);
static inline void process_column_array(struct qs_relation *, const qs_sm *);
static inline int qs_register_factor(qs_sheet *);
static inline void register_relation_kind_2(qs_sheet *, const qs_sm *, const cint *, const cint *);
static inline void register_relation_kind_1(qs_sheet *, const cint *, qs_sm *, const qs_sm *, qs_sm *, const qs_sm *);
static inline void register_relations(qs_sheet *, const cint *, const cint *, const cint *);
static inline void finalization_part_1(qs_sheet *, const uint64_t *);
static inline void finalization_part_2(qs_sheet *);
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
static inline void lanczos_build_array(qs_sheet *, uint64_t **, size_t, size_t);
static inline uint64_t *lanczos_block_worker(qs_sheet *);
static inline uint64_t * lanczos_block(qs_sheet *);
#endif //QS_HEADERS
