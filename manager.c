void print_help(const char *path) {
	printf("=== Factorization software using Quadratic Sieve ===\n");
	printf("\n");
	printf("This software is released \"as it is\" into the public domain, without any warranty, express or implied.\n");
	printf("\n");
	printf("DESCRIPTION:\n");
	printf("  This software supports factoring numbers through a Self-Initializing Quadratic Sieve (SIQS).\n");
	printf("  The Factorization Manager reads numbers to be factored from either a file or the command line.\n");
	printf("  It performs preliminary check before invoking more advanced algorithms like the Quadratic Sieve.\n");
	printf("\n");
	printf("USAGE:\n");
	printf("  %s [options] [numbers]\n", path);
	printf("\n");
	printf("OPTIONS:\n");
	printf("  -i, --input-file <FILE>        Factor all numbers from the specified input file.\n");
	printf("  -o, --output-file <FILE>       Write results to the specified output file.\n");
	printf("  -t, --timeout <SECONDS>        Set a timeout in seconds to interrupt the Quadratic Sieve after the specified duration.\n");
	printf("  -f, --force                    Override default limits (8191 digits for numbers and 220-bit for the Quadratic Sieve).\n");
	printf("  -v, --verbose                  Display detailed information, including Quadratic Sieve progress.\n");
	printf("  -h, --help                     Show this help message and exit.\n");
	printf("\n");
	printf("QUADRATIC SIEVE OPTIONS:\n");
	printf("  --qs-multiplier <value>\n");
	printf("  --qs-base-size <value>\n");
	printf("  --qs-large-prime <value>\n");
	printf("  --qs-alloc-mb <value>\n");
	printf("  --qs-sieve <value>\n");
	printf("  --qs-threshold <value>\n");
	printf("  --qs-error-bits <value>\n");
	printf("  --qs-laziness <value>\n\n");
	printf("  Navigate through the source code to see their default value and usage.\n");
	printf("\n");
	printf("EXAMPLES:\n");
	printf("  %s  -i input.txt -o output.txt --output-csv   # Factor all numbers from \"input.txt\" to \"output.txt\" in CSV.\n", path);
	printf("  %s 27333597444727959277 36190584594536893817  # Factor the numbers.\n", path);
	printf("\n");
	printf("EXIT STATUS:\n");
	printf("  0  All numbers were successfully fully factored.\n");
	printf("  1  At least one number among the results is not fully factored.\n");
	printf("  \n");
	printf("REPORTING BUGS:\n");
	printf("  You can read the full documentation and report issues to the \"github.com/michel-leonard/C-Quadratic-Sieve\" repository.\n");
	printf("\n");
	printf("TESTING:\n");
	printf("  %s -g -r 123         # provide a 'generated.txt' file that depends on the seed 123, suitable as input file.\n", path);
	printf("  %s -g 150            # provide a 'generated.txt' file containing a single 150-bit sample number.\n", path);
	printf("  %s -g 60 150         # provide a 'generated.txt' file containing sample numbers ranging from 60-bit to 150-bit.\n", path);
	printf("  %s -g 150 150 1000   # provide a 'generated.txt' file containing 1000 sample numbers of 150-bit.\n", path);
	printf("\n");
	printf("ADDITIONAL RESOURCES:\n");
	printf("  For online factorization tasks, consider tools such as:\n");
	printf("    - Dario Alpern's Integer Factorization Calculator: https://www.alpertron.com.ar/ECM.HTM\n");
	printf("    - Number Empire Factoring Calculator: https://www.numberempire.com/factoringcalculator.php\n");
	printf("\n");
	printf("  The last source code update by Michel was made on Tuesday, October 8, 2024.\n\n");
}

qs_md get_num(char *s) {
	char *end = 0;
	const qs_md res = strtoull(s + (*s == '-'), &end, 10);
	return end != s && !*end ? *s == '-' ? -res : res : 0;
}

int cli_param_match(const char * str, const char * long_name, const char * short_name){
	return (short_name && !strcmp(str, short_name)) || !strcmp(str, long_name);
}

int read_key_and_3_values(char **argv, state *state) {
#define READ(name_1, shortcut, name_2) \
    if (cli_param_match(key, "--" #name_1, "-" #shortcut) && \
        (n_1 = get_num(val_1)) && \
        (n_2 = get_num(val_2)) && \
        (n_3 = get_num(val_3))) { \
        state->params.name_2[0] = n_1; \
        state->params.name_2[1] = n_2; \
        state->params.name_2[2] = n_3; \
    }
	char *key = *argv, *val_1 = *(argv + 1), *val_2 = *(argv + 2), *val_3 = *(argv + 3);
	qs_md n_1, n_2, n_3;
	READ(generate, g, generate)
	else
		return 0;
	*argv = *(argv + 1) = *(argv + 2) = *(argv + 3) = 0;
	return 1;
}

int read_key_and_2_values(char **argv, state *state) {
#define FETCH ((n_1 = get_num(val_1)) && (n_2 = get_num(val_2)))
	char *key = *argv, *val_1 = *(argv + 1), *val_2 = *(argv + 2);
	qs_md n_1, n_2;
	if (cli_param_match(key, "--generate", "-g") && FETCH)
		state->params.generate[0] = n_1, state->params.generate[1] = n_2;
	else
		return 0;
	*argv = *(argv + 1) = *(argv + 2) = 0;
	return 1;
#undef FETCH
}

int read_key_value(char **argv, state *state) {
	char *key = *argv, *value = *(argv + 1);
	if (cli_param_match(key, "--verbose", "-v") && *value >= '0' && *value <= '9' && !*(value + 1))
		state->params.verbose = *value - '0' ;
	else if (cli_param_match(key, "--input-file", "-i"))
		state->params.input_file = value;
	else if (cli_param_match(key, "--output-file", "-o"))
		state->params.output_file = value;
	else if (cli_param_match(key, "--timeout", "-t"))
		state->params.timeout = get_num(value) ;
	else if (cli_param_match(key, "--rand-seed", "-r"))
		state->params.rand.seed ^= (state->params.rand.custom = get_num(value));
	else if (cli_param_match(key, "--generate", "-g") && get_num(value))
		state->params.generate[0] = get_num(value);
	// Quadratic Sieve specific parameters.
#define QS_PARAM(name_1, name_2) \
    else if (cli_param_match(key, "--qs-" #name_1, 0)) \
        state->params.qs_##name_2 = get_num(value);
	QS_PARAM(multiplier, multiplier)
	QS_PARAM(base-size, base_size)
	QS_PARAM(large-prime, large_prime)
	QS_PARAM(alloc-mb, alloc_mb)
	QS_PARAM(sieve, sieve)
	QS_PARAM(threshold, threshold)
	QS_PARAM(error-bits, error_bits)
	QS_PARAM(laziness, laziness)
	else
		return 0;
#undef QS_PARAM
	*argv = *(argv + 1) = 0;
	return 1;
}

int read_flags(char **argv, state *state) {
	char *key = *argv;
	if (cli_param_match(key, "--verbose", "-v"))
		state->params.verbose = 1;
	else if (cli_param_match(key, "--output-json", "-j"))
		state->params.output_format = 'J';
	else if (cli_param_match(key, "--output-json-compact", "-J"))
		state->params.output_format = 'j';
	else if (cli_param_match(key, "--output-csv", "-c"))
		state->params.output_format = 'c';
	else if (cli_param_match(key, "--output-sql", "-sql"))
		state->params.output_format = 's';
	else if (cli_param_match(key, "--force", "-f"))
		state->params.force = 1;
	else if (cli_param_match(key, "--help", "-h"))
		state->params.help = 1;
	else if (cli_param_match(key, "--generate", "-g"))
		state->params.generate[0] = -1;
	else
		return 0;
	*argv = 0;
	return 1;
}

void simple_rand(cint_sheet *sheet, uint64_t *seed, cint *arr, char * comment, int n_factors, int n_bits) {
	// The specified number of bits and prime factors will be assembled to provide a
	// non-trivial sample number for factoring, with all choices based on the seed.
	cint *res = arr, * calc = arr + 1, *tmp;
	arr += 2 ;
	int bits[n_factors], i;
	const int begin = n_bits + (n_factors >> 1), big_size = n_factors + 1, upto = n_bits / n_factors + 2;
	do {
		do {
			// Determinate the bit-size of each prime factor of N.
			for (bits[0] = begin, i = 1; i < n_factors; bits[0] -= bits[i], ++i)
				while (bits[i] = (int) xor_rand(seed, 0, upto - 1), bits[i] * big_size < n_bits);
		} while (bits[0] * big_size < n_bits);
		cint_reinit(res, 1);
		*comment = 0 ;
		for (i = 0; i < n_factors; ++i) {
			// Occasionally propose a square (N = P * Q^2), a cube (N = P * Q^3) or a quartic (N = P * Q^4).
			if (i == 2 && 2 < n_factors && bits[1] == bits[2] && !xor_rand(seed, 0, 4))
				sprintf(comment, " and a square"), cint_dup(arr + i, arr + i - 1); // 1 for 350
			else if(i == 3 && 3 < n_factors && !h_cint_compare(arr + i - 1, arr + i - 2))
				sprintf(comment, " and a cube"), cint_dup(arr + i, arr + i - 1); // 1 for 900
			else if(i == 4 && 4 < n_factors && !h_cint_compare(arr + i - 1, arr + i - 2) && !h_cint_compare(arr + i - 2, arr + i - 3))
				sprintf(comment, " and a quartic"), cint_dup(arr + i, arr + i - 1); // 1 for 4,500
			else
				do cint_random_bits(arr + i, bits[i], seed), *arr[i].mem |= 1;
				while (!cint_is_prime(sheet, arr + i, -1));
			cint_mul(res, arr + i, calc), tmp = res, res = calc, calc = tmp;
		}
		// Ensure the size of the resulting N correspond to the request.
	} while (cint_count_bits(res) != n_bits);
	if (res != arr - 2)
		cint_dup(calc, res);
}

void generate_input_file(state *state) {
	// Generate sample numbers based on the command line option -g <min-bits> <max-bits> <count>.
	// The "generated.txt" file is dependent on --rand-seed and consistent across all platforms.
	// For example "--generate 130 140 1000" propose 1000 non-trivial numbers between 130 and 140 bits.
	FILE *fp = fopen("generated.txt", "w");
	if (fp) {
		qs_md *p = state->params.generate;
		if (p[0] == -1)
			p[0] = 60, p[1] = 220, p[2] = 0;
		else if(p[1] == 0)
			p[1] = p[0], p[2] = 0 ;
		for(int i = 0; i < 2; ++i)
			p[i] = p[i] < 16 ? 16 : 512 < p[i] ? p[i] : p[i] ;
		qs_md seed = state->params.rand.seed, *r = &seed;
		int min_bits = (int)p[p[1] < p[0]] , max_bits = (int)p[p[0] < p[1]] ;
		int delta = max_bits - min_bits + 1 ;
		int limit = (int)(p[2] < delta ? (p[2] = delta, 1) : p[2] / delta) ;
		int start_bits = min_bits + (int)p[2] - limit * delta ;
		int count = limit * delta + start_bits - min_bits ;
		cint_sheet *sheet = cint_new_sheet(max_bits << 2);
		cint nums[7];
		int n_nums = sizeof(nums) / sizeof(*nums), max_len = cint_approx_digits_from_bits(max_bits, 10);
		char buf[max_len], title[127], comment[63], * _s = 1 < count ? "s" : "";
		for(int i = 0; i < n_nums; ++i)
			cint_init(nums + i, max_bits << 1, 1);
		sprintf(title, "# Generated %d sample number%s ", count, _s);
		if (min_bits == max_bits)
			sprintf(title + strlen(title), "of %d-bit", max_bits);
		else
			sprintf(title + strlen(title), "ranging from %d-bit to %d-bit", min_bits, max_bits);
		if (state->params.rand.custom)
			sprintf(title + strlen(title), " using seed %" PRIu64, state->params.rand.custom);
		fprintf(fp, "%s\n# Simply use \"--generate %d %d %d" , title, min_bits, max_bits, count);
		if (state->params.rand.custom)
			fprintf(fp, " --rand-seed %" PRIu64, state->params.rand.custom);
		fprintf(fp, "\" to retrieve this file\n\n");
		for(int b = min_bits, total = 0; b <= max_bits; ++b){
			for(int i = -(b < start_bits); i < limit; ++i){
				*r ^= *r << 11, *r ^= *r >> 27, *r ^= (1 + *r) << 26;
				int n_factors = xor_rand(r, 2, xor_rand(r, 2, 5));
				simple_rand(sheet, r, nums, comment, n_factors, b);
				cint_to_string_buffer(nums, buf, 10);
				fprintf(fp, "%-*s # %d bits with %d factors %s\n", max_len, buf, b, n_factors, comment);
				if (!(++total & 0XF))
					display_progress("Factorization file preparation", (double)total * 100.0 / (double)count);
			}
		}
		display_progress(0, 100);
		for(int i = 0; i < n_nums; ++i)
			free(nums[i].mem);
		cint_clear_sheet(sheet);
		fprintf(stdout, "%s in file 'generated.txt'.\n", title);
	} else
		perror("Factorization program generator");
}

void output_sql(state *state, int has_prev, int has_next) {
	if (has_prev == 0)
		fprintf(state->out, "INSERT INTO factorizations (number, factor, power, is_prime, duration_ms) VALUES\n");
	for (int i = 0; state->session.res[i].power; ++i) {
		const char *s = cint_to_string_buffer(&state->session.res[i].num, state->session.output_string, 10);
		const char * comma = has_next || state->session.res[i + 1].power ? ",\n" : "";
		fprintf(state->out, "('%s','%s',%d,%d,%" PRIu64 ")%s", state->session.input_string, s, state->session.res[i].power, state->session.res[i].prime, state->session.duration_ms, comma);
	}
	if (has_next == 0)
		fprintf(state->out, ";\n");
}

void output_json_pretty_print(state *state, int has_prev, int has_next) {
	fprintf(state->out, has_prev ? "    {" : "[\n    {");
	fprintf(state->out, "\n        \"input\": \"%s\",\n        \"factors\": [", state->session.input_string);
	for (int i = 0, pow; (pow = state->session.res[i].power); ++i) {
		const char *s = cint_to_string_buffer(&state->session.res[i].num, state->session.output_string, 10);
		const char *c = i ? "," : "", *p = state->session.res[i].prime ? "true" : "false";
		fprintf(state->out, "%s\n            {\n                \"num\": \"%s\",\n                \"power\": %d,\n                \"prime\": %s\n            }", c, s, pow, p);
	}
	fprintf(state->out, has_next ? "%s]%s},\n" : "%s]%s}\n]\n", "\n        ", "\n    ");
}

void output_json_compact(state *state, int has_prev, int has_next) {
	fprintf(state->out, has_prev ? "{" : "[\n{");
	fprintf(state->out, "\"input\":\"%s\",\"factors\":[", state->session.input_string);
	for (int i = 0, pow; (pow = state->session.res[i].power); ++i) {
		const char *s = cint_to_string_buffer(&state->session.res[i].num, state->session.output_string, 10);
		const char *c = i ? "," : "", *p = state->session.res[i].prime ? "true" : "false";
		fprintf(state->out, "%s{\"num\":\"%s\",\"power\":%d,\"prime\":%s}", c, s, pow, p);
	}
	fprintf(state->out, has_next ? "]},\n" : "]}\n]\n");
}

void output_csv(state *state, int has_prev, int has_next) {
	assert(has_next != -1);
	if (has_prev == 0)
		fprintf(state->out, "Input,Factors\r\n");
	fprintf(state->out, "%s,", state->session.input_string);
	for (int i = 0, pow; (pow = state->session.res[i].power); ++i) {
		const char *s = cint_to_string_buffer(&state->session.res[i].num, state->session.output_string, 10);
		for (int j = 0; j < pow; ++j)
			fprintf(state->out, !i && !j ? "%s" : ";%s", s);
	}
	fprintf(state->out, "\r\n");
}

void output_default(state *state, int has_prev, int has_next) {
	assert(has_prev != -1);
	fprintf(state->out, "Number: %s\nFactors: ", state->session.input_string);
	for (int i = 0, pow; (pow = state->session.res[i].power); ++i) {
		const char *s = cint_to_string_buffer(&state->session.res[i].num, state->session.output_string, 10);
		const char *c = i ? ", " : "";
		if (pow == 1)
			fprintf(state->out, "%s%s (%s)", c, s, state->session.res[i].prime ? "prime" : "not prime");
		else
			fprintf(state->out, "%s%s^%d", c, s, pow);
	}
	fprintf(state->out, has_next ? "\n\n" : "\n");
}

void display_progress(const char *name, double percentage) {
	static int chars = 0;
	if (percentage < 100.)
		// Print a progression percentage.
		chars = printf("\r- %s at %.02f %% ...\r", name, percentage);
	else
		// Clear the progression line.
		chars = !printf("\r%*s\r", chars, "");
	fflush(stdout);
}

void output(state *state) {
	int has_prev = state->scale.row_idx != 0, has_next = state->scale.row_idx + 1 != state->scale.total_rows ;
	display_progress(0, 100);
	switch (state->params.output_format) {
		case 'J' : output_json_pretty_print(state, has_prev, has_next); break;
		case 'j' : output_json_compact(state, has_prev, has_next); break;
		case 'c' : output_csv(state, has_prev, has_next); break;
		case 's' : output_sql(state, has_prev, has_next); break;
		default : output_default(state, has_prev, has_next); break;
	}
	++state->scale.row_idx; // Update the index after each factorization.
	if (1 < state->params.verbose)
		display_progress("Overall Progress", (double)state->scale.row_idx / (double) state->scale.total_rows * 100.0);
	fflush(state->out);
}

int validate_input_file(state *state) {
	FILE *fp = state->in;
	qs_md line = 0;
	while (!feof(fp)) {
		++line;
		char c = fgetc(fp);
		size_t digits = c >= '1' && c <= '9';
		if (digits || c == '-' || c == '+') {
			while (!feof(fp) && (c = fgetc(fp)) >= '0' && c <= '9')
				++digits;
			if (feof(fp) || c == ' ' || c == '\t' || c == '\r' || c == '\n') {
				++state->scale.total_rows;
				if (state->scale.max_digits < digits)
					state->scale.max_digits = digits;
			} else if (!feof(fp))
				return !fprintf(stderr, "Factorization program input: Unknown character '%c' (0x%x) at line %" PRIu64 "\n", c, c, line);
		}
		if (c != '\n')
			while (!feof(fp) && fgetc(fp) != '\n');
	}
	fseek(fp, 0, SEEK_SET);
	return 1;
}

size_t prepare_file_descriptors(state *state) {
	if (state->params.output_file) {
		state->out = fopen(state->params.output_file, "w");
		if (!state->out)
			return perror("Factorization program output"), 0;
	} else
		state->out = stdout; // Standard output.
	if (state->params.input_file) {
		state->in = fopen(state->params.input_file, "rb");
		if (!state->in) {
			perror("Factorization program input");
			if (state->out != stdout)
				fclose(state->out);
			return 0;
		} else if (!validate_input_file(state)) {
			if (state->out != stdout)
				fclose(state->out);
			fclose(state->in);
			return 0;
		}
	}
	// Set a limit for the input size.
	if (state->scale.max_digits >> 13 && !state->params.force)
		return !fprintf(stderr, "A number of %" PRIu64 " digits when \x1b[37;40;1moption -f\033[0m isn't set is too large for the %d limit.\n", state->scale.max_digits, 1 << 13);
	return state->scale.total_rows;
}

int validate_string_number(const char *s, state *state) {
	// Ensure that the string (number to factor) is well-formed, count the 
	// total of submitted numbers, and note the size (decimal digits) of the largest.
	size_t digits;
	s += *s == '-' || *s == '+';
	if (*s >= '1' && *s <= '9' && !s[digits = 1 + strspn(s + 1, "0123456789")]) {
		if (state->scale.max_digits < digits)
			state->scale.max_digits = digits;
		return ++state->scale.total_rows, 1;
	}
	return 0;
}

void debug_print(const state * state, int level, const char *format, ...) {
	if (level < state->params.verbose) {
		va_list args;
		va_start(args, format);
		vfprintf(stderr, format, args);
		va_end(args);
	}
}

char * simple_cint_string(state * state, const cint * N){
	char * s = cint_to_string_buffer(N, state->session.output_string, 10);
	if (0)	// Add a thousand separator for large number.
		for(int len = (int)strlen(s), i = len - 3, j = *s == '-'; j < i; i -= 3)
			memmove(s + i + 1, s + i, ++len - i), s[i] = ',';
	return s ;
}

// cint shortcuts
void simple_inline_cint(cint *N, const size_t size, void **mem) {
	// Fixed size cint is inlined, given mem is updated accordingly.
	N->mem = N->end = (h_cint_t *) *mem;
	*mem = N->mem + (N->size = size + 1);
}

void simple_dup_cint(cint *A, const cint *B, void **mem) {
	// Duplicates cint using the given memory, which is updated accordingly.
	// It uses the minimal size, the duplicate is not resizable.
	A->mem = A->end = (h_cint_t *) *mem;
	cint_dup(A, B);
	A->size = A->end - A->mem + 1;
	*mem = A->mem + A->size;
}

void simple_int_to_cint(cint *num, qs_md value) {
	// Pass the given 64-bit number into the given cint (positive only).
	for (cint_erase(num); value; *num->end++ = (h_cint_t) (value & cint_mask), value >>= cint_exponent);
}

qs_md simple_cint_to_int(const cint *num) {
	// Return the value of a cint as a 64-bit integer (sign is ignored).
	qs_md res = 0;
	for (h_cint_t *ptr = num->end; ptr > num->mem; res = (qs_md) (res * cint_base + *--ptr));
	return res;
}

// Avl
struct avl_node *avl_cint_inserter(void *args, const void *key_to_save) {
	// it expects as result a new node containing the given constant key.
	void *mem = *(void **) args;
	struct avl_node *res = mem;
	res->key = (cint *) (res + 1);
	mem = (cint *) res->key + 1;
	simple_dup_cint(res->key, key_to_save, &mem);
	assert(res->value == 0);
	*(void **) args = mem;
	return res;
}

// System
void *mem_aligned(void *ptr) {
	// Default alignment of the return value is 64.
	char *res __attribute__((aligned(64)));
	res = (char *) ptr + (64 - (uintptr_t) ptr) % 64;
	return res;
}

qs_md get_time_ms() {
	// returns the current Unix timestamp with milliseconds.
	struct timeval time;
	gettimeofday(&time, NULL);
	return (qs_md) time.tv_sec * 1000 + (qs_md) time.tv_usec / 1000;
}


void manager_add_factor(state *state, cint *num, int pow, int is_prime) {
	assert(pow);
	int i = 0;
	while (state->session.res[i].power && h_cint_compare(&state->session.res[i].num, num))
		++i;
	simple_inline_cint(&state->session.res[i].num, num->end - num->mem, &state->session.mem.now);
	cint_dup(&state->session.res[i].num, num);
	state->session.res[i].power = state->session.power * pow;
	state->session.res[i].prime = is_prime;
	//
}

void manager_add_simple_factor(state *state, qs_md num, int pow, int is_prime) {
	assert(pow);
	simple_int_to_cint(state->session.tmp, num);
	manager_add_factor(state, state->session.tmp, pow, is_prime);
}

void factorization_64_bits(state *state) {
	fac64_row res[16];
	qs_md num = simple_cint_to_int(&state->session.num);
	fac_64_worker(state, num, res);
	for (fac64_row *r = res; (*r).power; ++r)
		manager_add_simple_factor(state, (*r).prime, (*r).power, (*r).prime != 1);
	cint_reinit(&state->session.num, 1);
}

int factorization_trial_division(state *state, int stage, int bits) {
	assert(64 < bits);
	int calc = stage == 1 ? (1 << 20) - 23250 * bits + 127 * bits * bits : 4669914 ;
	const qs_md trial_upto = calc < 65522 ? 65522 : 4669914 < calc ? 4669914 : calc;
	cint *a = state->session.tmp, *b = a + 1, *c = a + 2, *tmp;
	cint *n = &state->session.num;
	cint_sheet *sheet = state->session.sheet;
	cint_reinit(a, 1);
	qs_md trial = state->session.trial_start;
	for (; trial < trial_upto; trial += 2)
		if (is_prime_4669913(trial)) {
			a->mem[0] = (h_cint_t) trial;
			cint_div(sheet, n, a, b, c);
			if (c->mem == c->end) {
				int pow = 0;
				do {
					++pow;
					tmp = n, n = b, b = tmp;
					cint_div(sheet, n, a, b, c);
				} while (c->mem == c->end);
				manager_add_simple_factor(state, trial, pow, 1);
				if (n != &state->session.num)
					cint_dup(&state->session.num, n);
				state->session.trial_start = trial + 2;
				return 1;
			}
		}
	state->session.trial_start = trial + 2;
	return 0;
}

int factorization_any_root_checker(state *state, const cint *n, cint *root, cint *rem) {
	int res = 0;
	cint_sheet *sheet = state->session.sheet;
	cint *max = state->session.tmp;
	cint_reinit(max, (h_cint_t) state->session.trial_start - 2);
	const int max_root = (int) cint_count_bits(n);
	for (int nth = 2; nth < max_root; nth += 2)
		if (is_prime_4669913((qs_md) nth)) {
			cint_nth_root_remainder(sheet, n, nth, root, rem);
			if (rem->mem == rem->end) {
				res = nth;
				break;
			}
			if (h_cint_compare(root, max) <= 0)
				break;
			nth -= !(nth & 1);
		}
	return res;
}

int factorization_perfect_power_checker(state *state, int bits) {
	assert(64 < bits);
	cint *root = state->session.tmp + 1, *rem = root + 1;
	int power = factorization_any_root_checker(state, &state->session.num, root, rem);
	if (power) {
		manager_add_factor(state, root, power, -1);
		cint_reinit(&state->session.num, 1);
	}
	return power;
}

int factorization_prime_number_checker(state *state, int bits) {
	assert(64 < bits);
	cint_sheet *sheet = state->session.sheet;
	int is_prime = cint_is_prime(sheet, &state->session.num, -1) != 0;
	if (is_prime) {
		manager_add_factor(state, &state->session.num, 1, 1);
		cint_reinit(&state->session.num, 1);
	}
	return is_prime;
}

int factorization_give_up(state *state, int bits) {
	assert(64 < bits);
	manager_add_factor(state, &state->session.num, state->session.power, 0);
	cint_reinit(&state->session.num, 1);
	return 1;
}

void factor(state *state) {
	state->session.duration_ms = get_time_ms();
	state->session.trial_start = 3;
	state->session.power = 1;
	int start_idx = 0, end_idx;
	cint_dup(state->session.tmp + 9, &state->session.num);
	if (state->session.num.nat < 0) {
		// Add -1 as factor for a negative number.
		cint_reinit(state->session.tmp, -1);
		manager_add_factor(state, state->session.tmp, 1, 0);
		state->session.num.nat = 1;
		++start_idx;
	}
	int bits = (int) cint_count_zeros(&state->session.num);
	if (bits) {
		// Remove the powers of two from the number.
		manager_add_simple_factor(state, 2, bits, 1);
		cint_right_shifti(&state->session.num, bits);
		// The number is odd.
		++start_idx;
	}
	start :;
	bits = (int) cint_count_bits(&state->session.num);
	display_progress(0, 100);
	if (bits < 65) {
		// 64-bit simple Pollard's Rho.
		if (1 < bits || start_idx == 0)
			factorization_64_bits(state);
	} else {
		int res = factorization_trial_division(state, 1, bits)
				  || factorization_perfect_power_checker(state, bits)
				  || factorization_prime_number_checker(state, bits)
				  || factorization_quadratic_sieve(state, bits)
				  || factorization_trial_division(state, 2, bits)
				  || factorization_give_up(state, bits);
		assert(res);
		if (cint_compare_char(&state->session.num, 1))
			manager_add_factor(state, &state->session.num, 1, -1);
	}
	end_idx = start_idx ;
	for (int i = (int)state->scale.max_factors - 1; start_idx <= i  ; --i)
		if (state->session.res[i].prime == -1) {
			cint_dup(&state->session.num, &state->session.res[i].num);
			state->session.power = state->session.res[i].power;
			cint_erase(&state->session.res[i].num);
			if (end_idx == start_idx)
				state->session.mem.now = state->session.res[i].num.mem ;
			state->session.res[i].power = state->session.res[i].prime = 0;
			goto start;
		} else if (end_idx == start_idx && state->session.res[i].power)
			end_idx = 1 + i;

	// Sort the results (they start with a cint) using the unsigned cint comparator.
	qsort(state->session.res + start_idx, end_idx - start_idx, sizeof(*state->session.res), (int (*)(const void *, const void *)) h_cint_compare);

	// Verify the correctness of all factorizations with a fatal error level.
	cint *A = state->session.tmp, *B = A + 1, *PRODUCT_OF_FACTORS = A + 2, *INPUT_NUMBER = A + 9, *TMP;
	cint_reinit(PRODUCT_OF_FACTORS, 1);
	for (int i = 0; state->session.res[i].power; ++i) {
		cint_reinit(B, state->session.res[i].power);
		cint_pow(state->session.sheet, &state->session.res[i].num, B, A);
		cint_mul(PRODUCT_OF_FACTORS, A, B), TMP = PRODUCT_OF_FACTORS, PRODUCT_OF_FACTORS = B, B = TMP;
	}
	// Exit with a non-zero status code if the product of all factors isn't equal to the input number.
	assert(cint_compare(INPUT_NUMBER, PRODUCT_OF_FACTORS) == 0);

	if(state->code == 0)
		for (int i = start_idx; state->session.res[i].power; ++i)
			// Exit code 0 will mean that the software has fully factored all inputs.
			state->code |= !state->session.res[i].prime && cint_count_bits(&state->session.res[i].num) != 1 ;

	state->session.duration_ms = get_time_ms() - state->session.duration_ms;
}

// Manager

void prepare_sessions(state *state) {
	// Prepare a session containing enough memory to handle the largest (in terms of digits) number.
	state->scale.max_bits = cint_approx_bits_from_digits(state->scale.max_digits, 10);
	for (int i = 3, bits = (int) state->scale.max_bits; 0 < bits; i += 2)
		if (is_prime_4669913(i))
			for (int j = (++state->scale.max_factors, i); j >>= 1; --bits);
	state->session.output_string = malloc(state->scale.max_bits + 32 - state->scale.max_bits % 16) ;
	assert(state->session.output_string); // String buffer to store any number represented in any base.
	const size_t bits = (2 + (state->scale.max_bits >> 4)) << 5;
	cint_init(&state->session.num, bits, 0); // The number to be factored.
	const size_t el_count = sizeof(state->session.tmp) / sizeof(*state->session.tmp);
	const size_t el_size = state->session.num.size * sizeof(*state->session.num.mem);
	cint *T = state->session.tmp;
	T[0].mem = T[0].end = calloc(el_count, el_size);
	assert(T[0].mem); // Temporary variables.
	T[0].size = state->session.num.size;
	for (size_t i = 1; i < el_count; ++i)
		T[i].mem = T[i].end = T[i - 1].mem + (T[i].size = T[i - 1].size);
	state->session.sheet = cint_new_sheet(bits);
	assert(state->session.sheet); // A computation sheet.
	state->session.mem.base = calloc(1, state->scale.max_factors * (sizeof(*state->session.res) + sizeof(*state->session.num.mem)) + (el_size << 1));
	assert(state->session.mem.base);
	state->session.res = state->session.mem.base ;
	state->session.mem.now = state->session.res + state->scale.max_factors ;
}

void erase_session(state *state) {
	// Erase the session (clear variables, number, and results).
	for (size_t i = 0; i < sizeof(state->session.tmp) / sizeof(*state->session.tmp); ++i)
		if (state->session.tmp[i].mem != state->session.tmp[i].end)
			cint_erase(&state->session.tmp[i]);
	cint_erase(&state->session.num);
	memset(state->session.mem.base, 0, (char*)state->session.mem.now - (char*)state->session.mem.base);
	state->session.mem.now = state->session.res + state->scale.max_factors ;
}

void clear_sessions(state *state) {
	// Free all memory, close all descriptors.
	free(state->session.output_string);
	free(state->session.tmp[0].mem);
	free(state->session.num.mem);
	cint_clear_sheet(state->session.sheet);
	free(state->session.mem.base);
	if (state->in)
		fclose(state->in);
	if (state->out != stdout)
		fclose(state->out);
}

/*
 How are performed memory allocations ?
 - A session is created for the largest number to be processed.
 - The session contains the memory (string buffer, numbers, results).
 - When processing multiple numbers in a row, the session is reused.
 - The Quadratic Sieve manage its memory independently.
 */

void process_single(state *state) {
	cint_reinit_by_string(&state->session.num, state->session.input_string, 10);
	factor(state);
	output(state);
	state->duration_ms += state->session.duration_ms ;
	if (state->scale.row_idx != state->scale.total_rows)
		// Get ready for the next request.
		erase_session(state);
}

void process_multi(int argc, char **argv, state *state) {
	prepare_sessions(state);
	if (state->params.input_file) {
		// Process request(s) incoming from a file.
		const size_t buf_size = state->scale.max_digits + 4096;
		char *s = state->session.input_string = malloc(buf_size);
		assert(s);
		while (fgets(s, buf_size, state->in))
			if ((*s >= '1' && *s <= '9') || ((*s == '+' || *s == '-') && *(s + 1) >= '1' && *(s + 1) <= '9')) {
				s[strspn(s, "+-0123456789")] = 0;
				process_single(state);
			}
		free(s);
	} else
		// Process request(s) incoming from the command line.
		for (int i = 1; i < argc; ++i)
			if (argv[i]) {
				state->session.input_string = argv[i];
				process_single(state);
			}
	if (1 < state->params.verbose)
		fprintf(stderr, "\nTook: %.02f s.\n", (double)state->duration_ms * 0.001);
	// Clear all memory used to process the requests.
	clear_sessions(state);
}
