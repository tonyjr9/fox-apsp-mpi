/*
 * Usage:  ./matgen N [count] [seed]
 *           N     matrix dimension (read from stdin if omitted)
 *           count number of matrices to generate (default 1)
 *           seed  RNG seed (default: current time)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

#define MAX_VAL 65536u /* 2^16, inclusive upper bound */

/* rand() may only provide 15 bits (e.g. on Windows), so combine two calls. */
static uint32_t rand30(void)
{
    return ((uint32_t)rand() << 15) ^ (uint32_t)rand();
}

/* Uniform-ish random value in [0, 2^16]. */
static uint32_t rand_val(void)
{
    return rand30() % (MAX_VAL + 1u);
}

static uint32_t *make_matrix(size_t n)
{
    uint32_t *m = malloc(n * n * sizeof *m);
    if (!m) {
        fprintf(stderr, "Out of memory for %zux%zu matrix\n", n, n);
        exit(EXIT_FAILURE);
    }
    for (size_t i = 0; i < n * n; i++)
        m[i] = rand_val();
    return m;
}

static void print_matrix(const uint32_t *m, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++)
            printf("%s%u", j ? " " : "", m[i * n + j]);
        putchar('\n');
    }
}

int main(int argc, char **argv)
{
    long n = 0, count = 1;
    unsigned seed = (unsigned)time(NULL);

    if (argc > 1) {
        n = strtol(argv[1], NULL, 10);
    } else {
        fprintf(stderr, "Enter N: "); /* stderr so it doesn't pollute redirected output */
        if (scanf("%ld", &n) != 1) {
            fprintf(stderr, "Invalid input\n");
            return EXIT_FAILURE;
        }
    }
    if (argc > 2) count = strtol(argv[2], NULL, 10);
    if (argc > 3) seed = (unsigned)strtoul(argv[3], NULL, 10);

    if (n <= 0 || count <= 0) {
        fprintf(stderr, "N and count must be positive integers\n");
        return EXIT_FAILURE;
    }

    srand(seed);

    for (long k = 0; k < count; k++) {
        uint32_t *m = make_matrix((size_t)n);
        printf("%ld\n", n); /* header: matrix dimension, as the readers expect */
        print_matrix(m, (size_t)n);
        free(m);
    }
    return EXIT_SUCCESS;
}
