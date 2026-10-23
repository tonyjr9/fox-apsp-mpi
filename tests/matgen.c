/*
 * Usage:  ./matgen N [count] [seed]
 *       N   matrix dimension (read from stdin if omitted)
 *       count number of matrices to generate (default 1)
 *       seed  RNG seed (default: current time)
 */
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "imports.h"

#define MAX_VAL   65536u /* 2^16, inclusive upper bound */
#define THRESHOLD 0.1

/* rand() may only provide 15 bits (e.g. on Windows), so combine two calls. */
static u32 rand30(void) {
  return ((u32)rand() << 15) ^ (u32)rand();
}

/* Uniform-ish random value in [0, 2^16]. */
static u32 rand_val(void) {
  return rand30() % (MAX_VAL + 1u);
}

static void make_random_matrix(Matrix m) {
  u64 N = (u64)(m.n) * (u64)(m.n);
  for (u64 i = 0; i < N; i++) {
    if (i % m.n == 0) {
      m.data[i] = 0;
      continue;
    }

    m.data[i] = (((f32)rand_val() / MAX_VAL) < THRESHOLD) ? rand_val() : INF;
  }
}

i32 main(const i32 argc, const char* argv[]) {
  u32 n = 0, count = 1;
  u32 seed = (u32)time(NULL);

  if (argc > 1) {
    n = strtol(argv[1], NULL, 10);

  } else {
    fprintf(stderr, "Enter N: "); /* stderr so it doesn't pollute redirected output */
    if (scanf("%" SCNu32, &n) != 1) {
      fprintf(stderr, "Invalid input\n");
      return EXIT_FAILURE;
    }
  }

  if (argc > 2) 
    count = strtol(argv[2], NULL, 10);
  if (argc > 3) 
    seed = (u32)strtoul(argv[3], NULL, 10);

  if (n <= 0 || count <= 0) {
    fprintf(stderr, "N and count must be positive i32egers\n");
    return EXIT_FAILURE;
  }

  srand(seed);

  u32* data = (u32*)malloc((u64)n * n * sizeof(u32));
  if (data == NULL) {
    fprintf(stderr, "Out of memory for %ux%u matrix\n", n, n);
    exit(EXIT_FAILURE);
  }
  Matrix m = (Matrix){ .n = n, .data = data };

  for (u32 k = 0; k < count; k++) {
    make_random_matrix(m);

    char filename[32];
    if (snprintf(filename, 32, "rand%u_input%u.bin", k, n) <= 0) {
      fprintf(stderr, "Could not create file name string");
      matrix_free(m);
      exit(EXIT_FAILURE);
    }

    FILE* out = fopen(filename, "wb");
    if (out == NULL) {
      fprintf(stderr, "Could not open matrix binary file to write new matrix, iteration %u\n", k);
      matrix_free(m);
      exit(EXIT_FAILURE);
    }

    matrix_write_bin(out, m);
    fclose(out);
  }

  matrix_free(m);

  return 0;
}
