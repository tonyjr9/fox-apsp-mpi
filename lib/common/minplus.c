#include "minplus.h"

void minplus_matrix_reset(Matrix c) {
  for (u64 i = 0; i < (u64)c.n * c.n; i++) {
    c.data[i] = INF;
  }
}

Matrix minplus_multiply(const Matrix a, const Matrix b) {
  if (a.n != b.n)
    return MATRIX_NULL;

  Matrix c = matrix_create(a.n, 0);
  minplus_multiply_into(c, a, b);
  return c;
}

void minplus_multiply_into(
  Matrix c,
  const Matrix a,
  const Matrix b) {

  minplus_matrix_reset(c);
  minplus_relax(c, a, b);
}

#define BLK 128

/*
 * Min-plus accumulate: c[i][j] = min(c[i][j], min over k of (a[i][k] + b[k][j])).
 * c keeps its previous values, so calling this once per block of k (as Fox
 * does) builds the full product. c must not share memory with a or b.
 */
void minplus_relax(Matrix c, const Matrix a, const Matrix b) {
  if (c.n != a.n || a.n != b.n)
    return;

  const u32 n = c.n;

  /*
   * Tiled in BLK x BLK blocks so the touched slices of a, b and c stay in
   * cache. Inside a tile the order is i-k-j, not i-j-k: the inner loop walks
   * along rows of b and c (contiguous memory) instead of down a column of b
   * (stride n), which lets the compiler vectorize it.
   */
  for (u32 ii = 0; ii < n; ii += BLK) {
    const u32 i_end = MIN(ii + BLK, n);

    for (u32 kk = 0; kk < n; kk += BLK) {
      const u32 k_end = MIN(kk + BLK, n);

      for (u32 jj = 0; jj < n; jj += BLK) {
        const u32 j_end = MIN(jj + BLK, n);

        for (u32 i = ii; i < i_end; i++) {
          // restrict: promises the compiler these pointers don't overlap, so it
          // can keep values in registers and vectorize. Breaking it is undefined behaviour.
          u32*       restrict out    = &c.data[i * n];   // row i of c
          const u32* restrict read_a = &a.data[i * n];   // row i of a

          for (u32 k = kk; k < k_end; k++) {
            const u32 left = read_a[k];
            if (left >= INF) // no path i -> k, row k of b can't help
              continue;

            const u32* restrict read_b = &b.data[k * n];  // row k of b

            // MIN compiles to a single SIMD min instruction per vector.
            for (u32 j = jj; j < j_end; j++) {
              const u32 candidate = left + read_b[j];
              out[j] = MIN(candidate, out[j]);
            }
          }
        }
      }
    }
  }
}

Matrix minplus_apsp(const Matrix m) {
  Matrix
    d    = matrix_copy(m),
    next = matrix_create(m.n, INF);

  if (d.data == NULL || next.data == NULL) {
    matrix_free(d);
    matrix_free(next);
    return MATRIX_NULL;
  }

  for (
    u32 covered = 1;
    covered + 1 < d.n; // d.n is unsigned and for n == 0, it would cause an underflow
    covered <<= 1
  ) {
    minplus_multiply_into(next, d, d);

    if (matrix_equal(next, d))
      break;

    u32* temp = next.data;
    next.data = d.data;
    d.data    = temp;
  }

  matrix_free(next);

  return d;
}
