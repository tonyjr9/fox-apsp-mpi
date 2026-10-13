#include "minplus.h"

/*
 * Min-plus product: c[i][j] = min over k of (a[i][k] + b[k][j]).
 * Overwrites c, which must not share memory with a or b.
 */
void minplus_accumulate(const Matrix c, const Matrix a, const Matrix b) {
  if (c.n != a.n || a.n != b.n)
    return;

  const u32 n = c.n;

  /*
   * Loop order is i-k-j, not i-j-k: the inner loop walks along rows of
   * b and c (contiguous memory) instead of down a column of b (stride n).
   * This avoids a cache miss on nearly every access and lets the
   * compiler vectorize the inner loop.
   */
  for (u32 i = 0; i < n; i++) {
    // restrict: promises the compiler these pointers don't overlap, so it
    // can keep values in registers and vectorize. Breaking it is undefined behaviour.
    u32*       restrict out    = &c.data[i * n];   // row i of c
    const u32* restrict read_a = &a.data[i * n];   // row i of a

    // Row i starts as unreachable
    //for (
    //  u32 j = 0;
    //  j < n;
    //  out[j++] = INF
    //);

    // Relax row i of c through each intermediate k:
    // c[i][j] = min(c[i][j], a[i][k] + b[k][j])
    for (u32 k = 0; k < n; k++) {
      const u32 left = read_a[k];
      if (left >= INF) // no path i -> k, row k of b can't help
        continue;

      const u32* restrict read_b = &b.data[k * n];  // row k of b

      for (u32 j = 0; j < n; j++) {
        const u32 right = read_b[j];
        if (right >= INF)
          continue;
        const u32 candidate = left + right;

        out[j] = MIN(candidate, out[j]);
        // MIN(a, b) is defined as (a <= b ? a : b)
        // this can be vectorized by the compiler into a single CPU
        // SIMD instruction helping a lot in the speed of the process
        // this can only be achieved if the compiler knows for sure
        // it can vectorize, and by giving it restrict in the declaration
        // of the pointers above, it can trust us
      }
    }
  }
}

Matrix minplus_apsp(const Matrix m) {
  Matrix 
    d    = matrix_copy(m),
    next = matrix_create(m.n, 0);

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
    minplus_accumulate(next, d, d);

    u32* temp = next.data;
    next.data = d.data;
    d.data    = temp;
  }

  matrix_free(next);

  return d;
}
