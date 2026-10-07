#include "minplus.h"

Matrix minplus_multiply(const Matrix a, const Matrix b) {
  if (a.n != b.n)
    return MATRIX_NULL;

  Matrix c = matrix_create(a.n, INF);
  minplus_accumulate(c, a, b);
  return c;
}

void minplus_accumulate(const Matrix c, const Matrix a, const Matrix b) {
  if (c.n != a.n || a.n != b.n)
    return;

  for (u32 i = 0; i < c.n; i++) {
    for (u32 j = 0; j < c.n; j++) {
      for (u32 k = 0; k < c.n; k++) {
        u32 left = MAT(a, i, k);
        u32 right = MAT(b, k, j);

        if (left == INF || right == INF)
          continue;

        u32 candidate = left + right;
        if (candidate < MAT(c, i, j))
          MAT(c, i, j) = candidate;
      }
    }
  }
}

Matrix minplus_apsp(const Matrix m) {
  u32 covered = 1;
  Matrix d = matrix_copy(m);

  while (covered < d.n - 1) {
    Matrix next = minplus_multiply(d, d);
    matrix_free(d);
    d = next;
    covered *= 2;
  }

  return d;
}
