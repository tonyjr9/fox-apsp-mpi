#include "apsp_seq.h"

Matrix apsp(Matrix m) {
  Matrix r = matrix_create(m.n, 0);

  u32 n = m.n;
  for (u32 k = 0; k < n; k++) {
    for (u32 i = 0; i < n; i++) {
      for (u32 j = 0; j < n; j++) {
        MAT(r,i,j) = MIN(MAT(m, i, j), (MAT(m, i, k) + MAT(m, k, j));
      }
    }
  }

  return r;
}
