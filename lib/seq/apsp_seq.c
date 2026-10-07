#include "apsp_seq.h"
#include "minplus.h"

Matrix floyd_warshall(Matrix m) {
  u32 n = m.n;
  Matrix r = matrix_copy(m);

  for (u32 k = 0; k < n; k++) {
    for (u32 i = 0; i < n; i++) {
      for (u32 j = 0; j < n; j++) {
        u32 tmp = MAT(r, i, k) + MAT(r, k, j);
        if (tmp < MAT(r, i, j)) MAT(r, i, j) = tmp;
      }
    }
  }

  return r;
}

Matrix apsp(Matrix m) { return minplus_apsp(m); }
