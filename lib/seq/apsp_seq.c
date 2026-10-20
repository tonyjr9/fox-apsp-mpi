#include "apsp_seq.h"
#include "minplus.h"

Matrix floyd_warshall(Matrix m) {
  const u32 n = m.n;
  Matrix r = matrix_copy(m);

  for (u32 k = 0; k < n; k++) {
    const u32* row_k = &r.data[k * n];

    for (u32 i = 0; i < n; i++) {
      const u32 ik = MAT(r, i, k);
      if (ik >= INF) // no path i -> k, row k can't help
        continue;

      // When i == k, row_i aliases row_k, but ik == 0 so nothing changes.
      u32* row_i = &r.data[i * n];
      for (u32 j = 0; j < n; j++)
        row_i[j] = MIN(ik + row_k[j], row_i[j]);
    }
  }

  return r;
}

Matrix apsp(Matrix m) { return minplus_apsp(m); }
