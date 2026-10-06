#include "matrix_io.h"

#include <inttypes.h>

Matrix matrix_read(FILE *in) {
  u32 n;
  if (fscanf(in, "%" SCNu32, &n) != 1 || n == 0) {
    fprintf(stderr, "matrix_read: missing or invalid dimension N\n");
    return NULL;
  }

  Matrix m = matrix_create(n, INF);

  for (u32 i = 0; i < n; i++) {
    for (u32 j = 0; j < n; j++) {
      u32 v;

      if (fscanf(in, "%" SCNu32, &v) != 1) {
        fprintf(
          stderr, 
          "matrix_read: expected %" PRIu32 "x%" PRIu32
          " values, input ended at (%" PRIu32 ",%" PRIu32 ")\n",
          n, n, i, j
        );
        matrix_free(m);
        return NULL;
      }

      if (i == j)
        MAT(m, i, j) = 0;
      else
        MAT(m, i, j) = (v == 0) ? INF : v;
    }
  }

  return m;
}

void matrix_print(FILE *out, const Matrix m) {
  for (u32 i = 0; i < m->n; i++) {
    for (u32 j = 0; j < m->n; j++) {
      u32 v = MAT(m, i, j);
      if (v == INF)
        v = 0;

      fprintf(
        out, 
        j + 1 < m->n ? "%" PRIu32 " " : "%" PRIu32, 
        v
      );
    }

    fputc('\n', out);
  }
}
