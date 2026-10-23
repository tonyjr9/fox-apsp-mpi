#include "matrix_io.h"

#include <inttypes.h>

Matrix matrix_read(FILE* in) {
  u32 n;
  if (fscanf(in, "%" SCNu32, &n) != 1 || n == 0) {
    fprintf(stderr, "matrix_read: missing or invalid dimension N\n");
    return MATRIX_NULL;
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
        return MATRIX_NULL;
      }

      if (i == j)
        MAT(m, i, j) = 0;
      else
        MAT(m, i, j) = (v == 0) ? INF : v;
    }
  }

  return m;
}


Matrix matrix_read_bin(FILE* in) {
  if (in == NULL) {
    fprintf(stderr, "matrix_read_bin: in file ptr is NULL\n");
    return MATRIX_NULL;
  }

  u32 n;
  if (fread(&n, sizeof(u32), 1, in) != 1 || n == 0) {
    fprintf(stderr, "matrix_read_bin: missing or invalid dimension N\n");
    return MATRIX_NULL;
  }

  Matrix m = matrix_create(n, INF);

  u64 n_read = fread(m.data, sizeof(u32), (u64)n * n, in);
  if (n_read != (u64)n * n) {
    fprintf(stderr, "matrix_read_bin: failed to read N * N values from the matrix file\n");
    matrix_free(m);
    return MATRIX_NULL;
  }

  for (u32 i = 0; i < n; i++) {
    if (m.data[i * n + i] != 0) {
      fprintf(stderr, "matrix_read_bin: input matrix file has (i, i) values different from 0\n");
      matrix_free(m);
      return MATRIX_NULL;
    }
  }

  return m;
}

void matrix_print(FILE* out, const Matrix m) {
  for (u32 i = 0; i < m.n; i++) {
    for (u32 j = 0; j < m.n; j++) {
      u32 v = MAT(m, i, j);
      if (v == INF)
        v = 0;

      fprintf(
        out, 
        j + 1 < m.n ? "%" PRIu32 " " : "%" PRIu32, 
        v
      );
    }

    fputc('\n', out);
  }
}

void matrix_write_bin(FILE* out, const Matrix m) {
  if (out == NULL) {
    fprintf(stderr, "matrix_write_bin: out file ptr is NULL\n");
    return;
  }

  if (fwrite(&m.n, sizeof(u32), 1, out) != 1 || m.n == 0) {
    fprintf(stderr, "matrix_write_bin: failed to write or invalid dimension N\n");
    return;
  }

  u64 n_write = fwrite(m.data, sizeof(u32), (u64)m.n * m.n, out);
  if (n_write != (u64)m.n * m.n) {
    fprintf(stderr, "matrix_write_bin: failed to write N * N values to the matrix file\n");
    return;
  }
}
