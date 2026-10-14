#include "matrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void die_oom(void) {
  fprintf(stderr, "matrix: out of memory\n");
  exit(EXIT_FAILURE);
}

Matrix matrix_create(const u32 n, const u32 fill) {
  Matrix m;
  u64 count = (u64)n * n;

  m.n = n;
  m.data = malloc(count * sizeof(u32));
  if (m.data == NULL && count > 0)
    die_oom();

  for (u64 i = 0; i < count; i++)
    m.data[i] = fill;

  return m;
}

void matrix_copy_into(Matrix dst, const Matrix src) {
  if (dst.n != src.n || dst.data == NULL || src.data == NULL)
    return;

  if (dst.data != src.data)
    memcpy(dst.data, src.data, (u64)dst.n * dst.n * sizeof(u32));
}

void matrix_free(Matrix m) {
  if (m.data != NULL)
    free(m.data);
}

Matrix matrix_copy(const Matrix m) {
  Matrix c = matrix_create(m.n, 0);
  memcpy(c.data, m.data, (u64)m.n * m.n * sizeof(u32));
  return c;
}

bool matrix_equal(const Matrix a, const Matrix b) {
  if (a.n != b.n)
    return false;
  return memcmp(a.data, b.data, (u64)a.n * a.n * sizeof(u32)) == 0;
}
