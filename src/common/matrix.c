#include "matrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void die_oom(void) {
    fprintf(stderr, "matrix: out of memory\n");
    exit(EXIT_FAILURE);
}

Matrix *matrix_create(u32 n, u32 fill) {
  Matrix *m = malloc(sizeof *m);

  if (m == NULL) 
    die_oom();

  size_t count = (size_t)n * n;
    m->n = n;
    m->data = malloc(count * sizeof *m->data);
    if (m->data == NULL && count > 0)
        die_oom();

    for (size_t k = 0; k < count; k++)
        m->data[k] = fill;

    return m;
}

void matrix_free(Matrix *m) {
  if (m == NULL)
    return;
  free(m->data);
  free(m);
}

Matrix *matrix_copy(const Matrix *m) {
  Matrix *c = matrix_create(m->n, 0);
  memcpy(c->data, m->data, (size_t)m->n * m->n * sizeof *m->data);
  return c;
}

bool matrix_equal(const Matrix *a, const Matrix *b) {
  if (a->n != b->n)
        return false;
  return memcmp(a->data, b->data,
          (size_t)a->n * a->n * sizeof *a->data) == 0;
}
