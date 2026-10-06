#include "matrix.h"
#include "minplus.h"
#include "test_util.h"

#define X INF

static const uint32_t W[16]  = { 0, 5, X, X,   X, 0, 1, X,
                               8, X, 0, 3,   2, X, X, 0 };
static const uint32_t D2[16] = { 0, 5, 6, X,   9, 0, 1, 4,
                               5,13, 0, 3,   2, 7, X, 0 };
static const uint32_t D4[16] = { 0, 5, 6, 9,   6, 0, 1, 4,
                               5,10, 0, 3,   2, 7, 8, 0 };

static void test_assignment_squaring(void) {
  Matrix *w  = mat_from(4, W);
  Matrix *e2 = mat_from(4, D2);
  Matrix *e4 = mat_from(4, D4);

  Matrix *d2 = minplus_multiply(w, w);
  CHECK(matrix_equal(d2, e2));

  Matrix *d4 = minplus_multiply(d2, d2);
  CHECK(matrix_equal(d4, e4));

  matrix_free(w);  matrix_free(e2);  matrix_free(e4);
  matrix_free(d2); matrix_free(d4);
}

static void test_identity(void) {
  const uint32_t I[16] = { 0, X, X, X,   X, 0, X, X,
                           X, X, 0, X,   X, X, X, 0 };
  Matrix *w = mat_from(4, W);
  Matrix *i = mat_from(4, I);

  Matrix *r1 = minplus_multiply(w, i);
  Matrix *r2 = minplus_multiply(i, w);
  CHECK(matrix_equal(r1, w));
  CHECK(matrix_equal(r2, w));

  matrix_free(w);  matrix_free(i);
  matrix_free(r1); matrix_free(r2);
}

static void test_all_inf_no_overflow(void) {
  Matrix *a = matrix_create(3, INF);
  Matrix *c = minplus_multiply(a, a);
  for (uint32_t k = 0; k < 9; k++)
      CHECK(c->data[k] == INF);
  matrix_free(a);
  matrix_free(c);
}

static void test_accumulate_keeps_existing_minimum(void) {
  Matrix *a = mat_from(4, W);
  Matrix *c = matrix_create(4, 1);

  minplus_accumulate(c, a, a);
  CHECK(MAT(c, 0, 0) == 0);
  CHECK(MAT(c, 0, 1) == 1);
  CHECK(MAT(c, 0, 3) == 1);

  matrix_free(a);
  matrix_free(c);
}

static void test_accumulate_over_two_halves(void) {
  Matrix *w = mat_from(4, W);
  Matrix *expected = minplus_multiply(w, w);

  Matrix *left = matrix_create(4, X), *right = matrix_create(4, X);
  for (uint32_t i = 0; i < 4; i++)
      for (uint32_t j = 0; j < 4; j++) {
          MAT(left,  i, j) = (j < 2) ? MAT(w, i, j) : X;
          MAT(right, i, j) = (j < 2) ? X : MAT(w, i, j);
      }
  Matrix *top = matrix_create(4, X), *bottom = matrix_create(4, X);
  for (uint32_t i = 0; i < 4; i++)
      for (uint32_t j = 0; j < 4; j++) {
          MAT(top,    i, j) = (i < 2) ? MAT(w, i, j) : X;
          MAT(bottom, i, j) = (i < 2) ? X : MAT(w, i, j);
      }

  Matrix *c = matrix_create(4, X);
  minplus_accumulate(c, left,  top);
  minplus_accumulate(c, right, bottom);
  CHECK(matrix_equal(c, expected));

  matrix_free(w); matrix_free(expected);
  matrix_free(left); matrix_free(right);
  matrix_free(top);  matrix_free(bottom);
  matrix_free(c);
}

int main(void) {
  test_assignment_squaring();
  test_identity();
  test_all_inf_no_overflow();
  test_accumulate_keeps_existing_minimum();
  test_accumulate_over_two_halves();
  return test_report("test_minplus");
}
