#include "matrix.h"
#include "minplus.h"
#include "apsp_seq.h"
#include "test_util.h"

#define X INF

static Matrix mat_from(const u32 n, const u32 *values) {
  Matrix m = matrix_create(n, 0);

  for (u32 i = 0; i < n * n; i++)
    m.data[i] = values[i];

  return m;
}

static const u32 W[16] = {
  0, 5, X, X,
  8, 0, 1, X,
  2, X, 0, 3,
  X, X, X, 0
};

static const u32 D2[16] = {
  0, 5, 6, X,
  3, 0, 1, 4,
  2, 7, 0, 3,
  X, X, X, 0
};

static const u32 D4[16] = {
  0, 5, 6, 9,
  3, 0, 1, 4,
  2, 7, 0, 3,
  X, X, X, 0
};

static void test_assignment_squaring(void) {
  Matrix w = mat_from(4, W);
  Matrix expected_d2 = mat_from(4, D2);
  Matrix expected_d4 = mat_from(4, D4);
  Matrix d2 = minplus_multiply(w, w);
  Matrix d4 = minplus_multiply(d2, d2);

  CHECK(matrix_equal(d2, expected_d2));
  CHECK(matrix_equal(d4, expected_d4));

  matrix_free(w);
  matrix_free(expected_d2);
  matrix_free(expected_d4);
  matrix_free(d2);
  matrix_free(d4);
}

static void test_identity(void) {
  const u32 identity[16] = {
    0, X, X, X,
    X, 0, X, X,
    X, X, 0, X,
    X, X, X, 0
  };
  Matrix w     = mat_from(4, W);
  Matrix id    = mat_from(4, identity);
  Matrix left  = minplus_multiply(w, id);
  Matrix right = minplus_multiply(id, w);

  CHECK(matrix_equal(left, w));
  CHECK(matrix_equal(right, w));

  matrix_free(w);
  matrix_free(id);
  matrix_free(left);
  matrix_free(right);
}

static void test_all_inf_no_overflow(void) {
  Matrix a = matrix_create(3, INF);
  Matrix c = minplus_multiply(a, a);

  for (u32 i = 0; i < 9; i++)
    CHECK(c.data[i] == INF);

  matrix_free(a);
  matrix_free(c);
}

static void test_relax_keeps_existing_minimum(void) {
  Matrix a = mat_from(4, W);
  Matrix c = matrix_create(4, 1);

  minplus_relax(c, a, a);

  CHECK(MAT(c, 0, 0) == 0);
  CHECK(MAT(c, 0, 1) == 1);
  CHECK(MAT(c, 0, 3) == 1);

  matrix_free(a);
  matrix_free(c);
}

static void test_relax_over_two_halves(void) {
  Matrix w = mat_from(4, W);
  Matrix expected = minplus_multiply(w, w);
  Matrix left = matrix_create(4, X);
  Matrix right = matrix_create(4, X);
  Matrix top = matrix_create(4, X);
  Matrix bottom = matrix_create(4, X);
  Matrix c = matrix_create(4, X);

  for (u32 i = 0; i < 4; i++) {
    for (u32 j = 0; j < 4; j++) {
      MAT(left, i, j) = (j < 2) ? MAT(w, i, j) : X;
      MAT(right, i, j) = (j < 2) ? X : MAT(w, i, j);
      MAT(top, i, j) = (i < 2) ? MAT(w, i, j) : X;
      MAT(bottom, i, j) = (i < 2) ? X : MAT(w, i, j);
    }
  }

  minplus_relax(c, left, top);
  minplus_relax(c, right, bottom);
  CHECK(matrix_equal(c, expected));

  matrix_free(w);
  matrix_free(expected);
  matrix_free(left);
  matrix_free(right);
  matrix_free(top);
  matrix_free(bottom);
  matrix_free(c);
}

/* Repeated squaring must agree with Floyd-Warshall on a random graph.
 * n = 150 is not a multiple of BLK, so partial tiles are exercised too. */
static void test_apsp_matches_floyd_warshall(void) {
  const u32 n = 150;
  Matrix w = matrix_create(n, INF);
  u32 seed = 12345;

  for (u32 i = 0; i < n; i++) {
    for (u32 j = 0; j < n; j++) {
      seed = seed * 1103515245u + 12345u;
      const u32 v = (seed >> 16) % 100;
      if (i == j)
        MAT(w, i, j) = 0;
      else if (v < 10) // ~10% density
        MAT(w, i, j) = 1 + (seed >> 8) % 1000;
    }
  }

  Matrix squared = apsp(w);
  Matrix fw      = floyd_warshall(w);

  CHECK(matrix_equal(squared, fw));

  matrix_free(w);
  matrix_free(squared);
  matrix_free(fw);
}

i32 main(void) {
  test_assignment_squaring();
  test_identity();
  test_all_inf_no_overflow();
  test_relax_keeps_existing_minimum();
  test_relax_over_two_halves();
  test_apsp_matches_floyd_warshall();
  return test_report("test_minplus");
}
