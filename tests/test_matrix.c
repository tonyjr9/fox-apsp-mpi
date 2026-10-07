#include "matrix.h"
#include "test_util.h"

static void test_create_and_access(void) {
  Matrix m = matrix_create(3, INF);
  CHECK(m.n == 3);
  for (uint32_t i = 0; i < 3; i++)
      for (uint32_t j = 0; j < 3; j++)
          CHECK(MAT(m, i, j) == INF);

  MAT(m, 1, 2) = 42;
  CHECK(MAT(m, 1, 2) == 42);
  CHECK(m.data[1 * 3 + 2] == 42);
  matrix_free(m);
}

static void test_copy_is_deep(void) {
  Matrix a = matrix_create(2, 7);
  Matrix b = matrix_copy(a);
  CHECK(b.data != a.data);
  CHECK(matrix_equal(a, b));

  MAT(b, 0, 1) = 9;
  CHECK(MAT(a, 0, 1) == 7);
  CHECK(!matrix_equal(a, b));
  matrix_free(a);
  matrix_free(b);
}

static void test_equal_size_mismatch(void) {
  Matrix a = matrix_create(2, 0);
  Matrix b = matrix_create(3, 0);
  CHECK(!matrix_equal(a, b));
  matrix_free(a);
  matrix_free(b);
}

static void test_free_null(void) {
  matrix_free(MATRIX_NULL);
  CHECK(1);
}

static void test_inf_does_not_overflow(void) {
  uint32_t sum = INF + INF;
  CHECK(sum >= INF);
}

int main(void) {
  test_create_and_access();
  test_copy_is_deep();
  test_equal_size_mismatch();
  test_free_null();
  test_inf_does_not_overflow();
  return test_report("test_matrix");
}
