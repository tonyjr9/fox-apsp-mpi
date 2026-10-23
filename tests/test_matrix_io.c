#include "matrix.h"
#include "matrix_io.h"
#include "test_util.h"

#include <string.h>

static Matrix read_str(const char *text) {
  FILE *f = tmpfile();
  fputs(text, f);
  rewind(f);
  Matrix m = matrix_read(f);
  fclose(f);
  return m;
}

static void print_str(const Matrix m, char *buf, size_t size) {
  FILE *f = tmpfile();
  matrix_print(f, m);
  rewind(f);
  size_t got = fread(buf, 1, size - 1, f);
  buf[got] = '\0';
  fclose(f);
}

static void test_assignment_example(void) {
  const char *in =
      "6\n0 2 0 5 0 0\n0 0 0 0 0 0\n0 2 0 0 0 5\n"
      "0 0 0 0 1 0\n3 9 3 0 0 0\n0 0 0 0 1 0\n";
  Matrix m = read_str(in);
  CHECK(m.n == 6);
  if (m.data == NULL)
      return;

  CHECK(MAT(m, 0, 1) == 2);
  CHECK(MAT(m, 0, 2) == INF);
  CHECK(MAT(m, 4, 1) == 9);
  for (uint32_t i = 0; i < 6; i++)
      CHECK(MAT(m, i, i) == 0);

  char out[256];
  print_str(m, out, sizeof out);
  CHECK(strcmp(out, in + 2) == 0);
  matrix_free(m);
}

static void test_tiny_sizes(void) {
  Matrix m1 = read_str("1\n0\n");
  CHECK(m1.n == 1 && MAT(m1, 0, 0) == 0);
  matrix_free(m1);

  Matrix m2 = read_str("2\n0 3\n4 0\n");
  CHECK(m2.n == 2 && MAT(m2, 0, 1) == 3 && MAT(m2, 1, 0) == 4);
  matrix_free(m2);
}

static void test_diagonal_forced_to_zero(void) {
  Matrix m = read_str("2\n5 1\n1 9\n");
  CHECK(m.data != NULL);
  if (m.data == NULL)
      return;
  CHECK(MAT(m, 0, 0) == 0 && MAT(m, 1, 1) == 0);
  matrix_free(m);
}

static void test_no_edges_prints_zeros(void) {
  Matrix m = read_str("3\n0 0 0\n0 0 0\n0 0 0\n");
  CHECK(m.data != NULL);
  if (m.data == NULL)
      return;
  CHECK(MAT(m, 0, 1) == INF);

  char out[64];
  print_str(m, out, sizeof out);
  CHECK(strcmp(out, "0 0 0\n0 0 0\n0 0 0\n") == 0);
  matrix_free(m);
}

static void test_print_inf_as_zero(void) {
  Matrix m = matrix_create(2, INF);
  MAT(m, 0, 0) = 0;
  MAT(m, 1, 1) = 0;
  MAT(m, 0, 1) = 5;
  char out[64];
  print_str(m, out, sizeof out);
  CHECK(strcmp(out, "0 5\n0 0\n") == 0);
  matrix_free(m);
}

static void test_malformed_input(void) {
  CHECK(read_str("").data == NULL);
  CHECK(read_str("abc\n").data == NULL);
  CHECK(read_str("0\n").data == NULL);
  CHECK(read_str("3\n0 1 2\n").data == NULL);
  CHECK(read_str("2\n0 1\n2").data == NULL);
}

static void test_bin_round_trip(void) {
  Matrix m = matrix_create(3, INF);
  MAT(m, 0, 0) = 0;
  MAT(m, 1, 1) = 0;
  MAT(m, 2, 2) = 0;
  MAT(m, 0, 1) = 7;
  MAT(m, 1, 2) = INF;
  MAT(m, 2, 0) = 65536;

  FILE *f = tmpfile();
  CHECK(f != NULL);
  if (f == NULL)
      return;
  matrix_write_bin(f, m);
  rewind(f);
  Matrix got = matrix_read_bin(f);
  CHECK(got.n == 3 && matrix_equal(got, m));
  matrix_free(got);
  fclose(f);
  matrix_free(m);
}

static void test_bin_malformed_input(void) {
  /* truncated: only the dimension */
  u32 n = 2;
  FILE *f = tmpfile();
  CHECK(f != NULL);
  if (f == NULL)
      return;
  fwrite(&n, sizeof(u32), 1, f);
  rewind(f);
  CHECK(matrix_read_bin(f).data == NULL);
  fclose(f);

  /* non-zero diagonal element */
  u32 bad[5] = { 2, 0, 3, 4, 5 };
  f = tmpfile();
  fwrite(bad, sizeof(u32), 5, f);
  rewind(f);
  CHECK(matrix_read_bin(f).data == NULL);
  fclose(f);
}

int main(void) {
  test_assignment_example();
  test_tiny_sizes();
  test_diagonal_forced_to_zero();
  test_no_edges_prints_zeros();
  test_print_inf_as_zero();
  test_malformed_input();
  test_bin_round_trip();
  test_bin_malformed_input();
  return test_report("test_matrix_io");
}
