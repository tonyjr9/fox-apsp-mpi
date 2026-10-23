/*
 * Convert a text matrix on stdin to the binary format on stdout.
 *   ./binconv < tests/examples/input/input6 > input6.bin
 * Same conventions as the text format: off-diagonal 0 becomes INF (no edge);
 * the diagonal is 0.
 */
#include "imports.h"

i32 main(void) {
  Matrix m = matrix_read(stdin);
  if (m.data == NULL)
    return 1;

  matrix_write_bin(stdout, m);

  matrix_free(m);
  return 0;
}