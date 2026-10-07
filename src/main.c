#include <stdlib.h>
#include "imports.h"
#include "apsp_seq.h"

i32 main(const i32 argc, const char *argv[]) {
  if (argc < 2) {
    fprintf(
        stderr,
        "[ERROR]: invalid command line arguments - ./common path/to/test/file");
    return 1;
  }

  FILE *file = fopen(argv[1], "r");
  if (file == NULL) {
    fprintf(stderr, "[ERROR]: could not open file %s", argv[1]);
    return 1;
  }

  Matrix m = matrix_read(file);
  Matrix r = apsp(m);

  matrix_print(stdout, m);
  printf("\n");
  matrix_print(stdout, r);

  matrix_free(r);
  matrix_free(m);

  return 0;
}
