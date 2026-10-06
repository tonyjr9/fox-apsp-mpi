#include <stdlib.h>
#include "minplus.h"

i32 main(const i32 argc, const char* argv[]) {
  if (argc < 2) {
    fprintf(stderr, "[ERROR]: invalid command line arguments - ./common path/to/test/file");
    return 1;
  }

  FILE* file = fopen(argv[1], "r");
  if (file == NULL) {
    fprintf(stderr, "[ERROR]: could not open file %s", argv[1]);
    return 1;
  }

  Matrix m = matrix_read(file);
  minplus_solve(m);
  matrix_free(m);
  
  return 0;
}
