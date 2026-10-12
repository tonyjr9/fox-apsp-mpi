#include "imports.h"
#include "fox.h"

int main(void) {
  Matrix m = matrix_read(stdin);
  if (m.data == NULL)
    return 1;

  Matrix r;
  matrix_print(stdout, r);

  matrix_free(r);
  matrix_free(m);

  return 0;
}
