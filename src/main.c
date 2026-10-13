#include "imports.h"

#ifdef USE_MPI

#include "fox.h"

int main(void) {
  Matrix m = matrix_read(stdin);
  if (m.data == NULL)
    return 1;

  GridData grid;
  fox_setup_grid(&g);

  Matrix r = matrix_copy(m);
  matrix_print(stdout, r);

  matrix_free(r);
  matrix_free(m);
  return 0;
}

#else // Sequencial 

#include "apsp_seq.h"
  
int main(void){
  Matrix m = matrix_read(stdin);
  if (m.data == NULL)
    return 1;

  Matrix r = floyd_warshall(m);

  matrix_print(stdout, r);
  
  matrix_free(r);
  matrix_free(m);
  return 0;
}

#endif
