#include "imports.h"

#ifdef USE_MPI

#include "fox.h"
Matrix calc_fox(Matrix m){
    MPI_Init(&argc, &argv);
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);  // total number of processes
    if (!fox_validate(m.n, size))              // Validate the size of the matrix
      return 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // this process's id (0..size-1)
    

    MPI_Finalize();
}

int main(void) {
  Matrix m = matrix_read(stdin);
  if (m.data == NULL)
    return 1;

  Matrix r = calc_fox(m);
  matrix_print(stdout, r);

  matrix_free(r);
  matrix_free(m);
  return 0;
}

#else // Sequencial 

#include "apsp_seq.h"
  
int main(){
  Matrix m = matrix_read(stdin);
  if (m.data == NULL)
    return 1;

  Matrix r = apsp(m);

  matrix_print(stdout, r);
  
  matrix_free(r);
  matrix_free(m);
  return 0;
}
