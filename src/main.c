#include "imports.h"
#include <string.h>

#ifdef USE_MPI
#include <mpi.h>

#include "fox.h"

int main(int argc, char **argv) {

  int n = 0, my_rank, nproc;

  Matrix input = MATRIX_NULL;
  
  FoxGrid grid;


  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  if (my_rank == 0) {
    input = matrix_read(stdin);
    n = (input.data == NULL) ? 0 : input.n;
  }
 
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (n <= 0 || !fox_validate((u32)n, (u16)nproc, my_rank)) {
    if (my_rank == 0) matrix_free(input);
    MPI_Finalize();
    return 1;
  }

  fox_setup_grid(&grid);

  // const u32 n_local = (u32)(n/grid.q);
  // const u32 count   = n_local * n_local;
  // Matrix A = matrix_create(n_local, 0);

  Matrix result = mpi_apsp(&grid, input, (u32)n);

  if (my_rank == 0)
    matrix_print(stdout, result);

  matrix_free(result);

  // Free memory
  fox_grid_free(&grid);
  MPI_Finalize();
  return 0;
}

#else // Sequencial 

#include "apsp_seq.h"
  
int main(void) {
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
