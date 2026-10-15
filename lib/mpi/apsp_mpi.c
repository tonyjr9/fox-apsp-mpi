#include "apsp_mpi.h"

#include "fox.h"

Matrix mpi_apsp(const FoxGrid* grid, Matrix input, u32 n) {
  const u32 n_local = (u32)(n/grid->q);
  Matrix A          = matrix_create(n_local, 0);
  Matrix B          = matrix_create(n_local, 0);
  Matrix C          = matrix_create(n_local, INF);

  FoxWorkspace workspace = { 
    .broadcast_A = matrix_create(n_local, 0) 
  };

  fox_scatter(grid, input, n, A);

  matrix_copy_into(B, A);

  for (u32 step = 1; step + 1 < n; step <<= 1) {
    fox_multiply((i32)n, grid, C, A, B, &workspace);
    matrix_copy_into(A, C);
    matrix_copy_into(B, A);
  }

  Matrix result = fox_gather(grid, input, n, A);

  matrix_free(A);
  matrix_free(B);
  matrix_free(C);
  matrix_free(workspace.broadcast_A);

  return result;
}

