#include "apsp_mpi.h"

#include "fox.h"
#include "fox_config.h"

Matrix mpi_apsp(const FoxGrid* grid, Matrix input, u32 n) {
  // Ident distribution computes on the full matrix everywhere: n_local == n.
  const u32 n_local = (FOX_DIST == FOX_DIST_IDENT) ? n : n / (u32)grid->q;
  Matrix A          = matrix_create(n_local, 0);
  Matrix B          = matrix_create(n_local, 0);
  Matrix C          = matrix_create(n_local, 0);

#if FOX_DIST == FOX_DIST_IDENT
  // No communication: fox_multiply is a plain local product, so the
  // receive buffers are never allocated (matrix_free is NULL-safe).
  FoxWorkspace workspace = {
    .A_buf  = { MATRIX_NULL, MATRIX_NULL },
    .B_next = MATRIX_NULL
  };
#else
  FoxWorkspace workspace = {
    .A_buf  = { matrix_create(n_local, 0), matrix_create(n_local, 0) },
    .B_next = matrix_create(n_local, 0)
  };
#endif

  fox_scatter(grid, input, A);

  matrix_copy_into(B, A);

  for (u32 covered = 1; covered + 1 < n; covered <<= 1) {
    fox_multiply(grid, C, A, B, &workspace);

    // if no matrix changed in all processes -> break
    int changed = !matrix_equal(A, C);
    MPI_Allreduce(MPI_IN_PLACE, &changed, 1, MPI_INT, MPI_LOR, grid->comm);
    if (!changed)
      break;

    // A <- C by swapping buffers; B needs its own copy because Fox uses it as a buffer.
    u32* temp = A.data;
    A.data    = C.data;
    C.data    = temp;
    matrix_copy_into(B, A);
  }

  Matrix result = fox_gather(grid, input, A);

  matrix_free(A);
  matrix_free(B);
  matrix_free(C);
  matrix_free(workspace.A_buf[0]);
  matrix_free(workspace.A_buf[1]);
  matrix_free(workspace.B_next);

  return result;
}
