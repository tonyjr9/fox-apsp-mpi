#include "apsp_mpi.h"

#include "fox.h"

#ifdef FOX_SHARED

// Shared variant (make fox_shared): see fox.h.
Matrix mpi_apsp(const FoxGrid* grid, Matrix input, u32 n) {
  const u32 n_local = n / (u32)grid->q;
  Matrix shared     = matrix_create(n, 0);      // every block, in rank order
  Matrix C          = matrix_create(n_local, 0);

  fox_share(grid, input, shared);

  for (u32 covered = 1; covered + 1 < n; covered <<= 1) {
    fox_multiply_shared(grid, C, shared);

    // if no matrix changed in all processes -> break
    const Matrix A = fox_block(grid, shared, grid->my_row, grid->my_col);
    i32 changed    = !matrix_equal(A, C);
    MPI_Allreduce(MPI_IN_PLACE, &changed, 1, MPI_INT, MPI_LOR, grid->comm);
    if (!changed)
      break;

    // A <- C everywhere: each process's new block goes into every copy.
    fox_allgather(grid, shared, C);
  }

  Matrix result = fox_unshare(grid, input, shared);

  matrix_free(shared);
  matrix_free(C);

  return result;
}

#else

Matrix mpi_apsp(const FoxGrid* grid, Matrix input, u32 n) {
  const u32 n_local = n / (u32)grid->q;
  Matrix A          = matrix_create(n_local, 0);
  Matrix B          = matrix_create(n_local, 0);
  Matrix C          = matrix_create(n_local, 0);

  FoxWorkspace workspace = {
    .A_buf  = { matrix_create(n_local, 0), matrix_create(n_local, 0) },
    .B_next = matrix_create(n_local, 0)
  };

  fox_scatter(grid, input, A);

  matrix_copy_into(B, A);

  for (u32 covered = 1; covered + 1 < n; covered <<= 1) {
    fox_multiply(grid, C, A, B, &workspace);

    // if no matrix changed in all processes -> break
    i32 changed = !matrix_equal(A, C);
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

#endif // FOX_SHARED
