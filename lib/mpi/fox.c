#include "fox.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

void fox_grid_free(FoxGrid* g){
  MPI_Comm_free(&g->row_comm);
  MPI_Comm_free(&g->col_comm);
  MPI_Comm_free(&g->comm);
}

bool fox_validate(u32 n, u16 p, int rank) {
  f32 q = sqrt(p);
  if (floorf(q) != q)
  {
    if (rank == 0) fprintf(stderr, "[ERROR]: number of processes must be a perfect square.\n");
    return false;
  }

  if (n%(u32)q != 0)
  {
    if (rank == 0) fprintf(stderr, "[ERROR]: matrix size isn't divisible by sqrt(#processes).\n");
    return false;
  }

  return true;
}

void fox_multiply(
  i32 n,
  const FoxGrid* grid,
  Matrix local_C,
  const Matrix local_A,
  Matrix local_B,
  FoxWorkspace* workspace) {

  i32 step, bcast_root, source, dest;
  const u32 n_local = n/grid->q;
  const int block_elements = (int)((u64)n_local * n_local);
  const int tag = 69;

  minplus_matrix_reset(local_C);
  /* Calculate ranks for the circular shift of B. */
  source = (grid->my_row + 1)           % grid->q;
  dest   = (grid->my_row + grid->q - 1) % grid->q;

  /* Reuse this storage for broadcasts from other columns. */

  for (step = 0; step < grid->q; step++) {
    bcast_root = (grid->my_row + step) % grid->q;
    memcpy(
        workspace->broadcast_A.data,
        local_A.data,
        block_elements * sizeof(u32)
    );
    MPI_Bcast(
      workspace->broadcast_A.data, 
      block_elements, 
      MPI_UINT32_T,
      bcast_root, 
      grid->row_comm);
    minplus_relax(local_C, workspace->broadcast_A, local_B);

    MPI_Sendrecv_replace(
      local_B.data, 
      block_elements, 
      MPI_UINT32_T,
      dest, tag, source, tag, 
      grid->col_comm,
      MPI_STATUS_IGNORE);
  }

}

/*
 * Distribute the N x N global matrix (root holds it) as contiguous
 * n_local x n_local blocks to every rank in MPI_COMM_WORLD.
 * Rank r receives grid block (r/q, r%q) packed row-major.
 */
void fox_scatter(
  const FoxGrid* grid,
  const Matrix global,
  u32 n,
  Matrix A) {

  const u32 n_local = (u32)(n/grid->q);
  const u32 count = n_local * n_local;
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  if (my_rank == 0) {
    u32* sendbuf = malloc((size_t)grid->p * count * sizeof(u32));
    if (sendbuf == NULL) {
      fprintf(stderr, "fox_scatter: out of memory\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    for (u32 r = 0; r < (u32)grid->p; r++) {
      const u32 row = r / grid->q;
      const u32 col = r % grid->q;
      u32* dst = &sendbuf[r * count];
      for (u32 i = 0; i < n_local; i++) {
        for (u32 j = 0; j < n_local; j++) {
          dst[i * n_local + j] =
              global.data[(row * n_local + i) * n + (col * n_local + j)];
        }
      }
    }
    MPI_Scatter(sendbuf, (int)count, MPI_UINT32_T,
                A.data, (int)count, MPI_UINT32_T,
                0, MPI_COMM_WORLD);
    free(sendbuf);
  } else {
    MPI_Scatter(NULL, 0, MPI_DATATYPE_NULL,
                A.data, (int)count, MPI_UINT32_T,
                0, MPI_COMM_WORLD);
  }
}

/*
 * Gather contiguous n_local x n_local blocks from every rank and write them
 * back into the N x N global matrix on root. Rank r holds grid block (r/q, r%q).
 */
Matrix fox_gather(
    const FoxGrid* grid,
    Matrix global,
    u32 n,
    Matrix A) {

  const u32 n_local = (u32)(n/grid->q);
  const u32 count   = n_local * n_local;
  int my_rank; 
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  if (my_rank == 0) {
    u32* recvbuf = malloc((size_t)grid->p * count * sizeof(u32));
    if (recvbuf == NULL) {
      fprintf(stderr, "fox_gather: out of memory\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Gather(A.data, (int)count, MPI_UINT32_T,
               recvbuf, (int)count, MPI_UINT32_T,
               0, MPI_COMM_WORLD);

    for (u32 r = 0; r < (u32)grid->p; r++) {
      const u32 row = r / grid->q;
      const u32 col = r % grid->q;
      const u32* src = &recvbuf[r * count];
      for (u32 i = 0; i < n_local; i++) {
        for (u32 j = 0; j < n_local; j++) {
          const u32 v = src[i * n_local + j];
          global.data[(row * n_local + i) * n + (col * n_local + j)] =
              (v == INF) ? INF : v;
        }
      }
    }
    free(recvbuf);
    return global;
  } else {
    MPI_Gather(A.data, (int)count, MPI_UINT32_T,
               NULL, 0, MPI_DATATYPE_NULL,
               0, MPI_COMM_WORLD);
  }
  return MATRIX_NULL;
}
