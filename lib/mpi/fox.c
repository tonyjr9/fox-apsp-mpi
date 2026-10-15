#include "fox.h"

#include <math.h>
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

void fox_scatter(
  const FoxGrid* grid,
  const Matrix global,
  u32 n,
  Matrix A) {

  const u32 n_local = (u32)(n/grid->q);
  const u32 count = n_local * n_local;
  minplus_matrix_reset(A);
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  if (my_rank == 0) {
    int dest = 0;
    Matrix tmp = matrix_create(n_local, 0);
    for (int row = 0; row != grid->q; ++row) {
      for (int col = 0; col != grid->q; ++col) {

        // Copy block (row, col) of the input matrix into tmp
        for (int i = 0; i != (int)n_local; ++i) {
          for (int j = 0; j != (int)n_local; ++j) {
            tmp.data[i * n_local + j] =
                global.data[(row * n_local + i) * n + (col * n_local + j)];
          }
        }

        if (dest == 0) {
          // MPI_Send to self would deadlock
          memcpy(A.data, tmp.data, (size_t)count * sizeof(u32));
        } else { 
          MPI_Send(tmp.data, (int)count, MPI_UINT32_T,
                   dest, 0, MPI_COMM_WORLD);
        }
        ++dest;
      }
    }
    matrix_free(tmp);
  } else {
    MPI_Recv(A.data, (int)count, MPI_UINT32_T,
             0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
}

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
    Matrix tmp = matrix_create(n_local, 0);
    int src = 0;

    for (int row = 0; row != grid->q; ++row) {
      for (int col = 0; col != grid->q; ++col) {

        if (src == 0) {
          memcpy(tmp.data, A.data, (size_t)count * sizeof(u32));
        } else {
          MPI_Recv(tmp.data, (int)count, MPI_UINT32_T,
                    src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        ++src;

        for (int i = 0; i != (int)n_local; ++i) {
          for (int j = 0; j != (int)n_local; ++j) {
            int dest = i * n_local + j;
            global.data[(row * n_local + i) * n + (col * n_local + j)] =
                  (tmp.data[dest] == INF) ? INF : tmp.data[dest];
          }
        }

      }
    }
    
    // matrix_print(stdout, global);

    matrix_free(tmp);
    return global;
  } else {
    MPI_Send(A.data, (int)count, MPI_UINT32_T,
              0, 0, MPI_COMM_WORLD);
  }
  return MATRIX_NULL;
}
 
#ifndef  
