#include "fox.h"

#include <math.h>
#include <mpi.h>

bool fox_validate(u32 n, u16 p, int rank)
{
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

// https://www.dcc.fc.up.pt/~ricroc/aulas/1011/ppd/apontamentos/guide_mpi.pdf
void fox_setup_grid(GridData* g)
{
  int old_rank;
  int dims[2]    = { 0, 0 };
  int periods[2] = { 1, 1 };      // wrap-around in both dimensions
  int coords[2];
  int free_coords[2];

  MPI_Comm_size(MPI_COMM_WORLD, &g->p);
  MPI_Comm_rank(MPI_COMM_WORLD, &old_rank);

  g->q    = (u16)sqrt((double)g->p);
  dims[0] = dims[1] = g->q;
 
  // q x q periodic grid
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &g->comm);
 
  MPI_Comm_rank(g->comm, &(g->my_rank));
  MPI_Cart_coords(g->comm, g->my_rank, 2, coords);
  g->my_row  = (u16)coords[0];
  g->my_col  = (u16)coords[1];
 
  // row communicator: the second coordinate varies
  free_coords[0] = 0; free_coords[1] = 1;
  MPI_Cart_sub(g->comm, free_coords, &g->row_comm);
 
  // column communicator: the first coordinate varies
  free_coords[0] = 1; free_coords[1] = 0;
  MPI_Cart_sub(g->comm, free_coords, &g->col_comm);
}

void fox(i32 n,
         const GridData* grid,
         Matrix local_C,
         const Matrix local_A,
         Matrix local_B)
{
  Matrix temp_A;
  i32 step, bcast_root, source, dest;
  const u32 n_local = n/grid->q;
  const int block_elements = (int)((u64)n_local * n_local);
  const int tag = 69;

  /* Calculate ranks for the circular shift of B. */
  source = (grid->my_row + 1)           % grid->q;
  dest   = (grid->my_row + grid->q - 1) % grid->q;

  /* Reuse this storage for broadcasts from other columns. */
  temp_A = matrix_create(n_local, 0);

  for (step = 0; step < grid->q; step++) {
    bcast_root = (grid->my_row + step) % grid->q;

    if (bcast_root == grid->my_col) {
      MPI_Bcast(
        local_A.data, 
        block_elements, 
        MPI_UINT32_T,
        bcast_root, 
        grid->row_comm);
      minplus_accumulate(local_C, local_A, local_B);
    } else {
      MPI_Bcast(
        temp_A.data, 
        block_elements, 
        MPI_UINT32_T,
        bcast_root, 
        grid->row_comm);
      minplus_accumulate(local_C, temp_A, local_B);
    }
    MPI_Sendrecv_replace(
      local_B.data, 
      block_elements, 
      MPI_UINT32_T,
      dest, tag, source, tag, 
      grid->col_comm,
      MPI_STATUS_IGNORE);
  }

  matrix_free(temp_A);
}
