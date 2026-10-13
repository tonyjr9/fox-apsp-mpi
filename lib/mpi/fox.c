#include "fox.h"

#include <math.h>
#include <mpi.h>

bool fox_validate(u32 n, u16 p)
{
  f32 q = sqrt(p);
  if (floorf(q) != q)
  {
    fprintf(stderr, "[ERROR]: number of processes must be a perfect square.\n");
    return false;
  }

  if (n%(u32)q != 0)
  {
    fprintf(stderr, "[ERROR]: matrix size isn't divisible by sqrt(#processes).\n");
    return false;
  }

  return true;
}

// https://www.dcc.fc.up.pt/~ricroc/aulas/1011/ppd/apontamentos/guide_mpi.pdf
void fox_setup_grid(GRID* g)
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
