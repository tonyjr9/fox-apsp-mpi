#include <math.h>
#include <stdio.h>

#include "foxgrid.h"

// q such that q * q == p, or 0 if p isn't a perfect square.
static int exact_sqrt(int p) {
  const int q = (int)lround(sqrt((double)p));
  return (q * q == p) ? q : 0;
}

bool fox_validate(u32 n, int p, int rank) {
  const int q = exact_sqrt(p);
  const char* reason = NULL;

  if (q == 0)
    reason = "number of processes must be a perfect square";
  else if (n % (u32)q != 0)
    reason = "matrix size isn't divisible by sqrt(#processes)";

  if (reason == NULL)
    return true;

  // The message on stdout matches tests/examples/output/output5 byte for byte.
  if (rank == 0) {
    fprintf(stderr, "[ERROR]: %s.\n", reason);
    printf("ERROR: Invalid configuration!");
  }
  return false;
}

// https://www.dcc.fc.up.pt/~ricroc/aulas/1011/ppd/apontamentos/guide_mpi.pdf
void fox_setup_grid(FoxGrid* g) {
  int dims[2]    = { 0, 0 };
  int periods[2] = { 1, 1 };      // wrap-around in both dimensions
  int coords[2];
  int free_coords[2];

  MPI_Comm_size(MPI_COMM_WORLD, &g->p);

  g->q    = exact_sqrt(g->p);
  dims[0] = dims[1] = g->q;

  // q x q periodic grid (reorder=0 keeps cartesian rank == MPI_COMM_WORLD rank,
  // so rank 0 of the grid is the process that read the input)
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &g->comm);

  MPI_Comm_rank(g->comm, &g->my_rank);
  MPI_Cart_coords(g->comm, g->my_rank, 2, coords);
  g->my_row = coords[0];
  g->my_col = coords[1];

  // row communicator: the second coordinate varies
  free_coords[0] = 0; free_coords[1] = 1;
  MPI_Cart_sub(g->comm, free_coords, &g->row_comm);

  // column communicator: the first coordinate varies
  free_coords[0] = 1; free_coords[1] = 0;
  MPI_Cart_sub(g->comm, free_coords, &g->col_comm);
}

void fox_grid_free(FoxGrid* g) {
  MPI_Comm_free(&g->row_comm);
  MPI_Comm_free(&g->col_comm);
  MPI_Comm_free(&g->comm);
}
