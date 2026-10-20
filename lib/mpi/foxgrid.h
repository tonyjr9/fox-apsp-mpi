#ifndef FOX_GRID_H
#define FOX_GRID_H

#include <mpi.h>

#include "imports.h"

typedef struct fox_grid FoxGrid;

struct fox_grid {
  MPI_Comm
    comm,     // q x q periodic cartesian comm
    row_comm, // my row; rank inside it == my_col
    col_comm; // my column; rank inside it == my_row

  int
    p,        // #processes
    q,        // q * q = processors
    my_rank,  // rank in comm
    my_row,
    my_col;
};

/*
 * Validate input according to Fox's rules:
 * - Number of processes must be a perfect square;
 * - Sqrt(#processes) must divide matrix size;
 *
 * Only the ROOT process prints: the reason on stderr and
 * "ERROR: Invalid configuration!" on stdout, as the assignment expects.
 */
bool fox_validate   (u32 n, int p, int rank);

void fox_setup_grid (FoxGrid* g);
void fox_grid_free  (FoxGrid* g);

#endif
