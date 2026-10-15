#ifndef FOX_H
#define FOX_H

#include <mpi.h>

#include "imports.h"

typedef struct fox_grid FoxGrid;

struct fox_grid {
  MPI_Comm
    comm,     // global comm
    row_comm,
    col_comm;

  i32
    p,        // #processes
    my_rank;  // global comm rank

  u16
    q,        // q * q = processors
    my_row,
    my_col;
};

typedef struct fox_workspace FoxWorkspace;

struct fox_workspace {
  Matrix broadcast_A;
};

void   fox_grid_free   (FoxGrid* g);
/*
 * Validate input according to Fox's rules:
 * - Number os processses must be a perfect square;
 * - Sqrt(#processes) must divide matrix size;
 *
 * Only prints error if run by ROOT process.
*/
bool   fox_validate    (u32 n, u16 p, int rank);
void   fox_setup_grid  (FoxGrid* g);

void   fox_scatter     (const FoxGrid* g, const Matrix global, u32 n, Matrix a);
Matrix fox_gather      (const FoxGrid* g, Matrix global, u32 n, Matrix a);
Matrix mpi_apsp        (const FoxGrid* g, Matrix input, u32 n);
void   fox_multiply    (i32 n, const FoxGrid* g, Matrix local_c,
                        const Matrix local_a, Matrix local_b,
                        FoxWorkspace* workspace);

#endif