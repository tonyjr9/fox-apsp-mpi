#ifndef FOX_H
#define FOX_H

#include <mpi.h>

#include "imports.h"
#include "foxgrid.h"

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

void   fox_scatter     (const FoxGrid* g, const Matrix global, u32 n, Matrix a);
Matrix fox_gather      (const FoxGrid* g, Matrix global, u32 n, Matrix a);
Matrix mpi_apsp        (const FoxGrid* g, Matrix input, u32 n);
void   fox_multiply    (i32 n, const FoxGrid* g, Matrix local_c,
                        const Matrix local_a, Matrix local_b,
                        FoxWorkspace* workspace);

#endif
