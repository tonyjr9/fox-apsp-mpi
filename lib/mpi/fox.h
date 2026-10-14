#ifndef FOX_H
#define FOX_H

#include <mpi.h>

#include "imports.h"

typedef struct {
  MPI_Comm 
    comm,     // global comm
    row_comm,
    col_comm;

  i32 
    p,       // #processes
    my_rank; // global comm rank

  u16 
    q,       // q * q = processors
    my_row,
    my_col;

} GridData;


typedef struct {
    Matrix broadcast_A;
} FoxWorkspace;

void fox_grid_free(GridData* g);
/*
 * Validate input according to Fox's rules:
 * - Number os processses must be a perfect square;
 * - Sqrt(#processes) must divide matrix size;
 *
 * Only prints error if run by ROOT process.
*/
bool   fox_validate(u32 n, u16 p, int rank);
void   fox_setup_grid(GridData* g);

void   fox_scatter(
  const GridData* grid, 
  const Matrix global, 
  u32 n,
  Matrix A
);

Matrix fox_gather(
  const GridData* grid, 
  Matrix global, 
  u32 n,
  Matrix A
);

Matrix mpi_apsp(
  const GridData* grid, 
  Matrix input,
  u32 n
);

void fox(
  i32 n,
  const GridData* g,
  Matrix local_C, 
  const Matrix local_A, 
  Matrix local_B,
  FoxWorkspace* workspace);

#endif
