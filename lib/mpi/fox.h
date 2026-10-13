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

} GRID;


/*
 * Validate input according to Fox's rules:
 * - Number os processses must be a perfect square;
 * - Sqrt(#processes) must divide matrix size
*/
bool fox_validate(u32 n, u16 p);

void fox_setup_grid(GRID* g);

void fox(const GRID* g,
         const Matrix local_A, Matrix local_B, Matrix local_C);

#endif
