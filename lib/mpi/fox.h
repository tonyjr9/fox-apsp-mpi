#ifndef FOX_H
#define FOX_H

#include <mpi.h>

#include "imports.h"

typedef struct
{
  int p;             // #processes
  MPI_Comm comm;     // global comm
  MPI_Comm row_comm;
  MPI_Comm col_comm;
  u16 q;             // q*q=#processors
  u16 my_row;
  u16 my_col;
  int my_rank;       // global comm rank
} GRID;


/*
 * Validate input according to Fox's rules:
 * - Number os processses must be a perfect square;
 * - Sqrt(#processes) must divide matrix size
*/
bool validate(u32 n, u16 p);

void setup(GRID* g);

void fox(const GRID* g,
         const Matrix local_A, Matrix local_B, Matrix local_C);
#endif
