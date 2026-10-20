#ifndef APSP_MPI_H
#define APSP_MPI_H

#include "imports.h"
#include "foxgrid.h"

/*
 * All-pairs shortest paths by repeated min-plus squaring, each product
 * computed with Fox's algorithm. Every rank of the grid must call it.
 * `input` is the full matrix on grid rank 0 (ignored elsewhere); the result is
 * written back into it and returned there, other ranks get MATRIX_NULL.
 */
Matrix mpi_apsp(const FoxGrid* grid, Matrix input, u32 n);

#endif
