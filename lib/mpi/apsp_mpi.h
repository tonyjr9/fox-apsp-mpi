#ifndef FOX_MPI_H
#define FOX_MPI_H

#include "imports.h"
#include "foxgrid.h"

Matrix mpi_apsp(const FoxGrid* grid, Matrix input, u32 n);

#endif
