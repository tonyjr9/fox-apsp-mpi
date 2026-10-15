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

void   fox_setup_grid  (FoxGrid* g);
