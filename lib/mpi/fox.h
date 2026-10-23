#ifndef FOX_H
#define FOX_H

#include "imports.h"
#include "foxgrid.h"

typedef struct fox_workspace FoxWorkspace;

struct fox_workspace {
  Matrix A_buf[2]; // A blocks broadcast along the row (current / next step)
  Matrix B_next;   // receives the next B block while the current one is in use
};

/*
 * Distribute / gather the N x N matrix (held by grid rank 0) across the
 * q x q grid, as directed by FOX_DIST (fox_config.h):
 * - row / colmajor: contiguous (n/q) x (n/q) blocks via MPI_Scatter /
 *   MPI_Gather; rank r owns block (r/q, r%q) or (r%q, r/q) respectively.
 * - ident: a single full-matrix broadcast; every process owns the whole
 *   matrix and nothing is ever communicated again.
 * `global` is only used on rank 0; fox_gather returns it there and
 * MATRIX_NULL on every other rank (ident returns a copy of the local
 * result on rank 0).
 */
void   fox_scatter  (const FoxGrid* g, const Matrix global, Matrix local);
Matrix fox_gather   (const FoxGrid* g, Matrix global, const Matrix local);

/*
 * Fox's algorithm: local_c = local_a (x) local_b over the whole grid, where
 * (x) is the min-plus product. Communication for step s+1 overlaps the
 * computation of step s. local_b is used as a buffer: its content is
 * undefined afterwards.
 * With FOX_DIST_IDENT this is a plain full-matrix min-plus product with no
 * communication: local_a is the whole matrix and workspace is unused.
 */
void   fox_multiply (const FoxGrid* g, Matrix local_c,
                     const Matrix local_a, Matrix local_b,
                     FoxWorkspace* workspace);

#endif
