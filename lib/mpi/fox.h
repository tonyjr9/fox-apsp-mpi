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
 * Distribute / gather the N x N matrix (held by grid rank 0) as contiguous
 * n_local x n_local blocks across the q x q grid using MPI_Scatter / MPI_Gather.
 * Grid rank r owns block (r/q, r%q). `global` is only used on rank 0;
 * fox_gather returns it there and MATRIX_NULL on every other rank.
 */
void   fox_scatter  (const FoxGrid* g, const Matrix global, Matrix local);
Matrix fox_gather   (const FoxGrid* g, Matrix global, const Matrix local);

/*
 * Fox's algorithm: local_c = local_a (x) local_b over the whole grid, where
 * (x) is the min-plus product. Communication for step s+1 overlaps the
 * computation of step s. local_b is used as a buffer: its content is
 * undefined afterwards.
 */
void   fox_multiply (const FoxGrid* g, Matrix local_c,
                     const Matrix local_a, Matrix local_b,
                     FoxWorkspace* workspace);

/*
 * Shared variant (the fox_shared binary). Every process keeps a full copy of
 * the matrix, stored block by block in rank order: block (r/q, r%q) is the
 * r-th contiguous (n/q) x (n/q) piece. Knowing that layout, each process reads
 * the A and B blocks of every Fox step straight from memory, so the multiply
 * itself communicates nothing; the only messages are the initial broadcast
 * and one all-gather of the new C blocks per squaring.
 *
 * fox_block           view (no copy) of block (row, col) of the shared copy
 * fox_share           ROOT packs global into the block layout, then broadcasts it
 * fox_unshare         ROOT unpacks the shared copy back into global and returns
 *                     it; MATRIX_NULL on every other rank
 * fox_multiply_shared local_c = (my row of blocks) (x) (my column of blocks)
 * fox_allgather       every rank's local block into its place in shared
 */
Matrix fox_block           (const FoxGrid* g, const Matrix shared, i32 row, i32 col);
void   fox_share           (const FoxGrid* g, const Matrix global, Matrix shared);
Matrix fox_unshare         (const FoxGrid* g, Matrix global, const Matrix shared);
void   fox_multiply_shared (const FoxGrid* g, Matrix local_c, const Matrix shared);
void   fox_allgather       (const FoxGrid* g, Matrix shared, const Matrix local);

#endif
