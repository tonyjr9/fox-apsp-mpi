#include "fox.h"

#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "fox_config.h"

#define ROOT 0

#if FOX_DIST == FOX_DIST_IDENT
/*
 * No distribution. The input only exists on the root, so the scatter is a
 * single full-matrix broadcast into the (already sized) local buffer; after
 * that, nothing is ever sent.
 */
void fox_scatter(const FoxGrid* grid, const Matrix global, Matrix local) {
  const u64 count = (u64)local.n * local.n;

  if (grid->my_rank == ROOT)
    memcpy(local.data, global.data, count * sizeof(u32));

  MPI_Bcast(local.data, (int)count, MPI_UINT32_T, ROOT, grid->comm);
}

Matrix fox_gather(const FoxGrid* grid, Matrix global, const Matrix local) {
  (void)global;

  if (grid->my_rank != ROOT)
    return MATRIX_NULL;

  return matrix_copy(local);
}

/* No communication: the full product is a plain min-plus multiply. */
void fox_multiply(
  const FoxGrid* grid,
  Matrix local_C,
  const Matrix local_A,
  Matrix local_B,
  FoxWorkspace* workspace) {

  (void)grid;
  (void)local_B;
  (void)workspace;

  minplus_multiply_into(local_C, local_A, local_A);
}

#else

static u32* alloc_or_abort(u64 count, const char* who) {
  u32* buf = malloc(count * sizeof(u32));
  if (buf == NULL) {
    fprintf(stderr, "%s: out of memory\n", who);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  return buf;
}

/* Top-left element of grid block (row, col) inside the global matrix. */
static u32* block_origin(const Matrix global, u32 n_local, int row, int col) {
  return &global.data[(u64)row * n_local * global.n + (u64)col * n_local];
}

void fox_multiply(
  const FoxGrid* grid,
  Matrix local_C,
  const Matrix local_A,
  Matrix local_B,
  FoxWorkspace* workspace) {

  const int block_elements = (int)((u64)local_A.n * local_A.n);
  const int tag = 69;

  /*
   * Which comm the A broadcast and the B circulation run over, and where in
   * it I sit. Row distribution (default) and column-major distribution
   * mirror each other: rank (r, c) holds block (r, c) or block (c, r), so
   * the stepped coordinate and the "am I the bcast root" position swap.
   */
#if FOX_DIST == FOX_DIST_COLMAJOR
  const int  a_pos  = grid->my_row;   // my rank in a_comm
  const int  a_base = grid->my_col;   // the coordinate that steps with the iterations
  MPI_Comm   a_comm = grid->col_comm;
  const int  b_src  = (grid->my_col + 1)           % grid->q;
  const int  b_dst  = (grid->my_col + grid->q - 1) % grid->q;
  MPI_Comm   b_comm = grid->row_comm;
#else // FOX_DIST_ROW
  const int  a_pos  = grid->my_col;
  const int  a_base = grid->my_row;
  MPI_Comm   a_comm = grid->row_comm;
  const int  b_src  = (grid->my_row + 1)           % grid->q;
  const int  b_dst  = (grid->my_row + grid->q - 1) % grid->q;
  MPI_Comm   b_comm = grid->col_comm;
#endif

  /* The process at a_base of this comm broadcasts its A block, straight from
   * local_A; the others receive into A_buf[step % 2]. */
  Matrix A_cur = (a_base == a_pos) ? local_A : workspace->A_buf[0];
  MPI_Bcast(A_cur.data, block_elements, MPI_UINT32_T, a_base, a_comm);

  /* Ranks in b_comm for the circular shift of B, one step per iteration. */

  minplus_matrix_reset(local_C);

  Matrix B_cur = local_B;
  Matrix B_nxt = workspace->B_next;

  for (int step = 0; step < grid->q; step++) {
    MPI_Request req[3];
    int nreq = 0;
    Matrix A_nxt = MATRIX_NULL;

    /* Start the transfers for step + 1; the last step needs none. */
    if (step + 1 < grid->q) {
      const int bcast_root = (a_base + step + 1) % grid->q;
      A_nxt = (bcast_root == a_pos) ? local_A : workspace->A_buf[(step + 1) % 2];

      MPI_Irecv(
        B_nxt.data,
        block_elements,
        MPI_UINT32_T,
        b_src, tag,
        b_comm,
        &req[nreq++]);

      MPI_Isend(
        B_cur.data,
        block_elements,
        MPI_UINT32_T,
        b_dst,   tag,
        b_comm,
        &req[nreq++]);

      MPI_Ibcast(
        A_nxt.data,
        block_elements,
        MPI_UINT32_T,
        bcast_root,
        a_comm,
        &req[nreq++]);
    }

    /* Runs while the transfers are in flight; they only read A_cur and B_cur. */
    minplus_relax(local_C, A_cur, B_cur);

    MPI_Waitall(nreq, req, MPI_STATUSES_IGNORE);

    Matrix temp = B_cur;
    B_cur = B_nxt;
    B_nxt = temp;
    A_cur = A_nxt;
  }
}

void fox_scatter(
  const FoxGrid* grid,
  const Matrix global,
  Matrix local) {

  const u32 n_local = local.n;
  const u64 count   = (u64)n_local * n_local;
  u32* sendbuf      = NULL;

  /* A block's rows are global.n apart, so pack each block contiguously, in rank order. */
  if (grid->my_rank == ROOT) {
    sendbuf = alloc_or_abort(grid->p * count, "fox_scatter");

    for (int r = 0; r < grid->p; r++) {
#if FOX_DIST == FOX_DIST_COLMAJOR
      const u32* src = block_origin(global, n_local, r % grid->q, r / grid->q);
#else // FOX_DIST_ROW
      const u32* src = block_origin(global, n_local, r / grid->q, r % grid->q);
#endif
      u32*       dst = &sendbuf[r * count];

      for (u32 i = 0; i < n_local; i++)
        memcpy(&dst[i * n_local], &src[(u64)i * global.n], n_local * sizeof(u32));
    }
  }

  /* Rank r receives piece r of sendbuf. The send arguments only matter on ROOT. */
  MPI_Scatter(sendbuf,    (int)count, MPI_UINT32_T,
              local.data, (int)count, MPI_UINT32_T,
              ROOT, grid->comm);

  free(sendbuf);
}

Matrix fox_gather(
  const FoxGrid* grid,
  Matrix global,
  const Matrix local) {

  const u32 n_local = local.n;
  const u64 count   = (u64)n_local * n_local;
  u32* recvbuf      = NULL;

  if (grid->my_rank == ROOT)
    recvbuf = alloc_or_abort(grid->p * count, "fox_gather");

  /* ROOT receives every rank's block, in rank order. The receive arguments only matter on ROOT. */
  MPI_Gather(local.data, (int)count, MPI_UINT32_T,
             recvbuf,    (int)count, MPI_UINT32_T,
             ROOT, grid->comm);

  if (grid->my_rank != ROOT)
    return MATRIX_NULL;

  for (int r = 0; r < grid->p; r++) {
#if FOX_DIST == FOX_DIST_COLMAJOR
    u32* dst = block_origin(global, n_local, r % grid->q, r / grid->q);
#else // FOX_DIST_ROW
    u32* dst = block_origin(global, n_local, r / grid->q, r % grid->q);
#endif
    const u32* src = &recvbuf[r * count];

    for (u32 i = 0; i < n_local; i++)
      memcpy(&dst[(u64)i * global.n], &src[i * n_local], n_local * sizeof(u32));
  }

  free(recvbuf);
  return global;
}

#endif // FOX_DIST != FOX_DIST_IDENT
