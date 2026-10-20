#include "fox.h"

#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#define ROOT 0

void fox_multiply(
  const FoxGrid* grid,
  Matrix local_C,
  const Matrix local_A,
  Matrix local_B,
  FoxWorkspace* workspace) {

  const int block_elements = (int)((u64)local_A.n * local_A.n);
  const int tag = 69;

  /* The process in column (my_row + step) of this row broadcasts its A block,
   * straight from local_A; the others receive into A_buf[step % 2]. */
  const int first_root = grid->my_row;
  Matrix A_cur = (first_root == grid->my_col) ? local_A : workspace->A_buf[0];
  MPI_Bcast(A_cur.data, block_elements, MPI_UINT32_T, first_root, grid->row_comm);

  /* Ranks in col_comm for the circular upward shift of B. */
  const int source = (grid->my_row + 1)           % grid->q;
  const int dest   = (grid->my_row + grid->q - 1) % grid->q;

  minplus_matrix_reset(local_C);

  Matrix B_cur = local_B;
  Matrix B_nxt = workspace->B_next;

  for (int step = 0; step < grid->q; step++) {
    MPI_Request req[3];
    int nreq = 0;
    Matrix A_nxt = MATRIX_NULL;

    /* Start the transfers for step + 1; the last step needs none. */
    if (step + 1 < grid->q) {
      const int bcast_root = (grid->my_row + step + 1) % grid->q;
      A_nxt = (bcast_root == grid->my_col) ? local_A : workspace->A_buf[(step + 1) % 2];

      MPI_Irecv(
        B_nxt.data,
        block_elements,
        MPI_UINT32_T,
        source, tag,
        grid->col_comm,
        &req[nreq++]);

      MPI_Isend(
        B_cur.data,
        block_elements,
        MPI_UINT32_T,
        dest,   tag,
        grid->col_comm,
        &req[nreq++]);

      MPI_Ibcast(
        A_nxt.data,
        block_elements,
        MPI_UINT32_T,
        bcast_root,
        grid->row_comm,
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
      const u32* src = block_origin(global, n_local, r / grid->q, r % grid->q);
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
    const u32* src = &recvbuf[r * count];
    u32*       dst = block_origin(global, n_local, r / grid->q, r % grid->q);

    for (u32 i = 0; i < n_local; i++)
      memcpy(&dst[(u64)i * global.n], &src[i * n_local], n_local * sizeof(u32));
  }

  free(recvbuf);
  return global;
}
