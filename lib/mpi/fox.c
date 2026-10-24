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

  const i32 block_elements = (i32)((u64)local_A.n * local_A.n);
  const i32 tag = 69;

  /* The process in column (my_row + step) of this row broadcasts its A block,
   * straight from local_A; the others receive into A_buf[step % 2]. */
  const i32 first_root = grid->my_row;
  Matrix A_cur = (first_root == grid->my_col) ? local_A : workspace->A_buf[0];
  MPI_Bcast(A_cur.data, block_elements, MPI_UINT32_T, first_root, grid->row_comm);

  /* Ranks in col_comm for the circular upward shift of B. */
  const i32 source = (grid->my_row + 1)           % grid->q;
  const i32 dest   = (grid->my_row + grid->q - 1) % grid->q;

  minplus_matrix_reset(local_C);

  Matrix B_cur = local_B;
  Matrix B_nxt = workspace->B_next;

  for (i32 step = 0; step < grid->q; step++) {
    MPI_Request req[3];
    i32 nreq = 0;
    Matrix A_nxt = MATRIX_NULL;

    /* Start the transfers for step + 1; the last step needs none. */
    if (step + 1 < grid->q) {
      const i32 bcast_root = (grid->my_row + step + 1) % grid->q;
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
static u32* block_origin(const Matrix global, u32 n_local, i32 row, i32 col) {
  return &global.data[(u64)row * n_local * global.n + (u64)col * n_local];
}

/*
 * Copy the global matrix into dst block by block, in rank order: block (r/q, r%q)
 * becomes the r-th n_local x n_local piece of dst. A block's rows are global.n
 * apart in the global matrix; in dst each block is contiguous.
 */
static void pack_blocks(const FoxGrid* grid, const Matrix global, u32 n_local, u32* dst) {
  const u64 count = (u64)n_local * n_local;

  for (i32 r = 0; r < grid->p; r++) {
    const u32* src   = block_origin(global, n_local, r / grid->q, r % grid->q);
    u32*       block = &dst[r * count];

    for (u32 i = 0; i < n_local; i++)
      memcpy(&block[i * n_local], &src[(u64)i * global.n], n_local * sizeof(u32));
  }
}

/* Inverse of pack_blocks: src holds the blocks in rank order. */
static void unpack_blocks(const FoxGrid* grid, const u32* src, u32 n_local, Matrix global) {
  const u64 count = (u64)n_local * n_local;

  for (i32 r = 0; r < grid->p; r++) {
    const u32* block = &src[r * count];
    u32*       dst   = block_origin(global, n_local, r / grid->q, r % grid->q);

    for (u32 i = 0; i < n_local; i++)
      memcpy(&dst[(u64)i * global.n], &block[i * n_local], n_local * sizeof(u32));
  }
}

void fox_scatter(
  const FoxGrid* grid,
  const Matrix global,
  Matrix local) {

  const u32 n_local = local.n;
  const u64 count   = (u64)n_local * n_local;
  u32* sendbuf      = NULL;

  if (grid->my_rank == ROOT) {
    sendbuf = alloc_or_abort(grid->p * count, "fox_scatter");
    pack_blocks(grid, global, n_local, sendbuf);
  }

  /* Rank r receives piece r of sendbuf. The send arguments only matter on ROOT. */
  MPI_Scatter(sendbuf,    (i32)count, MPI_UINT32_T,
              local.data, (i32)count, MPI_UINT32_T,
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
  MPI_Gather(local.data, (i32)count, MPI_UINT32_T,
             recvbuf,    (i32)count, MPI_UINT32_T,
             ROOT, grid->comm);

  if (grid->my_rank != ROOT)
    return MATRIX_NULL;

  unpack_blocks(grid, recvbuf, n_local, global);

  free(recvbuf);
  return global;
}

/* ---- Shared variant: every process holds every block (see fox.h). ---- */

Matrix fox_block(const FoxGrid* grid, const Matrix shared, i32 row, i32 col) {
  const u32 n_local = shared.n / (u32)grid->q;
  const u64 count   = (u64)n_local * n_local;

  return (Matrix){ .n = n_local, .data = &shared.data[(u64)(row * grid->q + col) * count] };
}

void fox_share(
  const FoxGrid* grid,
  const Matrix global,
  Matrix shared) {

  if (grid->my_rank == ROOT)
    pack_blocks(grid, global, shared.n / (u32)grid->q, shared.data);

  MPI_Bcast(shared.data, (i32)((u64)shared.n * shared.n), MPI_UINT32_T, ROOT, grid->comm);
}

Matrix fox_unshare(
  const FoxGrid* grid,
  Matrix global,
  const Matrix shared) {

  if (grid->my_rank != ROOT)
    return MATRIX_NULL;

  unpack_blocks(grid, shared.data, shared.n / (u32)grid->q, global);
  return global;
}

void fox_multiply_shared(
  const FoxGrid* grid,
  Matrix local_C,
  const Matrix shared) {

  minplus_matrix_reset(local_C);

  /* Same steps as fox_multiply, but A(my_row, k) and B(k, my_col) are read
   * straight from the shared copy instead of being broadcast / shifted in. */
  for (i32 step = 0; step < grid->q; step++) {
    const i32 k = (grid->my_row + step) % grid->q;

    minplus_relax(local_C,
                  fox_block(grid, shared, grid->my_row, k),
                  fox_block(grid, shared, k, grid->my_col));
  }
}

void fox_allgather(
  const FoxGrid* grid,
  Matrix shared,
  const Matrix local) {

  const i32 count = (i32)((u64)local.n * local.n);

  /* Rank r's block lands at piece r of shared, i.e. block (r/q, r%q). */
  MPI_Allgather(local.data,  count, MPI_UINT32_T,
                shared.data, count, MPI_UINT32_T,
                grid->comm);
}
