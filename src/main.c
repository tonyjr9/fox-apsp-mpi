#include "imports.h"

#ifdef USE_MPI
#include <mpi.h>

#include "apsp_mpi.h"

int main(int argc, char **argv) {

  int n = 0, my_rank, nproc;

  Matrix input = MATRIX_NULL;

  FoxGrid grid;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  // Only the root reads stdin; the others learn n from the broadcast.
  if (my_rank == 0) {
    input = matrix_read(stdin);
    n = (input.data == NULL) ? 0 : (int)input.n;
  }

  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (n <= 0 || !fox_validate((u32)n, nproc, my_rank)) {
    matrix_free(input);
    MPI_Finalize();
    return 1;
  }

  fox_setup_grid(&grid);

  // Start the clock together; I/O is left out of the measurement.
  MPI_Barrier(MPI_COMM_WORLD);
  const double start = MPI_Wtime();

  Matrix result = mpi_apsp(&grid, input, (u32)n);

  const double elapsed_ms = (MPI_Wtime() - start) * 1e3;

  if (my_rank == 0) {
    fprintf(stderr, "compute_ms=%.3f\n", elapsed_ms);
    matrix_print(stdout, result);
  }

  // On the root, result is input (written in place); elsewhere both are MATRIX_NULL.
  matrix_free(result);

  fox_grid_free(&grid);
  MPI_Finalize();
  return 0;
}

#else // Sequential

#include <time.h>

#include "apsp_seq.h"

static double now_ms(void) {
  struct timespec ts;
  timespec_get(&ts, TIME_UTC);
  return ts.tv_sec * 1e3 + ts.tv_nsec / 1e6;
}

int main(void) {
  Matrix m = matrix_read(stdin);
  if (m.data == NULL)
    return 1;

  const double start = now_ms();
  Matrix r = apsp(m);
  fprintf(stderr, "compute_ms=%.3f\n", now_ms() - start);

  matrix_print(stdout, r);

  matrix_free(r);
  matrix_free(m);
  return 0;
}

#endif
