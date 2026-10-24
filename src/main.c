#include "imports.h"

#ifdef USE_MPI
#include <mpi.h>

#include "apsp_mpi.h"

i32 main(const i32 argc, const char* argv[]) {
  if (argc < 3) {
    fprintf(stderr, "Expected file name path to matrix binary and output file name\n");
    return 1;
  }

  i32 n = 0, my_rank, nproc;

  Matrix input = MATRIX_NULL;

  FoxGrid grid;

  MPI_Init((i32*)&argc, (char***)&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  // Only the root reads the input file; the others learn n from the broadcast.
  if (my_rank == 0) {
    FILE* in = fopen(argv[1], "rb");
    if (in != NULL) {
      input = matrix_read_bin(in);
      fclose(in);
    } else {
      fprintf(stderr, "Could not open input file %s\n", argv[1]);
    }
    n = (input.data == NULL) ? 0 : (i32)input.n;
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
  const f64 start = MPI_Wtime();

  Matrix result = mpi_apsp(&grid, input, (u32)n);

  const f64 elapsed_ms = (MPI_Wtime() - start) * 1e3;

  if (my_rank == 0) {
    fprintf(stderr, "compute_ms=%.3f\n", elapsed_ms);

    FILE* out = fopen(argv[2], "wb");
    if (out == NULL) {
      fprintf(stderr, "Could not open output file %s\n", argv[2]);
      matrix_free(result);
      return 1;
    }

    matrix_write_bin(out, result);
    fclose(out);
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

static f64 now_ms(void) {
  struct timespec ts;
  timespec_get(&ts, TIME_UTC);
  return ts.tv_sec * 1e3 + ts.tv_nsec / 1e6;
}

i32 main(const i32 argc, const char* argv[]) {
  if (argc < 3) {
    fprintf(stderr, "Expected file name path to matrix binary and output file name\n");
    return 1;
  }

  FILE* in = fopen(argv[1], "rb");
  if (in == NULL) {
    fprintf(stderr, "Could not open input file %s\n", argv[1]);
    return 1;
  }

  Matrix m = matrix_read_bin(in);
  fclose(in);
  if (m.data == NULL)
    return 1;

  const f64 start = now_ms();
  Matrix r = apsp(m);
  fprintf(stderr, "compute_ms=%.3f\n", now_ms() - start);

  FILE* out = fopen(argv[2], "wb");
  if (out == NULL) {
    fprintf(stderr, "Could not open output file %s\n", argv[2]);
    matrix_free(r);
    matrix_free(m);
    return 1;
  }

  matrix_write_bin(out, r);
  fclose(out);

  matrix_free(r);
  matrix_free(m);

  return 0;
}

#endif
