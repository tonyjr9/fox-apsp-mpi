#include "imports.h"
#include <mpi.h>
#include <string.h>

#ifdef USE_MPI

#include "fox.h"

int main(int argc, char **argv) {

  int n = 0, my_rank, nproc;

  Matrix input = MATRIX_NULL;
  
  GridData grid;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  if (my_rank == 0) {
    input = matrix_read(stdin);
    n = (input.data == NULL) ? 0 : input.n;
  }
 
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (n <= 0 || !fox_validate((u32)n, (u16)nproc, my_rank)) {
    if (my_rank == 0) matrix_free(input);
    MPI_Finalize();
    return 1;
  }

  fox_setup_grid(&grid);

  const u32 n_local = (u32)(n/grid.q);
  const u32 count   = n_local * n_local;
  Matrix A = matrix_create(n_local, 0);

  // Shares work between processes
  if (my_rank == 0) {
    int dest = 0;
    Matrix tmp = matrix_create(n_local, 0);
    for (int row = 0; row != grid.q; ++row) {
      for (int col = 0; col != grid.q; ++col) {

        // Copy block (row, col) of the input matrix into tmp
        for (int i = 0; i != (int)n_local; ++i) {
          for (int j = 0; j != (int)n_local; ++j) {
            tmp.data[i * n_local + j] =
                input.data[(row * n_local + i) * n + (col * n_local + j)];
          }
        }

        if (dest == 0) {
          // MPI_Send to self would deadlock
          memcpy(A.data, tmp.data, (size_t)count * sizeof(u32));
        } else { 
          MPI_Send(tmp.data, (int)count, MPI_UINT32_T,
                   dest, 0, MPI_COMM_WORLD);
        }
        ++dest;
      }
    }
    matrix_free(tmp);
  } else {
    MPI_Recv(A.data, (int)count, MPI_UINT32_T,
             0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  Matrix B = matrix_copy(A);
  Matrix C = matrix_create(n_local, INF);

  // APSP by repeated squaring
  // If n <= 2 no product runs and A is already the answer.
  u32 step = 1;
  while (step < (u32)n - 1) {
    fox(n, &grid, C, A, B);
    memcpy(A.data, C.data, (size_t)count * sizeof(u32));
    memcpy(B.data, C.data, (size_t)count * sizeof(u32));

    step *= 2;
  }

  // Gather Results
  if (my_rank == 0) {
    Matrix tmp = matrix_create(n_local, 0);
    int src = 0;

    for (int row = 0; row != grid.q; ++row) {
      for (int col = 0; col != grid.q; ++col) {

        if (src == 0) {
          memcpy(tmp.data, A.data, (size_t)count * sizeof(u32));
        } else {
          MPI_Recv(tmp.data, (int)count, MPI_UINT32_T,
                    src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        ++src;

        for (int i = 0; i != (int)n_local; ++i) {
          for (int j = 0; j != (int)n_local; ++j) {
            int dest = i * n_local + j;
            input.data[(row * n_local + i) * n + (col * n_local + j)] =
                  (tmp.data[dest] == INF) ? 0 : tmp.data[dest];
          }
        }

      }
    }
    
    matrix_print(stdout, input);

    matrix_free(input);
    matrix_free(tmp);
  } else {
    MPI_Send(A.data, (int)count, MPI_UINT32_T,
              0, 0, MPI_COMM_WORLD);
  }

  matrix_free(A);
  matrix_free(B);
  matrix_free(C);

  MPI_Finalize();
  return 0;
}

#else // Sequencial 

#include "apsp_seq.h"
  
int main(void){
  Matrix m = matrix_read(stdin);
  if (m.data == NULL)
    return 1;

  Matrix r = apsp(m);

  matrix_print(stdout, r);
  
  matrix_free(r);
  matrix_free(m);
  return 0;
}

#endif
