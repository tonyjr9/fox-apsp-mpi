#include "../common/matrix.h"
#include "../common/matrix_io.h"

Matrix apsp(Matrix m) {
  Matrix r = matrix_create(m.n, 0);

  u32 n = m.n;
  for (u32 k = 0; k < n; k++) {
    for (u32 i = 0; i < n; i++) {
      for (u32 j = 0; j < n; j++) {
        MAT(r,i,j) = MIN(MAT(m, i, j), (MAT(m, i, k) + MAT(m, k, j));
      }
    }
  }

  return r;
}

int main(int argc, char **argv) { 

  Matrix m = matrix_read(stdin);

  return 0;
}
