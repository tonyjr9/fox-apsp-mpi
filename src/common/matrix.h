#ifndef MATRIX_H
#define MATRIX_H

#include <stdint.h>
#include <stdbool.h>

#define INF ((u32)(UINT32_MAX / 2))
#define MAT(m, i, j) ((m)->data[(i) * (m)->n + (j)]) 

typedef struct Matrix {
  u32 n;
  u32 *data;
} Matrix;

Matrix matrix_create (u32 n, u32 fill);
Matrix matrix_copy   (const Matrix m);
bool   matrix_equal  (const Matrix a, const Matrix b);
void   matrix_free   (Matrix m);

#endif
