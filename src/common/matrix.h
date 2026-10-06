#ifndef MATRIX_H
#define MATRIX_H

#include "utils.h"

#define INF ((u32)(UINT32_MAX / 2))
#define MAT(m, i, j) ((m).data[(i) * (m).n + (j)]) 

typedef struct matrix Matrix;

#define MATRIX_NULL { .n = 0, .data = NULL }

struct matrix {
  u32  n;
  u32* data;
};

Matrix matrix_create (u32 n, u32 fill);
Matrix matrix_copy   (const Matrix m);
bool   matrix_equal  (const Matrix a, const Matrix b);
void   matrix_free   (Matrix m);

#endif
