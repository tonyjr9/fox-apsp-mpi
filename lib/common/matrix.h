#ifndef MATRIX_H
#define MATRIX_H

#include <stddef.h>

#include "utils.h"

typedef struct matrix Matrix;

#define MATRIX_NULL (Matrix){ .n = 0, .data = NULL }

struct matrix {
  u32  n;
  u32* data;
};

Matrix matrix_create (u32 n, u32 fill);
Matrix matrix_copy   (const Matrix m);
bool   matrix_equal  (const Matrix a, const Matrix b);
void   matrix_free   (Matrix m);

#endif
