#ifndef MATRIX_H
#define MATRIX_H

#include <stdint.h>
#include <stdbool.h>

#define INF ((uint32_t)(UINT32_MAX / 2))
#define MAT(m, i, j) ((m)->data[(i) * (m)->n + (j)]) 
typedef uint32_t u32;

typedef struct Matrix {
  uint32_t n;
  uint32_t *data;
} Matrix;

Matrix *matrix_create(uint32_t n, uint32_t fill);
void    matrix_free(Matrix *m);
Matrix *matrix_copy(const Matrix *m);
bool    matrix_equal(const Matrix *a, const Matrix *b);

#endif
