#ifndef MINPLUS_H
#define MINPLUS_H

#include "matrix.h"
#include "utils.h"

void   minplus_matrix_reset  (Matrix c);
void   minplus_relax         (Matrix c, const Matrix a, const Matrix b);
Matrix minplus_multiply      (const Matrix a, const Matrix b);
void   minplus_multiply_into (Matrix c, const Matrix a, const Matrix b);
//void   minplus_accumulate    (Matrix c, const Matrix a, const Matrix b);
Matrix minplus_apsp          (Matrix m);

#endif
