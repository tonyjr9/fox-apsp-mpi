#ifndef MINPLUS_H
#define MINPLUS_H

#include "matrix.h"
#include "utils.h"

Matrix minplus_multiply(Matrix A, Matrix B);
void   minplus_accumulate(Matrix C, Matrix A, Matrix B);
Matrix minplus_apsp(Matrix M);

#endif
