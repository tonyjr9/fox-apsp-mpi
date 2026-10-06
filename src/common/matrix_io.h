#ifndef MATRIX_IO_H
#define MATRIX_IO_H

#include <stdio.h>
#include "matrix.h"

/* Reads "N" followed by N*N integers from `in`.
 * Off-diagonal 0 becomes INF (no edge); the diagonal is forced to 0.
 * Returns a new matrix, or NULL (with a message on stderr) on bad input. */
Matrix *matrix_read(FILE *in);

/* Prints one row per line, values separated by one space.
 * INF (unreachable) is printed as 0, as the assignment requires. */
void matrix_print(FILE *out, const Matrix *m);

#endif
