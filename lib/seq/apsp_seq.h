#ifndef APSP_SEQ_H
#define APSP_SEQ_H

#include "matrix_io.h"
#include "minplus.h"

/* Repeated min-plus squaring: the sequential baseline for Fox. */
Matrix apsp           (Matrix m);

/* O(n^3) reference, used by the tests to cross-check apsp(). */
Matrix floyd_warshall (Matrix m);

#endif /* !APSP_SEQ_H */
