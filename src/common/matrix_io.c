#include "matrix_io.h"

#include <inttypes.h>

Matrix *matrix_read(FILE *in)
{
    uint32_t n;
    if (fscanf(in, "%" SCNu32, &n) != 1 || n == 0) {
        fprintf(stderr, "matrix_read: missing or invalid dimension N\n");
        return NULL;
    }

    Matrix *m = matrix_create(n, INF);

    for (uint32_t i = 0; i < n; i++) {
        for (uint32_t j = 0; j < n; j++) {
            uint32_t v;
            if (fscanf(in, "%" SCNu32, &v) != 1) {
                fprintf(stderr, "matrix_read: expected %" PRIu32 "x%" PRIu32
                        " values, input ended at (%" PRIu32 ",%" PRIu32 ")\n",
                        n, n, i, j);
                matrix_free(m);
                return NULL;
            }
            if (i == j)
                MAT(m, i, j) = 0;
            else
                MAT(m, i, j) = (v == 0) ? INF : v;
        }
    }
    return m;
}

void matrix_print(FILE *out, const Matrix *m)
{
    for (uint32_t i = 0; i < m->n; i++) {
        for (uint32_t j = 0; j < m->n; j++) {
            uint32_t v = MAT(m, i, j);
            if (v == INF)
                v = 0;
            fprintf(out, j + 1 < m->n ? "%" PRIu32 " " : "%" PRIu32, v);
        }
        fputc('\n', out);
    }
}
