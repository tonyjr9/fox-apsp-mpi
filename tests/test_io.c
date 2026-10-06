#include <stdio.h>
#include <string.h>

#include "matrix.h"
#include "matrix_io.h"

static void print_raw(const Matrix *m)
{
    for (uint32_t i = 0; i < m->n; i++) {
        for (uint32_t j = 0; j < m->n; j++) {
            if (MAT(m, i, j) == INF)
                printf("INF ");
            else
                printf("%u ", (unsigned)MAT(m, i, j));
        }
        putchar('\n');
    }
}

int main(int argc, char **argv)
{
    Matrix *m = matrix_read(stdin);
    if (m == NULL)
        return 1;

    /* matrix_copy / matrix_equal sanity checks */
    Matrix *c = matrix_copy(m);
    if (!matrix_equal(m, c)) {
        fprintf(stderr, "FAIL: copy differs from original\n");
        return 2;
    }
    MAT(c, 0, m->n - 1) += 1;
    if (matrix_equal(m, c)) {
        fprintf(stderr, "FAIL: equal reports true after modification\n");
        return 3;
    }

    if (argc > 1 && strcmp(argv[1], "--raw") == 0)
        print_raw(m);
    else
        matrix_print(stdout, m);

    matrix_free(c);
    matrix_free(m);
    return 0;
}
