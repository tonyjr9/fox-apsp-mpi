#ifndef __APSP_SEQ_H__
#define __APSP_SEQ_H__

#include "matrix_io.h"
#include "minplus.h"

Matrix apsp           (Matrix m);
Matrix floyd_warshall (Matrix m);

#endif /* !__APSP_SEQ_H__ */
