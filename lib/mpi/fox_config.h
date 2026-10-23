#ifndef FOX_CONFIG_H
#define FOX_CONFIG_H

/*
 * How the grid distributes the matrix, and therefore how fox_multiply
 * communicates. Selected at compile time with -DFOX_DIST=<option>
 * (see the fox / fox_shared / fox-colmajor targets in the Makefile).
 *
 * - FOX_DIST_ROW      rank r holds block (r/q, r%q). A is broadcast along
 *                     the row, B circulates along the column. (default)
 * - FOX_DIST_COLMAJOR rank r holds block (r%q, r/q). A is broadcast along
 *                     the column, B circulates along the row.
 * - FOX_DIST_IDENT    no scatter, no communication: every process holds a
 *                     full copy of the matrix and computes a full min-plus
 *                     product on it at each step (the communication cost
 *                     of the distributed algorithm is zero, but so is the
 *                     speedup: it measures the compute cost alone).
 *                     Memory per process ~ 5 * n^2 * 4 bytes — small/medium
 *                     N only.
 */
#define FOX_DIST_ROW      1
#define FOX_DIST_COLMAJOR 2
#define FOX_DIST_IDENT    3

#ifndef FOX_DIST
#define FOX_DIST FOX_DIST_ROW
#endif

#endif // !FOX_CONFIG_H
