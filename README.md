# fox-apsp-mpi

All-pairs shortest path (APSP) solver for directed weighted graphs, parallelised
with MPI. Distances are computed by min-plus repeated squaring, and every matrix
product is distributed with Fox's algorithm on a Q x Q process grid.

Project Assignment I, Parallel Computing (2026/27), University of Porto.

## Authors

- Name 1 (up2026XXXXX)
- Name 2 (up2026XXXXX)

## Algorithm

1. Read the N x N matrix. Off the diagonal, 0 means "no edge" and is stored as INF;
   the diagonal is stored as 0.
2. Set D = W and m = 1. While m < N-1: D = D (x) D, m = 2m.
   (x) is the min-plus product: C[i][j] = min over k of (A[i][k] + B[k][j]).
3. Fox's algorithm computes each product. With P = Q*Q processes, each owns one
   (N/Q) x (N/Q) block. In stage s, the process in column (i+s) mod Q of grid
   row i broadcasts its A block along the row, every process multiplies it with
   its current B block (min-plus, accumulating into C), and B blocks shift up
   one row of the grid.
4. Print the result. INF (unreachable) is printed as 0.

Constraint: P = Q*Q and N mod Q = 0. Otherwise the program prints an error and
aborts without writing to stdout.

## Repository layout

```
.
├── Makefile
├── README.md
├── src/
│   └── main.c            # entry point (built twice: MPI and sequential)
├── lib/
│   ├── common/           # shared, no MPI
│   │   ├── matrix.[ch]       # Matrix struct, create/free/copy/equal
│   │   ├── matrix_io.[ch]    # read from stdin, print, INF <-> 0 mapping
│   │   ├── minplus.[ch]      # min-plus kernel (accumulate / multiply)
│   │   └── utils.h
│   ├── seq/              # sequential repeated squaring (baseline)
│   │   └── apsp_seq.[ch]
│   └── mpi/              # distributed version
│       └── fox.c             # grid setup, scatter/gather, Fox stages
├── tests/
│   ├── test_util.h       # CHECK() macro and helpers
│   ├── test_matrix.c
│   ├── test_matrix_io.c
│   ├── test_minplus.c
│   ├── run_tests.sh      # regression on tests/examples
│   └── examples/
│       ├── input/            # input5, input6, input300 ... input1200
│       └── output/           # expected outputs
├── bench/
│   ├── run_bench.sh
│   └── results.csv
└── report/
    └── report.typ
```

## Build

```
make            # builds ./fox (MPI) and ./seq (sequential baseline)
make clean
```

Requires `gcc` and an MPI implementation providing `mpicc` (Open MPI or MPICH).
`src/main.c` is compiled with `-DUSE_MPI` for `fox` and without it for `seq`.

## Run

Input is read from stdin and the result is written to stdout:

```
mpirun -np 4 ./fox < tests/examples/input/input6
./seq < tests/examples/input/input6
```

Invalid combinations (for example `input6` with 16 processes) print an error
on stderr and exit with a non-zero status.

## Tests

```
make test                     # unit tests: matrix, matrix_io, minplus
make check                    # build, then compare against tests/examples
MAXN=600 PROCS="1 4" tests/run_tests.sh      # quick subset
```

`tests/run_tests.sh` runs `seq` and `fox` on every input, checks the output byte
by byte against `tests/examples/output/`, and checks that invalid process counts
are rejected. Set `MPIRUN` to change the launcher (default
`mpirun --oversubscribe`).

## Benchmarks

```
bench/run_bench.sh            # P = 1, 4, 9, 16, 25 -> bench/results.csv
```

Times are in milliseconds and exclude I/O: `MPI_Wtime` around the compute phase
only, maximum over ranks, reported on stderr. Run on the LabCC machines, not on
an oversubscribed laptop, when producing numbers for the report.

## Report

`report/report.typ` (Typst). Compile with `typst compile report/report.typ`.

## Status

- [x] Matrix struct, I/O and unit tests
- [ ] Min-plus kernel
- [ ] Sequential baseline
- [ ] Input validation and process grid
- [ ] Scatter / gather
- [ ] Fox's algorithm
- [ ] All example tests passing
- [ ] Benchmarks
- [ ] Report
