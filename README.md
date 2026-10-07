# fox-apsp-mpi

All-pairs shortest path (APSP) solver for directed weighted graphs. The current
sequential baseline computes distances with min-plus repeated squaring; the MPI
version will distribute the products with Fox's algorithm.

Project Assignment I, Parallel Computing (2026/27), University of Porto.

## Authors

- Name 1 (up2026XXXXX)
- Name 2 (up2026XXXXX)

## Algorithm

1. Read the N x N matrix. Off the diagonal, 0 means "no edge" and is stored as INF;
   the diagonal is stored as 0.
2. Set D = W and m = 1. While m < N-1: D = D (x) D, m = 2m.
   (x) is the min-plus product: C[i][j] = min over k of (A[i][k] + B[k][j]).
3. The sequential version repeats the min-plus product until the covered path
   length is at least N-1.
4. The MPI version will apply Fox's algorithm to each product. INF
   (unreachable) is printed as 0.

Constraint: P = Q*Q and N mod Q = 0. Otherwise the program prints an error and
aborts without writing to stdout.

## Repository layout

```
.
├── Makefile
├── README.md
├── src/
│   └── main.c            # current sequential stdin/stdout entry point
├── lib/
│   ├── common/           # shared, no MPI
│   │   ├── matrix.[ch]       # Matrix struct, create/free/copy/equal
│   │   ├── matrix_io.[ch]    # read from stdin, print, INF <-> 0 mapping
│   │   ├── minplus.[ch]      # min-plus kernel (accumulate / multiply)
│   │   └── utils.h
│   ├── seq/              # sequential repeated squaring (baseline)
│   │   └── apsp_seq.[ch]
│   └── mpi/              # distributed version
│       └── fox.c             # MPI grid and Fox implementation
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
make            # currently builds ./seq
make seq        # builds ./seq
make clean
```

The current sequential build requires `gcc`. The MPI build will require an
implementation providing `mpicc` (Open MPI or MPICH).

## Run

Input is read from stdin and the result is written to stdout:

```
./seq < tests/examples/input/input6
```

MPI compatibility checks and invalid process-count handling will be added with
the Fox implementation.

## Tests

```
make test                     # unit tests: matrix, matrix_io, minplus
make check-seq                # compare ./seq with tests/examples
make check                    # unit tests plus sequential regression
make time-seq INPUT=tests/examples/input/input600
```

`tests/run_tests.sh` currently runs `seq` against the provided expected outputs
for `input6`, `input300`, `input600`, `input900`, and `input1200`. `input5` is
reserved for MPI compatibility testing because its dimension is incompatible
with the required square process grids. Outputs are compared byte by byte.

## Benchmarks

```
bench/run_test_time.sh        # current sequential timing helper
```

The final benchmark will measure computation-only times for sequential,
one-process MPI, and 4, 9, 16, and 25 MPI processes on the LabCC cluster.

## Report

`report/report.typ` (Typst). Compile with `typst compile report/report.typ`.

## Status

- [x] Matrix struct, I/O and unit tests
- [x] Min-plus kernel
- [x] Sequential baseline
- [x] Sequential regression against professor examples
- [ ] Input validation and process grid
- [ ] Scatter / gather
- [ ] Fox's algorithm
- [ ] MPI example and invalid-configuration tests
- [ ] Benchmarks
- [ ] Report
