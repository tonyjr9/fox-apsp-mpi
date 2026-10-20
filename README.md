# fox-apsp-mpi

All-pairs shortest path (APSP) solver for directed weighted graphs. Distances are
computed with min-plus repeated squaring; the MPI version distributes each
product with Fox's algorithm.

Project Assignment I, Parallel Computing (2026/27), University of Porto.

## Authors

- Name 1 (up2026XXXXX)
- Name 2 (up2026XXXXX)

## Algorithm

1. Read the N x N matrix. Off the diagonal, 0 means "no edge" and is stored as INF;
   the diagonal is stored as 0.
2. Set D = W and m = 1. While m < N-1: D = D (x) D, m = 2m.
   (x) is the min-plus product: C[i][j] = min over k of (A[i][k] + B[k][j]).
   Both versions stop early as soon as a squaring leaves D unchanged.
3. The sequential version (`seq`) runs the squarings with a tiled, vectorised
   min-plus kernel.
4. The MPI version (`fox`) scatters the matrix as Q x Q blocks, applies Fox's
   algorithm to each product and gathers the result on rank 0. The A broadcast
   and B shift for step s+1 (`MPI_Ibcast`, `MPI_Isend`/`MPI_Irecv`) overlap the
   computation of step s.
5. INF (unreachable) is printed as 0.

Constraint: P = Q*Q and N mod Q = 0. Otherwise rank 0 prints
`ERROR: Invalid configuration!` on stdout (the reason goes to stderr) and the
program exits with status 1.

Both binaries print the computation time (excluding I/O) to stderr as
`compute_ms=...`.

## Repository layout

```
.
├── Makefile
├── README.md
├── src/
│   └── main.c            # entry point; built as ./seq, or as ./fox with -DUSE_MPI
├── lib/
│   ├── common/           # shared, no MPI
│   │   ├── matrix.[ch]       # Matrix struct, create/free/copy/equal
│   │   ├── matrix_io.[ch]    # read from stdin, print, INF <-> 0 mapping
│   │   ├── minplus.[ch]      # min-plus kernel (accumulate / multiply)
│   │   └── utils.h
│   ├── seq/              # sequential baseline + Floyd-Warshall reference
│   │   └── apsp_seq.[ch]
│   └── mpi/              # distributed version
│       ├── foxgrid.[ch]      # validation, Q x Q grid, row/column communicators
│       ├── fox.[ch]          # Fox multiply, block scatter/gather
│       └── apsp_mpi.[ch]     # repeated squaring over Fox
├── tests/
│   ├── test_util.h       # CHECK() macro and helpers
│   ├── test_matrix.c
│   ├── test_matrix_io.c
│   ├── test_minplus.c
│   ├── matgen.c          # random N x N input generator
│   ├── run_tests.sh      # regression on tests/examples (seq or fox)
│   └── examples/
│       ├── input/            # input5, input6, input300 ... input1200
│       └── output/           # expected outputs
├── bench/
│   ├── run_bench.sh      # compute-time benchmark -> results.csv
│   ├── run_test_time.sh  # end-to-end seq timing per example
│   └── results.csv
└── report/
    ├── template.typ      # minimal Typst report template
    └── report.typ
```

## Build

```
make            # builds ./seq and ./fox
make seq        # builds ./seq only
make fox        # builds ./fox only
make clean
```

The MPI build requires an implementation providing `mpicc` (Open MPI or MPICH).
`-march=native` tunes for the build machine: on the cluster, compile on the
compute nodes.

## Run

Input is read from stdin and the result is written to stdout:

```
./seq < tests/examples/input/input6
mpirun -np 4 ./fox < tests/examples/input/input6
```

## Tests

```
make test                     # unit tests: matrix, matrix_io, minplus
make check-seq                # compare ./seq with tests/examples
make check-fox                # compare ./fox with tests/examples, np = 1 4 9
make check-fox NP="4 16"      # other process counts
make check                    # all of the above
make time-seq INPUT=tests/examples/input/input600
```

`tests/run_tests.sh` compares outputs byte by byte. For `fox`, inputs whose N
is not divisible by sqrt(NP) must be rejected with exactly the content of
`output5`. `input5` is skipped when its configuration is valid (seq, np=1),
since its reference output is that error message. `MPIRUN` defaults to
`mpirun --oversubscribe` (Open MPI); override it for MPICH.

## Benchmarks

```
make bench                                        # seq + fox, np = 1 4 9 16, 3 reps
NPS="1 4 9 16 25" REPS=5 make bench               # other process counts / reps
MPIRUN="mpirun --hostfile hosts" make bench       # several nodes
INPUTS="input300 input600" make bench             # subset of tests/examples
```

`bench/run_bench.sh` writes one CSV row per run to `bench/results.csv`
(`impl,np,input,n,rep,compute_ms`) and prints the best time and the speedup
over `seq` for each configuration on stderr. Configurations where sqrt(np)
does not divide N are skipped; failed runs are recorded as `NA`. When calling
the script directly on a machine with fewer cores than `np`, set
`MPIRUN="mpirun --oversubscribe"` (`make bench` does this by default).

The final benchmark measures computation-only times for sequential,
one-process MPI, and 4, 9, 16, and 25 MPI processes on the LabCC cluster.
`bench/run_test_time.sh` times the whole `./seq` run (I/O included) for each
example.

## Report

`report/report.typ` (Typst), styled by `report/template.typ`:

```
typst compile report/report.typ        # -> report/report.pdf
typst watch report/report.typ          # recompile on save
```

Unfinished parts are marked with `#todo[...]` and show up in red in the PDF.

## Status

- [x] Matrix struct, I/O and unit tests
- [x] Min-plus kernel
- [x] Sequential baseline
- [x] Sequential regression against professor examples
- [x] Input validation and process grid
- [x] Scatter / gather
- [x] Fox's algorithm
- [x] MPI example and invalid-configuration tests
- [ ] Benchmarks
- [ ] Report (structure and implementation sections written; results pending)
