# fox-apsp-mpi

All-pairs shortest path (APSP) solver for directed weighted graphs, parallelised
with MPI. It computes the distance matrix by min-plus repeated squaring, and each
matrix product is distributed with Fox's algorithm.

Project Assignment I, Parallel Computing (2026/27), University of Porto.

## Authors

- António Lanção (up202205416) 
- 
- 

## Algorithm

1. Read the N x N adjacency matrix. Treat 0 off the diagonal as "no edge" (INF).
2. Set D1 = W. While m < N-1, compute D(2m) = D(m) (x) D(m) and set m = 2m.
   (x) is the min-plus product: C[i][j] = min over k of (A[i][k] + B[k][j]).
3. Each product runs as Fox's algorithm on a Q x Q process grid (P = Q^2).
   Each process owns one (N/Q) x (N/Q) block.
4. Print the result. Unreachable pairs are printed as 0.

Constraints: P must be a perfect square (P = Q^2) and N mod Q = 0.
Otherwise the program prints an error and aborts.

## Repository layout

```
.
├── Makefile
├── README.md
├── src/
│   ├── fox.c          # MPI implementation
│   └── seq.c          # sequential baseline
├── tests/
│   ├── examples/      # input/expected output pairs
│   └── run_tests.sh   # compares output against expected files
├── bench/
│   ├── run_bench.sh   # runs P = 1, 4, 9, 16, 25
│   └── results.csv
└── report/
    ├── report.typ
    └── report.pdf
```

## Build

```
make          # builds ./fox (MPI) and ./seq (sequential)
make clean
```

Requires an MPI implementation (Open MPI or MPICH) and `mpicc`.

## Run

Input is read from stdin and output goes to stdout:

```
mpirun -np 4 ./fox < tests/examples/example.in
```

### Example

Input:

```
6
0 2 0 5 0 0
0 0 0 0 0 0
0 2 0 0 0 5
0 0 0 0 1 0
3 9 3 0 0 0
0 0 0 0 1 0
```

Expected output:

```
0 2 9 5 6 14
0 0 0 0 0 0
9 2 0 14 6 5
4 6 4 0 1 9
3 5 3 8 0 8
4 6 4 9 1 0
```

## Testing

```
./tests/run_tests.sh      # runs every example with P = 1, 4, 9, 16, 25 where valid
```

## Benchmarks

```
./bench/run_bench.sh      # writes bench/results.csv
```

Times are in milliseconds and exclude I/O (measured with MPI_Wtime around
the compute phase only, maximum over ranks).

## Status

- [ ] Sequential baseline
- [ ] Input parsing and validation
- [ ] Grid setup, scatter and gather
- [ ] Fox's algorithm
- [ ] Repeated squaring loop
- [ ] Tests passing
- [ ] Benchmarks
- [ ] Report

