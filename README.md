# fox-apsp-mpi

All-pairs shortest path (APSP) solver for directed weighted graphs. Distances are
computed with min-plus repeated squaring; the MPI version distributes each
product with Fox's algorithm.

Project Assignment I, Parallel Computing (2026/27), University of Porto.

## Authors

- António Lanção    (up202205416)
- David Sá          (up202303580)
- Henrique Teixeira (up202306640)

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
│   └── main.c            # entry point: ./seq; with -DUSE_MPI ./fox_bin, plus -DTEXT_IO ./fox
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
│   ├── matgen.c          # random N x N input generator (binary output)
│   ├── binconv.c         # text matrix -> binary matrix (stdin -> stdout)
│   ├── run_tests.sh      # regression on tests/examples (seq or fox)
│   └── examples/
│       ├── input/            # input5, input6, input300 ... input1200 (+ .bin)
│       ├── output/           # expected outputs (+ .bin)
│       └── to_bin.sh         # convert the text examples (via make examples-bin)
├── bench/
│   ├── run_bench.sh      # compute-time benchmark -> results.csv
│   ├── run_test_time.sh  # end-to-end seq timing per example
│   ├── run_rand.sh       # fox on the random matrices, whole lab (see below)
│   └── results.csv
└── report/
    ├── template.typ      # minimal Typst report template
    └── report.typ
```

## Build

```
make            # builds ./seq, ./fox_bin, ./fox_shared and ./fox
make fox        # builds ./fox only: text on stdin -> stdout (the assignment)
make fox_bin    # builds ./fox_bin only: binary files (tests, benchmarks)
make seq        # builds ./seq only
make gen        # random benchmark inputs: N = 2^6 .. 2^14 (seed 42), zipped
make examples-bin   # convert text examples to the binary format
make report     # graphs + report/report.pdf from the CSVs in bench/report/
make dist       # make report, then fox-apsp-mpi.zip for delivery
make clean
```

The MPI build requires an implementation providing `mpicc` (Open MPI or MPICH).
`-march=native` tunes for the build machine: on the cluster, compile on the
compute nodes.

## Run

### Assignment interface (text): `./fox`

As the assignment specifies, `./fox` reads the text matrix from stdin and
prints the result to stdout:

```
make
mpirun -np 4 ./fox < tests/examples/input/input6 > output6
```

### Binary files: `./fox_bin`

We also implemented reading and writing a binary format (a `u32` N header
followed by the N x N raw `u32` values; the reference conversion is
`tests/binconv.c`). It is faster than parsing and printing text and uses less
memory, so the tests and benchmarks use it. `./fox_bin` is the same program as
`./fox` (`src/main.c` without `-DTEXT_IO`), with input and output files:

```
./seq input6.bin output6.bin
mpirun -np 4 ./fox_bin input6.bin output6.bin
```

The text examples can be converted on the fly: `make examples-bin`
writes `input/*.bin` / `output/*.bin` next to the originals.

## Tests

```
make test                     # unit tests: matrix, matrix_io, minplus
make check-seq                # compare ./seq with tests/examples
make check-fox                # compare ./fox_bin with tests/examples, np = 1 4 9
make check-fox NP="4 16"      # other process counts
make check                    # all of the above
make time-seq INPUT=tests/examples/input/input600
```

`tests/run_tests.sh` compares outputs byte by byte. For `fox`, inputs whose N
is not divisible by sqrt(NP) must be rejected with exactly the content of
`output5`. `input5` is skipped when its configuration is valid (seq, np=1),
since its reference output is that error message. `MPIRUN` defaults to
`mpirun --oversubscribe` (Open MPI); override it for MPICH.

## Report data

`report/report.pdf` is built from the CSVs in `bench/report/` (one row per run,
`bench/run_bench.sh` format): the graphs by `bench/plot_lab.py` and
`bench/plot_power2.py`, the result tables by Typst itself (`report/bench.typ`
takes the minimum per configuration). To add results, append their rows to
the CSV and run `make report`; missing configurations show as `---`.

## Benchmarks

```
make bench                                        # seq + fox_bin, np = 1 4 9 16, 3 reps
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

## Random matrices on the lab cluster

`bench/run_rand.sh` runs `fox` on every random matrix in a zip on the DCC
Lab8 machines, picking the biggest square process grid each size allows
(`fox` requires P = Q*Q and N mod Q = 0, so the grid is the largest Q
with Q*Q <= total live slots and Q dividing N). Replace `upXXXXYYYYY`
with your student number and `l815` with any lab machine, as in
`setup.md`.

### 1. Generate the input matrices

`make gen` builds `tests/matgen` and writes `rand0_input<N>.bin` for
every N = 2^6 .. 2^14 (fixed seed 42, ~1.4 GB total), then zips them to
`rand0_input.zip` in the current directory. Three ways to run it:

- **On your own PC** (needs gcc; no MPI for this target):

  ```sh
  make gen
  ```

- **On a lab machine** (already has MPI): log in per `setup.md`,
  `cd ~/Paralela`, `make gen`.

- **Size it down** if disk or time are short — the generator takes
  `N [count] [seed]` directly, or override the Makefile variables:

  ```sh
  make gen MATRICES_2_TO_THE="8 12 14" GEN_SEED=7
  # or one matrix, no zip
  make build/tests/matgen && build/tests/matgen 1024 1 42
  ```

  A matrix of dimension N occupies N*N*4 + 4 bytes: N = 16384 is 1 GB.

**Do not commit the zip or the `.bin` files to the repo: the zip is
~145 MiB, above GitHub's 100 MiB file limit.**

### 2. Transfer the zip to the lab machine

The inputs and `bench/run_rand.sh` must sit in the same directory on
the machine you run from:

```sh
# from your PC, with the lab machines already set up per setup.md
scp -o ProxyJump=upXXXXYYYYY@ssh.alunos.dcc.fc.up.pt rand0_input.zip \
    upXXXXYYYYY@l815.alunos.dcc.fc.up.pt:~/Paralela/
```

(rsync works too: `rsync -avz -e "ssh -J upXXXXYYYYY@ssh.alunos.dcc.fc.up.pt"
rand0_input.zip upXXXXYYYYY@l815.alunos.dcc.fc.up.pt:~/Paralela/`; and
if you generated the matrices on a lab machine, no transfer is needed.)

### 3. Run

On the lab machine (logged in via `setup.md` section 2), from
`~/Paralela`:

```sh
export HWLOC_COMPONENTS=-gl
nohup sh bench/run_rand.sh > rand_job.log 2>&1 &
tail -f rand_job.log
```

The script:

- probes every machine in `hostfile.txt` and builds a temporary
  hostfile of the live ones, with the logged-in machine first (a hard
  mpirun requirement); set `SLOTS=` (default 8) if a machine has a
  different core count
- runs each `rand0_input<N>.bin` with its own grid, writing
  `rand0_output<N>.bin`, and zips the results into `rand0_output.zip`
- prints one `=== <file>: N=... grid=QxQ nproc=... ===` line per run,
  so the log shows which grid each size used

Note on scale: with N = 16384 the largest legal grid is 8x8 (64
processes), since 64 is the largest perfect square that divides 16384;
extra machines help the mid-size runs. Expect several hours for the
full set (the 16384 matrix alone is O(n^3) on ~1 GB of data). To run a
subset, zip fewer inputs, e.g. `make gen MATRICES_2_TO_THE="10 12 14"`.

When the lab is reserved (all machines down), re-run the script later
the same way; the log shows which machines were live.

### 4. Collect the results

On your own PC:

```sh
scp -o ProxyJump=upXXXXYYYYY@ssh.alunos.dcc.fc.up.pt \
    upXXXXYYYYY@l815.alunos.dcc.fc.up.pt:~/Paralela/rand0_output.zip .
```

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
