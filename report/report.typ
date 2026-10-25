#import "template.typ": report, todo, note-box, step-block

#show: report.with(
  title: "All-Pairs Shortest Paths with Fox's Algorithm",
  subtitle: "Min-plus repeated squaring distributed with MPI",
  authors: (
    (name: "Antonio Lanção",    id: "up202205416"),
    (name: "David Sá",          id: "up202303580"),
    (name: "Henrique Teixeira", id: "up202306640"),
  ),
  course: "Parallel Computing 2026/27",
  institution: "Faculty of Sciences · University of Porto",
  abstract: [
    We solve the all-pairs shortest path problem for directed weighted graphs by
    repeated squaring of the distance matrix under the min-plus semiring. Each
    product is distributed over a $Q times Q$ process grid with Fox's
    algorithm, implemented in C with MPI. The communication of each step
    overlaps the computation of the previous one, and the iteration stops as soon
    as the distances no longer change. #todo[one sentence with the headline result.]
  ],
)

#pagebreak()

= Introduction

Given a directed graph with $N$ nodes and non-negative edge weights, the
all-pairs shortest path (APSP) problem asks for the length of the shortest path
between every ordered pair of nodes. The graph is given as an $N times N$
adjacency matrix $W$, where $W_(i j)$ is the weight of the edge $i -> j$ and
$0$ off the diagonal means there is no edge.

This report describes a sequential baseline and a parallel implementation
using Fox's algorithm over MPI, how both were validated, and how they perform
on the LabCC cluster.

= Algorithm

== Min-plus product

The min-plus product of two $N times N$ matrices replaces the sum by $min$
and the product by $+$:

$ (A ⊗ B)_(i j) = min_(k = 1)^N (A_(i k) + B_(k j)) $ <minplus>

If $D^((m))$ holds the shortest distances using at most $m$ edges, then
$D^((2m)) = D^((m)) ⊗ D^((m))$. Starting from $D^((1)) = W$ (with
$infinity$ for missing edges and $0$ on the diagonal), $ceil(log_2 (N - 1))$
squarings are enough, since a shortest path has at most $N - 1$ edges. The
total cost is $O(N^3 log N)$, a $log N$ factor above Floyd--Warshall, but every
step is a matrix product, which parallelises well.

Both implementations also stop early: if a squaring leaves the matrix
unchanged, every later squaring would too, so the current matrix is
already the answer.

== Fox's algorithm <fox-alg>

The $P = Q^2$ processes form a $Q times Q$ grid, and process $(i, j)$ owns
block $(i, j)$ of $A$, $B$ and $C$, each of size $N slash Q times N slash Q$.
Fox's algorithm computes $C = A ⊗ B$ in $Q$ steps:

```
for s = 0 .. Q-1:
    k = (i + s) mod Q
    broadcast A[i][k] along row i          # process (i, k) is the root
    C[i][j] = min(C[i][j], A[i][k] (x) B[k][j])
    shift B one row up (circularly) along column j
```

After step $s$, process $(i, j)$ holds $B_((i + s + 1) mod Q, j)$, which is
exactly the block it needs for the next value of $k$. Each process performs
$Q dot (N slash Q)^3 = N^3 slash P$ operations per product and, per step,
receives one $A$ block through a row broadcast and one $B$ block from its
neighbour.

#pagebreak()

= Implementation

The code is split into a library without MPI (`lib/common`, `lib/seq`) and the
distributed part (`lib/mpi`). A single `src/main.c` is compiled three times:
as `seq`, as `fox` (with `-DUSE_MPI`) and as `fox_shared` (with
`-DUSE_MPI -DFOX_SHARED`). The three programs share the input reader, the
output writer and the min-plus kernel, and differ only in how the squarings
are organised.

== Data representation <data-repr>

Distances are stored as `uint32_t` in a row-major array. Missing edges are
stored as $"INF" = (2^32 - 1) slash 2$, so the sum of two values that are each
at most $"INF"$ cannot overflow. The programs read and write a binary format:
$N$ as a 4-byte unsigned integer, followed by the $N^2$ entries in row-major
order, with missing edges already stored as $"INF"$. The reader rejects a
matrix whose diagonal is not $0$. For the text format of the assignment,
reading maps off-diagonal zeros to $"INF"$, and printing maps $"INF"$ back
to $0$.

== Sequential version <seq-impl>

`seq` reads the matrix, runs `minplus_apsp` and writes the result. Only
`minplus_apsp` is timed. It keeps two $N times N$ buffers, the current
distances $D$ and the next product:

```
D = W
for covered = 1; covered + 1 < N; covered *= 2:
    Next = D (x) D            # minplus_multiply_into
    if Next == D: break       # early termination
    swap(D, Next)             # pointer swap, no copy
return D
```

`covered` is the largest number of edges that the paths in $D$ can have.
It doubles with each squaring, so the loop runs at most
$ceil(log_2 (N - 1))$ times. The product is written into the second buffer
because the min-plus product cannot be computed in place: later entries
of the product still need the old values of $D$.

=== Min-plus kernel

The min-plus kernel, `minplus_relax`, accumulates

$ C_(i j) <- min(C_(i j), A_(i k) + B_(k j)). $

`minplus_multiply_into` fills
$C$ with $"INF"$ and calls it once. Fox calls it once per step on the same
$C$, which builds the full product. Three choices make it fast:

- *Loop order $i$--$k$--$j$.* The inner loop walks along a row of $B$ and a row
  of $C$, both contiguous, so the compiler vectorises it with SIMD `min`
  instructions. With `-O3 -march=native` this is 8 lanes on AVX2.
- *Tiling.* The three loops are blocked in $128 times 128$ tiles, so the parts
  of $A$, $B$ and $C$ in use stay in cache.
- *Skipping unreachable rows.* If $A_(i k) = "INF"$, row $k$ of $B$ cannot
  improve row $i$ of $C$ and is skipped. Sparse inputs benefit most.

The `restrict` qualifiers tell the compiler that the three matrices do not
alias. `floyd_warshall` is kept only as an independent reference for the unit
tests.

== Fox with scatter, broadcast and shift <fox-impl>

`fox` follows the textbook algorithm in three stages. Rank 0 distributes the
blocks once, every product moves $B$ by point-to-point shifts and $A$ by row
broadcasts, and rank 0 collects the result at the end:

#step-block(1, [Scatter])[
  Rank 0 reads $W$ and packs its blocks in rank order. One `MPI_Scatter`
  gives process $(i, j)$ its block $A = W_(i j)$, and $B$ starts as a copy
  of $A$.
]
#step-block(2, [Squarings --- $Q$ steps of `Irecv` $B$, `Isend` $B$, `Ibcast` $A$])[
  Each squaring is a Fox product of $Q$ steps. In every step the process
  receives the next $B$ from the process below, sends its current $B$ to the
  process above, takes part in the row broadcast of the next $A$, and
  computes with the blocks it already has. An `MPI_Allreduce` then decides
  whether to stop; otherwise $A$ takes the new block and $B$ a copy of it,
  since the product rotates $B$.
]
#step-block(3, [Gather])[
  One `MPI_Gather` collects the blocks in rank order, and rank 0 unpacks
  them into $W$.
]

#block(breakable: false)[
In code (`mpi_apsp`):

```
rank 0: read W, pack its blocks in rank order
MPI_Scatter                   # process (i, j) receives A = W[i][j]
B = A
repeat at most ceil(log2(N - 1)) times:
    C = fox_multiply(A, B)
    changed = MPI_Allreduce(C != A, LOR)
    if not changed: break
    swap(A, C);  B = A        # B is rotated during the product, so it gets a copy
MPI_Gather                    # rank 0 unpacks the blocks into W
```
]

Each product (`fox_multiply`) runs the $Q$ steps of @fox-alg, with the
transfers for step $s + 1$ started before the computation of step $s$:

```
MPI_Bcast A[i][i] along row i                          # step 0
C = INF
for s = 0 .. Q-1:
    if s + 1 < Q:
        MPI_Irecv  B_next from (i + 1) mod Q along column j
        MPI_Isend  B_cur  to   (i - 1) mod Q along column j
        MPI_Ibcast A[i][(i + s + 1) mod Q] along row i  -> A_next
    C = min(C, A_cur (x) B_cur)                        # minplus_relax
    MPI_Waitall
    swap(B_cur, B_next);  A_cur = A_next
```

=== Process grid and distribution

Rank 0 reads the matrix and broadcasts $N$. Every rank then checks that $P$ is a
perfect square and that $Q$ divides $N$. Otherwise rank 0 prints
`ERROR: Invalid configuration!` and all ranks exit. The grid is a periodic
two-dimensional Cartesian communicator (`MPI_Cart_create`, created without
reordering so that grid rank 0 is the rank that read the input), and
`MPI_Cart_sub` derives one communicator per row and one per column. The
broadcasts of $A$ use the row communicator, and the shifts of $B$ use the
column communicator.

Blocks are not contiguous in the global matrix, so rank 0 packs them in rank
order into a buffer: block $(i, j)$ becomes the $(i Q + j)$-th contiguous
$N slash Q times N slash Q$ piece. A single `MPI_Scatter` then gives each
process its block. At the end, `MPI_Gather` collects the blocks in the same
order, and rank 0 unpacks them in place into the input matrix.

=== Overlapping communication and computation

The transfers of step $s + 1$ run while step $s$ computes:

- the next $A$ block is broadcast with a non-blocking `MPI_Ibcast` into one of
  two alternating buffers, while the current one is in use;
- the $B$ shift is an `MPI_Isend` of the current block to the process above
  and an `MPI_Irecv` of the next one from the process below, into a spare
  buffer. The two buffers are swapped after `MPI_Waitall`.

The computation only reads the buffers being sent, which MPI allows, so the
step costs roughly $max(t_"comp", t_"comm")$ instead of their sum. Only the
broadcast of step 0 is blocking, since nothing can be computed before it
arrives. The process that is the broadcast root sends directly from its own
$A$ block, without a copy. All buffers are allocated once, in `mpi_apsp`, and
reused across squarings. Each process holds six blocks ($A$, $B$, $C$, two
$A$ receive buffers and the spare $B$), i.e. $6 (N slash Q)^2 dot 4$ bytes.

=== Early termination

After each product, every process compares its new block with the previous one,
and an `MPI_Allreduce` with `MPI_LOR` decides whether any block changed. The
cost is one small collective per squaring, which can save several full
products. `fox_shared` uses the same test.

== Shared-copy variant with `MPI_Allgather` <shared-impl>

`fox_shared` removes the communication from inside the product. Every
process keeps a full copy $M$ of the matrix, stored block by block in the same
rank order as the scatter buffer of `fox`. Because every process knows where
each block is, step $s$ of the product reads $A_(i k)$ and $B_(k j)$, with
$k = (i + s) mod Q$, straight from memory (`fox_block` returns a view, not a
copy), instead of receiving them by row broadcast and column shift. The
communication moves to the end of each squaring:

#step-block(1, [Broadcast])[
  Rank 0 reads $W$ and packs its blocks in rank order into $M$. One
  `MPI_Bcast` gives every process the whole matrix.
]
#step-block(2, [Squarings --- local product, then `MPI_Allgather`])[
  Each process computes its block $C$ in $Q$ steps with no communication.
  After the `MPI_Allreduce` stop test, one `MPI_Allgather` puts every
  process's new block into every copy of $M$.
]
#step-block(3, [Unpack])[
  Rank 0 already holds the whole result and unpacks $M$ into $W$; no
  gather is needed.
]

#block(breakable: false)[
In code:

```
rank 0: read W, pack its blocks in rank order -> M
MPI_Bcast M                   # every process gets the whole matrix
repeat at most ceil(log2(N - 1)) times:
    C = INF
    for s = 0 .. Q-1:
        k = (i + s) mod Q
        C = min(C, M[i][k] (x) M[k][j])
    changed = MPI_Allreduce(C != M[i][j], LOR)
    if not changed: break
    MPI_Allgather C -> M      # block of rank r lands at piece r
rank 0: unpack M into W
```
]

The block layout is what makes the all-gather work without packing: the
piece that rank $r$ contributes is exactly where block $(r slash Q, r mod Q)$
belongs.

The trade-off is volume and memory, summarised in @variants. Per
squaring, `fox` delivers $2 Q$ blocks to each process, $2 N^2 slash Q$
values, in $2 Q$ messages that overlap the computation. `fox_shared`
delivers $N^2$ values to each process in one collective that does not
overlap anything, $Q slash 2$ times more data, and each process holds
$N^2 + (N slash Q)^2$ values instead of $6 (N slash Q)^2$.

#figure(
  caption: [Communication and memory of the two parallel versions, per process.],
  table(
    columns: 3,
    align: (left, left, left),
    table.header([], [`fox`], [`fox_shared`]),
    [start],        [`MPI_Scatter`, $N^2 slash P$ values], [`MPI_Bcast`, $N^2$ values],
    [each step],    [`MPI_Isend`/`MPI_Irecv` of $B$, `MPI_Ibcast` of $A$, overlapped], [none],
    [each squaring], [$2 Q$ blocks, $2 N^2 slash Q$ values], [`MPI_Allgather`, $N^2$ values, not overlapped],
    [end],          [`MPI_Gather`, $N^2 slash P$ values], [none],
    [memory],       [$6 (N slash Q)^2$ values], [$N^2 + (N slash Q)^2$ values],
  ),
) <variants>

= Validation

- *Unit tests* (`make test`): matrix allocation and copy, input parsing
  (including malformed input), and the min-plus kernel. The kernel is checked
  against the assignment's example, the identity matrix, all-INF matrices
  (overflow), split accumulation, and Floyd--Warshall on a random graph with
  $N = 150$, which is not a multiple of the tile size.
- *Regression* (`make check`): `seq` and `fox` are compared byte for byte with
  the reference outputs for $N in {6, 300, 600, 900, 1200, 3000}$, with
  $P in {1, 4, 9, 16}$. Configurations where $Q$ does not divide $N$ must
  be rejected with exactly the expected error message.

All tests pass.

#pagebreak()

= Random inputs <matgen>

The example inputs of the assignment stop at $N = 3000$ and are not powers
of two, so a $Q times Q$ grid only fits some of them. For the benchmarks we
generate our own inputs with `tests/matgen`.

== How `matgen` works

`matgen N [count] [seed]` writes `count` matrices of size $N times N$,
named `rand`$k$`_input`$N$`.bin`, in the binary format of
@data-repr. It seeds the C generator with `seed` (the current time
if omitted) and fills the matrix in row-major order:

- the diagonal is $0$;
- every other entry is an edge with probability 10%, otherwise $"INF"$;
- an edge gets a weight drawn uniformly from $[0, 2^16]$.

`rand()` may only produce 15 random bits, so each draw combines two calls
into a 30-bit value before reducing it to $[0, 2^16]$. The weight bound
keeps every distance far from $"INF"$: a shortest path has at most
$N - 1$ edges, so for $N = 2^14$ it costs at most
$(2^14 - 1) dot 2^16 approx 1.07 dot 10^9$, below
$"INF" approx 2.15 dot 10^9$. `make gen` runs `matgen` for
$N = 2^6, ..., 2^14$ with seed 42 and zips the files, so every machine
benchmarks exactly the same matrices.

== Why 90% of the entries are $"INF"$ <density>

The density decides what the input looks like as a graph, so we measured
its effect before choosing it. @density-table runs the sequential
`minplus_apsp` on matrices generated like `matgen`'s, with seed 42, at
several densities.

#figure(
  caption: [Effect of edge density on the sequential algorithm: squarings
    performed (including the last, which only confirms convergence) out of
    $ceil(log_2 (N - 1))$, pairs left unreachable, and time on a laptop
    (AMD Ryzen 7 7840HS).],
  table(
    columns: 5,
    align: (right, right, right, right, right),
    table.header([$N$], [edges], [squarings], [unreachable], [time (ms)]),
    [64],   [100%], [5 / 6],  [0%],    [0.2],
    [64],   [10%],  [5 / 6],  [0%],    [0.1],
    [64],   [1%],   [4 / 6],  [95.9%], [0.0],
    [256],  [100%], [5 / 8],  [0%],    [8.9],
    [256],  [10%],  [5 / 8],  [0%],    [8.9],
    [256],  [1%],   [6 / 8],  [18.1%], [6.2],
    [2048], [100%], [6 / 11], [0%],    [4996],
    [2048], [10%],  [6 / 11], [0%],    [4572],
    [2048], [1%],   [6 / 11], [0%],    [3636],
    [2048], [0.1%], [6 / 11], [33.1%], [1477],
  ),
) <density-table>

#note-box(label: "Why 10% edges")[
  - *Every pair is still connected.* At 10% no pair is unreachable for any
    $N >= 64$, so the output is a full distance matrix. At 1% a graph with
    64 nodes leaves 96% of the pairs unreachable, and its "result" is
    mostly $"INF"$.
  - *The work is not made easier.* A 10% matrix needs as many squarings as
    a complete graph and, at $N = 2048$, takes 92% of its time. Sparser
    inputs finish faster only because the kernel skips $"INF"$ entries in
    the first squarings, which would flatter the measurements.
  - *Missing edges are exercised.* A complete graph never contains
    $"INF"$, so the overflow guard and the skip in `minplus_relax` would
    never run. With 90% $"INF"$ the first squarings meet them everywhere.
]

Two consequences follow. The number of squarings barely depends on
density, between 5 and 6 here, so early termination saves about half of
the $ceil(log_2 (N - 1))$ products on these inputs. And since the
weights can be $0$, `matgen` may produce zero-weight edges; the binary
format keeps them as edges, unlike the text format, where $0$ means no
edge.

#pagebreak()

= Performance evaluation

== Methodology

#todo[LabCC node specification: CPU, cores per node, memory, network,
  compiler and MPI versions, and how processes were mapped to nodes.]

Times are wall-clock times of the computation only (scatter, squarings and
gather), measured on rank 0 with `MPI_Wtime`. The clock starts after an
`MPI_Barrier` and stops once the gather returns.
Input and output are excluded. Each configuration was run
#todo[R] times and the minimum is reported. The benchmark is driven by
`bench/run_bench.sh`. Speedup and efficiency are computed relative to the
sequential program:

$ S(P) = T_"seq" / T(P), quad E(P) = S(P) / P $

== Results

#figure(
  caption: [Execution time in milliseconds.],
  table(
    columns: 7,
    align: (left, ..range(6).map(_ => right)),
    [$N$], [seq], [$P = 1$], [$P = 4$], [$P = 9$], [$P = 16$], [$P = 25$],
    [300],  [], [], [], [], [], [],
    [600],  [], [], [], [], [], [],
    [900],  [], [], [], [], [], [],
    [1200], [], [], [], [], [], [],
    [3000], [], [], [], [], [], [],
  ),
) <times>

#figure(
  caption: [Speedup relative to the sequential version.],
  table(
    columns: 6,
    align: (left, ..range(5).map(_ => right)),
    [$N$], [$P = 1$], [$P = 4$], [$P = 9$], [$P = 16$], [$P = 25$],
    [300],  [], [], [], [], [],
    [600],  [], [], [], [], [],
    [900],  [], [], [], [], [],
    [1200], [], [], [], [], [],
    [3000], [], [], [], [], [],
  ),
) <speedup>

== Random matrices on LabCC

@lab-rand0-plot and @lab-rand0 show `seq` and `fox` on the random matrices
`rand0_input`$N$ of @matgen,
$N = 2^6, ..., 2^14$, run on LabCC with `bench/run_bench.sh` and a hostfile
of 29 nodes with 8 slots each. Each configuration ran twice and the
minimum is reported. Of the requested
$P in {4, 9, 25, 36, 64, 81, 100, 144, 196}$, only $P = 4$ ($Q = 2$) and
$P = 64$ ($Q = 8$) have a $Q$ that divides a power of two; the others were
skipped. Processes fill the hostfile slot by slot, so $P = 4$ runs on a single
node and $P = 64$ on eight full nodes.

#figure(
  image("lab_rand0_speedup.png", width: 100%),
  caption: [Left: minimum time per $N$ on a log scale, `seq` dashed. Right:
    speedup against `seq`; below the dashed line, `fox` is slower than
    `seq`. Speedups under 0.5 are left unlabelled; they are in @lab-rand0.],
) <lab-rand0-plot>

#figure(
  caption: [Random matrices on LabCC: minimum of 2 runs in milliseconds,
    speedup $S$ and efficiency $E$ relative to `seq`.],
  table(
    columns: 8,
    align: (left, ..range(7).map(_ => right)),
    stroke: (_, y) => (
      top: if y == 0 { 0.8pt } else if y == 2 { 0.4pt } else { 0pt },
    ),
    table.header(
      [], [], table.cell(colspan: 3, align: center)[$P = 4$],
      table.cell(colspan: 3, align: center)[$P = 64$],
      [$N$], [seq], [time], [$S$], [$E$], [time], [$S$], [$E$],
    ),
    [64],    [0.137],     [0.197],    [0.70], [0.17], [12.9],     [0.01], [0.000],
    [128],   [0.732],     [0.527],    [1.39], [0.35], [36.2],     [0.02], [0.000],
    [256],   [5.24],      [2.54],     [2.06], [0.52], [131.5],    [0.04], [0.001],
    [512],   [43.1],      [15.7],     [2.74], [0.69], [465.1],    [0.09], [0.001],
    [1024],  [464.6],     [129.6],    [3.58], [0.90], [2869],     [0.16], [0.003],
    [2048],  [4547],      [1027],     [4.43], [1.11], [9174],     [0.50], [0.008],
    [4096],  [39508],     [11157],    [3.54], [0.89], [39133],    [1.01], [0.016],
    [8192],  [349261],    [92488],    [3.78], [0.94], [153850],   [2.27], [0.035],
    [16384], [2964332],   [738863],   [4.01], [1.00], [592689],   [5.00], [0.078],
  ),
) <lab-rand0>

With four processes on one node, `fox` is slower than `seq` only for
$N = 64$, and from $N = 1024$ on its efficiency stays between 0.89 and 1.11.
The values slightly above 1 ($N = 2048$ and $N = 16384$) come from each
process working on an $N slash 2 times N slash 2$ block, which uses the
cache better than the full matrix in `seq`.

With 64 processes on eight nodes the picture is reversed. Up to $N = 2048$,
`fox` is slower than `seq`, by up to two orders of magnitude for small
$N$: each product has 8 steps, each with a row broadcast and a column shift
that now cross the network, and with $N slash 8 <= 256$ there is too little
computation per step to hide them. It breaks even at $N = 4096$, and only at
$N = 16384$ (blocks of $2048 times 2048$) does it beat $P = 4$, with a
speedup of 5.0 against 4.0. Even then the efficiency is 0.08, so a 16-fold
increase in processes buys a 1.25-fold reduction in time. $P = 64$ only pays
off for matrices far larger than those in the assignment.

The two runs of each configuration agree within 4% for $N >= 2048$; for
smaller $N$, where times are a few milliseconds, they differ by up to 32%.

== Fox against the shared-copy variant

@fox-vs-shared and @fox-vs-shared-times compare `fox` and `fox_shared` on
random matrices with $N = 2^6, ..., 2^10$ (`matgen`, seed 42), measured
with `bench/run_power2.sh` on a laptop with an AMD Ryzen 7 7840HS (8 cores,
16 hardware threads). Values are the median of 3 runs (6 for $P in {4, 16}$).
$P = 64$ oversubscribes the machine, four processes per hardware thread, so
those points show contention, not scaling.

#figure(
  image("power2_fox_vs_shared.png", width: 100%),
  caption: [Top: median time per $P$, one panel per $N$, each with its own
    scale. Bottom: `fox_shared` time divided by `fox` time; above 1,
    `fox_shared` is slower. Shaded: $P$ above the 16 hardware threads.],
) <fox-vs-shared>

#figure(
  caption: [Median execution time in milliseconds, `fox` against `fox_shared`.],
  table(
    columns: 9,
    align: (left, ..range(8).map(_ => right)),
    stroke: (_, y) => (
      top: if y == 0 { 0.8pt } else if y == 2 { 0.4pt } else { 0pt },
    ),
    table.header(
      [], table.cell(colspan: 2, align: center)[$P = 1$],
      table.cell(colspan: 2, align: center)[$P = 4$],
      table.cell(colspan: 2, align: center)[$P = 16$],
      table.cell(colspan: 2, align: center)[$P = 64$],
      [$N$], [fox], [shared], [fox], [shared], [fox], [shared], [fox], [shared],
    ),
    [64],   [0.18],  [0.13],  [0.19],  [0.15],  [0.47], [0.39],  [7.92], [2.59],
    [128],  [0.76],  [0.60],  [0.45],  [0.39],  [0.76], [0.62],  [4.67], [5.88],
    [256],  [5.78],  [5.37],  [1.93],  [1.96],  [2.23], [1.99],  [8.44], [14.6],
    [512],  [46.0],  [50.0],  [14.1],  [13.0],  [7.42], [11.3],  [18.5], [39.4],
    [1024], [541.3], [529.4], [128.6], [126.3], [81.9], [107.1], [94.6], [204.6],
  ),
) <fox-vs-shared-times>

For $N <= 256$, `fox_shared` is slightly faster (0.8--0.9 of the time of
`fox` at $P <= 16$). The blocks are tiny there, so the cost of `fox` is the
latency of its $2 Q$ messages per product rather than their size, and one
all-gather per squaring is cheaper. As $N$ grows the volume takes over: at
$P = 16$, `fox_shared` takes 1.52 times as long as `fox` for $N = 512$ and
1.31 times for $N = 1024$, which is the $Q slash 2 = 2$ times larger
transfer, not hidden behind computation. At $P = 4$ the two are within a few
percent, since $Q slash 2 = 1$: both move the same amount of data. At $P = 1$
both run the same single local product.

== Discussion

#todo[Discuss @times and @speedup: overhead of $P = 1$ against `seq`;
  scaling with $P$ for each $N$; where communication starts to dominate
  (small $N$, large $P$); effect of using several nodes; effect of early
  termination.]

#pagebreak()

= Conclusions

#todo[Summary of results and possible improvements.]
