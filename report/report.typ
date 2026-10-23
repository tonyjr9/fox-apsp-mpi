#import "template.typ": report, todo

#show: report.with(
  title: "All-Pairs Shortest Paths with Fox's Algorithm",
  subtitle: "Min-plus repeated squaring distributed with MPI",
  authors: (
    (name: "Antonio Lanção", id: "up202205416"),
    (name: "David Sá", id: "up202303580"),
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

== Fox's algorithm

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

= Implementation

The code is split into a library without MPI (`lib/common`, `lib/seq`) and the
distributed part (`lib/mpi`). A single `src/main.c` is compiled twice, as
`seq` and, with `-DUSE_MPI`, as `fox`.

== Data representation

Distances are stored as `uint32_t` in a row-major array. Missing edges are
stored as $"INF" = (2^32 - 1) slash 2$, so the sum of two values that are each
at most $"INF"$ cannot overflow. On output, $"INF"$ is printed as $0$, as the
assignment requires. Reading the input forces the diagonal to $0$ and maps
off-diagonal zeros to $"INF"$.

== Sequential kernel

The min-plus kernel (`minplus_relax`) accumulates
$C_(i j) <- min(C_(i j), A_(i k) + B_(k j))$, so calling it once per Fox step
builds the full product. Three choices make it fast:

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

== Process grid and distribution

Rank 0 reads the matrix and broadcasts $N$. Every rank then checks that $P$ is a
perfect square and that $Q$ divides $N$. Otherwise rank 0 prints
`ERROR: Invalid configuration!` and all ranks exit. The grid is a periodic
two-dimensional Cartesian communicator (`MPI_Cart_create`, created without
reordering so that grid rank 0 is the rank that read the input), and
`MPI_Cart_sub` derives one communicator per row and one per column.

Blocks are not contiguous in the global matrix, so rank 0 packs them in rank
order into a buffer and distributes them with a single `MPI_Scatter`. The result
is collected with `MPI_Gather` and unpacked in place into the input matrix.

== Overlapping communication and computation

The transfers for step $s + 1$ start before the computation of step $s$
(`fox_multiply`):

- the next $A$ block is broadcast with a non-blocking `MPI_Ibcast` into one of
  two alternating buffers, while the current one is in use;
- the $B$ shift is an `MPI_Isend` of the current block and an `MPI_Irecv` into
  a spare buffer, and the two buffers are swapped after `MPI_Waitall`.

The computation only reads the buffers being sent, which MPI allows, so the
step costs roughly $max(t_"comp", t_"comm")$ instead of their sum. The
process that is the broadcast root sends directly from its own $A$ block,
without a copy. All buffers are allocated once, in `mpi_apsp`, and reused
across squarings. Each process holds six blocks, i.e.
$6 (N slash Q)^2 dot 4$ bytes.

== Early termination

After each product, every process compares its new block with the previous one,
and an `MPI_Allreduce` with `MPI_LOR` decides whether any block changed. The
cost is one small collective per squaring, which can save several full
products. #todo[report how many squarings each input actually needs.]

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

== Discussion

#todo[Discuss @times and @speedup: overhead of $P = 1$ against `seq`;
  scaling with $P$ for each $N$; where communication starts to dominate
  (small $N$, large $P$); effect of using several nodes; effect of early
  termination.]

= Conclusions

#todo[Summary of results and possible improvements.]
