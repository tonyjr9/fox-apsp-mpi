#import "template.typ": report, todo, note-box

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
    We solve the all-pairs shortest path problem by repeated squaring of the
    distance matrix under the min-plus semiring, with each product distributed
    over a $Q times Q$ process grid by Fox's algorithm in C and MPI. On LabCC,
    four processes on one node reach a speedup of 3.9 on the assignment inputs
    and 4.0 on a random matrix with $N = 16384$. Once the grid spans several
    nodes, the 1 Gb/s network dominates, and more processes only help on the
    largest random matrices.
  ],
)

= Algorithm

The graph is an $N times N$ matrix $W$, with $"INF"$ for missing edges and
$0$ on the diagonal. The min-plus product replaces the sum by $min$ and the
product by $+$:

$ (A ⊗ B)_(i j) = min_k (A_(i k) + B_(k j)) $

If $D$ holds the shortest distances over paths of at most $m$ edges, $D ⊗ D$
holds them for $2m$ edges, so $ceil(log_2 (N - 1))$ squarings of $W$ give
the answer, in $O(N^3 log N)$. All versions stop early once a squaring leaves
the matrix unchanged.

= Implementation

One `src/main.c` is compiled as `seq`, `fox` (`-DUSE_MPI`) and `fox_shared`
(`-DUSE_MPI -DFOX_SHARED`); the three share the I/O and the min-plus kernel.
Distances are `uint32_t` in row-major order, with
$"INF" = (2^32 - 1) slash 2$ so that adding two of them cannot overflow.
Inputs are binary: $N$, then the $N^2$ entries.

== Sequential version

`seq` times only `minplus_apsp`, which squares into a second buffer (the
product cannot be computed in place) and swaps pointers:

```
D = W
for covered = 1; covered + 1 < N; covered *= 2:
    Next = D (x) D            # minplus_multiply_into
    if Next == D: break       # early termination
    swap(D, Next)
```

The kernel `minplus_relax` accumulates
$C_(i j) <- min(C_(i j), A_(i k) + B_(k j))$, so Fox can call it once per
step on the same $C$. It runs the loops in $i$--$k$--$j$ order, so the inner
loop walks contiguous rows and vectorises into SIMD `min`; tiles them in
$128 times 128$ blocks to stay in cache; and skips row $k$ of $B$ when
$A_(i k) = "INF"$.

== Fox: scatter, $("Irecv" B, "Isend" B, "Ibcast" A)^*$, gather

Rank 0 checks that $P = Q^2$ and that $Q$ divides $N$, builds a periodic
Cartesian grid, and derives one communicator per row and per column. It packs
the blocks of $W$ in rank order, so one `MPI_Scatter` hands out the blocks
and one `MPI_Gather` collects them at the end. Within each product, the
transfers for step $s + 1$ are started before step $s$ computes, so a step
costs about $max(t_"comp", t_"comm")$:

```
MPI_Scatter W -> A;  B = A
repeat:                                         # one squaring
    MPI_Bcast A[i][i] along row i
    C = INF
    for s = 0 .. Q-1:
        if s + 1 < Q:
            MPI_Irecv  B_next from (i + 1) mod Q along column j
            MPI_Isend  B      to   (i - 1) mod Q along column j
            MPI_Ibcast A[i][(i + s + 1) mod Q] along row i -> A_next
        C = min(C, A_cur (x) B)                 # minplus_relax
        MPI_Waitall;  swap(B, B_next);  A_cur = A_next
    if not MPI_Allreduce(C != A, LOR): break    # early termination
    swap(A, C);  B = A
MPI_Gather A -> W
```

The computation only reads the buffers being sent, and the broadcast root
sends straight from its own block. Each process holds six blocks: $A$, $B$,
$C$, two $A$ receive buffers and a spare $B$.

== Shared copy: broadcast, (local product, `MPI_Allgather`)#super[$*$]

`fox_shared` takes the communication out of the product. Every process keeps
the whole matrix $M$, stored block by block in rank order, and reads
$A_(i k)$ and $B_(k j)$ straight from it. After each squaring, one
`MPI_Allgather` writes every new block into every copy; the rank-order layout
puts rank $r$'s block exactly where it belongs, with no packing:

```
MPI_Bcast M
repeat:
    C = INF
    for s = 0 .. Q-1:  k = (i + s) mod Q;  C = min(C, M[i][k] (x) M[k][j])
    if not MPI_Allreduce(C != M[i][j], LOR): break
    MPI_Allgather C -> M
```

Per squaring, `fox` moves $2 N^2 slash Q$ values to each process and overlaps
them with computation. `fox_shared` moves $N^2$, $Q slash 2$ times more,
without overlap, and needs $N^2$ values of memory per process instead of
$6 (N slash Q)^2$.

== Validation

Unit tests cover I/O and the kernel, including all-$"INF"$ matrices and a
comparison with Floyd--Warshall. `seq` and `fox` match the reference outputs
byte for byte for $N in {6, 300, ..., 3000}$ with $P in {1, 4, 9, 16}$, and
invalid grids are rejected with the expected message.

= Random inputs

The assignment inputs are not powers of two, so we generate matrices with
`matgen N [count] [seed]`. The diagonal is $0$, and every other entry is an
edge with probability 10%, weighted uniformly in $[0, 2^16]$, and $"INF"$
otherwise. Two `rand()` calls give 30 random bits per draw. A path of at most
$N - 1$ such edges stays below $1.07 dot 10^9 < "INF"$ for $N <= 2^14$.
`make gen` writes $N = 2^6, ..., 2^14$ with seed 42.

#grid(
  columns: (55%, 1fr),
  column-gutter: 14pt,
  align: (top, horizon),
  [#figure(
    caption: [Edge density against squarings done (of
      $ceil(log_2 (N - 1))$), unreachable pairs and `seq` time, on an AMD
      Ryzen 7 7840HS.],
    table(
      columns: 5,
      align: (right, right, right, right, right),
      table.header([$N$], [edges], [squar.], [unreach.], [ms]),
      [64],   [100%], [5 / 6],  [0%],    [0.2],
      [64],   [10%],  [5 / 6],  [0%],    [0.1],
      [64],   [1%],   [4 / 6],  [95.9%], [0.0],
      [2048], [100%], [6 / 11], [0%],    [4996],
      [2048], [10%],  [6 / 11], [0%],    [4572],
      [2048], [0.1%], [6 / 11], [33.1%], [1477],
    ),
  ) <density>],
  note-box(label: "Why 90% INF")[
    At 10% edges every pair is still reachable for $N >= 64$, while at 1% a
    64-node graph leaves 96% of the pairs unconnected. The work is not made
    easier: 10% needs as many squarings as a complete graph and 92% of its
    time at $N = 2048$. And unlike a complete graph, it exercises $"INF"$
    handling and the kernel's skip (@density).
  ],
)

= Performance evaluation

Times cover the computation only (scatter, squarings, gather), measured on
rank 0 with `MPI_Wtime` after a barrier. We report the minimum of 3 runs
on the assignment inputs and of 2 on the random matrices, with speedups $S = T_"seq" slash T(P)$
against `seq` and $S_1 = T(1) slash T(P)$ against `fox` with one process,
and efficiency $E = S slash P$. Each LabCC node has an AMD Ryzen 7 5700G (8 cores, 16
threads, 16 MB L3) and 14 GB of memory, runs Debian 12, and is connected by
1 Gb/s Ethernet; we used Open MPI 4.1.4 and GCC 12.2.

== Assignment inputs on LabCC

All five assignment inputs fit $P in {1, 4, 9, 16, 25}$. The hostfile lists
29 nodes with 8 slots each, filled in order, so $P = 4$ runs on one node,
$P = 9$ and $P = 16$ on two, and $P = 25$ on four
(`bench/report/lab_examples.csv`). #todo[$P = 25$ at $N = 3000$ is still running.]

#figure(
  image("lab_examples_speedup.png", width: 74%),
  caption: [Time (log scale) and speedup on the assignment inputs on LabCC.
    Speedups below 0.5 are left unlabelled; they are in @lab-speedup.],
) <lab-examples-plot>

#figure(
  caption: [Assignment inputs on LabCC: minimum of 3 runs, in milliseconds.],
  table(
    columns: 7,
    align: (left, ..range(6).map(_ => right)),
    table.header([$N$], [seq], [$P = 1$], [$P = 4$], [$P = 9$], [$P = 16$], [$P = 25$]),
    [300],  [14.7],   [16.0],   [6.9],   [170.2],  [257.5],  [4628],
    [600],  [145.5],  [148.5],  [44.3],  [734.5],  [1159.8], [24145],
    [900],  [430.6],  [433.0],  [139.5], [1820.5], [3151.9], [91880],
    [1200], [1238.9], [1221.0], [370.6], [3729.3], [5563.7], [225073],
    [3000], [10249.0], [10334.4], [2663.4], [13346.9], [18320.4], [---],
  ),
) <lab-times>

#figure(
  caption: [Speedup $S$ against `seq` and $S_1$ against `fox` with $P = 1$,
    from @lab-times.],
  table(
    columns: 9,
    align: (left, ..range(8).map(_ => right)),
    stroke: (_, y) => (
      top: if y == 0 { 0.8pt } else if y == 2 { 0.4pt } else { 0pt },
    ),
    table.header(
      [], table.cell(colspan: 2, align: center)[$P = 4$],
      table.cell(colspan: 2, align: center)[$P = 9$],
      table.cell(colspan: 2, align: center)[$P = 16$],
      table.cell(colspan: 2, align: center)[$P = 25$],
      [$N$], [$S$], [$S_1$], [$S$], [$S_1$], [$S$], [$S_1$], [$S$], [$S_1$],
    ),
    [300],  [2.12], [2.31], [0.086], [0.094], [0.057], [0.062], [0.003], [0.003],
    [600],  [3.29], [3.35], [0.20],  [0.20],  [0.13],  [0.13],  [0.006], [0.006],
    [900],  [3.09], [3.10], [0.24],  [0.24],  [0.14],  [0.14],  [0.005], [0.005],
    [1200], [3.34], [3.30], [0.33],  [0.33],  [0.22],  [0.22],  [0.006], [0.005],
    [3000], [3.85], [3.88], [0.77],  [0.77],  [0.56],  [0.56],  [---],   [---],
  ),
) <lab-speedup>

With $P = 1$, `fox` is within 2% of `seq` from $N = 600$; at $N = 300$, a
15 ms run, it is 9% slower, and $S_1$ is correspondingly higher than $S$.
On one node, four processes reach a speedup of 3.1--3.9 from $N = 600$,
up to 3.85 at $N = 3000$ (efficiency 0.96). From $P = 9$ the grid spans two nodes, and the
picture changes completely: with a single process on the second node,
$P = 9$ is 5--25 times slower than $P = 4$ and slower than `seq` for every
input, though the gap narrows as $N$ grows (0.77 at $N = 3000$). Each Fox step synchronises a row and a column, so the whole grid
waits for the transfers over the 1 Gb/s link. $P = 16$ is slower still.
$P = 25$ takes minutes and varies up to 3.7 times between runs (24 to 90 s
at $N = 600$), which points to other load on the shared nodes during the
run, so we draw no conclusion from it. As the random matrices below show,
two nodes need $N >= 4096$ to break even with `seq`, beyond the assignment
inputs.

== Assignment inputs on an Apple M1 Max

For comparison we ran the same inputs on a MacBook with an Apple M1 Max
(10 cores: 8 performance, 2 efficiency) and 32 GB of memory
(`bench/report/mac_examples.csv`). With only 10 cores it can run $P = 1$, 4
and 9, but not 16 or 25.

#figure(
  image("mac_examples_speedup.png", width: 74%),
  caption: [Time (log scale) and speedup on the assignment inputs on the
    M1 Max.],
) <mac-plot>

#figure(
  caption: [M1 Max: minimum of 3 runs in ms, speedup $S$ against `seq` and
    $S_1$ against `fox` with $P = 1$.],
  table(
    columns: 9,
    align: (left, ..range(8).map(_ => right)),
    stroke: (_, y) => (
      top: if y == 0 { 0.8pt } else if y == 2 { 0.4pt } else { 0pt },
    ),
    table.header(
      [], [], [], table.cell(colspan: 3, align: center)[$P = 4$],
      table.cell(colspan: 3, align: center)[$P = 9$],
      [$N$], [seq], [$P = 1$], [time], [$S$], [$S_1$], [time], [$S$], [$S_1$],
    ),
    [300],  [24.0],    [23.7],    [34.2],   [0.70], [0.69], [37.9],   [0.63], [0.63],
    [600],  [154.8],   [154.6],   [77.9],   [1.99], [1.98], [127.9],  [1.21], [1.21],
    [900],  [551.9],   [545.0],   [186.5],  [2.96], [2.92], [229.1],  [2.41], [2.38],
    [1200], [1363.2],  [1373.9],  [377.5],  [3.61], [3.64], [386.3],  [3.53], [3.56],
    [3000], [11658.2], [11749.5], [3114.7], [3.74], [3.77], [1849.1], [6.30], [6.35],
  ),
) <mac-table>

All processes share one machine here, so $P = 9$ does not suffer the
network penalty it has on LabCC: it is slower than $P = 4$ up to
$N = 1200$, with three communication steps per product instead of two and
one process on an efficiency core, but reaches 6.30 at $N = 3000$.

== Random matrices on LabCC

Here $P = 4$ runs on one node, $P = 16$ on two and $P = 64$ on eight; these
are the process counts whose $Q$ divides a power of two
(`bench/report/lab_rand0.csv`).

#figure(
  image("lab_rand0_speedup.png", width: 74%),
  caption: [Time (log scale) and speedup on the random matrices on LabCC.
    Labels close to one already shown are omitted.],
) <lab-plot>

With four processes, efficiency stays between 0.89 and 1.11 from
$N = 1024$, slightly above 1 where the smaller blocks fit the cache better.
Across nodes, $P = 16$ and $P = 64$ take about the same time (within 15%
from $N = 128$ to 8192): the network, not the process count, sets the time. Both are slower
than `seq` up to $N = 2048$ and break even at $N = 4096$. At $N = 16384$,
$P = 16$ reaches 3.88, still below $P = 4$ (4.01), and $P = 64$ reaches
5.00, with efficiencies of 0.24 and 0.08.

== Fox against the shared copy

#figure(
  image("power2_fox_vs_shared.png", width: 80%),
  caption: [`fox` against `fox_shared` on an AMD Ryzen 7 7840HS (16 hardware
    threads; shaded $P$ oversubscribes it), median of 3 to 6 runs. Bottom:
    `fox_shared` time over `fox` time.],
) <shared-plot>

For $N <= 256$, `fox_shared` is mostly slightly faster, down to 0.8 of the
time of `fox`: the blocks are tiny, so one all-gather beats $2 Q$
latency-bound messages. From $N = 512$ the volume
dominates: at $P = 16$, `fox_shared` takes 1.3--1.5 times as long, the
$Q slash 2 = 2$ times larger transfer that it cannot overlap. At $P = 4$, where
$Q slash 2 = 1$, the two are within a few percent.

= Conclusions

Fox scales well within a node: four processes reach a speedup of up to 3.9
on the assignment inputs. Across nodes the 1 Gb/s network dominates: several
nodes beat `seq` clearly only from $N = 8192$, and four processes only at 16384. Overlapping the
transfers keeps `fox` ahead of `fox_shared` once the matrices are large.
