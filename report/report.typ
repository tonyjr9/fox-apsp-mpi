#import "template.typ": report, todo, note-box
#import "bench.typ": load-best, fox-time, seq-time, fmt-ms, fmt-speedup, ratio, time-rows, speedup-rows, combined-rows

#let lab = load-best("/bench/report/lab_examples.csv")
#let mac = load-best("/bench/report/mac_examples.csv")
#let p9-8plus1 = load-best("/bench/report/lab_p9_8plus1.csv")
#let p9-3x3 = load-best("/bench/report/lab_p9_3x3.csv")
#let p9-5plus4 = load-best("/bench/report/lab_p9_5plus4.csv")
#let p9-2x4plus1 = load-best("/bench/report/lab_p9_2x4plus1.csv")
#let sizes = (300, 600, 900, 1200, 3000)

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
    nodes, the 1 Gb/s network costs more than the extra processes gain, except
    on the largest random matrices.
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
(`bench/report/lab_examples.csv`).
#if fox-time(lab, 3000, 25) == none [#todo[$P = 25$ at $N = 3000$ is still running.]]

#figure(
  image("lab_examples_speedup.png", width: 68%),
  caption: [Time (log scale) and speedup on the assignment inputs on LabCC.
    Speedups below 0.5 are left unlabelled; they are in @lab-speedup.],
) <lab-examples-plot>

#figure(
  caption: [Assignment inputs on LabCC: minimum of 3 runs in milliseconds.],
  table(
    columns: 7,
    align: (left, ..range(6).map(_ => right)),
    table.header([$N$], [seq], [$P = 1$], [$P = 4$], [$P = 9$], [$P = 16$], [$P = 25$]),
    ..time-rows(lab, sizes, (1, 4, 9, 16, 25)),
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
    ..speedup-rows(lab, sizes, (4, 9, 16, 25)),
  ),
) <lab-speedup>

With $P = 1$, `fox` is within 2% of `seq` from $N = 600$ (9% slower on the
15 ms run at $N = 300$). On one node, four processes reach 3.1--3.9 from
$N = 600$, with an efficiency of 0.96 at $N = 3000$. With 9, 16 and 25
processes, on two to four nodes, this run is slower than `seq` for every
input, and $P = 25$ varies up to 3.7 times between runs. A rerun of $P = 9$
(@placement) shows these multi-node times were inflated by other load on the
shared nodes; we read them as upper bounds.

== Placement of nine processes <placement>

Nine processes do not fit one 8-slot node. We reran $P = 9$ with the default
placement (8 + 1) and three others (`bench/report/lab_p9_*.csv`, minimum of
3 runs; 2 for 2 + 2 + 2 + 2 + 1 at $N = 3000$).

#figure(
  caption: [$P = 9$ on LabCC by placement: time in ms and speedup $S$ against
    `seq`. "First run" is $P = 9$ from @lab-times.],
  table(
    columns: 12,
    align: (left, ..range(11).map(_ => right)),
    stroke: (_, y) => (
      top: if y == 0 { 0.8pt } else if y == 2 { 0.4pt } else { 0pt },
    ),
    table.header(
      [], [], [], [], table.cell(colspan: 2, align: center)[8 + 1],
      table.cell(colspan: 2, align: center)[5 + 4],
      table.cell(colspan: 2, align: center)[3 + 3 + 3],
      table.cell(colspan: 2, align: center)[2+2+2+2+1],
      [$N$], [seq], [$P = 4$], [first run], [time], [$S$], [time], [$S$], [time], [$S$], [time], [$S$],
    ),
    ..sizes.map(n => {
      let s = seq-time(lab, n)
      let a = fox-time(p9-8plus1, n, 9)
      let b = fox-time(p9-5plus4, n, 9)
      let c = fox-time(p9-3x3, n, 9)
      let d = fox-time(p9-2x4plus1, n, 9)
      (
        [#n], fmt-ms(s), fmt-ms(fox-time(lab, n, 4)), fmt-ms(fox-time(lab, n, 9)),
        fmt-ms(a), fmt-speedup(ratio(s, a)), fmt-ms(b), fmt-speedup(ratio(s, b)),
        fmt-ms(c), fmt-speedup(ratio(s, c)), fmt-ms(d), fmt-speedup(ratio(s, d)),
      )
    }).flatten(),
  ),
) <placement-table>

With the same 8 + 1 placement, the rerun is 4--8 times faster than the first
run and beats `seq` from $N = 600$ (3.29 at $N = 3000$), but never beats four
processes on one node. Ranks fill the grid row by row, so with 8 + 1 only the
row and column of the remote process cross the 1 Gb/s link; with 3 + 3 + 3
every column shift of $B$ does (8--47% slower); with 5 + 4 every column shift
and one row broadcast do (21--90% slower); with 2 + 2 + 2 + 2 + 1 every row
and column spans two nodes, so all transfers cross it (8--17 times slower).

== Assignment inputs on an Apple M1 Max

For comparison we ran the same inputs on a MacBook with an Apple M1 Max
(10 cores: 8 performance, 2 efficiency) and 32 GB of memory
(`bench/report/mac_examples.csv`). With only 10 cores it can run $P = 1$, 4
and 9, but not 16 or 25.

#grid(
  columns: (46%, 1fr),
  column-gutter: 10pt,
  align: horizon,
  [#figure(
    image("mac_examples_speedup.png", width: 100%),
    caption: [Time (log scale) and speedup on the M1 Max.],
  ) <mac-plot>],
  [#figure(
    caption: [M1 Max: minimum of 3 runs in ms, $S$ against `seq`, $S_1$
      against $P = 1$.],
    {
      set text(8pt)
      table(
        columns: 9,
        inset: (x: 0.35em, y: 0.4em),
        align: (left, ..range(8).map(_ => right)),
        stroke: (_, y) => (
          top: if y == 0 { 0.8pt } else if y == 2 { 0.4pt } else { 0pt },
        ),
        table.header(
          [], [], [], table.cell(colspan: 3, align: center)[$P = 4$],
          table.cell(colspan: 3, align: center)[$P = 9$],
          [$N$], [seq], [$P = 1$], [time], [$S$], [$S_1$], [time], [$S$], [$S_1$],
        ),
        ..combined-rows(mac, sizes, (4, 9)),
      )
    },
  ) <mac-table>],
)

On one machine $P = 9$ avoids the network: it trails $P = 4$ up to
$N = 1200$ (three steps per product, one process on an efficiency core) but
reaches 6.30 at $N = 3000$.

== Random matrices on LabCC

Here $P = 4$ runs on one node, $P = 16$ on two and $P = 64$ on eight; these
are the process counts whose $Q$ divides a power of two
(`bench/report/lab_rand0.csv`).

#figure(
  image("lab_rand0_speedup.png", width: 68%),
  caption: [Time (log scale) and speedup on the random matrices on LabCC.
    Labels close to one already shown are omitted.],
) <lab-plot>

With four processes, efficiency stays between 0.89 and 1.11 from $N = 1024$.
$P = 16$ and $P = 64$ take about the same time (within 15% from $N = 128$ to
8192), so the network, not the process count, sets it. They break even with
`seq` at $N = 4096$ and reach 3.88 and 5.00 at $N = 16384$ ($P = 4$: 4.01).

== Fox against the shared copy


#figure(
  image("power2_fox_vs_shared.png", width: 72%),
  caption: [`fox` against `fox_shared`, median of 3 to 6 runs. Bottom:
    `fox_shared` time over `fox` time; shaded $P$ oversubscribes the 16
    hardware threads.],
) <shared-plot>

On an AMD Ryzen 7 7840HS laptop (`bench/report/results_power2.csv`, random
matrices with $N = 2^6, ..., 2^10$), `fox_shared` is mostly faster, by up to
28%, for $N <= 256$, where one all-gather beats $2 Q$ latency-bound messages. From
$N = 512$, at $P = 16$, it takes 1.3--1.5 times as long: it moves
$Q slash 2 = 2$ times more data, without overlap. At $P = 4$ ($Q slash 2 = 1$)
the two are within a few percent.

= Conclusions

Fox scales well within a node: four processes reach a speedup of up to 3.9.
Across nodes the 1 Gb/s link costs more than the extra processes gain, and
only 64 processes beat four, at $N = 16384$. Overlapping the transfers keeps
`fox` ahead of `fox_shared` on large matrices.
