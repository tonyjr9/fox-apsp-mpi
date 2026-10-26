// Result tables computed from the benchmark CSVs (bench/run_bench.sh format:
// impl,np,input,n,rep,compute_ms). Compile with --root at the repository
// root (make report) so the CSVs under /bench/report/ are readable.

// Minimum compute_ms per configuration, keyed "n/impl/np"; NA rows are skipped.
#let load-best(path) = {
  let best = (:)
  for r in csv(path, row-type: dictionary) {
    if r.compute_ms == "NA" { continue }
    let key = r.n + "/" + r.impl + "/" + r.np
    let v = float(r.compute_ms)
    if key not in best or v < best.at(key) { best.insert(key, v) }
  }
  best
}

#let lookup(best, n, impl, p) = best.at(str(n) + "/" + impl + "/" + str(p), default: none)
#let seq-time(best, n) = lookup(best, n, "seq", 1)
#let fox-time(best, n, p) = lookup(best, n, "fox", p)

// v with exactly d decimals ("15.0", not "15").
#let fixed(v, d) = {
  if d == 0 { return str(int(calc.round(v))) }
  let parts = str(calc.round(v, digits: d)).split(".")
  let frac = if parts.len() > 1 { parts.at(1) } else { "" }
  parts.at(0) + "." + frac + "0" * (d - frac.len())
}

#let missing = [---]

// Times: one decimal, whole milliseconds from 100 s up.
#let fmt-ms(v) = if v == none { missing } else { [#fixed(v, if v >= 100000 { 0 } else { 1 })] }

// Speedups: two decimals, three below 0.1 so small values don't read as 0.00.
#let fmt-speedup(v) = if v == none { missing } else { [#fixed(v, if v < 0.1 { 3 } else { 2 })] }

#let ratio(a, b) = if a == none or b == none { none } else { a / b }

// N, seq, then the time for each P.
#let time-rows(best, ns, ps) = ns.map(n => (
  [#n],
  fmt-ms(seq-time(best, n)),
  ..ps.map(p => fmt-ms(fox-time(best, n, p))),
)).flatten()

// N, then S (against seq) and S1 (against fox with P = 1) for each P.
#let speedup-rows(best, ns, ps) = ns.map(n => (
  [#n],
  ..ps.map(p => {
    let t = fox-time(best, n, p)
    (fmt-speedup(ratio(seq-time(best, n), t)), fmt-speedup(ratio(fox-time(best, n, 1), t)))
  }).flatten(),
)).flatten()

// N, seq, P = 1, then time, S and S1 for each P.
#let combined-rows(best, ns, ps) = ns.map(n => (
  [#n],
  fmt-ms(seq-time(best, n)),
  fmt-ms(fox-time(best, n, 1)),
  ..ps.map(p => {
    let t = fox-time(best, n, p)
    (
      fmt-ms(t),
      fmt-speedup(ratio(seq-time(best, n), t)),
      fmt-speedup(ratio(fox-time(best, n, 1), t)),
    )
  }).flatten(),
)).flatten()
