#!/bin/sh
#
# Compare the fox distribution variants on power-of-two matrices.
#   ./bench/run_power2.sh
#   NPS="1 4 9 16" REPS=3 VARIANTS="fox fox_shared" ./bench/run_power2.sh
#
# Variants:
#   fox          default row distribution (baseline)
#   fox-colmajor mirrored scatter: rank r holds block (r%q, r/q)
#   fox_shared   no communication: every process holds the full matrix
#
# Inputs are rand0_input<N>.bin in the repo root, generated with a fixed
# seed (matgen, seed 42) when missing. Times are the compute_ms the
# binaries print on stderr (I/O excluded). Configurations the grid can't
# handle (N % sqrt(np) != 0) are skipped.

set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"

INPUTS=${INPUTS:-"64 128 256 512 1024"}
NPS=${NPS:-"1 4 9 16"}
REPS=${REPS:-3}
MPIRUN=${MPIRUN:-mpirun --oversubscribe}
VARIANTS=${VARIANTS:-"fox fox-colmajor fox_shared"}
OUT=${OUT:-$ROOT/bench/results_power2.csv}

# q with q*q == $1, or 0 if $1 isn't a perfect square
grid_q() {
  q=1; while [ $((q * q)) -lt "$1" ]; do q=$((q + 1)); done
  [ $((q * q)) -eq "$1" ] && echo "$q" || echo 0
}

# Generate the inputs if missing (fixed seed: reproducible across machines)
make -s build/tests/matgen
for n in $INPUTS; do
  if [ ! -f "rand0_input$n.bin" ]; then
    echo "generating rand0_input$n.bin (seed 42)" >&2
    ./build/tests/matgen "$n" 1 42
  fi
done

for v in $VARIANTS; do
  [ -x "$v" ] || { echo "error: $v not found, run make first" >&2; exit 1; }
done

ERR=$(mktemp)
trap 'rm -f "$ERR"' EXIT

# run <variant> <np> <n> <command...>: one CSV row per repetition
run() {
  variant=$1; np=$2; n=$3; shift 3
  rep=1
  while [ "$rep" -le "$REPS" ]; do
    ms=
    if "$@" "rand0_input$n.bin" /dev/null 2> "$ERR"; then
      ms=$(sed -n 's/^compute_ms=//p' "$ERR")
    fi
    echo "$n,$variant,$np,$rep,${ms:-NA}"
    rep=$((rep + 1))
  done
}

{
  echo "n,variant,np,rep,compute_ms"
  for n in $INPUTS; do
    for np in $NPS; do
      q=$(grid_q "$np")
      [ "$q" -eq 0 ] || [ $((n % q)) -ne 0 ] && { echo "skip: N=$n with np=$np" >&2; continue; }
      for v in $VARIANTS; do
        if [ "$np" -eq 1 ]; then
          run "$v" 1 "$n" "./$v"
        else
          run "$v" "$np" "$n" "$MPIRUN" -np "$np" "./$v"
        fi
      done
    done
  done
} | tee "$OUT"

# Summary on stderr: best of REPS per (N, variant, np), and how the variant
# compares to the default fox at the same N and np (ratio > 1: fox is slower,
# i.e. that's the communication the no-comm variant avoids)
echo >&2
awk -F, 'NR > 1 && $5 != "NA" {
  key = $1 "," $2 "," $3
  if (!(key in best) || $5 < best[key]) best[key] = $5
  if (!(key in seen)) { seen[key] = 1; order[++k] = key }
}
END {
  printf "%6s %-14s %4s %12s %10s\n", "N", "variant", "np", "best_ms", "fox/this"
  for (i = 1; i <= k; i++) {
    split(order[i], f, ",")
    ratio = "-"
    foxms = (f[1] ",fox," f[3]) in best ? best[f[1] ",fox," f[3]] : 0
    if (foxms > 0) ratio = foxms / best[order[i]]
    printf "%6s %-14s %4s %12.1f %10s\n", f[1], f[2], f[3], best[order[i]], ratio
  }
}' "$OUT" >&2
