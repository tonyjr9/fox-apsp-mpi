#!/bin/sh
#
# Benchmark seq and fox on tests/examples: one CSV row per run, then a summary.
#   ./bench/run_bench.sh
#   NPS="1 4 9 16 25" REPS=5 MPIRUN="mpirun --hostfile hosts" ./bench/run_bench.sh
#
# Times are the compute_ms the binaries print on stderr (I/O excluded).
# Configurations the grid can't handle (N % sqrt(np) != 0) are skipped.

set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
INPUT_DIR="$ROOT/tests/examples/input"
BIN_SEQ=${BIN_SEQ:-$ROOT/seq}
BIN_FOX=${BIN_FOX:-$ROOT/fox}
INPUTS=${INPUTS:-"input300 input600 input900 input1200 input3000"}
NPS=${NPS:-"1 4 9 16"}
REPS=${REPS:-3}
MPIRUN=${MPIRUN:-mpirun}
OUT=${OUT:-$ROOT/bench/results.csv}

for bin in "$BIN_SEQ" "$BIN_FOX"; do
  [ -x "$bin" ] || { echo "error: $bin not found, run make first" >&2; exit 1; }
done

ERR=$(mktemp)
trap 'rm -f "$ERR"' EXIT

# q with q*q == $1, or 0 if $1 isn't a perfect square
grid_q() {
  q=1; while [ $((q * q)) -lt "$1" ]; do q=$((q + 1)); done
  [ $((q * q)) -eq "$1" ] && echo "$q" || echo 0
}

# run <impl> <np> <input> <command...>: one CSV row per repetition
run() {
  impl=$1; np=$2; input=$3; shift 3
  n=$(head -n 1 "$INPUT_DIR/$input")
  rep=1
  while [ "$rep" -le "$REPS" ]; do
    ms=
    if "$@" < "$INPUT_DIR/$input" > /dev/null 2> "$ERR"; then
      ms=$(sed -n 's/^compute_ms=//p' "$ERR")
    fi
    echo "$impl,$np,$input,$n,$rep,${ms:-NA}"
    rep=$((rep + 1))
  done
}

{
  echo "impl,np,input,n,rep,compute_ms"
  for input in $INPUTS; do
    n=$(head -n 1 "$INPUT_DIR/$input")
    run seq 1 "$input" "$BIN_SEQ"

    for np in $NPS; do
      q=$(grid_q "$np")
      if [ "$q" -eq 0 ] || [ $((n % q)) -ne 0 ]; then
        echo "skip: $input (N=$n) with np=$np" >&2
        continue
      fi
      # $MPIRUN is split on purpose so it can carry flags (--hostfile, ...)
      run fox "$np" "$input" $MPIRUN -np "$np" "$BIN_FOX"
    done
  done
} | tee "$OUT"

# Summary on stderr: best of REPS per configuration, speedup against seq
echo >&2
awk -F, 'NR > 1 && $6 != "NA" {
  key = $3 "," $1 "," $2
  if (!(key in best) || $6 < best[key]) best[key] = $6
  if (!(key in seen)) { seen[key] = 1; order[++k] = key }
}
END {
  printf "%-10s %-4s %4s %12s %8s\n", "input", "impl", "np", "best_ms", "speedup"
  for (i = 1; i <= k; i++) {
    split(order[i], f, ",")
    base = best[f[1] ",seq,1"]
    printf "%-10s %-4s %4s %12.1f %8.2f\n", f[1], f[2], f[3], best[order[i]], base / best[order[i]]
  }
}' "$OUT" >&2
