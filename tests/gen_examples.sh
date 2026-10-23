#!/bin/sh
#
# Generate random inputs and their reference outputs.
#   ./tests/gen_examples.sh                       # default sizes into tests/large
#   ./tests/gen_examples.sh 5040 10080            # chosen sizes
#   OUT_DIR=tests/examples ./tests/gen_examples.sh 840
#   REF="mpirun -np 16 ./fox" ./tests/gen_examples.sh 10080   # faster reference
#
# Files are written as $OUT_DIR/input/inputN and $OUT_DIR/output/outputN, so
# the result can be checked with EXAMPLES_DIR=$OUT_DIR ./tests/run_tests.sh.
# The seed defaults to N, so the same N always gives the same matrix.
#
# Pick N divisible by every Q you want to run with (P = Q*Q):
#   2520  -> Q = 1..10, 12, 14, 15        (up to 225 processes)
#   5040  -> Q = 1..10, 12, 14, 15, 16    (up to 256 processes)
#   10080 -> as 5040, plus 18, 20, 24 ...

set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
OUT_DIR=${OUT_DIR:-$ROOT/tests/large}
CC=${CC:-gcc}
REF=${REF:-$ROOT/seq}
SIZES=${*:-5040}

mkdir -p "$OUT_DIR/input" "$OUT_DIR/output"

# Build the generator (and ./seq, if it is the reference) when missing or stale
if [ ! -x "$ROOT/matgen" ] || [ "$ROOT/tests/matgen.c" -nt "$ROOT/matgen" ]; then
  $CC -O2 -Wall -Wextra -std=c11 -o "$ROOT/matgen" "$ROOT/tests/matgen.c"
fi
[ "$REF" = "$ROOT/seq" ] && make -C "$ROOT" --no-print-directory seq >/dev/null

for n in $SIZES; do
  input="$OUT_DIR/input/input$n"
  output="$OUT_DIR/output/output$n"

  if [ -s "$input" ]; then
    echo "input$n: exists, keeping it"
  else
    echo "input$n: generating (seed $n)"
    "$ROOT/matgen" "$n" 1 "$n" > "$input.tmp" && mv "$input.tmp" "$input"
  fi

  if [ -s "$output" ] && [ "$output" -nt "$input" ]; then
    echo "output$n: exists, keeping it"
  else
    echo "output$n: computing with $REF"
    $REF < "$input" > "$output.tmp" && mv "$output.tmp" "$output"
  fi
done
