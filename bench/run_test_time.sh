#!/bin/sh

set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
INPUT_DIR="$ROOT/tests/examples/input"
printf 'input,n,time_ms\n'

for input in "$INPUT_DIR"/input*; do
  name=$(basename "$input")
  n=$(sed -n '1p' "$input")

  start=$(date +%s%N)
  "$ROOT/seq" < "$input" > /dev/null
  end=$(date +%s%N)

  awk -v name="$name" -v n="$n" -v start="$start" -v end="$end" \
    'BEGIN { printf "%s,%s,%.3f\n", name, n, (end - start) / 1000000 }'
done
