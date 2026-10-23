#!/bin/sh

set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
INPUT_DIR="$ROOT/tests/examples/input"
BIN_SUFFIX=${BIN_SUFFIX:-.bin}
printf 'input,n,time_ms\n'

for input in "$INPUT_DIR"/*.bin; do
  name=$(basename "$input" "$BIN_SUFFIX")
  n=$(head -c 4 "$input" | od -An -tu4 -j 0 | tr -d ' ')

  start=$(date +%s%N)
  "$ROOT/seq" "$input" /dev/null
  end=$(date +%s%N)

  awk -v name="$name" -v n="$n" -v start="$start" -v end="$end" \
    'BEGIN { printf "%s,%s,%.3f\n", name, n, (end - start) / 1000000 }'
done