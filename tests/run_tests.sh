#!/bin/sh

set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
INPUT_DIR="$ROOT/tests/examples/input"
OUTPUT_DIR="$ROOT/tests/examples/output"
TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

for input in "$INPUT_DIR"/input*; do
  name=$(basename "$input")
  expected="$OUTPUT_DIR/output${name#input}"

  # input5 checks MPI compatibility only: N=5 is not divisible by
  # the Q x Q process grids required by the assignment.
  if [ "$name" = input5 ]; then
    continue
  fi

  actual="$TMP_DIR/$name"
  "$ROOT/seq" < "$input" > "$actual"

  if ! diff -u "$expected" "$actual"; then
    echo "FAIL: $name" >&2
    exit 1
  fi

  echo "PASS: $name"
done
