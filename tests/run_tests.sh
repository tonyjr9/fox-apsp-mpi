#!/bin/sh

set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
INPUT_DIR="$ROOT/tests/examples/input"
OUTPUT_DIR="$ROOT/tests/examples/output"
BIN="$ROOT/seq"
TIME_LIMIT=${TIME_LIMIT:-10}   # seconds per test, override with TIME_LIMIT=30 ./test.sh

TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

# Colors only when writing to a terminal
if [ -t 1 ]; then
  GREEN=$(printf '\033[32m'); RED=$(printf '\033[31m')
  YELLOW=$(printf '\033[33m'); DIM=$(printf '\033[2m'); RESET=$(printf '\033[0m')
else
  GREEN=; RED=; YELLOW=; DIM=; RESET=
fi

# Milliseconds since epoch (GNU date); empty if unsupported (e.g. macOS)
now_ms() {
  t=$(date +%s%N 2>/dev/null) || return 0
  case $t in *N) ;; *) echo "${t%??????}" ;; esac
}

# Use `timeout` when available so a hang fails the test instead of the run
if command -v timeout >/dev/null 2>&1; then
  RUN="timeout $TIME_LIMIT"
else
  RUN=
fi

[ -x "$BIN" ] || { echo "${RED}error:${RESET} $BIN not found, run make first" >&2; exit 1; }

pass=0; fail=0; skip=0

for input in "$INPUT_DIR"/input*; do
  name=$(basename "$input")
  expected="$OUTPUT_DIR/output${name#input}"

  # input5 checks MPI compatibility only: N=5 is not divisible by
  # the Q x Q process grids required by the assignment.
  if [ "$name" = input5 ]; then
    printf '%sSKIP%s  %s %s(MPI only)%s\n' "$YELLOW" "$RESET" "$name" "$DIM" "$RESET"
    skip=$((skip + 1))
    continue
  fi

  actual="$TMP_DIR/$name"
  start=$(now_ms)
  if $RUN "$BIN" < "$input" > "$actual" 2> "$TMP_DIR/$name.err"; then status=0; else status=$?; fi
  end=$(now_ms)

  elapsed=
  [ -n "$start" ] && [ -n "$end" ] && elapsed=" ${DIM}($((end - start)) ms)${RESET}"

  if [ "$status" -eq 124 ]; then
    printf '%sTIME%s  %s %s(> %ss)%s\n' "$RED" "$RESET" "$name" "$DIM" "$TIME_LIMIT" "$RESET"
    fail=$((fail + 1))
  elif [ "$status" -ne 0 ]; then
    printf '%sCRASH%s %s%s exit %s\n' "$RED" "$RESET" "$name" "$elapsed" "$status"
    sed 's/^/        /' "$TMP_DIR/$name.err" | head -n 10
    fail=$((fail + 1))
  elif diff -u "$expected" "$actual" > "$TMP_DIR/$name.diff"; then
    printf '%sPASS%s  %s%s\n' "$GREEN" "$RESET" "$name" "$elapsed"
    pass=$((pass + 1))
  else
    printf '%sFAIL%s  %s%s\n' "$RED" "$RESET" "$name" "$elapsed"
    sed 's/^/        /' "$TMP_DIR/$name.diff" | head -n 20
    fail=$((fail + 1))
  fi
done

echo
printf '%s%d passed%s, %s%d failed%s, %d skipped\n' \
  "$GREEN" "$pass" "$RESET" "$RED" "$fail" "$RESET" "$skip"

[ "$fail" -eq 0 ]
