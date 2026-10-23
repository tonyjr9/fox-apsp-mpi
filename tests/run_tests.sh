#!/bin/sh
#
# Regression on tests/examples.
#   ./tests/run_tests.sh                    # sequential: ./seq
#   BIN=./fox NP=4 ./tests/run_tests.sh     # MPI: $MPIRUN -np 4 ./fox
#
# With NP set, inputs whose N is not divisible by sqrt(NP) (or any input, when
# NP is not a perfect square) must be rejected: non-zero exit and exactly
# "$INVALID_MSG" on stdout (the content of output5).
# An input whose reference output is that message (input5) is skipped when the
# configuration is valid, since there is no matrix to compare against.

set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
EXAMPLES_DIR=${EXAMPLES_DIR:-$ROOT/tests/examples}  # e.g. EXAMPLES_DIR=tests/large
INPUT_DIR="$EXAMPLES_DIR/input"
OUTPUT_DIR="$EXAMPLES_DIR/output"
BIN=${BIN:-$ROOT/seq}
NP=${NP:-}
MPIRUN=${MPIRUN:-mpirun}
INVALID_MSG='ERROR: Invalid configuration!'
TIME_LIMIT=${TIME_LIMIT:-120}  # seconds per test, override with TIME_LIMIT=300 ./tests/run_tests.sh

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
  TIMEOUT="timeout $TIME_LIMIT"
else
  TIMEOUT=
fi

[ -x "$BIN" ] || { echo "${RED}error:${RESET} $BIN not found, run make first" >&2; exit 1; }

LAUNCH=
Q=0
if [ -n "$NP" ]; then
  LAUNCH="$MPIRUN -np $NP"
  Q=1; while [ $((Q * Q)) -lt "$NP" ]; do Q=$((Q + 1)); done
  [ $((Q * Q)) -eq "$NP" ] || Q=0   # not a perfect square: every input is invalid
fi

echo "${DIM}== $(basename "$BIN")${NP:+ np=$NP}${RESET}"

pass=0; fail=0; skip=0

for input in "$INPUT_DIR"/input*; do
  name=$(basename "$input")
  expected="$OUTPUT_DIR/output${name#input}"
  n=$(head -n 1 "$input")

  valid=1
  [ -n "$NP" ] && { [ "$Q" -eq 0 ] || [ $((n % Q)) -ne 0 ]; } && valid=0

  if [ "$valid" -eq 1 ] && printf '%s' "$INVALID_MSG" | cmp -s - "$expected"; then
    printf '%sSKIP%s  %s %s(reference is for an invalid grid)%s\n' "$YELLOW" "$RESET" "$name" "$DIM" "$RESET"
    skip=$((skip + 1))
    continue
  fi

  actual="$TMP_DIR/$name"
  start=$(now_ms)
  if $TIMEOUT $LAUNCH "$BIN" < "$input" > "$actual" 2> "$TMP_DIR/$name.err"; then status=0; else status=$?; fi
  end=$(now_ms)

  elapsed=
  [ -n "$start" ] && [ -n "$end" ] && elapsed=" ${DIM}($((end - start)) ms)${RESET}"

  if [ "$status" -eq 124 ]; then
    printf '%sTIME%s  %s %s(> %ss)%s\n' "$RED" "$RESET" "$name" "$DIM" "$TIME_LIMIT" "$RESET"
    fail=$((fail + 1))
  elif [ "$valid" -eq 0 ]; then
    if [ "$status" -ne 0 ] && printf '%s' "$INVALID_MSG" | cmp -s - "$actual"; then
      printf '%sPASS%s  %s %s(rejected)%s\n' "$GREEN" "$RESET" "$name" "$DIM" "$RESET"
      pass=$((pass + 1))
    else
      printf '%sFAIL%s  %s should have been rejected (exit %s)\n' "$RED" "$RESET" "$name" "$status"
      fail=$((fail + 1))
    fi
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
