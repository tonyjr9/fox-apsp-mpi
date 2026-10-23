#!/bin/sh
#
# Convert the text matrices in tests/examples to the binary format used by
# seq/fox, next to the originals (input6 -> input6.bin).
#   ./tests/examples/to_bin.sh
# Input files carry a dimension header; output files do not, so for outputs
# the dimension is borrowed from the matching input file.
# Files that are already converted are skipped unless the text changed.
# Files that are not matrices (output5 is an error message) are skipped
# with a notice.

set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
CONV=${CONV:-$ROOT/build/tests/binconv}

[ -x "$CONV" ] || { echo "error: $CONV not found, run 'make binconv' first" >&2; exit 1; }

for dir in input output; do
  for f in "$ROOT/tests/examples/$dir"/*; do
    name=$(basename "$f")
    case $name in *.bin) continue ;; esac

    # Skip if the binary exists, is non-empty, and is not older than the text.
    if [ -s "$f.bin" ] && [ ! "$f" -nt "$f.bin" ]; then
      continue
    fi

    if "$CONV" < "$f" > "$f.bin" 2> /dev/null; then
      echo "converted $name -> $name.bin"
      continue
    fi
    rm -f "$f.bin"

    # Output files have no dimension header: borrow it from the matching input.
    if [ "$dir" = "output" ]; then
      src="$ROOT/tests/examples/input/input${name#output}"
      if [ -f "$src" ] && { sed -n '1p' "$src"; cat "$f"; } | "$CONV" > "$f.bin" 2> /dev/null; then
        echo "converted $name -> $name.bin (dimension from input${name#output})"
        continue
      fi
      rm -f "$f.bin"
    fi

    echo "skipped $name (not a matrix)"
  done
done