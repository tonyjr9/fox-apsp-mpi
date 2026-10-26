#!/bin/sh
# Run fox on every rand0_input*.bin in $IN_ZIP, one machine grid per size.
#
# The grid uses the largest q with q*q processes (q*q must divide N and fit
# the live slots), so bigger matrices automatically get more machines.
# Must be run from the machine that is first in the hostfile (the one you
# are logged in to); probe it with `export HWLOC_COMPONENTS=-gl` first.

set -eu

SLOTS=${SLOTS:-8}                    # cores per lab machine
HOSTFILE=${HOSTFILE:-hostfile.txt}   # full lab hostfile, all machines
IN_ZIP=${IN_ZIP:-rand0_input.zip}
OUT_ZIP=${OUT_ZIP:-rand0_output.zip}
WORK=${WORK:-rand}

probe() {
  echo "$(hostname)"
  for h in $(awk '{print $1}' "$HOSTFILE"); do
    [ "$h" = "$(hostname | tr 'A-Z' 'a-z')" ] && continue
    if ssh -o BatchMode=yes -o ConnectTimeout=3 "$h" true 2> /dev/null; then
      echo "$h"
    fi
  done
}

LIVE=$(probe)
N_MACH=$(printf '%s\n' "$LIVE" | wc -l)
SLOTS_TOTAL=$(( N_MACH * SLOTS ))
{
  for h in $LIVE; do
    printf '%s slots=%s\n' "$h" "$SLOTS"
  done
} > "$WORK.hostfile"
echo "live machines: $N_MACH (slots $SLOTS_TOTAL)"

mkdir -p "$WORK"
unzip -o "$IN_ZIP" -d "$WORK"

fail=0
for f in "$WORK"/rand0_input*.bin; do
  name=$(basename "$f")
  N=${name#rand0_input}
  N=${N%.bin}

  c=$(awk -v s="$SLOTS_TOTAL" 'BEGIN { c = int(sqrt(s)); while (c * c > s) c--; print c }')
  [ "$c" -gt "$N" ] && c=$N
  while [ "$c" -gt 1 ] && [ $(( N % c )) -ne 0 ]; do
    c=$(( c - 1 ))
  done
  p=$(( c * c ))

  echo "=== $name: N=$N grid=${c}x${c} nproc=$p ==="
  if mpirun -n "$p" --hostfile "$WORK.hostfile" \
      fox_bin "$f" "$WORK/rand0_output$N.bin"; then
    :
  else
    echo "FAIL $name" >&2
    fail=1
  fi
done

[ "$fail" -eq 0 ] || exit 1

zip -q "$OUT_ZIP" "$WORK"/rand0_output*.bin
echo DONE