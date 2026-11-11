#!/usr/bin/env bash
# Modular computations for N = 3 behind the proposition on integral solutions (prop:integral) of the paper.
# Requires msolve.  Run from the repository root.
set -euo pipefail
C=certificates/N3
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
for p in 1000000007 1073741827; do
  sed "2s/.*/$p/" "$C/R3_0_1.ms" > "$tmp/R3.ms"
  msolve -v 2 -f "$tmp/R3.ms" -o "$tmp/R3.out" -t 1 > "$tmp/log" 2>&1
  d=$(grep -o "Dimension of quotient: [0-9]*" "$tmp/log" | sort -u)
  echo "slice, p=$p: $d"
  python3 code/jacobian_modp.py "$C/R3_0_1.ms" "$C/r3q" "$p" | tail -1
done
for f in "$C"/modp/ch3_*.ms; do
  for p in 65521 1000000007 1073741827; do
    sed "2s/.*/$p/" "$f" > "$tmp/c.ms"
    msolve -f "$tmp/c.ms" -t 1 > "$tmp/c.out"
    echo "$(basename "$f") p=$p: $(cut -c1-5 "$tmp/c.out" | head -n 1)"
  done
done
