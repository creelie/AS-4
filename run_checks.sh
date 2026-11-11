#!/usr/bin/env bash
# Reruns the exact checks of the paper.  Run from the repository root.
set -euo pipefail
root=$(cd "$(dirname "$0")" && pwd)
cd "$root/code"
python3 verify_fterms.py
python3 expected_counts.py
C="$root/certificates"
python3 certify.py "$C/N1/N1_0_1.ms" "$C/N1/o1"
python3 certify.py "$C/N2/N2_0_1.ms" "$C/N2/o2q"
python3 certify.py "$C/N3/R3_0_1.ms" "$C/N3/r3q"
python3 certify_modp.py "$C/N4/M4.ms" "$C/N4/m4"
python3 jacobian_modp.py "$C/N1/N1_0_1.ms" "$C/N1/o1"
python3 jacobian_modp.py "$C/N2/N2_0_1.ms" "$C/N2/o2q"
python3 jacobian_modp.py "$C/N3/R3_0_1.ms" "$C/N3/r3q"
if command -v Singular >/dev/null; then
  for f in "$C/N1/param1.sing" "$C/N1/N1_0_1.sing" "$C/N2/R2_0_1.sing" "$C/N2/ex2_par.sing"; do
    echo "== $f"; Singular -q "$f"
  done
  echo "The exact N = 3 Gröbner basis (R3_std.sing) and the twelve charts ch3_*.sing over Q did not finish in our runs."
fi
if command -v msolve >/dev/null; then
  (cd "$root" && bash code/modp_n3.sh)
fi
