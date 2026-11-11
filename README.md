# Zero angular momentum of D2-D2-D2-D6 microstates

Code and certificates for the paper *Black Hole Microstates of Two D6-Branes Carry No Angular Momentum* by Deep Bhattacharjee.

- `paper/` LaTeX source and PDF.
- `code/` builds the superpotential and the gauge-fixed equations (`gen.py`, `reduce.py`, `strata.py`, `strata2.py`), checks the F-term list (`verify_fterms.py`), checks rational univariate representations exactly (`certify.py`, and `certify_modp.py` for N = 4), checks that every certified solution is nonsingular (`jacobian_modp.py`), runs the modular N = 3 computations (`modp_n3.sh`), prints the predicted counts (`expected_counts.py`).
- `certificates/N1` to `N4` hold the input systems (`.ms`, `.sing`), the msolve outputs and the Singular logs. `run_checks.sh` reruns the checks.
- `figures/` draws the figures.

Requirements: Python 3 with sympy, numpy, matplotlib and python-flint; Singular 4; msolve 0.7 or later.

DOI (all versions): see the Zenodo record linked from the paper.
