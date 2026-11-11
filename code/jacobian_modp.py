"""Check that every solution of a rational univariate representation is a
nonsingular (isolated, multiplicity one) solution of the system.

Usage: python3 jacobian_modp.py system.ms rur_output [prime]

The representation x_i = -num_i(t)/(c_i den(t)), chi(t) = 0 is read from
msolve output over Q.  Everything is reduced modulo a prime p; the script
checks that chi stays squarefree of the same degree modulo p and that den and
the c_i are units.  For a random p-matrix R it computes D(t) = det(R J(x(t)))
with J the Jacobian matrix of the system, and checks gcd(D, chi) = 1 modulo p.
Then Res(D_Q, chi) is nonzero modulo p, hence nonzero, so at every root of chi
the Jacobian has full rank over Q-bar.
"""
import sys, random, re, ast
from flint import nmod_poly, nmod
import sympy as sp
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from certify import parse

def main(msfile, outfile, p=2305843009213693951):
    vars_, lin, chi, den, nums = parse(outfile)
    L = open(msfile).read().split("\n", 2)
    xs = [sp.Symbol(v) for v in L[0].split(",")]
    polys = [sp.Poly(sp.sympify(t.replace("^", "**")), *xs) for t in L[2].strip().split(",\n")]
    n = len(xs)
    def red(f):  # fmpq_poly with integer coefficients -> nmod_poly
        cs = [int(sp.Rational(str(f[i])) % p) if sp.Rational(str(f[i])).q == 1 else
              int(sp.Rational(str(f[i])).p * pow(sp.Rational(str(f[i])).q, -1, p) % p) for i in range(f.degree() + 1)]
        return nmod_poly(cs, p)
    chip = red(chi)
    assert chip.degree() == chi.degree(), "leading coefficient of chi vanishes mod p"
    assert chip.gcd(chip.derivative()).degree() == 0, "chi not squarefree mod p"
    denp = red(den) % chip
    g = denp.gcd(chip); assert g.degree() == 0, "den not a unit mod (p, chi)"
    # inverse of den modulo chi by extended Euclid
    r0, r1 = chip, denp; s0, s1 = nmod_poly([0], p), nmod_poly([1], p)
    while r1 != 0:
        q = r0 // r1; r0, r1 = r1, r0 - q * r1; s0, s1 = s1, s0 - q * s1
    dinv = (s0 * nmod_poly([int(pow(int(r0[0]), -1, p))], p)) % chip
    vals = {}
    for (num, cc), v in zip(nums, vars_[:-1]):
        assert cc % p != 0
        vals[v] = (red(num) * dinv * nmod_poly([(-pow(cc, -1, p)) % p], p)) % chip
    vals[vars_[-1]] = nmod_poly([0, 1], p)
    powcache = {}
    def pw(v, e):
        if (v, e) not in powcache: powcache[(v, e)] = _pw(vals[v], e, chip)
        return powcache[(v, e)]
    def ev(P):
        acc = nmod_poly([0], p)
        for mon, co in P.terms():
            c = int(sp.Rational(co).p * pow(sp.Rational(co).q, -1, p) % p)
            term = nmod_poly([c], p)
            for v, e in zip(xs, mon):
                if e: term = (term * pw(str(v), e)) % chip
            acc = acc + term
        return acc % chip
    J = [[ev(sp.Poly(P.diff(x), *xs)) for x in xs] for P in polys]
    for attempt in range(3):
        rng = random.Random(attempt)
        R = [[rng.randrange(p) for _ in polys] for _ in range(n)]
        A = [[sum((J[k][j] * R[i][k] for k in range(len(polys))), nmod_poly([0], p)) % chip for j in range(n)] for i in range(n)]
        D = bareiss(A, p)
        gd = D.gcd(chip).degree()
        print(f"attempt {attempt}: deg chi {chip.degree()}, deg gcd(det, chi) = {gd}")
        if gd == 0:
            print(f"all {chip.degree()} solutions are nonsingular"); return True
    return False

def _pw(a, e, m):
    r = nmod_poly([1], a.modulus()); b = a
    while e:
        if e & 1: r = (r * b) % m
        b = (b * b) % m; e >>= 1
    return r

def bareiss(A, p):
    n = len(A); M = [row[:] for row in A]; prev = nmod_poly([1], p); sign = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            for i in range(k + 1, n):
                if M[i][k] != 0:
                    M[k], M[i] = M[i], M[k]; sign = -sign; break
            else:
                return nmod_poly([0], p)
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[k][k] * M[i][j] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return M[n - 1][n - 1] * sign

if __name__ == "__main__":
    p = int(sys.argv[3]) if len(sys.argv) > 3 else 2305843009213693951
    ok = main(sys.argv[1], sys.argv[2], p)
    sys.exit(0 if ok else 1)
