"""Exact upper bound for the number of solutions over Q from one prime p.
If the top-degree forms of the integer generators have only the trivial common
zero modulo p, then dim_Q Q[x]/I <= dim_Fp Fp[x]/(I mod p) (closed fibre of a
finite scheme over Z_(p)).  Writes a Singular script that checks both."""
import sys, sympy as sp
src, p, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
txt = open(src).read().split("\n", 2)
names = txt[0].strip().split(",")
X = sp.symbols(names)
body = txt[2].replace("\n", "").rstrip(",").strip()
eqs = [e for e in body.split(",") if e.strip()]
loc = {n: s for n, s in zip(names, X)}
tops, full, dropped = [], [], 0
for e in eqs:
    P = sp.Poly(sp.sympify(e, locals=loc), *X)
    P = P.clear_denoms()[1]
    if P.is_zero: continue
    d = P.total_degree()
    T = [(m, c) for m, c in P.terms() if sum(m) == d]
    Tp = [(m, int(c) % p) for m, c in T if int(c) % p]
    if not Tp: dropped += 1
    else: tops.append(Tp)
    full.append([(m, int(c) % p) for m, c in P.terms() if int(c) % p])
def fmt(terms):
    s = []
    for m, c in terms:
        mon = "*".join(f"x({i+1})^{k}" if k > 1 else f"x({i+1})" for i, k in enumerate(m) if k)
        s.append(f"{c}*{mon}" if mon else f"{c}")
    return "+".join(s) if s else "0"
with open(out, "w") as f:
    f.write(f"ring r={p},(x(1..{len(X)})),dp;\noption(redSB);\n")
    f.write("ideal T=" + ",\n".join(fmt(t) for t in tops) + ";\n")
    f.write("ideal I=" + ",\n".join(fmt(t) for t in full if t) + ";\n")
    f.write('ideal GT=std(T); print("top forms: dim="+string(dim(GT)));\n')
    f.write('ideal GI=std(I); print("mod p: dim="+string(dim(GI))+" vdim="+string(vdim(GI)));\nquit;\n')
print(src, "eqs", len(eqs), "vars", len(X), "top forms dropped", dropped)
