import sys, random, sympy as sp
from gen import system
from reduce import linelim
N=int(sys.argv[1]); char=int(sys.argv[2]); seed=int(sys.argv[3])
eqs,vs,c,W=system(N)
random.seed(seed)
cv={s: sp.Rational(random.choice([-1,1])*random.randint(1,30), random.randint(1,30)) for s in c.values()}
eqs=[sp.expand(e.subs(cv)) for e in eqs]
eqs,vs,subs=linelim(eqs,vs)
print(len(eqs),len(vs),file=sys.stderr)
eqs=[sp.expand(e*sp.ilcm(1,1,*[sp.fraction(t)[1] for t in sp.Poly(e,*vs).coeffs()])) for e in eqs]
ren={v: sp.Symbol(f"x{i+1}") for i,v in enumerate(vs)}
lines=[",".join(f"x{i+1}" for i in range(len(vs))), str(char)]
lines.append(",\n".join(str(e.subs(ren)).replace("**","^") for e in eqs))
open(f"R{N}_{char}_{seed}.ms","w").write("\n".join(lines)+"\n")
