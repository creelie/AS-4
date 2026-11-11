import sys, random, sympy as sp
from gen import system
N=int(sys.argv[1]); char=sys.argv[2]; seed=int(sys.argv[3])
eqs,vs,c,W=system(N)
random.seed(seed)
cv={s: sp.Rational(random.choice([-1,1])*random.randint(1,30), random.randint(1,30)) for s in c.values()}
eqs=[sp.expand(e.subs(cv)) for e in eqs]
names=[str(v) for v in vs]
ren={v: sp.Symbol(f"x({i+1})") for i,v in enumerate(vs)}
out=["LIB \"primdec.lib\";",f"ring r={char},(x(1..{len(vs)})),dp;", "option(redSB);","ideal I="]
out.append(",\n".join(str(e.subs(ren)).replace("**","^") for e in eqs)+";")
out+=["ideal G=std(I);","int d=dim(G); int v=vdim(G);", "print(\"dim=\"+string(d)+\" vdim=\"+string(v));",
      "ideal R=radical(I); ideal GR=std(R); print(\"radical vdim=\"+string(vdim(GR)));","quit;"]
open(f"N{N}_{char}_{seed}.sing","w").write("\n".join(out))
print(cv)
