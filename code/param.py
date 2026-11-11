import sys, sympy as sp
from gen import system
N=int(sys.argv[1]); eqs,vs,c,W=system(N)
cn=list(c.values())
# scaling: set c12=1
eqs=[sp.expand(e.subs({c[(1,2)]:1})) for e in eqs]
pars=[str(s) for s in cn if s!=c[(1,2)]]
ren={v: sp.Symbol(f"x({i+1})") for i,v in enumerate(vs)}
out=[f"ring r=(0,{','.join(pars)}),(x(1..{len(vs)})),dp;","option(redSB);","ideal I="+",\n".join(str(e.subs(ren)).replace("**","^") for e in eqs)+";","ideal G=std(I);","print(\"dim=\"+string(dim(G))+\" vdim=\"+string(vdim(G)));","quit;"]
open(f"param{N}.sing","w").write("\n".join(out))
