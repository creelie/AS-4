import sys, random, sympy as sp
from gen import system
from reduce import linelim
def build(N, rowsub, extra_vars=(), seed=1, name="s"):
    eqs,vs,c,W=system(N,gauge=False)
    sub={}
    for (k,m) in [(1,1),(1,2),(1,3),(2,1),(2,2),(3,3)]: sub[sp.Symbol(f"p{k}{m}_00")]=0
    sub[sp.Symbol("z12_00")]=1; sub[sp.Symbol("z23_00")]=1
    sub.update(rowsub)
    random.seed(seed)
    cv={s: sp.Rational(random.choice([-1,1])*random.randint(1,30), random.randint(1,30)) for s in c.values()}
    eqs=[sp.expand(e.subs(sub).subs(cv)) for e in eqs]; eqs=[e for e in eqs if e!=0]
    vs=[v for v in vs if v not in sub]+list(extra_vars)
    vs=[v for v in vs if any(v in e.free_symbols for e in eqs)]
    eqs,vs,_=linelim(eqs,vs)
    eqs=[sp.expand(e*sp.ilcm(1,1,*[sp.fraction(t)[1] for t in sp.Poly(e,*vs).coeffs()])) for e in eqs]
    ren={v: sp.Symbol(f"x({i+1})") for i,v in enumerate(vs)}
    out=[f"ring r=0,(x(1..{len(vs)})),dp;","ideal I="+",\n".join(str(e.subs(ren)).replace("**","^") for e in eqs)+";",
         "ideal G=std(I);","print(\"%s: \"+string(size(G))+\" gens, dim=\"+string(dim(G))+\" first=\"+string(G[1]));"%name,"quit;"]
    open(f"{name}.sing","w").write("\n".join(out))
    return len(eqs),len(vs)
S=sp.Symbol
z=lambda k,l,i,j: S(f"z{k}{l}_{i}{j}")
a,b,lam=S('a'),S('b'),S('lam')
# N=2 exceptional: r1 || r2 : r1=e1, r2=lam*e1
print(build(2,{z(1,4,0,0):1,z(1,4,0,1):0,z(2,4,0,0):lam,z(2,4,0,1):0},[lam],name="ex2_par"))
# N=3: rank R <= 2 strata (kernel vector e3: third components zero)
base3={z(k,4,0,2):0 for k in (1,2,3)}
d=dict(base3); d.update({z(1,4,0,0):1,z(1,4,0,1):0,z(2,4,0,0):0,z(2,4,0,1):1}); print(build(3,d,name="ex3_12"))
d=dict(base3); d.update({z(1,4,0,0):1,z(1,4,0,1):0,z(3,4,0,0):0,z(3,4,0,1):1}); print(build(3,d,name="ex3_13"))
d=dict(base3); d.update({z(2,4,0,0):1,z(2,4,0,1):0,z(3,4,0,0):0,z(3,4,0,1):1}); print(build(3,d,name="ex3_23"))
d={z(k,4,0,j):0 for k in (1,2,3) for j in (1,2)}; d.update({z(1,4,0,0):1}); print(build(3,d,name="ex3_rank1"))
