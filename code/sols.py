import sys, random, pickle, sympy as sp, mpmath as mp
from flint import fmpz_poly, fmpq_poly, acb
from gen import system
from reduce import linelim
from certify import parse
N=int(sys.argv[1]); out=sys.argv[2]; seed=1
eqs0,vs0,c,W=system(N)
random.seed(seed)
cv={s: sp.Rational(random.choice([-1,1])*random.randint(1,30), random.randint(1,30)) for s in c.values()}
eqs=[sp.expand(e.subs(cv)) for e in eqs0]
eqs,vs,subs=linelim(eqs,vs0)
vars_,lin,chi,den,nums=parse(out)
if len(vars_)==len(vs0): vs,subs=vs0,{}
# chi in Q[t]; numeric roots
num=fmpz_poly([int(x) for x in (chi*chi.denom()).numer().coeffs()]) if hasattr(chi,'denom') else None
q=fmpq_poly(chi); dd=q.denom(); zp=fmpz_poly([int(q[i]*dd) for i in range(q.degree()+1)])
roots=[r for r,m in zp.complex_roots()]
assert len(roots)==chi.degree()
mp.mp.dps=60
def ev(p,t):
    s=mp.mpc(0)
    for i in range(p.degree(),-1,-1): s=s*t+mp.mpf(int(p[i].p))/int(p[i].q)
    return s
sols=[]
idx={f"x{i+1}":v for i,v in enumerate(vs)}
for r in roots:
    t=mp.mpc(mp.mpf(r.real.mid().str(50,radius=False)),mp.mpf(r.imag.mid().str(50,radius=False)))
    val={}
    for (nm,cc),v in zip(nums,vars_[:-1]):
        val[idx[v]]=-ev(nm,t)/(cc*ev(den,t))
    val[idx[vars_[-1]]]=t
    base={s:sp.Float(str(x.real),50)+sp.I*sp.Float(str(x.imag),50) for s,x in val.items()}
    for k,e in subs.items():
        val[k]=complex(sp.N(e.subs(base),50))
    val={k:complex(v) if not isinstance(v,mp.mpc) else complex(v) for k,v in val.items()}
    sols.append(val)
# residuals on the full gauge-fixed system
fl=[sp.lambdify(vs0,e,"math") for e in [sp.expand(e.subs(cv)) for e in eqs0]]
import gen
Wf=sp.lambdify(vs0,W.subs(gen.LASTSUB).subs(cv),"math")
import cmath
res=0; Ws=[]
gaugevals={}
for s in sols:
    args=[s[v] for v in vs0]
    res=max(res,max(abs(f(*args)) for f in fl))
    Ws.append(Wf(*[s.get(v,0) for v in vs0]))
print("N",N,"solutions",len(sols),"max residual",res)
pickle.dump((sols,Ws,{str(k):float(v) for k,v in cv.items()}),open(f"sols{N}.pkl","wb"))
