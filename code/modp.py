import sys, sympy as sp
src,p,out=sys.argv[1],int(sys.argv[2]),sys.argv[3]
L=open(src).read().split("\n",2)
xs=[sp.Symbol(v) for v in L[0].split(",")]
polys=[sp.Poly(sp.sympify(t.replace("^","**")),*xs) for t in L[2].strip().split(",\n")]
res=[]
for P in polys:
    terms=[]
    for mon,co in P.terms():
        co=int(co)%p
        if co==0: continue
        m="*".join(f"{x}^{e}" if e>1 else str(x) for x,e in zip(xs,mon) if e)
        terms.append(f"{co}*{m}" if m else str(co))
    res.append("+".join(terms))
open(out,"w").write(L[0]+"\n"+str(p)+"\n"+",\n".join(res)+"\n")
