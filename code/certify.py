import re,ast,sys
from flint import fmpq_poly, fmpq, fmpz
import sympy as sp
def parse(fn):
    s=open(fn).read().strip().rstrip(':')
    s2=re.sub(r'(-?\d+) / 2\^(\d+)', r'("\1",\2)', s).replace("'", '"')
    T=ast.literal_eval(s2)
    P=T[1]; nv,deg,vars_,lin=P[1],P[2],P[3],P[4]
    el,den,nums=P[5][1]
    poly=lambda d: fmpq_poly([int(a) for a in d[1]])
    return vars_,lin,poly(el),poly(den),[(poly(n[0]),int(n[1])) for n in nums]
def certify(msfile, outfile):
    vars_,lin,chi,den,nums=parse(outfile)
    L=open(msfile).read().split("\n",2)
    xs=[sp.Symbol(v) for v in L[0].split(",")]
    polys=[sp.Poly(sp.sympify(t.replace("^","**")),*xs) for t in L[2].strip().split(",\n")]
    # separating element: linear form sum lin_i x_i ; msolve: param variable t = that form
    assert sum(1 for a in lin if a)==1 and lin[-1]==1, lin
    # variable order in output may be permuted
    order=vars_
    sep=order[-1]
    # x_i = -num_i(t)/(c_i * den(t)) for order[:-1]
    inv_den=None
    # compute den inverse mod chi
    g,s_,t_=chi.xgcd(den) if hasattr(chi,'xgcd') else (None,None,None)
    vals={}
    # work modulo chi with fractions: represent each x as (numerator poly) * inv(den) mod chi
    # inverse of den mod chi via extended gcd
    def inv_mod(a,m):
        # extended Euclid over Q[t]
        r0,r1=m,a%m; s0,s1=fmpq_poly(0),fmpq_poly(1)
        while r1!=0:
            q=r0//r1; r0,r1=r1,r0-q*r1; s0,s1=s1,s0-q*s1
        assert r0.degree()==0, "not coprime"
        return (s0*fmpq_poly([1/r0[0]]))%m if False else (s0* fmpq_poly([fmpq(1)/r0[0]]))%m
    di=inv_mod(den,chi)
    for (num,cc),v in zip(nums,order[:-1]):
        vals[v]=(-num*di*fmpq_poly([fmpq(1,cc)]))%chi
    vals[sep]=fmpq_poly([0,1])
    # squarefree check
    g=chi.gcd(chi.derivative())
    print("chi degree",chi.degree(),"gcd(chi,chi') degree",g.degree())
    worst=0
    for P in polys:
        acc=fmpq_poly(0)
        for mon,co in P.terms():
            term=fmpq_poly([fmpq(int(sp.numer(co)),int(sp.denom(co)))])
            for v,e in zip(xs,mon):
                if e: term=(term*pow(vals[str(v)],e))%chi
            acc=(acc+term)%chi
        if acc!=0: worst+=1
    print("equations not vanishing:",worst,"of",len(polys))
if __name__=="__main__": certify(sys.argv[1],sys.argv[2])
