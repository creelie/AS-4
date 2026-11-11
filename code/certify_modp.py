import ast, sys, re
from flint import nmod_poly
def run(msfile, outfile):
    T=ast.literal_eval(open(outfile).read().strip().rstrip(':'))
    P=T[1]; p,nv,deg,vars_,lin=P[0],P[1],P[2],P[3],P[4]
    el,den,params=P[5][1]
    poly=lambda c: nmod_poly([x%p for x in c],p)
    chi=poly(el[1]); d=poly(den[1])
    assert sum(1 for a in lin if a)==1 and lin[-1]==1
    g=chi.gcd(chi.derivative()); print("p",p,"deg chi",chi.degree(),"gcd(chi,chi') deg",g.degree(),"den",den)
    L=open(msfile).read().split("\n",2)
    names=L[0].strip().split(",")
    eqs=[e.strip() for e in L[2].replace("\n","").split(",") if e.strip()]
    for sign in (-1,1):
        vals={}
        for ent,v in zip(params,vars_[:-1]):
            num,c=(ent[0],ent[1]) if len(ent)==2 and not isinstance(ent[1],list) else (ent[0],1)
            inv=pow(c% p,-1,p)
            vals[v]=(poly(num[1])*nmod_poly([sign*inv%p],p))%chi
        vals[vars_[-1]]=nmod_poly([0,1],p)
        bad=0
        for e in eqs:
            acc=nmod_poly([0],p)
            for term in re.split(r'(?=[+-])',e.replace(' ','')):
                if not term: continue
                parts=term.split('*'); coef=1; mon=nmod_poly([1],p)
                for f in parts:
                    if f in ('+','-',''): continue
                    neg=f.startswith('-'); f=f.lstrip('+-')
                    if neg: coef=-coef
                    if '^' in f: v,k=f.split('^'); k=int(k)
                    else: v,k=f,1
                    if v[0].isdigit(): coef*=int(v)**k
                    else: mon=(mon*pow(vals[v],k))%chi if k<4 else (mon*pow(vals[v],k,chi))%chi
                acc=(acc+mon*nmod_poly([coef%p],p))%chi
            if acc!=0: bad+=1
        print("sign",sign,"equations not vanishing:",bad,"of",len(eqs))
run(sys.argv[1],sys.argv[2])
