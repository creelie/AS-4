import sympy as sp, sys, itertools
def system(N, gauge=True):
    eps=lambda a,b,c: sp.LeviCivita(a,b,c)
    Z={}; 
    for k in range(1,5):
        for l in range(1,5):
            if k==l: continue
            r = N if k==4 else 1; s = N if l==4 else 1
            Z[k,l]=sp.Matrix(r,s,lambda i,j: sp.Symbol(f"z{k}{l}_{i}{j}"))
    Phi={}
    for k in range(1,5):
        for m in range(1,4):
            r = N if k==4 else 1
            Phi[k,m]=sp.Matrix(r,r,lambda i,j: sp.Symbol(f"p{k}{m}_{i}{j}"))
    c={(k,l):sp.Symbol(f"c{k}{l}") for k in range(1,5) for l in range(k+1,5)}
    cc=lambda k,l: c[(min(k,l),max(k,l))]
    tr=lambda M: M.trace()
    W=0
    for k,l,m in itertools.permutations([1,2,3]):
        W+= eps(k,l,m)*(Phi[k,m]*Z[k,l]*Z[l,k])[0,0]
    for k in range(1,4):
        W+= tr(Z[4,k]*Phi[k,k][0,0]*Z[k,4]) - tr(Phi[4,k]*Z[4,k]*Z[k,4])
    T=[(3,1,2,1),(1,3,2,1),(1,2,4,1),(4,2,1,1),(1,3,4,-1),(3,1,4,1),(3,4,2,1),(4,3,2,1)]
    # terms: Z31 Z12 Z23 ; Z13 Z32 Z21 ; Z12 Z24 Z41 ; Z42 Z21 Z14 ; -Z13 Z34 Z41 ; Z31 Z14 Z43 ; Z34 Z42 Z23 ; Z43 Z32 Z24
    W+= (Z[3,1]*Z[1,2]*Z[2,3])[0,0] + (Z[1,3]*Z[3,2]*Z[2,1])[0,0] + tr(Z[1,2]*Z[2,4]*Z[4,1]) + tr(Z[4,2]*Z[2,1]*Z[1,4]) \
        - tr(Z[1,3]*Z[3,4]*Z[4,1]) + tr(Z[3,1]*Z[1,4]*Z[4,3]) + tr(Z[3,4]*Z[4,2]*Z[2,3]) + tr(Z[4,3]*Z[3,2]*Z[2,4])
    for k,l,m in itertools.permutations([1,2,3]):
        W+= cc(k,l)*eps(k,l,m)*Phi[k,m][0,0]
    for k in range(1,4):
        W+= cc(k,4)*(N*Phi[k,k][0,0] - tr(Phi[4,k]))
    P1,P2,P3=Phi[4,1],Phi[4,2],Phi[4,3]
    W+= -tr(P1*(P2*P3-P3*P2))
    W=sp.expand(W)
    allvars=[s for M in list(Z.values())+list(Phi.values()) for s in M]
    eqs=[sp.expand(sp.diff(W,v)) for v in allvars]
    sub={}
    if gauge:
        for (k,m) in [(1,1),(1,2),(1,3),(2,1),(2,2),(3,3)]: sub[Phi[k,m][0,0]]=0
        sub[Z[1,2][0,0]]=1; sub[Z[2,3][0,0]]=1
        rows=[Z[1,4],Z[2,4],Z[3,4]]
        for a in range(min(N,3)):
            for j in range(N): sub[rows[a][0,j]]= 1 if j==a else 0
        if N==4:
            for j in range(N): sub[Phi[4,2][0,j]]= 1 if j==3 else 0
    eqs=[sp.expand(e.subs(sub)) for e in eqs]
    eqs=[e for e in eqs if e!=0]
    vars_=[v for v in allvars if v not in sub]
    global LASTSUB; LASTSUB=sub
    return eqs, vars_, c, W
if __name__=="__main__":
    N=int(sys.argv[1]); eqs,vs,c,W=system(N)
    print(len(eqs),len(vs))
