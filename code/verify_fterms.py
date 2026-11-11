"""Check the F-term equations listed in Appendix A of the paper against symbolic
differentiation of the superpotential, for N = 1, 2, 3, 4 (no gauge fixing)."""
import sympy as sp
from gen import system
def check(N):
    eqs,vs,c,W=system(N,gauge=False)
    Z=lambda k,l: sp.Matrix(N if k==4 else 1, N if l==4 else 1, lambda i,j: sp.Symbol(f"z{k}{l}_{i}{j}"))
    P=lambda k,m: sp.Matrix(N if k==4 else 1, N if k==4 else 1, lambda i,j: sp.Symbol(f"p{k}{m}_{i}{j}"))
    C=lambda k,l: c[(min(k,l),max(k,l))]
    I=sp.eye(N); s=lambda k,m: P(k,m)[0,0]
    # paper's list: field -> matrix expression that must vanish (shape = transpose of field)
    L={
     (1,2): (s(1,3)-s(2,3))*Z(2,1)+Z(2,3)*Z(3,1)+Z(2,4)*Z(4,1),
     (2,1): (s(1,3)-s(2,3))*Z(1,2)+Z(1,3)*Z(3,2)+Z(1,4)*Z(4,2),
     (1,3): (s(3,2)-s(1,2))*Z(3,1)+Z(3,2)*Z(2,1)-Z(3,4)*Z(4,1),
     (3,1): (s(3,2)-s(1,2))*Z(1,3)+Z(1,2)*Z(2,3)+Z(1,4)*Z(4,3),
     (2,3): (s(2,1)-s(3,1))*Z(3,2)+Z(3,1)*Z(1,2)+Z(3,4)*Z(4,2),
     (3,2): (s(2,1)-s(3,1))*Z(2,3)+Z(2,1)*Z(1,3)+Z(2,4)*Z(4,3),
     (1,4): (s(1,1)*I-P(4,1))*Z(4,1)+Z(4,2)*Z(2,1)+Z(4,3)*Z(3,1),
     (2,4): (s(2,2)*I-P(4,2))*Z(4,2)+Z(4,1)*Z(1,2)+Z(4,3)*Z(3,2),
     (3,4): (s(3,3)*I-P(4,3))*Z(4,3)-Z(4,1)*Z(1,3)+Z(4,2)*Z(2,3),
     (4,1): Z(1,4)*(s(1,1)*I-P(4,1))+Z(1,2)*Z(2,4)-Z(1,3)*Z(3,4),
     (4,2): Z(2,4)*(s(2,2)*I-P(4,2))+Z(2,1)*Z(1,4)+Z(2,3)*Z(3,4),
     (4,3): Z(3,4)*(s(3,3)*I-P(4,3))+Z(3,1)*Z(1,4)+Z(3,2)*Z(2,4)}
    bad=0
    for (k,l),E in L.items():
        F=Z(k,l)
        for i in range(F.rows):
            for j in range(F.cols):
                d=sp.expand(sp.diff(W,F[i,j]))
                if sp.expand(d-E[j,i])!=0: bad+=1
    # adjoint equations (F1)-(F3); eps sign absorbed: compare up to sign
    for (k,m,l) in [(1,3,2),(2,3,1),(1,2,3),(3,2,1),(2,1,3),(3,1,2)]:
        d=sp.expand(sp.diff(W,P(k,m)[0,0])); e=sp.expand(Z(k,l)[0,0]*Z(l,k)[0,0]+C(k,l))
        if sp.expand(d-e)!=0 and sp.expand(d+e)!=0: bad+=1
    for k in (1,2,3):
        d=sp.expand(sp.diff(W,P(k,k)[0,0])); e=sp.expand((Z(k,4)*Z(4,k))[0,0]+N*C(k,4))
        if sp.expand(d-e)!=0: bad+=1
        a,b=[(k%3)+1,((k+1)%3)+1]
        E=-(Z(4,k)*Z(k,4)+P(4,a)*P(4,b)-P(4,b)*P(4,a)+C(k,4)*I)
        for i in range(N):
            for j in range(N):
                if sp.expand(sp.diff(W,P(4,k)[i,j])-E[j,i])!=0: bad+=1
    return bad
for N in (1,2,3,4): print("N =",N,"mismatches:",check(N))
