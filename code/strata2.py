"""N=3 rank-two strata with the residual gauge fixed completely.
Stratum (a,b): r_a=e1, r_b=e2, third row in span(e1,e2).  Some e_i Phi^(4)_j
(i=1,2) leaves span(e1,e2) (Lemma rowspan); chart sets that row equal to e3."""
import sympy as sp
exec(open("strata.py").read().split("S=sp.Symbol")[0])
S=sp.Symbol
z=lambda k,l,i,j: S(f"z{k}{l}_{i}{j}")
p=lambda m,i,j: S(f"p4{m}_{i}{j}")
for (a,b) in [(1,2),(1,3),(2,3)]:
    c=({1,2,3}-{a,b}).pop()
    base={z(k,4,0,2):0 for k in (1,2,3)}
    base.update({z(a,4,0,0):1,z(a,4,0,1):0,z(b,4,0,0):0,z(b,4,0,1):1})
    # rows e1 (=r_a) Phi_a and e2 (=r_b) Phi_b already lie in R
    for (i,j) in [(0,m) for m in (1,2,3) if m!=a]+[(1,m) for m in (1,2,3) if m!=b]:
        d=dict(base); d.update({p(j,i,0):0,p(j,i,1):0,p(j,i,2):1})
        name=f"ch3_{a}{b}_e{i+1}P{j}"
        print(name, build(3,d,name=name), flush=True)
