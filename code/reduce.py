import sympy as sp
def linelim(eqs, vs, protect=()):
    eqs=[sp.expand(e) for e in eqs if sp.expand(e)!=0]
    subs={}
    changed=True
    while changed:
        changed=False
        best=None
        for e in eqs:
            P=sp.Poly(e,*vs)
            for v in vs:
                if v in protect: continue
                d=P.degree(v)
                if d!=1: continue
                # coefficient of v must be a nonzero constant and v must not appear multiplied
                q,r=sp.div(e,v)
                q=sp.expand(q)
                if q.free_symbols & set(vs): continue
                if q==0: continue
                cand=(len(sp.Add.make_args(r)), e, v, sp.expand(-r/q))
                if best is None or cand[0]<best[0]: best=cand
        if best:
            _,e,v,val=best
            subs={k:sp.expand(x.subs(v,val)) for k,x in subs.items()}; subs[v]=val
            eqs=[sp.expand(x.subs(v,val)) for x in eqs]
            eqs=[x for x in eqs if x!=0]
            vs=[x for x in vs if x!=v]
            changed=True
    # dedupe
    seen=set(); out=[]
    for e in eqs:
        k=sp.srepr(sp.expand(e)); k2=sp.srepr(sp.expand(-e))
        if k in seen or k2 in seen: continue
        seen.add(k); out.append(e)
    return out, vs, subs
