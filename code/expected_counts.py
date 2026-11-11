"""Coefficients of -theta_1(z|tau)^2/eta(tau)^6 at y^0: the predicted number of vacua -c(4N)."""
import sympy as sp
q,y=sp.symbols('q y')
M=10
P=(y-2+1/y)
for n in range(1,M+1):
    P=sp.expand(sp.series(P*(1-q**n*y)**2*(1-q**n/y)**2*sum(sp.binomial(k+3,3)*q**(n*k) for k in range(M//n+1)),q,0,M+1).removeO())
P=sp.expand(P)
vals=[-P.coeff(q,N).coeff(y,0) for N in range(M+1)]
print(vals)
