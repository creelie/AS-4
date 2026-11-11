"""Figures 2 and 3: the vacua for N=2 and N=3 at c0, the coordinate u and the critical values."""
import pickle, sys, numpy as np, sympy as sp, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
for N in (2,3):
    sols,Ws,cv=pickle.load(open(f"data/sols{N}.pkl","rb"))
    u=np.array([s[sp.Symbol('z31_00')] for s in sols])   # Z12=Z23=1 in the slice, so u=Z31
    W=np.array(Ws,dtype=complex)
    fig,axs=plt.subplots(1,2,figsize=(11,4.8))
    for ax,z,name in ((axs[0],u,r"$u=Z^{(12)}Z^{(23)}Z^{(31)}$"),(axs[1],W,r"critical value $f_{c_0}$")):
        done=set()
        for i,a in enumerate(z):
            j=int(np.argmin(abs(z+a)))
            if (j,i) in done: continue
            done.add((i,j))
            ax.plot([a.real,-a.real],[a.imag,-a.imag],color="#bbbbbb",lw=0.5,zorder=1)
        ax.scatter(z.real,z.imag,s=14 if N==2 else 8,c=np.angle(z),cmap="twilight",edgecolors="k",linewidths=0.3,zorder=3)
        ax.axhline(0,color="k",lw=0.4); ax.axvline(0,color="k",lw=0.4)
        ax.set_xlabel("real part"); ax.set_ylabel("imaginary part"); ax.set_title(name+f"   ({len(z)} vacua, N={N})",fontsize=10)
        if N==3: ax.set_xscale("symlog",linthresh=1); ax.set_yscale("symlog",linthresh=1)
        ax.set_aspect("equal" if N==2 else "auto")
    plt.tight_layout(); plt.savefig(f"fig_N{N}.png",dpi=200,bbox_inches="tight"); plt.close()
