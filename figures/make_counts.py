"""Figure: the predicted number of vacua -c(4N) and the values established here."""
import sys, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
pred=[2,12,56,208,684,2032,5616,14592,36088,85500,195312]
n4=sys.argv[1] if len(sys.argv)>1 else "none"   # "mod" if N=4 was counted modulo a prime
pred=pred[1:]; N=np.arange(1,len(pred)+1)
fig,ax=plt.subplots(figsize=(7,4.4))
ax.semilogy(N,pred,"-",color="#999999",lw=1,zorder=1)
k=4 if n4=="mod" else 3
ax.scatter(N[k:],pred[k:],s=30,facecolors="white",edgecolors="#555555",zorder=2,label=r"predicted $-\hat c(4N)$, not counted")
ax.scatter([1],[12],s=70,marker="D",color="#1f77b4",zorder=3,label=r"$N=1$: exact for generic $c$")
ax.scatter([2,3],[56,208],s=70,marker="o",color="#c0392b",zorder=3,label=r"$N=2,3$: exact at $c_0$, lower bound for generic $c$")
if n4=="mod":
    ax.scatter([4],[684],s=70,marker="s",color="#e67e22",zorder=3,label=r"$N=4$: count modulo a large prime")
for n,v in zip(N[:4],pred[:4]): ax.annotate(str(v),(n,v),textcoords="offset points",xytext=(6,-12),fontsize=8)
ax.set_xlabel("number of D6-branes $N$"); ax.set_ylabel("number of vacua")
ax.set_ylim(6,4e5); ax.set_xticks(N); ax.legend(fontsize=8,loc="upper left",frameon=False); ax.grid(alpha=0.25,which="both")
plt.tight_layout(); plt.savefig("fig_counts.png",dpi=200)
