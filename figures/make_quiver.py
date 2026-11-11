"""Figure 1: the quiver of the D2-D2-D2-D6 system drawn on a tetrahedron."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa

V={1:np.array([1.0,0.0,-0.35]),2:np.array([-0.5,0.866,-0.35]),3:np.array([-0.5,-0.866,-0.35]),4:np.array([0,0,1.05])}
col={1:"#3b6fb6",2:"#3b6fb6",3:"#3b6fb6",4:"#c0392b"}
fig=plt.figure(figsize=(6.4,5.6)); ax=fig.add_subplot(111,projection="3d",computed_zorder=False)
def arrow(a,b,off,color):
    p,q=V[a],V[b]; d=q-p; n=np.cross(d,np.array([0.1,0.2,1.0])); n/=np.linalg.norm(n)
    p2=p+0.13*d/np.linalg.norm(d)+off*n; q2=q-0.13*d/np.linalg.norm(d)+off*n
    ax.quiver(*p2,*(q2-p2),color=color,arrow_length_ratio=0.12,linewidth=1.4)
for a in range(1,5):
    for b in range(a+1,5):
        c="#555555" if b<4 else "#8e44ad"
        arrow(a,b,0.045,c); arrow(b,a,0.045,c)
for k,p in V.items():
    ax.scatter(*p,s=900 if k==4 else 600,color=col[k],depthshade=False,edgecolors="k",zorder=5)
    ax.text(*(p+np.array([0,0,0.0])),f"{k}",color="white",fontsize=13,ha="center",va="center",zorder=10,weight="bold")
    t=np.linspace(0,2*np.pi,60); r=0.11
    for j in range(3):
        cen=p+np.array([0,0,0.16+0.06*j]) if k==4 else p*np.array([1.25,1.25,1])+np.array([0,0,0.05*j])
        ax.plot(cen[0]+r*np.cos(t),cen[1]+r*np.sin(t),cen[2]+0*t,color=col[k],lw=0.9,alpha=0.8)
lab={1:r"U(1)  D2 on $T^2_{45}$",2:r"U(1)  D2 on $T^2_{67}$",3:r"U(1)  D2 on $T^2_{89}$",4:r"U(N)  D6 on $T^6$"}
LP={1:(1.45,0.0,-0.7),2:(-0.9,1.3,-0.55),3:(0.05,-1.45,-0.62),4:(0,0,1.5)}
for k,p in V.items():
    ax.text(*LP[k],lab[k],fontsize=9.5,ha="center")
ax.set_axis_off(); ax.view_init(elev=24,azim=-118); ax.set_box_aspect((1,1,0.95))
plt.tight_layout(); plt.savefig("fig_quiver.png",dpi=200,bbox_inches="tight")
