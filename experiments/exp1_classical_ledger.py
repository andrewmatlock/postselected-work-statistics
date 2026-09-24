"""Classical: symmetric closure preserves the detailed FT; twisted retro-constraints
amputate the entropy-production spectrum (an instance of absolute irreversibility)."""
import numpy as np, itertools
f=[0.55,0.35,0.45]; r=[0.05,0.20,0.10]
M=np.zeros((3,3))
for i in range(3):
    M[i,(i+1)%3]=f[i]; M[i,(i-1)%3]=r[i]; M[i,i]=1-f[i]-r[i]
ev,evec=np.linalg.eig(M.T); pi=np.real(evec[:,np.argmin(np.abs(ev-1))]); pi/=pi.sum()
T=7
trajs=list(itertools.product(range(3),repeat=T+1))
prob=np.empty(len(trajs)); sig=np.empty(len(trajs))
for k,tr in enumerate(trajs):
    p=pi[tr[0]]; s=np.log(pi[tr[0]])-np.log(pi[tr[-1]])
    for t in range(T):
        p*=M[tr[t],tr[t+1]]
        if M[tr[t],tr[t+1]]>0: s+=np.log(M[tr[t],tr[t+1]]/M[tr[t+1],tr[t]])
    prob[k]=p; sig[k]=s
sel=prob>0; prob,sig=prob[sel],sig[sel]; trajs=[t for t,m in zip(trajs,sel) if m]
def dft(w,label):
    w=w/w.sum(); key=np.round(sig,8); out=[]
    for s in sorted(set(key[key>1e-8])):
        Ps=w[np.isclose(key,s)].sum(); Pms=w[np.isclose(key,-s)].sum()
        if Ps>1e-12 and Pms>1e-12: out.append(np.log(Ps/Pms)-s)
    mx=max(abs(o) for o in out) if out else float('nan')
    print(f"[{label}] paired points={len(out)}  max|Crooks deviation|={mx:.2e}")
print("=== Classical consistency ledger ==="); print(f"pi = {np.round(pi,4)}")
dft(prob.copy(),"free ensemble          ")
key=np.round(sig,6)
for twist,name in [(0,"symmetric  x0=xT   "),(1,"twisted    x0=xT+1 "),(2,"twisted    x0=xT+2 ")]:
    mask=np.array([tr[0]==(tr[-1]+twist)%3 for tr in trajs]); w=prob*mask
    dft(w,name); wn=w/w.sum()
    print(f"     support: P(s<0)={wn[key<-1e-8].sum():.4f}  P(0)={wn[np.abs(key)<=1e-8].sum():.4f}  P(s>0)={wn[key>1e-8].sum():.4f}")
