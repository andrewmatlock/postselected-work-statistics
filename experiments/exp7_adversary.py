"""The forgery class: generalized detailed balance w.r.t. any reference distribution
reproduces the factorization exactly. A bath at the wrong temperature suffices."""
import numpy as np
ES=np.array([-1.5,-0.5,0.5,1.5])
def sep_residual(D):
    rows=[];vals=[]
    for i in range(4):
        for j in range(4):
            if not np.isnan(D[i,j]):
                r=np.zeros(8); r[i]=-1; r[4+j]=1; rows.append(r); vals.append(D[i,j])
    A=np.array(rows); b=np.array(vals); x,*_=np.linalg.lstsq(A,b,rcond=None)
    return float(np.abs(A@x-b).max()), float(np.abs(b).mean())
print("=== Exact forgery by a wrong-temperature bath ===")
for bp in [0.2,0.4,0.8]:
    P=np.zeros((4,4))
    for i in range(4):
        for j in range(4):
            if i!=j: P[i,j]=(1/3)*min(1,np.exp(-bp*(ES[j]-ES[i])))
        P[i,i]=1-P[i].sum()
    D=np.array([[np.log(P[i,j]/P[j,i]) if i!=j else 0.0 for j in range(4)] for i in range(4)])
    res,mean=sep_residual(D)
    anti=np.linalg.norm(D+D.T)
    print(f"  beta'={bp}: residual={res:.2e}  mean|D|={mean:.3f}  antisymmetry ||D+D^T||={anti:.1e}"
          f"  max|D_ii|={max(abs(D[k,k]) for k in range(4)):.1e}")
print("\n=> single potential (antisymmetric, D_ii=0) marks reference equilibrium;")
print("   two potentials (D_ii != 0) marks post-selection. Taxonomy discriminator.")
