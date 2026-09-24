"""D-CTC deviations are non-separable; P-CTC deviations factorize exactly.
A work-statistics discriminator between the two standard CTC models."""
import numpy as np
rng=np.random.default_rng(17)
def ptrace(M,keep):
    T=M.reshape(4,2,4,2)
    return np.einsum('abcb->ac',T) if keep=='S' else np.einsum('abad->bd',T)
def fixed_point(U,rho):
    tau=np.eye(2,dtype=complex)/2
    for _ in range(600):
        new=ptrace(U@np.kron(rho,tau)@U.conj().T,'ctc')
        if np.linalg.norm(new-tau)<1e-13: return new
        tau=new
    return tau
def sep(D):
    rows=[];vals=[]
    for i in range(4):
        for j in range(4):
            if not np.isnan(D[i,j]):
                r=np.zeros(8); r[i]=-1; r[4+j]=1; rows.append(r); vals.append(D[i,j])
    A=np.array(rows); b=np.array(vals); x,*_=np.linalg.lstsq(A,b,rcond=None)
    return float(np.abs(A@x-b).max()), float(np.abs(b).mean())
print("=== D-CTC vs P-CTC deviation structure (identical unitaries) ===")
print("            D-CTC resid  mean|D|  ratio |  P-CTC resid")
for k in range(4):
    A=rng.normal(size=(8,8))+1j*rng.normal(size=(8,8)); U,_=np.linalg.qr(A)
    PF=np.array([np.real(np.diag(ptrace(U@np.kron(np.diag(np.eye(4)[i]).astype(complex),
        fixed_point(U,np.diag(np.eye(4)[i]).astype(complex)))@U.conj().T,'S'))) for i in range(4)])
    Ud=U.conj().T
    PR=np.array([np.real(np.diag(ptrace(Ud@np.kron(np.diag(np.eye(4)[j]).astype(complex),
        fixed_point(Ud,np.diag(np.eye(4)[j]).astype(complex)))@Ud.conj().T,'S'))) for j in range(4)])
    D=np.full((4,4),np.nan)
    for i in range(4):
        for j in range(4):
            if PF[i,j]>1e-12 and PR[j,i]>1e-12: D[i,j]=np.log(PF[i,j]/PR[j,i])
    rD,mD=sep(D)
    C=np.einsum('akbk->ab',U.reshape(4,2,4,2))
    u=np.log(np.diag(C@C.conj().T).real); v=np.log(np.diag(C.conj().T@C).real)
    rL,_=sep(np.array([[u[j]-v[i] for j in range(4)] for i in range(4)]))
    print(f"  haar #{k}:  {rD:9.4f} {mD:8.3f} {rD/mD:7.2f} |  {rL:.1e}")
