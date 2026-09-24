import numpy as np
rng=np.random.default_rng(17)
ES=np.array([-1.5,-0.5,0.5,1.5]); beta=1.0
pth=np.exp(-beta*ES); pth/=pth.sum()

def ptrace(M,keep):
    T=M.reshape(4,2,4,2)
    return np.einsum('abcb->ac',T) if keep=='S' else np.einsum('abad->bd',T)
def fixed_point(U,rho):
    tau=np.eye(2,dtype=complex)/2
    for _ in range(800):
        new=ptrace(U@np.kron(rho,tau)@U.conj().T,'ctc')
        if np.linalg.norm(new-tau)<1e-13: return new
        tau=new
    return tau
def sep_res(D):
    rows=[];vals=[]
    for i in range(4):
        for j in range(4):
            if not np.isnan(D[i,j]):
                r=np.zeros(8); r[i]=-1; r[4+j]=1; rows.append(r); vals.append(D[i,j])
    A=np.array(rows); b=np.array(vals); x,*_=np.linalg.lstsq(A,b,rcond=None)
    return float(np.abs(A@x-b).max()), float(np.abs(b).mean())

print("=== D-CTC non-separability under BOTH mixture conventions ===")
print("  (A) per-branch tau  |  (B) single tau from the full thermal mixture")
print("        conv A: resid  mean|D| | conv B: resid  mean|D| | P-CTC resid")
for k in range(5):
    A_=rng.normal(size=(8,8))+1j*rng.normal(size=(8,8)); U,_=np.linalg.qr(A_)
    rows={}
    for conv in ('A','B'):
        if conv=='B':
            rhoF=np.diag(pth).astype(complex); tauF=fixed_point(U,rhoF)
            rhoR=np.diag(pth).astype(complex); tauR=fixed_point(U.conj().T,rhoR)
        D=np.full((4,4),np.nan)
        PF=np.zeros((4,4)); PR=np.zeros((4,4))
        for i in range(4):
            rho=np.zeros((4,4),dtype=complex); rho[i,i]=1
            tau=fixed_point(U,rho) if conv=='A' else tauF
            PF[i]=np.real(np.diag(ptrace(U@np.kron(rho,tau)@U.conj().T,'S')))
        for j in range(4):
            rho=np.zeros((4,4),dtype=complex); rho[j,j]=1
            Ud=U.conj().T
            tau=fixed_point(Ud,rho) if conv=='A' else tauR
            PR[j]=np.real(np.diag(ptrace(Ud@np.kron(rho,tau)@Ud.conj().T,'S')))
        for i in range(4):
            for j in range(4):
                if PF[i,j]>1e-12 and PR[j,i]>1e-12: D[i,j]=np.log(PF[i,j]/PR[j,i])
        rows[conv]=sep_res(D)
    C=np.einsum('akbk->ab',U.reshape(4,2,4,2))
    u=np.log(np.diag(C@C.conj().T).real); v=np.log(np.diag(C.conj().T@C).real)
    DL=np.array([[u[j]-v[i] for j in range(4)] for i in range(4)])
    rL,_=sep_res(DL)
    print(f"  #{k}:   {rows['A'][0]:8.4f} {rows['A'][1]:7.3f} |"
          f" {rows['B'][0]:8.4f} {rows['B'][1]:7.3f} | {rL:.1e}")
