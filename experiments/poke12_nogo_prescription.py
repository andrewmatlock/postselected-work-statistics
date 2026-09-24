"""The attribution no-go, constructive: the Halmos dilation of any P-CTC operator gives
an ordinary heralded device with identical accepted-run statistics."""
import numpy as np
from scipy.linalg import sqrtm
rng=np.random.default_rng(9)
ES=np.array([-1.5,-0.5,0.5,1.5]); pth=np.exp(-ES); pth/=pth.sum()
A=rng.normal(size=(8,8))+1j*rng.normal(size=(8,8)); U,_=np.linalg.qr(A)
C=np.einsum('akbk->ab',U.reshape(4,2,4,2))
K=C/np.linalg.svd(C,compute_uv=False).max()
V=np.block([[K,sqrtm(np.eye(4)-K@K.conj().T)],[sqrtm(np.eye(4)-K.conj().T@K),-K.conj().T]])
print("=== Attribution no-go (Halmos dilation) ===")
print(f"dilation unitarity ||V+V - I|| = {np.linalg.norm(V.conj().T@V-np.eye(8)):.2e}")
mx=0
for i in range(4):
    a=np.abs(C[:,i])**2/np.vdot(C[:,i],C[:,i]).real
    b=np.abs(V[:4,:4][:,i])**2/np.vdot(V[:4,:4][:,i],V[:4,:4][:,i]).real
    mx=max(mx,np.abs(a-b).max())
print(f"forward conditionals |P-CTC - heralded| max = {mx:.2e}")
m1=pth*np.diag(C.conj().T@C).real; m1/=m1.sum()
m2=pth*np.diag(K.conj().T@K).real; m2/=m2.sum()
print(f"accepted-run initial marginals differ by {np.abs(m1-m2).max():.2e}")
print(f"the biased past: |marginal - thermal| = {np.abs(m1-pth).max():.4f}")
print(f"Jarzynski efficacy <e^-bW> = "
      f"{float(sum(pth[i]*abs(C[j,i])**2*np.exp(-(ES[j]-ES[i])) for i in range(4) for j in range(4))/sum(pth[k]*np.vdot(C[:,k],C[:,k]).real for k in range(4))):.4f}")
print("\n=> only the DISCARD COUNT distinguishes them: attribution needs the denominator.")
