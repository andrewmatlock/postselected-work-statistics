"""
analyze_shallow.py — analysis for the shallow-circuit experiment, tolerant of
partial data (missing reverse inputs simply drop the corresponding D column).

Usage:  python analyze_shallow.py ctc_shallow_salvaged.json
"""
import json, sys
import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator

path = sys.argv[1] if len(sys.argv) > 1 else "ctc_shallow_salvaged.json"
d = json.load(open(path))
P = np.array(d.get("ansatz_params",
    [2.774141, 1.438552, 5.364134, 1.634255, 2.682839, 4.893008,
     1.171051, 1.459015, 1.474151, 0.031240, 3.189096, 2.607931]))

qc = QuantumCircuit(3)
qc.ry(P[0],0); qc.ry(P[1],1); qc.ry(P[2],2); qc.cx(1,2)
qc.rz(P[3],1); qc.ry(P[4],2); qc.cx(0,2)
qc.rz(P[5],0); qc.ry(P[6],2); qc.cx(1,0)
qc.ry(P[7],0); qc.ry(P[8],1); qc.cx(0,2)
qc.ry(P[9],0); qc.rz(P[10],1); qc.ry(P[11],2)
M = Operator(qc).data[:4,:4]
u = np.log(np.diag(M@M.conj().T).real); v = np.log(np.diag(M.conj().T@M).real)
D_th = np.array([[u[j]-v[i] for j in range(4)] for i in range(4)])
acc_th = np.exp(v)
ES = np.array([-1.5,-0.5,0.5,1.5]); pth = np.exp(-ES); pth /= pth.sum()

def parse(block):
    P_ = np.full((4,4), np.nan); N = np.zeros((4,4)); acc = np.full(4, np.nan)
    for k, ci in block.items():
        i = int(k); tot = sum(ci.values()); kept = 0
        for bits, n in ci.items():
            if bits[0]=='0':
                N[i, int(bits[1:],2)] += n; kept += n
        acc[i] = kept/tot
        P_[i] = N[i]/max(kept,1)
    return P_, N, acc

PF, NF, accF = parse(d["forward"]); PR, NR, accR = parse(d["reverse"])
print(f"{path}  backend={d['backend']}")
print(f"forward inputs: {sorted(int(k) for k in d['forward'])}  "
      f"reverse inputs: {sorted(int(k) for k in d['reverse'])}")

ok = ~np.isnan(accF)
print(f"\n[acceptance]  measured {np.round(accF,3)}   predicted {np.round(acc_th,3)}")
if ok.all():
    marg = pth*accF; marg /= marg.sum()
    Ntot = sum(sum(c.values()) for c in d["forward"].values())
    Na = float((pth*accF).sum()*Ntot)
    z = np.abs(marg-pth)/np.sqrt(pth*(1-pth)/Na)
    marg_th = pth*acc_th; marg_th /= marg_th.sum()
    print(f"[biased past] deviation {np.abs(marg-pth).max():.4f} "
          f"(predicted {np.abs(marg_th-pth).max():.4f})  max z = {z.max():.1f} sigma")

D = np.full((4,4), np.nan); S = np.full((4,4), np.nan)
for i in range(4):
    for j in range(4):
        if np.isnan(PF[i,j]) or np.isnan(PR[j,i]): continue
        if NF[i,j] < 200 or NR[j,i] < 200: continue
        D[i,j] = np.log(PF[i,j]/PR[j,i])
        S[i,j] = np.sqrt((1-PF[i,j])/NF[i,j] + (1-PR[j,i])/NR[j,i])
used = ~np.isnan(D)
print(f"\n[deviation matrix] cells usable: {used.sum()}/16")

rows, vals, sig = [], [], []
for i in range(4):
    for j in range(4):
        if used[i,j]:
            r = np.zeros(8); r[i] = -1; r[4+j] = 1
            rows.append(r); vals.append(D[i,j]); sig.append(S[i,j])
Am = np.array(rows); b = np.array(vals); s = np.array(sig)
x,*_ = np.linalg.lstsq(Am/s[:,None], b/s, rcond=None)
resid = Am@x - b
print(f"   separability residual = {np.abs(resid).max():.4f}   mean|D| = {np.abs(b).mean():.3f}")
print(f"   RATIO residual/mean|D| = {np.abs(resid).max()/np.abs(b).mean():.3f}"
      f"   <-- generic driven channels ~1.0, exact theory 0.0")
print(f"   chi2/dof = {float(((resid/s)**2).sum())/max(len(b)-7,1):.1f}")
print(f"   max|D_meas - D_th| = {np.nanmax(np.abs(D-D_th)):.4f}")
print("\n measured D:\n", np.round(D,3))
print(" theory D:\n", np.round(D_th,3))
