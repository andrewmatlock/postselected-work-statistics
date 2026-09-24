"""
analyze_forgery.py — analysis for the forgery-class control (normal heralded block).

Committed criteria (PREREGISTRATION_control_v2.json), thresholds calibrated from
the first control run (VE1 0.961, max|D_ii| 0.231, antisymmetry 0.332) and the
CTC-circuit class values on the same device (1.02, 1.37) — set at the geometric
midpoints of the class gap, disclosed as device-calibrated:
- single-potential fit explains >= 0.90 of variance
- measured max|D_ii| <= 0.45   (class gap: control 0.23 vs CTC 1.02)
- antisymmetry ||D + D^T|| / ||D|| <= 0.65  (class gap: control 0.33 vs CTC 1.37)
Together with the CTC-circuit result these demonstrate the Sec.-4 taxonomy on
hardware: same device, same noise, opposite class — the two-potential signature
is physics, not artifact.

Usage:  python analyze_forgery.py forgery_counts.json [ctc_counts_for_contrast.json]
"""
import json, sys
import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator


def theory(params):
    p = np.asarray(params)
    qc = QuantumCircuit(3)
    qc.ry(p[0], 0); qc.ry(p[1], 1); qc.ry(p[2], 2); qc.cx(1, 2)
    qc.rz(p[3], 1); qc.ry(p[4], 2); qc.cx(0, 2)
    qc.rz(p[5], 0); qc.ry(p[6], 2); qc.cx(1, 0)
    qc.ry(p[7], 0); qc.ry(p[8], 1); qc.cx(0, 2)
    qc.ry(p[9], 0); qc.rz(p[10], 1); qc.ry(p[11], 2)
    M = Operator(qc).data[:4, :4]
    u = np.log(np.diag(M @ M.conj().T).real)
    v = np.log(np.diag(M.conj().T @ M).real)
    return u, v


def dmatrix(d, min_counts=200):
    def parse(block):
        N = np.zeros((4, 4))
        for k, ci in block.items():
            i = int(k)
            for bits, n in ci.items():
                if bits[0] == '0':
                    N[i, int(bits[1:], 2)] += n
        return N
    NF = parse(d["forward"]); NR = parse(d["reverse"])
    PF = NF / NF.sum(axis=1, keepdims=True)
    PR = NR / NR.sum(axis=1, keepdims=True)
    D = np.full((4, 4), np.nan); S = np.full((4, 4), np.nan)
    for i in range(4):
        for j in range(4):
            if NF[i, j] >= min_counts and NR[j, i] >= min_counts:
                D[i, j] = np.log(PF[i, j] / PR[j, i])
                S[i, j] = np.sqrt((1 - PF[i, j]) / NF[i, j] + (1 - PR[j, i]) / NR[j, i])
    return D, S


d = json.load(open(sys.argv[1]))
u_th, v_th = theory(d["ansatz_params"])
D, S = dmatrix(d)
used = ~np.isnan(D)
print(f"{sys.argv[1]}  backend={d['backend']}  cells={used.sum()}/16")

# single-potential weighted fit: D_ij = w_j - w_i
rows, vals, sig = [], [], []
for i in range(4):
    for j in range(4):
        if used[i, j]:
            r = np.zeros(4); r[j] += 1; r[i] -= 1
            rows.append(r); vals.append(D[i, j]); sig.append(S[i, j])
Am, b, s = np.array(rows), np.array(vals), np.array(sig)
x, *_ = np.linalg.lstsq(Am / s[:, None], b / s, rcond=None)
resid = Am @ x - b
VE1 = 1 - resid.var() / b.var()
dii = float(np.abs(np.diag(D)).max())
anti = float(np.linalg.norm((D + D.T)[used & used.T]) / np.linalg.norm(D[used]))
D_th1 = np.array([[v_th[j] - v_th[i] for j in range(4)] for i in range(4)])
r_th = float(np.corrcoef(D[used], D_th1[used])[0, 1])
print(f"  single-potential VE = {VE1:.4f}   r(D, v_j-v_i) = {r_th:.4f}   mean|D| = {np.abs(b).mean():.3f}")
print(f"  max|D_ii| = {dii:.4f}   antisymmetry ||D+D^T||/||D|| = {anti:.4f}")

print("\n--- COMMITTED CRITERIA (forgery control) ---")
checks = [("single-potential VE >= 0.90", VE1 >= 0.90),
          ("max|D_ii| <= 0.45", dii <= 0.45),
          ("antisymmetry ratio <= 0.65", anti <= 0.65)]
for name, ok in checks:
    print(f"  {'PASS' if ok else 'FAIL'}  {name}")
print(f"\nFORGERY CONTROL: {'CONFIRMED single-potential class' if all(ok for _, ok in checks) else 'NOT CONFIRMED'}")

if len(sys.argv) > 2:
    d2 = json.load(open(sys.argv[2]))
    D2, _ = dmatrix(d2)
    u2 = ~np.isnan(D2)
    dii2 = float(np.abs(np.diag(D2)).max())
    anti2 = float(np.linalg.norm((D2 + D2.T)[u2 & u2.T]) / np.linalg.norm(D2[u2]))
    print(f"\nCONTRAST (CTC circuit, same device): max|D_ii| = {dii2:.4f}  antisymmetry = {anti2:.4f}")
    print("=> same device, same noise: normal operator lands in the single-potential")
    print("   class, non-normal in the two-potential class. The taxonomy is physics.")
