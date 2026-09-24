"""
analyze_paired.py — analysis for PREREGISTRATION_paired.json.

Primary (drift-immune, paired within one job):
  R_D  = max|D_ii|(ctc) / max|D_ii|(control)      gate: >= 3.0
  R_A  = antisym(ctc)   / antisym(control)        gate: >= 1.5
CTC-arm structure (same gates as v2):
  VE >= 0.90 and r(D, D_th) >= 0.95 on 16/16 cells
Control secondary (reported, non-gating): single-potential VE, max|D_ii|,
antisymmetry, r against predicted v_j - v_i.

Usage:  python analyze_paired.py paired_counts.json
"""
import json, sys
import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator

GATES = dict(RD=3.0, RA=1.5, ve=0.90, r=0.95)


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


def dmatrix(block, min_counts=200):
    def parse(b):
        N = np.zeros((4, 4))
        for k, ci in b.items():
            for bits, n in ci.items():
                if bits[0] == '0':
                    N[int(k), int(bits[1:], 2)] += n
        return N
    NF, NR = parse(block["forward"]), parse(block["reverse"])
    PF = NF / NF.sum(axis=1, keepdims=True)
    PR = NR / NR.sum(axis=1, keepdims=True)
    D = np.full((4, 4), np.nan); S = np.full((4, 4), np.nan)
    for i in range(4):
        for j in range(4):
            if NF[i, j] >= min_counts and NR[j, i] >= min_counts:
                D[i, j] = np.log(PF[i, j] / PR[j, i])
                S[i, j] = np.sqrt((1 - PF[i, j]) / NF[i, j] + (1 - PR[j, i]) / NR[j, i])
    return D, S


def wfit(D, S, npar):
    """Weighted lstsq to D_ij = u_j - v_i (npar=8) or w_j - w_i (npar=4)."""
    rows, vals, sig = [], [], []
    for i in range(4):
        for j in range(4):
            if not np.isnan(D[i, j]):
                if npar == 8:
                    r = np.zeros(8); r[i] = -1; r[4 + j] = 1
                else:
                    r = np.zeros(4); r[j] += 1; r[i] -= 1
                rows.append(r); vals.append(D[i, j]); sig.append(S[i, j])
    Am, b, s = np.array(rows), np.array(vals), np.array(sig)
    x, *_ = np.linalg.lstsq(Am / s[:, None], b / s, rcond=None)
    resid = Am @ x - b
    return 1 - resid.var() / b.var(), x


def antisym(D):
    m = ~np.isnan(D) & ~np.isnan(D.T)
    return float(np.linalg.norm((D + D.T)[m]) / np.linalg.norm(D[~np.isnan(D)]))


d = json.load(open(sys.argv[1]))
u_c, v_c = theory(d["ctc_params"])
u_g, v_g = theory(d["ctrl_params"])
D_th_ctc = np.array([[u_c[j] - v_c[i] for j in range(4)] for i in range(4)])
D_th_ctrl = np.array([[v_g[j] - v_g[i] for j in range(4)] for i in range(4)])

Dc, Sc = dmatrix(d["ctc"])
Dg, Sg = dmatrix(d["control"])
uc = ~np.isnan(Dc); ug = ~np.isnan(Dg)

VEc, _ = wfit(Dc, Sc, 8)
rc = float(np.corrcoef(Dc[uc], D_th_ctc[uc])[0, 1])
VEg1, _ = wfit(Dg, Sg, 4)
rg = float(np.corrcoef(Dg[ug], D_th_ctrl[ug])[0, 1])

dii_c = float(np.abs(np.diag(Dc)).max())
dii_g = float(np.abs(np.diag(Dg)).max())
a_c, a_g = antisym(Dc), antisym(Dg)
RD = dii_c / dii_g if dii_g > 0 else np.inf
RA = a_c / a_g if a_g > 0 else np.inf

print(f"{sys.argv[1]}  backend={d['backend']}  layout={d.get('layout')}  "
      f"cells: ctc {uc.sum()}/16, control {ug.sum()}/16")
print(f"  CTC arm:     VE={VEc:.4f}  r(D,D_th)={rc:.4f}  max|D_ii|={dii_c:.4f}  antisym={a_c:.4f}")
print(f"  control arm: VE1={VEg1:.4f} r(D,v_j-v_i)={rg:.4f}  max|D_ii|={dii_g:.4f}  antisym={a_g:.4f}")
print(f"  PAIRED RATIOS:  R_D = {RD:.2f}   R_A = {RA:.2f}")

print("\n--- PRE-REGISTERED CRITERIA (paired) ---")
checks = [
    (f"cells 16/16 both arms", uc.sum() == 16 and ug.sum() == 16),
    (f"R_D >= {GATES['RD']}", RD >= GATES['RD']),
    (f"R_A >= {GATES['RA']}", RA >= GATES['RA']),
    (f"CTC VE >= {GATES['ve']}", VEc >= GATES['ve']),
    (f"CTC r >= {GATES['r']}", rc >= GATES['r']),
]
for name, ok in checks:
    print(f"  {'PASS' if ok else 'FAIL'}  {name}")
verdict = all(ok for _, ok in checks)
print(f"\nPAIRED TAXONOMY TEST: {'CONFIRMED' if verdict else 'NOT CONFIRMED'}")
print("(control-arm class metrics above are reported as secondary, non-gating)")
