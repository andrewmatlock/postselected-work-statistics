"""
analyze_witness_v2.py — analysis for PREREGISTRATION_witness_v2.json.

Identical ingestion and VE machinery to analyze_witness.py; the stochastic
gate is replaced by the scale-aware D-variance ratio K = Var(D_arm)/Var(D_B),
because run 1 on both devices showed VE (a scale-invariant ratio) is blind to
depolarizing compression once the baseline residual is non-negligible, while
K matched the exact parameter-free theory on both devices (theory 0.255;
measured 0.251 fez, 0.280 marrakesh).

The coherent gate is likewise scale-aware: R_C2 = MSresid(C2)/MSresid(B)
(ratio of weighted mean-square residuals, drift-immune within the job),
because the v1 VE-drop gate shares VE's ratio blindness and passed marrakesh
by only 0.011; R_C2 measured 3.34 (marrakesh) and 3.15 (fez) in v1 vs a
floor of 2. The VE drop is reported but not gated.

Registered gates (v2):
  G0 sanity      : VE_B >= 0.90, 16/16 cells all arms
  G1 covariant   : VE_C1 >= VE_B - 0.02  AND  r(D_C1, D_th(C')) > r(D_C1, D_th(C))
  G2 fresh coh.  : R_C2 >= 2  AND  K_C2 >= 0.8
  G3 stochastic  : K_I <= 0.5  AND  |K_I - 0.255| <= 0.10
Usage:  python analyze_witness_v2.py <counts.json>
"""
import json, sys
import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator

d = json.load(open(sys.argv[1]))
P0 = d["params"]["P0"]
EPS = d["params"]["eps_cov"]
K_I_THEORY = 0.255


def theory_D(p):
    qc = QuantumCircuit(3)
    qc.ry(p[0], 0); qc.ry(p[1], 1); qc.ry(p[2], 2); qc.cx(1, 2)
    qc.rz(p[3], 1); qc.ry(p[4], 2); qc.cx(0, 2)
    qc.rz(p[5], 0); qc.ry(p[6], 2); qc.cx(1, 0)
    qc.ry(p[7], 0); qc.ry(p[8], 1); qc.cx(0, 2)
    qc.ry(p[9], 0); qc.rz(p[10], 1); qc.ry(p[11], 2)
    M = Operator(qc).data[:4, :4]
    u = np.log(np.diag(M @ M.conj().T).real)
    v = np.log(np.diag(M.conj().T @ M).real)
    return np.array([[u[j] - v[i] for j in range(4)] for i in range(4)])


D_th0 = theory_D(P0)
P1 = list(P0); P1[4] += EPS
D_th1 = theory_D(P1)

N = {}
for row, c in zip(d["plan"], d["counts"]):
    key = (row["arm"], row["dir"])
    N.setdefault(key, np.zeros((4, 4)))
    for bits, n in c.items():
        if bits[0] == '0':
            N[key][row["input"], int(bits[1:], 2)] += n


def arm_stats(NF, NR):
    PF = NF / NF.sum(axis=1, keepdims=True)
    PR = NR / NR.sum(axis=1, keepdims=True)
    D = np.log(PF / PR.T)
    S = np.sqrt((1 - PF) / NF + (1 - PR.T) / NR.T)
    rows, vals, sig = [], [], []
    for i in range(4):
        for j in range(4):
            r = np.zeros(8); r[i] = -1; r[4 + j] = 1
            rows.append(r); vals.append(D[i, j]); sig.append(S[i, j])
    Am, b, s = np.array(rows), np.array(vals), np.array(sig)
    x, *_ = np.linalg.lstsq(Am / s[:, None], b / s, rcond=None)
    resid = Am @ x - b
    cells = int((NF >= 200).sum() == 16 and (NR >= 200).sum() == 16)
    return D, 1 - resid.var() / b.var(), cells, (resid ** 2).mean()


DB, VEB, okB, MSB = arm_stats(N[("B", "F")], N[("B", "R")])
DC1, VEC1, okC1, _ = arm_stats(N[("C1", "F")], N[("C1", "R")])
DC2, VEC2, okC2, MSC2 = arm_stats(N[("C2", "F")], N[("B", "R")])
DI, VEI, okI, _ = arm_stats(N[("I", "F")], N[("B", "R")])

RC2 = MSC2 / MSB

KC1 = DC1.var() / DB.var()
KC2 = DC2.var() / DB.var()
KI = DI.var() / DB.var()

r_shift = np.corrcoef(DC1.ravel(), D_th1.ravel())[0, 1]
r_base = np.corrcoef(DC1.ravel(), D_th0.ravel())[0, 1]

print(f"{d['backend']}  layout={d.get('layout')}")
print(f"  B  : VE={VEB:.4f}")
print(f"  C1 : VE={VEC1:.4f}   K={KC1:.3f}   r(D,D_th(C'))={r_shift:.4f} vs r(D,D_th(C))={r_base:.4f}")
print(f"  C2 : VE={VEC2:.4f}   K={KC2:.3f}   R_C2={RC2:.3f}   (VE drop {VEB-VEC2:+.4f}, reported not gated)")
print(f"  I  : VE={VEI:.4f}   K={KI:.3f}   (theory K {K_I_THEORY})")

print("\n--- PRE-REGISTERED GATES (witness v2) ---")
checks = [
    ("G0 cells 16/16 all arms", okB and okC1 and okC2 and okI),
    ("G0 VE_B >= 0.90", VEB >= 0.90),
    ("G1 VE_C1 >= VE_B - 0.02", VEC1 >= VEB - 0.02),
    ("G1 r(D_C1,D_th(C')) > r(D_C1,D_th(C))", r_shift > r_base),
    ("G2 R_C2 >= 2", RC2 >= 2),
    ("G2 K_C2 >= 0.8", KC2 >= 0.8),
    ("G3 K_I <= 0.5", KI <= 0.5),
    ("G3 |K_I - 0.255| <= 0.10", abs(KI - K_I_THEORY) <= 0.10),
]
for name, ok in checks:
    print(f"  {'PASS' if ok else 'FAIL'}  {name}")
print(f"\nNOISE-SPECTROMETER TEST v2: "
      f"{'CONFIRMED' if all(ok for _, ok in checks) else 'NOT CONFIRMED'}")
