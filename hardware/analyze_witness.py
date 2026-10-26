"""
analyze_witness.py — analysis for PREREGISTRATION_witness.json.

Per arm: merge counts (the I arm's Pauli variants merge by their shot weights
automatically, since counts are raw), herald on ancilla=0, form D against the
appropriate reverse block (arms C2 and I reuse arm B's reverse), weighted
least-squares fit of D_ij = u_j - v_i, report VE.

Registered gates:
  G0 sanity      : VE_B >= 0.90, 16/16 cells all arms
  G1 covariant   : VE_C1 >= VE_B - 0.02  AND  r(D_C1, D_th(C')) > r(D_C1, D_th(C))
  G2 fresh coh.  : VE_B - VE_C2 >= 0.05
  G3 stochastic  : VE_B - VE_I  >= 0.02
Usage:  python analyze_witness.py witness_counts.json
"""
import json, sys
import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator

d = json.load(sys.argv[1:] and open(sys.argv[1]) or open("witness_counts.json"))
P0 = d["params"]["P0"]
EPS = d["params"]["eps_cov"]


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

N = {}   # (arm, dir) -> 4x4 heralded counts
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
    return D, 1 - resid.var() / b.var(), cells


DB, VEB, okB = arm_stats(N[("B", "F")], N[("B", "R")])
DC1, VEC1, okC1 = arm_stats(N[("C1", "F")], N[("C1", "R")])
DC2, VEC2, okC2 = arm_stats(N[("C2", "F")], N[("B", "R")])
DI, VEI, okI = arm_stats(N[("I", "F")], N[("B", "R")])

r_shift = np.corrcoef(DC1.ravel(), D_th1.ravel())[0, 1]
r_base = np.corrcoef(DC1.ravel(), D_th0.ravel())[0, 1]

print(f"{d['backend']}  layout={d.get('layout')}")
print(f"  B  : VE={VEB:.4f}")
print(f"  C1 : VE={VEC1:.4f}   r(D, D_th(C'))={r_shift:.4f} vs r(D, D_th(C))={r_base:.4f}"
      f"   max|D-D_th(C')|={np.abs(DC1-D_th1).max():.3f}")
print(f"  C2 : VE={VEC2:.4f}   (drop {VEB-VEC2:+.4f})")
print(f"  I  : VE={VEI:.4f}   (drop {VEB-VEI:+.4f})")

print("\n--- PRE-REGISTERED GATES (witness) ---")
checks = [
    ("G0 cells 16/16 all arms", okB and okC1 and okC2 and okI),
    ("G0 VE_B >= 0.90", VEB >= 0.90),
    ("G1 VE_C1 >= VE_B - 0.02", VEC1 >= VEB - 0.02),
    ("G1 r(D_C1,D_th(C')) > r(D_C1,D_th(C))", r_shift > r_base),
    ("G2 VE_B - VE_C2 >= 0.05", VEB - VEC2 >= 0.05),
    ("G3 VE_B - VE_I >= 0.02", VEB - VEI >= 0.02),
]
for name, ok in checks:
    print(f"  {'PASS' if ok else 'FAIL'}  {name}")
print(f"\nNOISE-SPECTROMETER TEST: "
      f"{'CONFIRMED' if all(ok for _, ok in checks) else 'NOT CONFIRMED'}")
