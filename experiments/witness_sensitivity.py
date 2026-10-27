"""
witness_sensitivity.py — exact (infinite-shot) sensitivity curves for the
noise-spectrometer witness (spin-off A), computed AFTER registration and
BEFORE hardware data arrive; curves are theory, not gates, and do not touch
the registered criteria.

Three curves, same conventions as run_witness.py / analyze_witness.py
(herald on ancilla qubit 2 = 0; system outcome from qubits 0,1):

  1. covariant blindness : VE_C1(eps) for eps in [0, 1.5]   -> expected == 1
  2. fresh coherent      : VE drop vs eps_fresh in [0, 0.7]
  3. stochastic          : VE drop vs depolarizing q in [0, 0.8]

Output: witness_sensitivity.json next to this script.
"""
import json
import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector

P0 = [2.774141, 1.438552, 5.364134, 1.634255, 2.682839, 4.893008,
      1.171051, 1.459015, 1.474151, 0.031240, 3.189096, 2.607931]


def base(p, fresh_ry=None, pauli=None):
    qc = QuantumCircuit(3)
    qc.ry(p[0], 0); qc.ry(p[1], 1); qc.ry(p[2], 2)
    qc.cx(1, 2)
    qc.rz(p[3], 1); qc.ry(p[4], 2)
    qc.cx(0, 2)
    qc.rz(p[5], 0); qc.ry(p[6], 2)
    qc.cx(1, 0)
    qc.ry(p[7], 0); qc.ry(p[8], 1)
    if fresh_ry is not None:
        qc.ry(fresh_ry, 2)
    if pauli == 'x': qc.x(2)
    elif pauli == 'y': qc.y(2)
    elif pauli == 'z': qc.z(2)
    qc.cx(0, 2)
    qc.ry(p[9], 0); qc.rz(p[10], 1); qc.ry(p[11], 2)
    return qc


def heralded_P(bodies_with_weights):
    """4x4 heralded conditionals P[i,j] for a weighted mixture of circuit bodies."""
    P = np.zeros((4, 4))
    for i in range(4):
        probs = np.zeros(8)
        for body, w in bodies_with_weights:
            qc = QuantumCircuit(3)
            if i & 2: qc.x(1)
            if i & 1: qc.x(0)
            qc.compose(body, inplace=True)
            probs += w * np.abs(Statevector(qc).data) ** 2
        # statevector index q2 q1 q0; herald q2 = 0 -> indices 0..3 = (q1 q0)
        accepted = probs[:4]
        P[i] = accepted / accepted.sum()
    return P


def VE(PF, PR):
    D = np.log(PF / PR.T)
    rows, vals = [], []
    for i in range(4):
        for j in range(4):
            r = np.zeros(8); r[i] = -1; r[4 + j] = 1
            rows.append(r); vals.append(D[i, j])
    A, b = np.array(rows), np.array(vals)
    x, *_ = np.linalg.lstsq(A, b, rcond=None)
    resid = A @ x - b
    return 1 - resid.var() / b.var()


PR_B = heralded_P([(base(P0).inverse(), 1.0)])
PF_B = heralded_P([(base(P0), 1.0)])
VE_B = VE(PF_B, PR_B)

out = {"P0": P0, "VE_B_ideal": VE_B, "curves": {}}

# 1. covariant blindness: shift param 4 in BOTH directions, own reverse
cov = []
for eps in np.arange(0.0, 1.51, 0.1):
    P1 = list(P0); P1[4] += eps
    v = VE(heralded_P([(base(P1), 1.0)]),
           heralded_P([(base(P1).inverse(), 1.0)]))
    cov.append({"eps": round(float(eps), 2), "VE": float(v)})
out["curves"]["covariant"] = cov

# 2. fresh coherent error, forward only, reverse = arm B's
fresh = []
for eps in np.arange(0.0, 0.71, 0.05):
    v = VE(heralded_P([(base(P0, fresh_ry=float(eps)), 1.0)]), PR_B)
    fresh.append({"eps": round(float(eps), 2), "VE": float(v),
                  "drop": float(VE_B - v)})
out["curves"]["fresh_coherent"] = fresh

# 3. depolarizing via exact Pauli mixture, forward only, reverse = arm B's
stoch = []
for q in np.arange(0.0, 0.81, 0.05):
    mix = [(base(P0), 1 - 3 * q / 4)] + \
          [(base(P0, pauli=p), q / 4) for p in 'xyz']
    v = VE(heralded_P(mix), PR_B)
    stoch.append({"q": round(float(q), 2), "VE": float(v),
                  "drop": float(VE_B - v)})
out["curves"]["depolarizing"] = stoch

json.dump(out, open(__file__.replace(".py", ".json"), "w"), indent=1)

print(f"VE_B ideal = {VE_B:.6f}")
print("\ncovariant (expect VE == 1 at every eps):")
print("  max deviation from 1:",
      max(abs(1 - c["VE"]) for c in cov))
print("\nfresh coherent eps -> VE drop:")
for c in fresh:
    print(f"  {c['eps']:4.2f}  {c['drop']:+.4f}")
print("\ndepolarizing q -> VE drop:")
for c in stoch:
    print(f"  {c['q']:4.2f}  {c['drop']:+.4f}")
# registered points must reproduce the registration's ideal_math numbers
f35 = [c for c in fresh if c["eps"] == 0.35][0]
q40 = [c for c in stoch if c["q"] == 0.4][0]
print(f"\nregistered points: C2 drop {f35['drop']:.3f} (registered 0.105), "
      f"I drop {q40['drop']:.3f} (registered 0.048)")
