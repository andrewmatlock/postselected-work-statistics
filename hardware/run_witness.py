"""
run_witness.py — spin-off A: the pair-residual as a noise spectrometer
(PREREGISTRATION_witness.json). Four arms, ONE job, pinned layout, fixed seed:

  B  : baseline forward+reverse pair (8 pubs x 20k)
  C1 : inversion-COVARIANT miscalibration (ansatz param p4 += 0.3 in BOTH
       directions; the inverse circuit inverts the miscalibrated gate).
       Prediction: ledger intact, D tracks the shifted prediction D_th(C').
  C2 : FRESH coherent error (extra ry(0.35) on the ancilla before the 4th CX,
       forward only; reverse reuses arm B's). Prediction: ledger degrades.
  I  : FRESH stochastic error (exact depolarizing q=0.4 at the same site via
       a Pauli mixture {I,X,Y,Z} with shot weights {0.7,0.1,0.1,0.1}, forward
       only). Prediction: ledger degrades by the precomputed amount.

Usage:  python run_witness.py [hardware [backend]]
Output: witness_counts.json
"""
import json, sys, time
import numpy as np
from qiskit import QuantumCircuit, transpile

MODE = sys.argv[1] if len(sys.argv) > 1 else "sim"
BACKEND_NAME = sys.argv[2] if len(sys.argv) > 2 else None
LAYOUT = [53, 54, 55]
SEED = 42
SHOTS = 20000
EPS_COV = 0.3
EPS_FRESH = 0.35
Q_DEPOL = 0.4
OUT = "witness_counts.json"

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


def wrap(body, i):
    qc = QuantumCircuit(3, 3)
    if i & 2: qc.x(1)
    if i & 1: qc.x(0)
    qc.compose(body, inplace=True)
    qc.measure([0, 1, 2], [0, 1, 2])
    return qc


pubs, plan = [], []   # plan rows: (arm, direction, input, pauli_or_None, shots)

def add(body, arm, direction, i, shots, pauli=None):
    pubs.append(wrap(body, i))
    plan.append({"arm": arm, "dir": direction, "input": i,
                 "pauli": pauli, "shots": shots})

for i in range(4):   # B
    add(base(P0), "B", "F", i, SHOTS)
for i in range(4):
    add(base(P0).inverse(), "B", "R", i, SHOTS)
P1 = list(P0); P1[4] += EPS_COV
for i in range(4):   # C1
    add(base(P1), "C1", "F", i, SHOTS)
for i in range(4):
    add(base(P1).inverse(), "C1", "R", i, SHOTS)
for i in range(4):   # C2 (forward only; reverse = arm B's)
    add(base(P0, fresh_ry=EPS_FRESH), "C2", "F", i, SHOTS)
w0 = 1 - 3 * Q_DEPOL / 4
weights = {"i": w0, "x": Q_DEPOL/4, "y": Q_DEPOL/4, "z": Q_DEPOL/4}
for i in range(4):   # I (forward only)
    for pname, w in weights.items():
        pl = None if pname == "i" else pname
        add(base(P0, pauli=pl), "I", "F", i, int(round(SHOTS * w)), pauli=pname)

if MODE == "hardware":
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
    service = QiskitRuntimeService()
    backend = (service.backend(BACKEND_NAME) if BACKEND_NAME else
               service.backend("ibm_marrakesh"))
    name = backend.name
    isa = transpile(pubs, backend=backend, optimization_level=3,
                    initial_layout=LAYOUT, seed_transpiler=SEED)
    sampler = SamplerV2(mode=backend)
    sampler.options.twirling.enable_measure = True
    sampler.options.dynamical_decoupling.enable = True
    t0 = time.time()
    job = sampler.run([(c, None, row["shots"]) for c, row in zip(isa, plan)])
    print(f"job id: {job.job_id()}  ({len(pubs)} pubs, layout {LAYOUT}) — record this ID")
    result = job.result()
    print(f"completed in {time.time()-t0:.0f}s wall")
    counts = [r.data.c.get_counts() for r in result]
else:
    from qiskit_aer import AerSimulator
    backend = AerSimulator(); name = "aer_simulator"
    isa = transpile(pubs, backend=backend, seed_transpiler=SEED)
    counts = [backend.run(c, shots=row["shots"]).result().get_counts()
              for c, row in zip(isa, plan)]
    print(f"backend: {name} (dry run, {len(pubs)} pubs)")

json.dump({"plan": plan, "counts": counts, "backend": name,
           "layout": LAYOUT, "seed_transpiler": SEED,
           "params": {"P0": P0, "eps_cov": EPS_COV, "eps_fresh": EPS_FRESH,
                      "q_depol": Q_DEPOL}},
          open(OUT, "w"), indent=1)
print(f"wrote {OUT} — analyze with:  python analyze_witness.py {OUT}")
