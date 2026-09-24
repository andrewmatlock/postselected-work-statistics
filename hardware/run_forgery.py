"""
run_forgery.py — the forgery-class control experiment (note Sec. 4 taxonomy).

Same 4-CX ansatz family and protocol as the CTC circuit, but with parameters
chosen so the heralded block M is NORMAL (M M† = M† M to ~1e-12). A normal
operator sits in the single-potential (generalized-detailed-balance / forgery)
class: D_ij = v_j − v_i, antisymmetric, D_ii = 0 — despite violations as large
as the CTC circuit's. On the same device this is the negative control for the
two-potential signature: if hardware noise generically produced D_ii ≠ 0, it
would show here too.

Usage:
    python run_forgery.py                  # Aer dry run
    python run_forgery.py hardware [backend]
Output: forgery_counts.json  (analyze with analyze_forgery.py)
"""
import json, sys, time
import numpy as np
from qiskit import QuantumCircuit, transpile

MODE = sys.argv[1] if len(sys.argv) > 1 else "sim"
BACKEND_NAME = sys.argv[2] if len(sys.argv) > 2 else None
SHOTS = 20000
OUT = "forgery_counts.json"

PARAMS = np.array([0.664318, 4.528296, 2.612652, 1.285941, 2.099737, 3.013701,
                   1.345209, 2.715569, 4.528296, 6.044891, 3.141593, 4.712389])


def ansatz():
    p = PARAMS
    qc = QuantumCircuit(3)
    qc.ry(p[0], 0); qc.ry(p[1], 1); qc.ry(p[2], 2)
    qc.cx(1, 2)
    qc.rz(p[3], 1); qc.ry(p[4], 2)
    qc.cx(0, 2)
    qc.rz(p[5], 0); qc.ry(p[6], 2)
    qc.cx(1, 0)
    qc.ry(p[7], 0); qc.ry(p[8], 1)
    qc.cx(0, 2)
    qc.ry(p[9], 0); qc.rz(p[10], 1); qc.ry(p[11], 2)
    return qc


def circuits(inverse):
    base = ansatz()
    body = base.inverse() if inverse else base
    out = []
    for i in range(4):
        qc = QuantumCircuit(3, 3)
        if i & 2: qc.x(1)
        if i & 1: qc.x(0)
        qc.compose(body, inplace=True)
        qc.measure([0, 1, 2], [0, 1, 2])
        out.append(qc)
    return out


all_circuits = circuits(False) + circuits(True)

if MODE == "hardware":
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
    service = QiskitRuntimeService()
    backend = (service.backend(BACKEND_NAME) if BACKEND_NAME else
               service.least_busy(operational=True, simulator=False, min_num_qubits=3))
    name = backend.name
    isa = transpile(all_circuits, backend=backend, optimization_level=3)
    d0 = isa[0]
    print(f"backend: {name}   transpiled depth={d0.depth()}  2q gates="
          f"{sum(v for k, v in d0.count_ops().items() if k in ('cz', 'cx', 'ecr'))}")
    sampler = SamplerV2(mode=backend)
    sampler.options.twirling.enable_measure = True
    sampler.options.dynamical_decoupling.enable = True
    t0 = time.time()
    job = sampler.run(isa, shots=SHOTS)
    print(f"job id: {job.job_id()}  (8 pubs x {SHOTS} shots) — record this ID")
    result = job.result()
    print(f"completed in {time.time()-t0:.0f}s wall")
    counts = [r.data.c.get_counts() for r in result]
else:
    from qiskit_aer import AerSimulator
    backend = AerSimulator(); name = "aer_simulator"
    isa = transpile(all_circuits, backend=backend, optimization_level=3)
    counts = [backend.run(c, shots=SHOTS).result().get_counts() for c in isa]
    print(f"backend: {name} (dry run)")

json.dump({"forward": {str(i): counts[i] for i in range(4)},
           "reverse": {str(i): counts[4 + i] for i in range(4)},
           "backend": name, "shots_per_input": SHOTS,
           "ansatz_params": PARAMS.tolist(), "circuit": "forgery_normal_4cx"},
          open(OUT, "w"), indent=2)
print(f"wrote {OUT} — analyze with:  python analyze_forgery.py {OUT}")
