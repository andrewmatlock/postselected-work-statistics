"""
run_prereg.py — executes the pre-registered experiment in notes/PREREGISTRATION.json:
shallow 4-CX circuit, full 16-cell deviation matrix, 20k shots/input,
ALL EIGHT CIRCUITS SUBMITTED AS ONE JOB (per-job overhead dominates the
Open Plan budget, so one job with 8 pubs is ~1-2 QPU minutes total).

Analysis is done by analyze_shallow.py, unmodified, per the pre-registration.

Usage:
    python run_prereg.py                  # noiseless simulator dry run
    python run_prereg.py hardware         # least-busy Heron backend
    python run_prereg.py hardware ibm_x   # pin a backend
Output: ctc_prereg_counts.json  (same schema as ctc_shallow_counts.json)
"""
import json, sys, time
import numpy as np
from qiskit import QuantumCircuit, transpile

MODE = sys.argv[1] if len(sys.argv) > 1 else "sim"
BACKEND_NAME = sys.argv[2] if len(sys.argv) > 2 else None
SHOTS = 20000                      # fixed by PREREGISTRATION.json
OUT = "ctc_prereg_counts.json"

PARAMS = np.array(json.load(open(__file__.rsplit("/", 2)[0] +
                                 "/notes/PREREGISTRATION.json"))["ansatz_params"])


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


all_circuits = circuits(False) + circuits(True)   # F0..F3, R0..R3

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
    job = sampler.run(isa, shots=SHOTS)             # ONE job, 8 pubs
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

fw = {str(i): counts[i] for i in range(4)}
rv = {str(i): counts[4 + i] for i in range(4)}
json.dump({"forward": fw, "reverse": rv, "backend": name,
           "shots_per_input": SHOTS, "ansatz_params": PARAMS.tolist(),
           "circuit": "shallow_4cx", "preregistered": True},
          open(OUT, "w"), indent=2)
print(f"wrote {OUT} — now run:  python analyze_shallow.py {OUT}")
