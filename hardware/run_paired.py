"""
run_paired.py — the paired taxonomy experiment (PREREGISTRATION_paired.json).

Both arms in ONE job: 8 circuits of the CTC (non-normal) ansatz + 8 of the
high-contrast normal control, transpiled with the same pinned initial layout
and the same transpiler seed. One job = one calibration snapshot: common
drift cancels in the paired ratio statistics, which is what killed the two
absolute-gated control attempts.

Usage:
    python run_paired.py                    # Aer dry run
    python run_paired.py hardware [backend] [q0 q1 q2]
Output: paired_counts.json  {"ctc": {...}, "control": {...}}
"""
import json, sys, time
import numpy as np
from qiskit import QuantumCircuit, transpile

MODE = sys.argv[1] if len(sys.argv) > 1 else "sim"
BACKEND_NAME = sys.argv[2] if len(sys.argv) > 2 else None
LAYOUT = [int(x) for x in sys.argv[3:6]] if len(sys.argv) > 5 else [53, 54, 55]
SHOTS = 20000
SEED = 42
OUT = "paired_counts.json"

CTC_PARAMS = [2.774141, 1.438552, 5.364134, 1.634255, 2.682839, 4.893008,
              1.171051, 1.459015, 1.474151, 0.031240, 3.189096, 2.607931]
CTRL_PARAMS = [4.904035, 0.646883, 0.516469, 1.555524, 4.679263, 5.049346,
               3.141593, 5.493009, 0.646883, 3.647830, 3.141593, 3.633187]


def ansatz(p):
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


def circuits(params):
    out = []
    base = ansatz(params)
    for inverse in (False, True):
        body = base.inverse() if inverse else base
        for i in range(4):
            qc = QuantumCircuit(3, 3)
            if i & 2: qc.x(1)
            if i & 1: qc.x(0)
            qc.compose(body, inplace=True)
            qc.measure([0, 1, 2], [0, 1, 2])
            out.append(qc)
    return out


all_circuits = circuits(CTC_PARAMS) + circuits(CTRL_PARAMS)   # 16 pubs

if MODE == "hardware":
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
    service = QiskitRuntimeService()
    backend = (service.backend(BACKEND_NAME) if BACKEND_NAME else
               service.least_busy(operational=True, simulator=False, min_num_qubits=3))
    name = backend.name
    isa = transpile(all_circuits, backend=backend, optimization_level=3,
                    initial_layout=LAYOUT, seed_transpiler=SEED)
    for tag, k in (("ctc", 0), ("ctrl", 8)):
        d = isa[k]
        print(f"{tag}: transpiled depth={d.depth()}  2q="
              f"{sum(v for g, v in d.count_ops().items() if g in ('cz', 'cx', 'ecr'))}")
    sampler = SamplerV2(mode=backend)
    sampler.options.twirling.enable_measure = True
    sampler.options.dynamical_decoupling.enable = True
    t0 = time.time()
    job = sampler.run(isa, shots=SHOTS)
    print(f"job id: {job.job_id()}  (16 pubs x {SHOTS}, layout {LAYOUT}) — record this ID")
    result = job.result()
    print(f"completed in {time.time()-t0:.0f}s wall")
    counts = [r.data.c.get_counts() for r in result]
else:
    from qiskit_aer import AerSimulator
    backend = AerSimulator(); name = "aer_simulator"
    isa = transpile(all_circuits, backend=backend, seed_transpiler=SEED)
    counts = [backend.run(c, shots=SHOTS).result().get_counts() for c in isa]
    print(f"backend: {name} (dry run)")

def pack(offset):
    return {"forward": {str(i): counts[offset + i] for i in range(4)},
            "reverse": {str(i): counts[offset + 4 + i] for i in range(4)}}

json.dump({"ctc": pack(0), "control": pack(8), "backend": name,
           "shots_per_input": SHOTS, "layout": LAYOUT, "seed_transpiler": SEED,
           "ctc_params": CTC_PARAMS, "ctrl_params": CTRL_PARAMS},
          open(OUT, "w"), indent=1)
print(f"wrote {OUT} — analyze with:  python analyze_paired.py {OUT}")
