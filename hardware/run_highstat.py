"""
run_highstat.py — high-statistics hardware run for the CTC fluctuation experiment.

Batches shots across multiple jobs to stay inside per-job limits, accumulates
counts, and writes ctc_highstat_counts.json in the same format as before.

Usage:
    python run_highstat.py                  # 100k shots/input on least-busy device
    python run_highstat.py 200000           # custom total shots per input
    python run_highstat.py 100000 ibm_torino  # pin a specific backend

Runtime: expect queue waits. 8 circuits x (total/batch) jobs.
"""
import json, sys, time
import numpy as np
from scipy.linalg import sqrtm
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import UnitaryGate
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2

TOTAL_SHOTS = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
BACKEND_NAME = sys.argv[2] if len(sys.argv) > 2 else None
BATCH = 20000                      # shots per job
SEED = 0
OUT = "ctc_highstat_counts.json"

# --- operator (identical to previous runs via fixed seed) ---
rng = np.random.default_rng(SEED)
A = rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8))
U, _ = np.linalg.qr(A)
C = np.einsum('akbk->ab', U.reshape(4, 2, 4, 2))
K = C / np.linalg.svd(C, compute_uv=False).max()
DK = sqrtm(np.eye(4) - K.conj().T @ K)
DKd = sqrtm(np.eye(4) - K @ K.conj().T)
V = np.block([[K, DKd], [DK, -K.conj().T]])
assert np.allclose(V.conj().T @ V, np.eye(8), atol=1e-8)

service = QiskitRuntimeService()
if BACKEND_NAME:
    backend = service.backend(BACKEND_NAME)
else:
    backend = service.least_busy(operational=True, simulator=False, min_num_qubits=3)
print(f"backend: {backend.name} ({backend.num_qubits} qubits)")
print(f"target: {TOTAL_SHOTS} shots/input  ({-(-TOTAL_SHOTS//BATCH)} batches x 8 circuits)")

sampler = SamplerV2(mode=backend)
sampler.options.twirling.enable_measure = True
sampler.options.dynamical_decoupling.enable = True     # idle-qubit error suppression


def build(W, i):
    qc = QuantumCircuit(3, 3)
    if i & 2: qc.x(1)
    if i & 1: qc.x(0)
    qc.append(UnitaryGate(W, label='V'), [0, 1, 2])
    qc.measure([0, 1, 2], [0, 1, 2])
    return transpile(qc, backend=backend, optimization_level=3)


def accumulate(W, tag):
    out = {}
    for i in range(4):
        isa = build(W, i)
        acc = {}
        remaining = TOTAL_SHOTS
        while remaining > 0:
            n = min(BATCH, remaining)
            job = sampler.run([isa], shots=n)
            print(f"   {tag}{i}: job {job.job_id()} ({n} shots) ...", flush=True)
            counts = job.result()[0].data.c.get_counts()
            for k, v in counts.items():
                acc[k] = acc.get(k, 0) + int(v)
            remaining -= n
        out[i] = acc
        print(f"   {tag}{i}: total {sum(acc.values())} shots")
    return out


t0 = time.time()
print("\n--- forward ---")
fw = accumulate(V, "F")
print("--- reverse ---")
rv = accumulate(V.conj().T, "R")

json.dump({"forward": fw, "reverse": rv, "backend": backend.name,
           "shots_per_input": TOTAL_SHOTS, "seed": SEED},
          open(OUT, "w"), indent=2)
print(f"\nwrote {OUT}  ({time.time()-t0:.0f}s elapsed)")
print(f"analyze with:  python analyze_weighted.py {OUT}")
