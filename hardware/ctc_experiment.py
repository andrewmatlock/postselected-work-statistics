"""
CTC fluctuation-statistics experiment — full protocol.
Runs on Aer simulator by default; set USE_HARDWARE=True for IBM Quantum.

Measures:  (1) forward + reverse conditionals p(j|i), p~(i|j)
           (2) the biased past  (accepted-run initial marginal vs thermal)
           (3) separability residual of the Crooks deviation matrix
Exports raw counts to ctc_hardware_counts.json for analysis.
"""
import json
import numpy as np
from scipy.linalg import sqrtm
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import UnitaryGate

# ---------------- configuration ----------------
USE_HARDWARE = False          # flip to True for real qubits
SHOTS_PER_INPUT = 20000
SEED = 0
# ------------------------------------------------

rng = np.random.default_rng(SEED)
ES = np.array([-1.5, -0.5, 0.5, 1.5]); beta = 1.0
pth = np.exp(-beta * ES); pth /= pth.sum()

# fixed random CTC interaction (same seed => same operator everywhere)
A = rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8))
U, _ = np.linalg.qr(A)
C = np.einsum('akbk->ab', U.reshape(4, 2, 4, 2))
K = C / np.linalg.svd(C, compute_uv=False).max()
DK = sqrtm(np.eye(4) - K.conj().T @ K)
DKd = sqrtm(np.eye(4) - K @ K.conj().T)
V = np.block([[K, DKd], [DK, -K.conj().T]])          # Halmos dilation, unitary
assert np.allclose(V.conj().T @ V, np.eye(8), atol=1e-8), "dilation not unitary"

# ---------------- backend ----------------
if USE_HARDWARE:
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
    service = QiskitRuntimeService()
    backend = service.least_busy(operational=True, simulator=False, min_num_qubits=3)
    sampler = SamplerV2(mode=backend)
    sampler.options.twirling.enable_measure = True    # readout mitigation
    backend_name = backend.name
    print(f"backend: {backend_name} ({backend.num_qubits} qubits)")
else:
    from qiskit_aer import AerSimulator
    backend = AerSimulator(); backend_name = "aer_simulator"
    print("backend: Aer simulator")


def run_direction(W, tag):
    """W: 8x8 unitary. Returns {input i: counts dict} with ancilla = qubit 2."""
    out = {}
    for i in range(4):
        qc = QuantumCircuit(3, 3)
        if i & 2: qc.x(1)
        if i & 1: qc.x(0)
        qc.append(UnitaryGate(W, label='V'), [0, 1, 2])
        qc.measure([0, 1, 2], [0, 1, 2])
        if USE_HARDWARE:
            isa = transpile(qc, backend=backend, optimization_level=3)
            job = sampler.run([isa], shots=SHOTS_PER_INPUT)
            counts = job.result()[0].data.c.get_counts()
        else:
            counts = backend.run(transpile(qc, backend), shots=SHOTS_PER_INPUT).result().get_counts()
        out[i] = {k: int(v) for k, v in counts.items()}
        print(f"  {tag} input {i}: {sum(out[i].values())} shots")
    return out


def conditionals(counts):
    """Herald on ancilla == 0 (leftmost bit in qiskit strings). Returns P, acceptance."""
    P = np.zeros((4, 4)); acc = np.zeros(4)
    for i in range(4):
        tot = sum(counts[i].values()); kept = 0
        for bits, n in counts[i].items():
            if bits[0] == '0':
                P[i, int(bits[1:], 2)] += n; kept += n
        acc[i] = kept / tot
        P[i] /= max(kept, 1)
    return P, acc


def sep_residual(D):
    rows, vals = [], []
    for i in range(4):
        for j in range(4):
            if not np.isnan(D[i, j]):
                r = np.zeros(8); r[i] = -1; r[4 + j] = 1
                rows.append(r); vals.append(D[i, j])
    A_ = np.array(rows); b = np.array(vals)
    x, *_ = np.linalg.lstsq(A_, b, rcond=None)
    return float(np.abs(A_ @ x - b).max()), float(np.abs(b).mean())


print("\n--- forward ---")
fw = run_direction(V, "F")
print("--- reverse ---")
rv = run_direction(V.conj().T, "R")

PF, accF = conditionals(fw)
PR, accR = conditionals(rv)

# theory
PF_th = np.array([[abs(C[j, i])**2 / np.vdot(C[:, i], C[:, i]).real for j in range(4)] for i in range(4)])
Cd = C.conj().T
PR_th = np.array([[abs(Cd[i, j])**2 / np.vdot(Cd[:, j], Cd[:, j]).real for i in range(4)] for j in range(4)])

print(f"\n[1] theory check: max|measured - predicted| forward = {np.abs(PF - PF_th).max():.4f}")
print(f"                                        reverse = {np.abs(PR - PR_th).max():.4f}")

# the biased past
marg = pth * accF; marg /= marg.sum()
Na = SHOTS_PER_INPUT * float(pth @ accF) * 4
z = np.abs(marg - pth) / np.sqrt(pth * (1 - pth) / Na)
print(f"\n[2] biased past: max deviation {np.abs(marg - pth).max():.4f}   max z = {z.max():.1f} sigma")
print(f"    thermal  {np.round(pth,4)}")
print(f"    accepted {np.round(marg,4)}")

# separability
D = np.full((4, 4), np.nan)
for i in range(4):
    for j in range(4):
        if PF[i, j] > 1e-6 and PR[j, i] > 1e-6:
            D[i, j] = np.log(PF[i, j] / PR[j, i])
res, mean_dev = sep_residual(D)
print(f"\n[3] separability residual = {res:.4f}   (mean|deviation| = {mean_dev:.3f})")
print("    noiseless theory ~1e-16;  generic driven channels ~0.5-1.5")

json.dump({"forward": fw, "reverse": rv, "backend": backend_name,
           "shots_per_input": SHOTS_PER_INPUT, "seed": SEED},
          open("ctc_hardware_counts.json", "w"), indent=2)
print("\nwrote ctc_hardware_counts.json")
