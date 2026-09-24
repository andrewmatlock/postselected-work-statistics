"""
ctc_shallow.py — noise-robust redesign of the CTC fluctuation experiment.

The theorem D_ij = u_j - v_i holds for ANY operator, so instead of choosing a
random operator and dilating it into a deep circuit (~20+ CX, ~15% depolarizing
on current hardware), we choose a SHALLOW circuit (4 CX) and read the effective
operator off its heralded sub-block. Same physics, ~25x better signal-to-noise.

The ansatz was selected for: high herald acceptance (~0.4-0.5), no small
transition probabilities (min ~0.06, since depolarizing noise damages small
probabilities most), and strong two-potential structure (max|D_ii| ~ 1.2, the
signature that excludes every reference-equilibrium explanation).

Usage:
    python ctc_shallow.py                    # simulator, 20k shots/input
    python ctc_shallow.py 100000 hardware    # 100k shots/input on IBM hardware
    python ctc_shallow.py 100000 hardware ibm_torino   # pin a backend
"""
import json, sys, time
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import Operator

SHOTS = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
MODE = sys.argv[2] if len(sys.argv) > 2 else "sim"
BACKEND_NAME = sys.argv[3] if len(sys.argv) > 3 else None
BATCH = 20000
OUT = "ctc_shallow_counts.json"

PARAMS = np.array([2.774141, 1.438552, 5.364134, 1.634255, 2.682839, 4.893008,
                   1.171051, 1.459015, 1.474151, 0.031240, 3.189096, 2.607931])


def ansatz():
    """4-CX shallow circuit. Qubit 2 = ancilla (herald)."""
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


base = ansatz()
Umat = Operator(base).data
M = Umat[:4, :4]                     # ancilla-heralded block (qubit 2 = 0 in and out)
u = np.log(np.diag(M @ M.conj().T).real)      # ||M†|j>||²
v = np.log(np.diag(M.conj().T @ M).real)      # ||M|i>||²
D_th = np.array([[u[j] - v[i] for j in range(4)] for i in range(4)])
PF_th = np.array([[abs(M[j, i])**2 / np.exp(v[i]) for j in range(4)] for i in range(4)])
PR_th = np.array([[abs(M[a, b])**2 / np.exp(u[a]) for b in range(4)] for a in range(4)])
acc_th = np.exp(v)

ES = np.array([-1.5, -0.5, 0.5, 1.5]); beta = 1.0
pth = np.exp(-beta * ES); pth /= pth.sum()

print(f"shallow ansatz: {base.count_ops()}")
print(f"mean|D| = {np.abs(D_th).mean():.3f}   max|D_ii| = {max(abs(D_th[k,k]) for k in range(4)):.3f}")
print(f"acceptance = {np.round(acc_th,3)}   min cell prob = {PF_th.min():.3f}")

# ---------------- backend ----------------
if MODE == "hardware":
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
    service = QiskitRuntimeService()
    backend = (service.backend(BACKEND_NAME) if BACKEND_NAME else
               service.least_busy(operational=True, simulator=False, min_num_qubits=3))
    sampler = SamplerV2(mode=backend)
    sampler.options.twirling.enable_measure = True
    sampler.options.dynamical_decoupling.enable = True
    name = backend.name
else:
    from qiskit_aer import AerSimulator
    backend = AerSimulator(); name = "aer_simulator"
print(f"backend: {name}\n")


def circuits(inverse):
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


def run(inverse, tag):
    res = {}
    for i, qc in enumerate(circuits(inverse)):
        acc = {}; remaining = SHOTS
        isa = transpile(qc, backend=backend, optimization_level=3)
        if i == 0:
            print(f"  transpiled depth={isa.depth()}  2q gates="
                  f"{sum(v for k,v in isa.count_ops().items() if k in ('cz','cx','ecr'))}")
        while remaining > 0:
            n = min(BATCH, remaining)
            if MODE == "hardware":
                job = sampler.run([isa], shots=n)
                print(f"   {tag}{i}: {job.job_id()} ({n})", flush=True)
                c = job.result()[0].data.c.get_counts()
            else:
                c = backend.run(isa, shots=n).result().get_counts()
            for k, val in c.items():
                acc[k] = acc.get(k, 0) + int(val)
            remaining -= n
        res[i] = acc
    return res


t0 = time.time()
print("--- forward ---");  fw = run(False, "F")
print("--- reverse ---");  rv = run(True, "R")


def parse(block):
    P = np.zeros((4, 4)); N = np.zeros((4, 4)); tot = np.zeros(4)
    for i in range(4):
        b = block[i]; tot[i] = sum(b.values())
        for bits, n in b.items():
            if bits[0] == '0':
                N[i, int(bits[1:], 2)] += n
    P = N / N.sum(axis=1, keepdims=True)
    return P, N, N.sum(axis=1) / tot


PF, NF, accF = parse(fw)
PR, NR, accR = parse(rv)

print(f"\n[1] conditionals vs theory: fwd {np.abs(PF-PF_th).max():.4f}  rev {np.abs(PR-PR_th).max():.4f}")

marg = pth * accF; marg /= marg.sum()
marg_th = pth * acc_th; marg_th /= marg_th.sum()
Na = float((pth * accF).sum() * 4 * SHOTS)
z = np.abs(marg - pth) / np.sqrt(pth * (1 - pth) / Na)
print(f"\n[2] biased past: deviation {np.abs(marg-pth).max():.4f} (predicted {np.abs(marg_th-pth).max():.4f})"
      f"   max z = {z.max():.1f} sigma")

D = np.array([[np.log(PF[i, j] / PR[j, i]) if PF[i,j]>0 and PR[j,i]>0 else np.nan
               for j in range(4)] for i in range(4)])
rows, vals = [], []
for i in range(4):
    for j in range(4):
        if not np.isnan(D[i, j]):
            r = np.zeros(8); r[i] = -1; r[4+j] = 1
            rows.append(r); vals.append(D[i, j])
Am = np.array(rows); b = np.array(vals)
x, *_ = np.linalg.lstsq(Am, b, rcond=None)
print(f"\n[3] separability residual = {np.abs(Am@x-b).max():.4f}   mean|D| = {np.abs(b).mean():.3f}")
print(f"    ratio residual/mean|D| = {np.abs(Am@x-b).max()/np.abs(b).mean():.3f}")
print(f"    max|D_meas - D_th| = {np.nanmax(np.abs(D-D_th)):.4f}")
print(f"    reference: generic driven channels give ratio ~1.0")

json.dump({"forward": {str(k): v for k, v in fw.items()},
           "reverse": {str(k): v for k, v in rv.items()},
           "backend": name, "shots_per_input": SHOTS,
           "ansatz_params": PARAMS.tolist(), "circuit": "shallow_4cx"},
          open(OUT, "w"), indent=2)
print(f"\nwrote {OUT}  ({time.time()-t0:.0f}s)")
