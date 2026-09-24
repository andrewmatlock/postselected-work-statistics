"""
analyze_counts.py — re-analyze ctc_hardware_counts.json (no new shots needed).
Usage:  python analyze_counts.py [path_to_json]
Default path: ctc_hardware_counts.json in the current directory.
"""
import json, sys
import numpy as np

path = sys.argv[1] if len(sys.argv) > 1 else "ctc_hardware_counts.json"
d = json.load(open(path))
print(f"file: {path}   backend: {d['backend']}   shots/input: {d['shots_per_input']}")

# rebuild the same operator from the stored seed
rng = np.random.default_rng(d["seed"])
A = rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8))
U, _ = np.linalg.qr(A)
C = np.einsum('akbk->ab', U.reshape(4, 2, 4, 2))
ES = np.array([-1.5, -0.5, 0.5, 1.5]); beta = 1.0
pth = np.exp(-beta * ES); pth /= pth.sum()


def cond(counts):
    """counts[str(i)] -> dict. Returns P[input, output] and acceptance per input."""
    P = np.zeros((4, 4)); acc = np.zeros(4)
    for i in range(4):
        ci = counts[str(i)] if str(i) in counts else counts[i]
        tot = sum(ci.values()); kept = 0
        for bits, n in ci.items():
            if bits[0] == '0':                      # herald: ancilla clbit = 0
                P[i, int(bits[1:], 2)] += n; kept += n
        acc[i] = kept / tot
        P[i] /= max(kept, 1)
    return P, acc


PF, accF = cond(d["forward"])
PR, accR = cond(d["reverse"])

# --- [1] theory check on conditionals ---
PF_th = np.array([[abs(C[j, i])**2 / np.vdot(C[:, i], C[:, i]).real
                   for j in range(4)] for i in range(4)])
Cd = C.conj().T
# PR[a,b] = p~(b|a): input a evolves under V-dagger (operator C-dagger)
PR_th = np.array([[abs(Cd[b, a])**2 / np.vdot(Cd[:, a], Cd[:, a]).real
                   for b in range(4)] for a in range(4)])
print(f"\n[1] conditionals vs theory:  forward max|diff| = {np.abs(PF-PF_th).max():.4f}"
      f"   reverse max|diff| = {np.abs(PR-PR_th).max():.4f}")

# --- [2] the biased past ---
marg = pth * accF; marg /= marg.sum()
Na = d["shots_per_input"] * 4 * float(pth @ accF)
z = np.abs(marg - pth) / np.sqrt(pth * (1 - pth) / Na)
print(f"\n[2] biased past: max deviation {np.abs(marg-pth).max():.4f}   max z = {z.max():.1f} sigma")
print(f"    thermal  {np.round(pth,4)}")
print(f"    accepted {np.round(marg,4)}")

# --- [3] deviation matrix, corrected indexing ---
# D_ij = ln[ p(j|i) / p~(i|j) ];  p~(i|j) = PR[j, i]  (input j, output i)
D = np.full((4, 4), np.nan)
for i in range(4):
    for j in range(4):
        if PF[i, j] > 1e-6 and PR[j, i] > 1e-6:
            D[i, j] = np.log(PF[i, j] / PR[j, i])

u = np.log(np.diag(C @ C.conj().T).real)      # ||C†|j>||²
v = np.log(np.diag(C.conj().T @ C).real)      # ||C|i>||²
D_th = np.array([[u[j] - v[i] for j in range(4)] for i in range(4)])
print(f"\n[3] deviation matrix vs theorem D_ij = u_j - v_i:")
print(f"    max |D_measured - D_theory| = {np.nanmax(np.abs(D - D_th)):.4f}")

rows, vals = [], []
for i in range(4):
    for j in range(4):
        if not np.isnan(D[i, j]):
            r = np.zeros(8); r[i] = -1; r[4 + j] = 1
            rows.append(r); vals.append(D[i, j])
Am = np.array(rows); b = np.array(vals)
x, *_ = np.linalg.lstsq(Am, b, rcond=None)
res = float(np.abs(Am @ x - b).max())
print(f"    separability residual = {res:.4f}   (mean|deviation| = {np.abs(b).mean():.3f})")
print(f"    reference: exact theory ~1e-16 | generic driven channels 0.5-1.5")

print("\nmeasured D:\n", np.round(D, 3))
print("theory   D:\n", np.round(D_th, 3))
