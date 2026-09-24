"""
analyze_weighted.py — statistics-aware analysis of CTC counts.

Fixes the naive analysis by (a) propagating shot-noise error bars into every
cell of the deviation matrix, (b) excluding cells with too few counts to
support a log-ratio, and (c) reporting a chi-square goodness-of-fit against
the factorization theorem D_ij = u_j - v_i rather than a raw max-residual.

Usage:  python analyze_weighted.py [file1.json file2.json ...]
        (default: ctc_hardware_counts.json)
"""
import json, sys
import numpy as np

MIN_COUNTS = 200          # cells below this in either direction are excluded

paths = sys.argv[1:] if len(sys.argv) > 1 else ["ctc_hardware_counts.json"]


def load_raw(counts_block):
    """Return R[input, output] = raw heralded counts, and totals per input."""
    R = np.zeros((4, 4)); tot = np.zeros(4)
    for i in range(4):
        ci = counts_block[str(i)] if str(i) in counts_block else counts_block[i]
        tot[i] = sum(ci.values())
        for bits, n in ci.items():
            if bits[0] == '0':
                R[i, int(bits[1:], 2)] += n
    return R, tot


for path in paths:
    d = json.load(open(path))
    rng = np.random.default_rng(d["seed"])
    A = rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8))
    U, _ = np.linalg.qr(A)
    C = np.einsum('akbk->ab', U.reshape(4, 2, 4, 2))
    u = np.log(np.diag(C @ C.conj().T).real)
    v = np.log(np.diag(C.conj().T @ C).real)
    D_th = np.array([[u[j] - v[i] for j in range(4)] for i in range(4)])

    ES = np.array([-1.5, -0.5, 0.5, 1.5]); beta = 1.0
    pth = np.exp(-beta * ES); pth /= pth.sum()

    NF, totF = load_raw(d["forward"])
    NR, totR = load_raw(d["reverse"])
    keptF = NF.sum(axis=1); keptR = NR.sum(axis=1)
    PF = NF / keptF[:, None]; PR = NR / keptR[:, None]

    print("=" * 68)
    print(f"{path}   backend={d['backend']}   shots/input={d['shots_per_input']}")

    # ---- biased past ----
    accF = keptF / totF
    marg = pth * accF; marg /= marg.sum()
    Na = float((pth * accF).sum() * totF.sum())
    z = np.abs(marg - pth) / np.sqrt(pth * (1 - pth) / Na)
    print(f"\n[biased past]  max deviation {np.abs(marg-pth).max():.4f}   max z = {z.max():.1f} sigma")
    print(f"   thermal  {np.round(pth,4)}")
    print(f"   accepted {np.round(marg,4)}")

    # ---- deviation matrix with error bars ----
    D = np.full((4, 4), np.nan); S = np.full((4, 4), np.nan)
    for i in range(4):
        for j in range(4):
            nf, nr = NF[i, j], NR[j, i]
            if nf >= MIN_COUNTS and nr >= MIN_COUNTS:
                D[i, j] = np.log(PF[i, j] / PR[j, i])
                # multinomial variance of ln p: (1-p)/n  per direction
                S[i, j] = np.sqrt((1 - PF[i, j]) / nf + (1 - PR[j, i]) / nr)
    used = ~np.isnan(D)
    print(f"\n[deviation matrix]  cells used: {used.sum()}/16 "
          f"(excluded: counts < {MIN_COUNTS} in either direction)")

    # weighted least-squares fit to D_ij = u_j - v_i
    rows, vals, sig = [], [], []
    for i in range(4):
        for j in range(4):
            if used[i, j]:
                r = np.zeros(8); r[i] = -1; r[4 + j] = 1
                rows.append(r); vals.append(D[i, j]); sig.append(S[i, j])
    Am = np.array(rows); b = np.array(vals); s = np.array(sig)
    Aw = Am / s[:, None]; bw = b / s
    x, *_ = np.linalg.lstsq(Aw, bw, rcond=None)
    resid = Am @ x - b
    pulls = resid / s
    dof = max(len(b) - 7, 1)          # 8 params, 1 gauge freedom
    chi2 = float((pulls ** 2).sum())
    print(f"   weighted fit to u_j - v_i:")
    print(f"      max |residual|      = {np.abs(resid).max():.4f}")
    print(f"      max pull (sigmas)   = {np.abs(pulls).max():.2f}")
    print(f"      chi2/dof            = {chi2/dof:.2f}   (dof={dof})")

    # direct comparison to the analytic theorem
    diff = np.abs(D - D_th)[used]
    pull_th = (np.abs(D - D_th) / S)[used]
    print(f"   vs analytic theorem D_ij = u_j - v_i:")
    print(f"      max |D_meas - D_th| = {diff.max():.4f}")
    print(f"      max pull            = {pull_th.max():.2f} sigma")
    print(f"      mean |D_meas|       = {np.abs(D[used]).mean():.3f}")
    print(f"   reference: generic driven channels give residuals 0.5-1.5")

    print("\n   measured D (excluded cells = nan):")
    print("  ", str(np.round(D, 3)).replace("\n", "\n   "))
    print("   theory D:")
    print("  ", str(np.round(D_th, 3)).replace("\n", "\n   "))
    print("   counts (forward heralded):")
    print("  ", str(NF.astype(int)).replace("\n", "\n   "))
