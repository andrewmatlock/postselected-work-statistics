"""
analyze_v2.py — pre-registered analysis for PREREGISTRATION_v2.json.

Fixes relative to analyze_shallow.py (2026-09 audit):
- Biased past: proper delta-method error propagation from binomial acceptance
  errors (the v1 formula modeled the reweighted marginal as a multinomial draw,
  understating significance and gating on the lowest-power component), plus the
  model-free acceptance-uniformity chi-square and the attenuation coefficient of
  the acceptance contrast against theory.
- Factorization: gates on GLOBAL statistics (variance explained, correlation with
  the predicted D matrix, fitted-vs-predicted potentials); the v1 max-residual
  ratio is reported but does not gate (max-statistics are fragile and sit on a
  device-systematic floor ~0.35-0.40 at any shot count).
- Potentials appearing in no usable cell are excluded from the comparison
  (a missing reverse input leaves its u unconstrained).

Usage:  python analyze_v2.py primary.json [replication.json]
Exit: prints per-criterion PASS/FAIL and the pre-registered verdict.
"""
import json, sys
import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator

# ---- pre-registered thresholds (PREREGISTRATION_v2.json is the authority) ----
T = dict(
    ve_confirm=0.90, r_confirm=0.95, u_max=0.25, v_max=0.15,
    ve_refute=0.60, r_refute=0.60,
    bp_chi2=30.0, bp_att_min=0.35, bp_att_sig=5.0,
    repl_r=0.95,
)

ES = np.array([-1.5, -0.5, 0.5, 1.5]); beta = 1.0
pth = np.exp(-beta * ES); pth /= pth.sum()


def theory(params):
    p = np.asarray(params)
    qc = QuantumCircuit(3)
    qc.ry(p[0], 0); qc.ry(p[1], 1); qc.ry(p[2], 2); qc.cx(1, 2)
    qc.rz(p[3], 1); qc.ry(p[4], 2); qc.cx(0, 2)
    qc.rz(p[5], 0); qc.ry(p[6], 2); qc.cx(1, 0)
    qc.ry(p[7], 0); qc.ry(p[8], 1); qc.cx(0, 2)
    qc.ry(p[9], 0); qc.rz(p[10], 1); qc.ry(p[11], 2)
    M = Operator(qc).data[:4, :4]
    u = np.log(np.diag(M @ M.conj().T).real)
    v = np.log(np.diag(M.conj().T @ M).real)
    D_th = np.array([[u[j] - v[i] for j in range(4)] for i in range(4)])
    return u, v, D_th, np.exp(v)


def parse(block):
    N = np.zeros((4, 4)); shots = np.zeros(4)
    for k, ci in block.items():
        i = int(k); shots[i] = sum(ci.values())
        for bits, n in ci.items():
            if bits[0] == '0':
                N[i, int(bits[1:], 2)] += n
    return N, shots


def deviation_matrix(NF, NR, min_counts=200):
    with np.errstate(divide='ignore', invalid='ignore'):
        PF = NF / NF.sum(axis=1, keepdims=True)
        PR = NR / NR.sum(axis=1, keepdims=True)
    D = np.full((4, 4), np.nan); S = np.full((4, 4), np.nan)
    for i in range(4):
        for j in range(4):
            if NF[i, j] >= min_counts and NR[j, i] >= min_counts:
                D[i, j] = np.log(PF[i, j] / PR[j, i])
                S[i, j] = np.sqrt((1 - PF[i, j]) / NF[i, j] + (1 - PR[j, i]) / NR[j, i])
    return D, S


def analyze(path):
    d = json.load(open(path))
    u_th, v_th, D_th, acc_th = theory(d.get("ansatz_params",
        [2.774141, 1.438552, 5.364134, 1.634255, 2.682839, 4.893008,
         1.171051, 1.459015, 1.474151, 0.031240, 3.189096, 2.607931]))
    NF, shF = parse(d["forward"]); NR, shR = parse(d["reverse"])
    out = {"path": path, "backend": d.get("backend", "?")}

    # ---- biased past, proper statistics ----
    a = NF.sum(axis=1) / shF
    sa2 = a * (1 - a) / shF
    k_acc = NF.sum(axis=1)
    p_pool = k_acc.sum() / shF.sum()
    out["bp_chi2"] = float((((k_acc - shF * p_pool) ** 2) / (shF * p_pool * (1 - p_pool))).sum())
    ct = acc_th - acc_th.mean()
    att = float(ct @ (a - a.mean()) / (ct @ ct))
    att_sig = np.sqrt(float((ct ** 2) @ sa2) / (ct @ ct) ** 2)
    out["att"], out["att_z"] = att, att / att_sig
    norm = (pth * a).sum(); marg = pth * a / norm
    J = np.array([[pth[i] * ((i == k) / norm - a[i] * pth[k] / norm ** 2)
                   for k in range(4)] for i in range(4)])
    cov = J @ np.diag(sa2) @ J.T
    out["bp_dev"] = float(np.abs(marg - pth).max())
    out["bp_zmax"] = float((np.abs(marg - pth) / np.sqrt(np.diag(cov))).max())

    # ---- factorization, global statistics ----
    D, S = deviation_matrix(NF, NR)
    used = ~np.isnan(D)
    out["cells"] = int(used.sum())
    rows, vals, sig = [], [], []
    for i in range(4):
        for j in range(4):
            if used[i, j]:
                r_ = np.zeros(8); r_[i] = -1; r_[4 + j] = 1
                rows.append(r_); vals.append(D[i, j]); sig.append(S[i, j])
    Am, b, s = np.array(rows), np.array(vals), np.array(sig)
    x, *_ = np.linalg.lstsq(Am / s[:, None], b / s, rcond=None)
    resid = Am @ x - b
    out["VE"] = float(1 - resid.var() / b.var())
    out["r_th"] = float(np.corrcoef(D[used], D_th[used])[0, 1])
    out["ratio"] = float(np.abs(resid).max() / np.abs(b).mean())
    out["chi2_dof"] = float(((resid / s) ** 2).sum() / max(len(b) - 7, 1))
    rows_used = used.any(axis=1); cols_used = used.any(axis=0)
    vf = x[:4]; uf = x[4:]
    gauge = (uf[cols_used] - u_th[cols_used]).mean()
    out["u_err"] = float(np.abs((uf - gauge)[cols_used] - u_th[cols_used]).max())
    out["v_err"] = float(np.abs((vf - gauge)[rows_used] - v_th[rows_used]).max())
    out["D"], out["used"] = D, used
    return out


def verdicts(o):
    checks = [
        ("cells == 16",                 o["cells"] == 16),
        (f"VE >= {T['ve_confirm']}",    o["VE"] >= T['ve_confirm']),
        (f"r(D,D_th) >= {T['r_confirm']}", o["r_th"] >= T['r_confirm']),
        (f"max|u-u_th| <= {T['u_max']}", o["u_err"] <= T['u_max']),
        (f"max|v-v_th| <= {T['v_max']}", o["v_err"] <= T['v_max']),
    ]
    confirmed = all(ok for _, ok in checks)
    refuted = (o["VE"] <= T['ve_refute']) or (o["r_th"] <= T['r_refute'])
    bp = [
        (f"acceptance-uniformity chi2 >= {T['bp_chi2']}", o["bp_chi2"] >= T['bp_chi2']),
        (f"attenuation >= {T['bp_att_min']}",             o["att"] >= T['bp_att_min']),
        (f"attenuation significance >= {T['bp_att_sig']} sigma", o["att_z"] >= T['bp_att_sig']),
    ]
    return checks, confirmed, refuted, bp, all(ok for _, ok in bp)


primary = analyze(sys.argv[1])
print(f"PRIMARY: {primary['path']}  backend={primary['backend']}  cells={primary['cells']}/16")
print(f"  VE={primary['VE']:.4f}  r(D,D_th)={primary['r_th']:.4f}  "
      f"u_err={primary['u_err']:.4f}  v_err={primary['v_err']:.4f}")
print(f"  [reported, non-gating] ratio={primary['ratio']:.4f}  chi2/dof={primary['chi2_dof']:.1f}")
print(f"  biased past: max|dev|={primary['bp_dev']:.4f} at z={primary['bp_zmax']:.1f} (proper); "
      f"chi2(3)={primary['bp_chi2']:.0f}; attenuation={primary['att']:.3f} ({primary['att_z']:.1f} sigma)")

checks, confirmed, refuted, bp, bp_ok = verdicts(primary)
print("\n--- PRE-REGISTERED CRITERIA (v2) ---")
for name, ok in checks:
    print(f"  factorization  {'PASS' if ok else 'FAIL'}  {name}")
for name, ok in bp:
    print(f"  biased past    {'PASS' if ok else 'FAIL'}  {name}")
fv = "CONFIRMED" if confirmed else ("REFUTED" if refuted else "INCONCLUSIVE")
print(f"\nFACTORIZATION: {fv}")
print(f"BIASED PAST:   {'CONFIRMED' if bp_ok else 'NOT CONFIRMED'}")

if len(sys.argv) > 2:
    rep = analyze(sys.argv[2])
    both = primary["used"] & rep["used"]
    r_ab = float(np.corrcoef(primary["D"][both], rep["D"][both])[0, 1])
    print(f"\nREPLICATION: {rep['path']}  backend={rep['backend']}  cells={rep['cells']}/16")
    print(f"  VE={rep['VE']:.4f}  r(D,D_th)={rep['r_th']:.4f}  cross-device r(D_A,D_B)={r_ab:.4f}")
    _, rc, rr, _, rbp = verdicts(rep)
    print(f"  device-independence (secondary): factorization "
          f"{'CONFIRMED' if rc else ('REFUTED' if rr else 'INCONCLUSIVE')}; "
          f"cross-device r >= {T['repl_r']}: {'PASS' if r_ab >= T['repl_r'] else 'FAIL'}")
