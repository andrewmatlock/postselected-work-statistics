"""Generate all manuscript figures from the committed raw data in ../data/.
Run from paper/:  ../venv/bin/python make_figures.py
"""
import json, sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "hardware"))
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator

DATA = os.path.join(os.path.dirname(__file__), "..", "data")
plt.rcParams.update({"font.size": 8.5, "axes.titlesize": 9, "figure.dpi": 200,
                     "savefig.bbox": "tight", "font.family": "serif"})


def block(params):
    p = np.asarray(params)
    qc = QuantumCircuit(3)
    qc.ry(p[0], 0); qc.ry(p[1], 1); qc.ry(p[2], 2); qc.cx(1, 2)
    qc.rz(p[3], 1); qc.ry(p[4], 2); qc.cx(0, 2)
    qc.rz(p[5], 0); qc.ry(p[6], 2); qc.cx(1, 0)
    qc.ry(p[7], 0); qc.ry(p[8], 1); qc.cx(0, 2)
    qc.ry(p[9], 0); qc.rz(p[10], 1); qc.ry(p[11], 2)
    return Operator(qc).data[:4, :4]


def theory_D(params):
    M = block(params)
    u = np.log(np.diag(M @ M.conj().T).real)
    v = np.log(np.diag(M.conj().T @ M).real)
    return np.array([[u[j] - v[i] for j in range(4)] for i in range(4)]), u, v


def measured(path):
    d = json.load(open(os.path.join(DATA, path)))
    def parse(b):
        N = np.zeros((4, 4))
        for k, ci in b.items():
            for bits, n in ci.items():
                if bits[0] == '0':
                    N[int(k), int(bits[1:], 2)] += n
        return N
    NF, NR = parse(d["forward"]), parse(d["reverse"])
    PF = NF / NF.sum(axis=1, keepdims=True)
    PR = NR / NR.sum(axis=1, keepdims=True)
    D = np.full((4, 4), np.nan); S = np.full((4, 4), np.nan)
    for i in range(4):
        for j in range(4):
            if NF[i, j] >= 200 and NR[j, i] >= 200:
                D[i, j] = np.log(PF[i, j] / PR[j, i])
                S[i, j] = np.sqrt((1 - PF[i, j]) / NF[i, j] + (1 - PR[j, i]) / NR[j, i])
    acc = NF.sum(axis=1) / np.array([sum(d["forward"][str(i)].values()) for i in range(4)])
    sacc = np.sqrt(acc * (1 - acc) / np.array([sum(d["forward"][str(i)].values()) for i in range(4)]))
    return D, S, acc, sacc, d


def fit_potentials(D, S):
    rows, vals, sig = [], [], []
    for i in range(4):
        for j in range(4):
            if not np.isnan(D[i, j]):
                r = np.zeros(8); r[i] = -1; r[4 + j] = 1
                rows.append(r); vals.append(D[i, j]); sig.append(S[i, j])
    Am, b, s = np.array(rows), np.array(vals), np.array(sig)
    x, *_ = np.linalg.lstsq(Am / s[:, None], b / s, rcond=None)
    return x[:4], x[4:]


CTC_PARAMS = [2.774141, 1.438552, 5.364134, 1.634255, 2.682839, 4.893008,
              1.171051, 1.459015, 1.474151, 0.031240, 3.189096, 2.607931]
D_th, u_th, v_th = theory_D(CTC_PARAMS)
Dm, Sm, accm, saccm, dm = measured("ctc_prereg2_marrakesh.json")
Df, Sf, accf, saccf, df = measured("ctc_prereg2_fez.json")

# ---------------- Fig 1: the deviation-matrix ledger ----------------
fig, axes = plt.subplots(1, 4, figsize=(7.0, 1.9),
                         gridspec_kw={"width_ratios": [1, 1, 1, 1.35]})
vmax = 1.45
for ax, M, title in zip(axes[:3], [D_th, Dm, Df],
                        ["theory $u_j-v_i$", "ibm_marrakesh", "ibm_fez"]):
    im = ax.imshow(M, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_title(title)
    ax.set_xticks(range(4)); ax.set_yticks(range(4))
    ax.set_xlabel("$j$")
    if ax is axes[0]:
        ax.set_ylabel("$i$")
ax = axes[3]
for D_, S_, c, lab in [(Dm, Sm, "#c0392b", "marrakesh"), (Df, Sf, "#2471a3", "fez")]:
    m = ~np.isnan(D_)
    ax.errorbar(D_th[m], D_[m], yerr=S_[m], fmt="o", ms=3, lw=0.8, color=c,
                label=lab, alpha=0.85)
lim = [-1.55, 1.05]
ax.plot(lim, lim, "k-", lw=0.7)
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("predicted $D_{ij}$"); ax.set_ylabel("measured $D_{ij}$")
ax.legend(frameon=False, loc="upper left", fontsize=7.5)
fig.colorbar(im, ax=axes[:3], shrink=0.85, pad=0.02)
fig.savefig("fig1_ledger.pdf")
plt.close(fig)

# ---------------- Fig 2: potentials and the biased past ----------------
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 2.1))
xs = np.arange(4)
for D_, S_, c, lab, off in [(Dm, Sm, "#c0392b", "marrakesh", -0.07),
                            (Df, Sf, "#2471a3", "fez", 0.07)]:
    vf, uf = fit_potentials(D_, S_)
    g = (uf - u_th).mean()
    a1.plot(xs + off, uf - g, "o", ms=4, color=c, label=f"$u_j$ {lab}")
    a1.plot(xs + off, vf - g, "s", ms=4, color=c, mfc="none", label=f"$v_i$ {lab}")
a1.plot(xs, u_th, "k^", ms=6, mfc="none", label="theory $u$")
a1.plot(xs, v_th, "kv", ms=6, mfc="none", label="theory $v$")
a1.set_xticks(xs); a1.set_xlabel("state index"); a1.set_ylabel("potential")
a1.legend(frameon=False, fontsize=6.5, ncol=2)
a1.set_title("fitted vs.\\ predicted potentials")

acc_th = np.exp(v_th)
a2.plot(xs, acc_th, "k^", ms=6, mfc="none", label="theory")
a2.errorbar(xs - 0.05, accm, yerr=saccm, fmt="o", ms=4, color="#c0392b", label="marrakesh")
a2.errorbar(xs + 0.05, accf, yerr=saccf, fmt="o", ms=4, color="#2471a3", label="fez")
a2.axhline(float((accm).mean()), color="#c0392b", lw=0.6, ls=":")
a2.set_xticks(xs); a2.set_xlabel("input state $i$"); a2.set_ylabel("acceptance $\\|C|i\\rangle\\|^2$")
a2.legend(frameon=False, fontsize=7)
a2.set_title("input-dependent acceptance (the biased initial marginal)")
fig.tight_layout()
fig.savefig("fig2_potentials.pdf")
plt.close(fig)

# ---------------- Fig 3: taxonomy, paired same-job test ----------------
pp = os.path.join(DATA, "paired_counts.json")
if os.path.exists(pp):
    pd_ = json.load(open(pp))

    def measured_block(block):
        def parse(b):
            N = np.zeros((4, 4))
            for k, ci in b.items():
                for bits, n in ci.items():
                    if bits[0] == '0':
                        N[int(k), int(bits[1:], 2)] += n
            return N
        NF, NR = parse(block["forward"]), parse(block["reverse"])
        PF = NF / NF.sum(axis=1, keepdims=True)
        PR = NR / NR.sum(axis=1, keepdims=True)
        D = np.log(PF / PR.T)
        S = np.sqrt((1 - PF) / NF + (1 - PR.T) / NR.T)
        return D, S

    Dp, Sp = measured_block(pd_["ctc"])
    Dq, Sq = measured_block(pd_["control"])
    _, ug2, vg2 = theory_D(pd_["ctrl_params"])
    Dth_ctrl = np.array([[vg2[j] - vg2[i] for j in range(4)] for i in range(4)])

    fig, (b1, b2) = plt.subplots(1, 2, figsize=(7.0, 2.0))
    w = 0.35
    b1.bar(xs - w / 2, np.diag(Dp), w, color="#c0392b",
           label="non-normal $C$ (CTC circuit)")
    b1.bar(xs + w / 2, np.diag(Dq), w, color="#7f8c8d",
           label="normal $C$ (control)")
    b1.axhline(0, color="k", lw=0.6)
    b1.set_xticks(xs); b1.set_xlabel("state $i$"); b1.set_ylabel("$D_{ii}$")
    b1.set_title("zero-work deviations, same job, same qubits")
    b1.legend(frameon=False, fontsize=7)

    m = ~np.isnan(Dq)
    b2.errorbar(Dth_ctrl[m], Dq[m], yerr=Sq[m], fmt="o", ms=3, lw=0.8,
                color="#7f8c8d", label="control arm")
    lim2 = [-1.55, 1.55]
    b2.plot(lim2, lim2, "k-", lw=0.7)
    b2.set_xlim(lim2); b2.set_ylim(lim2)
    b2.set_xlabel("predicted $v_j - v_i$ (single potential)")
    b2.set_ylabel("measured $D_{ij}$")
    b2.set_title("control arm: single-potential class")
    fig.tight_layout()
    fig.savefig("fig3_control.pdf")
    plt.close(fig)
    print("fig3 (paired) written")

print("figures written")
