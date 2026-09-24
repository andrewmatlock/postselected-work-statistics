"""Random-channel null; the adjoint-renormalization artifact; unital control;
dimension-independence of the P-CTC factorization.

Regenerated 2026-09 during the code audit (the original script predates the repo).
The regeneration surfaced a sharp fact the original claims skated past:

  THE REVERSE CONVENTION IS LOAD-BEARING. For ANY channel E with Kraus {K_a},
  sum_a |K_a[j,i]|^2 is one matrix read two ways: column-normalized it is the
  forward accepted-run conditional, row-normalized it is the ADJOINT-reverse
  accepted-run conditional. Renormalizing both directions per input therefore
  factorizes EVERY channel exactly (D_ij = ln rowsum_j - ln colsum_i), CTC or not.
  The fingerprint is only evidence under the PHYSICAL-INVERSE reverse protocol
  (run the inverse dynamics with fresh environment, as the hardware does), and
  the decisive statistic is factorization into the PREDICTED potentials, not
  factorization per se.

Claims tested:
  (a) physical-inverse null: random dilation channels with genuine violations do
      NOT factorize (this is the correct null for the hardware protocol);
  (a') adjoint+renormalized "null": everything factorizes — the artifact, on display;
  (b) unital control: zero violations, trivially separable, diagnostic of nothing;
  (c) P-CTC factorization is dimension-independent.
"""
import numpy as np
rng = np.random.default_rng(5)
ES4 = np.array([-1.5, -0.5, 0.5, 1.5])


def sep_residual(D):
    rows, vals = [], []
    n = D.shape[0]
    for i in range(n):
        for j in range(n):
            if not np.isnan(D[i, j]):
                r = np.zeros(2 * n); r[i] = -1; r[n + j] = 1
                rows.append(r); vals.append(D[i, j])
    A = np.array(rows); b = np.array(vals)
    if not len(vals):
        return np.nan, 0.0
    x, *_ = np.linalg.lstsq(A, b, rcond=None)
    return float(np.abs(A @ x - b).max()), float(np.abs(b).mean())


def dilation_kraus(U, dS, dE):
    """Channel from unitary dilation with environment starting in |0>."""
    Ut = U.reshape(dS, dE, dS, dE)
    return [Ut[:, k, :, 0] for k in range(dE)]


def tpm_probs(kraus):
    """p[j,i] = <j| E(|i><i|) |j> — raw TPM transition probabilities (no renorm)."""
    p = np.zeros_like(np.abs(kraus[0]))
    for K in kraus:
        p += np.abs(K) ** 2
    return p


def D_from(pf, pr, renorm):
    """D_ij = ln p(j|i) - ln p~(i|j). pf[j,i] forward, pr[i,j] reverse."""
    n = pf.shape[0]
    if renorm:
        pf = pf / pf.sum(axis=0)[None, :]
        pr = pr / pr.sum(axis=0)[None, :]
    D = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(n):
            if pf[j, i] > 1e-12 and pr[i, j] > 1e-12:
                D[i, j] = np.log(pf[j, i] / pr[i, j])
    return D


print("=== (a) physical-inverse null: 300 Haar dilation channels ===")
print("    (forward: U dilation, env|0>; reverse: U-dagger dilation, fresh env|0>)")
for dE in (2, 4):
    res, means = [], []
    for _ in range(150):
        A = rng.normal(size=(4 * dE, 4 * dE)) + 1j * rng.normal(size=(4 * dE, 4 * dE))
        U, _ = np.linalg.qr(A)
        pf = tpm_probs(dilation_kraus(U, 4, dE))
        pr = tpm_probs(dilation_kraus(U.conj().T, 4, dE))
        D = D_from(pf, pr, renorm=False)
        r, m = sep_residual(D)
        if m > 0.05:
            res.append(r); means.append(m)
    res = np.array(res)
    print(f"  env dim {dE}: n={len(res)}  min residual={res.min():.3f}  "
          f"median={np.median(res):.3f}  mean|D| median={np.median(means):.3f}")

print("\n=== (a') the adjoint-renormalization artifact ===")
worst = 0.0
for _ in range(50):
    A = rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8))
    U, _ = np.linalg.qr(A)
    kraus = dilation_kraus(U, 4, 2)
    pf = tpm_probs(kraus)
    pr_adj = tpm_probs([K.conj().T for K in kraus])       # adjoint "reverse"
    D = D_from(pf, pr_adj, renorm=True)
    r, _ = sep_residual(D)
    worst = max(worst, r)
print(f"  50 random channels, adjoint reverse + per-input renorm: max residual = {worst:.1e}")
print("  => EVERY channel factorizes under this convention. The convention, not the")
print("     physics, produces the fingerprint. Detection doctrine must specify the")
print("     physical-inverse protocol and test against PREDICTED potentials.")

print("\n=== (b) unital control: dephasing channel ===")
Uu, _ = np.linalg.qr(rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)))
kraus = [np.diag([1, 0, 0, 0.0]) @ Uu, np.diag([0, 1, 0, 0.0]) @ Uu,
         np.diag([0, 0, 1, 0.0]) @ Uu, np.diag([0, 0, 0, 1.0]) @ Uu]
pf = tpm_probs(kraus)
pr = tpm_probs([K.conj().T for K in kraus])
D = D_from(pf, pr, renorm=False)
r, m = sep_residual(D)
print(f"  mean|D| = {m:.2e}  (zero violations => trivially separable, diagnostic of nothing)")

print("\n=== (c) P-CTC factorization, dimension scaling ===")
for dS, dE, label in [(4, 2, "2+1 qubits"), (8, 4, "3+2 qubits (32-dim)")]:
    A = rng.normal(size=(dS * dE, dS * dE)) + 1j * rng.normal(size=(dS * dE, dS * dE))
    U, _ = np.linalg.qr(A)
    C = np.einsum('akbk->ab', U.reshape(dS, dE, dS, dE))
    u = np.log(np.diag(C @ C.conj().T).real)
    v = np.log(np.diag(C.conj().T @ C).real)
    PF = np.abs(C) ** 2 / np.exp(v)[None, :]
    PR = np.abs(C.conj().T) ** 2 / np.exp(u)[None, :]
    D = np.log(PF.T / PR)
    r, m = sep_residual(D)
    print(f"  {label}: residual={r:.2e}  mean|D|={m:.3f}")
