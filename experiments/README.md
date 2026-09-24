# Experiments

Each script is self-contained (numpy/scipy), fixed-seed, runs in seconds, and prints the claim it supports.

| Script | Claim |
|---|---|
| `exp1_classical_ledger.py` | Reversal-symmetric closure preserves the detailed FT; twisted retro-constraints amputate the entropy-production spectrum (= absolute irreversibility) |
| `exp4_quantum_factorization.py` | No quantum amputation (zero-pattern identity); P-CTC Crooks deviations factorize exactly (residual ~1e-16) |
| `exp5_forgery_scaling.py` | Random-channel null (physical-inverse convention); the adjoint-renormalization artifact; unital control; dimension-independence. Regenerated 2026-09 — see `notes/audit-2026-09.md` |
| `exp7_adversary.py` | Exact forgery by generalized detailed balance (wrong-temperature bath); adversarial near-forgery |
| `poke12_nogo_prescription.py` | Halmos dilation reproduces all accepted-run statistics; prescription reconciliation; the biased past; Jarzynski efficacy |
| `deutsch_discriminator.py` | D-CTC deviations non-separable vs P-CTC separable — a work-statistics discriminator between CTC models |
| `deutsch_mixture_check.py` | The discriminator survives both mixture conventions numerically |
| `exp9_deutsch_analytic.py` | **Analytic separation (exact):** for an explicit circuit, Deutsch cycle invariant Σ = ln(5/3) under per-branch and ≈0.530 under mixture convention; P-CTC Σ = 0 identically. Sympy, exact arithmetic, asserted |
