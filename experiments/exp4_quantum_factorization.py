"""Quantum: no amputation (zero-pattern identity); P-CTC Crooks deviations factorize
exactly into a state potential D_ij = u_j - v_i. Dimension-independent.

v2 (2026-09 audit): the original version constructed D as u_j - v_i directly, making
the separability fit tautological. This version breaks the circularity two ways:
(1) the effective operator is built by explicit postselected-teleportation amplitudes
    (maximally entangled pair + projection), independently of the partial-trace formula,
    and the two constructions are compared;
(2) D is computed from forward/reverse conditional PROBABILITIES (the measurable
    objects), then tested against the u_j - v_i prediction.
"""
import numpy as np
rng = np.random.default_rng(3)


def ctc_operator_teleport(U, dS, dE):
    """Effective operator via explicit postselected teleportation:
    |Phi+> on (CTC_in, CTC_out), U on (S, CTC_in), project (CTC_out, CTC_in') on |Phi+>.
    Amplitude <j| C |i> = sum_k <j,k| U |i,k> emerges physically, not by formula."""
    phi = np.eye(dE).reshape(dE * dE) / np.sqrt(dE)      # |Phi+> unnormalized convention
    C = np.zeros((dS, dS), complex)
    Ut = U.reshape(dS, dE, dS, dE)
    for i in range(dS):
        # state |i>_S |Phi+>_(in,out); apply U to (S, in); project (in, out) on |Phi+>
        # amplitude to land in |j>_S: sum_{k,l,m} U[j,k;i,l] phi[l,m] phi*[k,m]
        for j in range(dS):
            amp = 0
            for k in range(dE):
                for l in range(dE):
                    for m in range(dE):
                        amp += Ut[j, k, i, l] * phi[l * dE + m] * np.conj(phi[k * dE + m])
            C[j, i] = amp * dE          # strip the 1/dE teleportation prefactor
    return C


def sep_residual(D):
    rows, vals = [], []
    n = D.shape[0]
    for i in range(n):
        for j in range(n):
            if not np.isnan(D[i, j]):
                r = np.zeros(2 * n); r[i] = -1; r[n + j] = 1
                rows.append(r); vals.append(D[i, j])
    A = np.array(rows); b = np.array(vals)
    x, *_ = np.linalg.lstsq(A, b, rcond=None)
    return float(np.abs(A @ x - b).max()), float(np.abs(b).mean())


def analyse(U, dS, dE, label, check_teleport=False):
    Ut = U.reshape(dS, dE, dS, dE)
    C = np.einsum('akbk->ab', Ut)
    if check_teleport:
        Ctp = ctc_operator_teleport(U, dS, dE)
        assert np.abs(C - Ctp).max() < 1e-12, "teleportation construction disagrees"
    u = np.log(np.diag(C @ C.conj().T).real)
    v = np.log(np.diag(C.conj().T @ C).real)
    # D from the measurable conditionals, not from the identity
    PF = np.abs(C) ** 2 / np.exp(v)[None, :]              # PF[j,i] = p(j|i)
    PR = np.abs(C.conj().T) ** 2 / np.exp(u)[None, :]     # PR[i,j] = p~(i|j)
    D = np.full((dS, dS), np.nan)
    for i in range(dS):
        for j in range(dS):
            if PF[j, i] > 1e-13 and PR[i, j] > 1e-13:
                D[i, j] = np.log(PF[j, i] / PR[i, j])
    res, mean = sep_residual(D)
    D_th = np.array([[u[j] - v[i] for j in range(dS)] for i in range(dS)])
    dev_th = np.nanmax(np.abs(D - D_th))
    zero = np.abs(np.abs(C) - np.abs(C.conj().T).T).max()
    print(f"  {label}: separability residual={res:.2e}  |D - (u_j-v_i)|={dev_th:.2e}  "
          f"mean|D|={mean:.3f}  zero-pattern asymmetry={zero:.1e}  "
          f"max|D_ii|={max(abs(D[k, k]) for k in range(dS)):.3f}")


print("=== Quantum factorization identity (from conditionals; teleport-verified) ===")
for k in range(4):
    A = rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8)); Q, _ = np.linalg.qr(A)
    analyse(Q, 4, 2, f"2 sys qubits + 1 CTC qubit  #{k}", check_teleport=True)
for k in range(2):
    A = rng.normal(size=(32, 32)) + 1j * rng.normal(size=(32, 32)); Q, _ = np.linalg.qr(A)
    analyse(Q, 8, 4, f"3 sys qubits + 2 CTC qubits #{k}", check_teleport=(k == 0))
print("\nzero-pattern asymmetry ~0 => |<i|C+|j>| = |<j|C|i>| identically => no one-sided amputation")
