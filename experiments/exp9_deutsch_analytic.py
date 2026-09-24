"""Analytic separation of the Deutsch and Lloyd CTC models by work statistics.

For a two-level system the deviation matrix D is separable (D_ij = u_j - v_i)
iff the cycle invariant Sigma = D00 + D11 - D01 - D10 vanishes. Every P-CTC has
Sigma = 0 identically (Theorem 1). This script evaluates Sigma for the Deutsch
model on the explicit circuit

    U0 = exp(-i pi/4 SWAP) . (Ry(pi/2) x I) . CRy(pi/2)

in EXACT symbolic arithmetic (entries in Q(i, sqrt2)), under both prescriptions
for mixed inputs (the BLSS ambiguity):
  - per-branch: each pure TPM input gets its own Deutsch fixed point
    (unique here, so no maximum-entropy tie-breaking is invoked);
  - mixture (convention B): one shared fixed point from the thermal input.

Result:  Sigma_per-branch = ln(5/3) exactly;  Sigma_mixture = ln R_B > 0.
Hence the separation is proven for this circuit under either convention.
Genericity across circuits remains numerical (deutsch_discriminator.py).
"""
import sympy as sp

i_ = sp.I; s2 = sp.sqrt(2); I2 = sp.eye(2)
SWAP = sp.Matrix([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])
Ry = sp.Matrix([[1,-1],[1,1]])/s2
CRY = sp.Matrix(sp.BlockDiagMatrix(I2, Ry))
PSWAP = (sp.eye(4) - i_*SWAP)/s2


def kron(A, B):
    return sp.Matrix(sp.KroneckerProduct(A, B).doit())


U = sp.simplify(PSWAP * kron(Ry, I2) * CRY)


def ptrace_sys(M):
    return sp.Matrix(2, 2, lambda a, b: M[a, b] + M[2 + a, 2 + b])


def ptrace_ctc(M):
    return sp.Matrix(2, 2, lambda a, b: M[2 * a, 2 * b] + M[2 * a + 1, 2 * b + 1])


def fixed_points(Um, rho_s):
    a, b, c, d = sp.symbols('a b c d', complex=True)
    tau = sp.Matrix([[a, b], [c, d]])
    E = sp.expand(ptrace_sys(Um * kron(rho_s, tau) * Um.H) - tau)
    eqs = [sp.Eq(E[k, l], 0) for k in range(2) for l in range(2)] + [sp.Eq(a + d, 1)]
    return [sp.Matrix([[s[a], s[b]], [s[c], s[d]]])
            for s in sp.solve(eqs, [a, b, c, d], dict=True)]


def conditionals(Um, rho, tau=None):
    if tau is None:
        taus = fixed_points(Um, rho)
        assert len(taus) == 1, f"fixed point not unique ({len(taus)})"
        tau = taus[0]
    ro = ptrace_ctc(Um * kron(rho, tau) * Um.H)
    return sp.simplify(ro[0, 0]), sp.simplify(ro[1, 1])


E0, E1 = sp.diag(1, 0), sp.diag(0, 1)

print("=== per-branch prescription (unique fixed point per branch) ===")
P = {}
for tag, Um in [("F", U), ("R", U.H)]:
    for i, rho in [(0, E0), (1, E1)]:
        p0, p1 = conditionals(Um, rho)
        P[(tag, i, 0)], P[(tag, i, 1)] = p0, p1
        print(f"  {tag}{i}: p(0|{i})={sp.nsimplify(p0)}  p(1|{i})={sp.nsimplify(p1)}")
R = sp.radsimp(sp.simplify(
    P[("F",0,0)]*P[("F",1,1)]*P[("R",1,0)]*P[("R",0,1)] /
   (P[("F",0,1)]*P[("F",1,0)]*P[("R",0,0)]*P[("R",1,1)])))
print(f"  cross-ratio R = {sp.nsimplify(R)}   Sigma = ln R = ln(5/3) = {sp.N(sp.log(R),10)}")
assert sp.simplify(R - sp.Rational(5, 3)) == 0, "expected R = 5/3 exactly"

print("\n=== mixture prescription (shared fixed point from thermal input) ===")
w = sp.exp(sp.Rational(1, 2)); Z = w + 1/w
rho_mix = sp.diag(w/Z, (1/w)/Z)
Q = {}
for tag, Um in [("F", U), ("R", U.H)]:
    taus = fixed_points(Um, rho_mix)
    assert len(taus) == 1
    for i, rho in [(0, E0), (1, E1)]:
        p0, p1 = conditionals(Um, rho, tau=taus[0])
        Q[(tag, i, 0)], Q[(tag, i, 1)] = p0, p1
RB = sp.simplify(
    Q[("F",0,0)]*Q[("F",1,1)]*Q[("R",1,0)]*Q[("R",0,1)] /
   (Q[("F",0,1)]*Q[("F",1,0)]*Q[("R",0,0)]*Q[("R",1,1)]))
print(f"  R_B = {sp.N(RB,15)}   Sigma_B = {sp.N(sp.log(RB),10)}")
assert sp.simplify(RB - 1) != 0

print("\n=> Deutsch: Sigma != 0 under BOTH conventions (exact).")
print("   P-CTC:   Sigma = 0 identically for every circuit (Theorem 1).")
print("   The two CTC models are separated by an exact work-statistics invariant.")
