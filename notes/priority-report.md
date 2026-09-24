# Priority Report — Systematic Prior-Art Check
**Date: 2026-09-20.** Scope: the 19 assessed papers in `Research Papers/` plus two
external leads flagged during the scan (Zeng–Yong 2017; Albash–Lidar–Marvian–Zanardi 2013),
checked against the note's five claims. This discharges the "gating task" of README §Open
question 3, within the scope below. *Not yet a full database sweep; see Residual risk.*

## Claims checked
- **C1** — factorization D_ij = u_j − v_i with u_j = ln‖C†|j⟩‖², v_i = ln‖C|i⟩‖²; D_ii ≠ 0 iff C non-normal
- **C2** — the biased past: P(i) = p_th(i)‖C|i⟩‖²/N_F with zero discarded runs
- **C3** — Jarzynski functional = constant efficacy N_R/N_F
- **C4** — attribution no-go via Halmos dilation (heralded device ≡ CTC on accepted runs)
- **C5** — D-CTC vs P-CTC work-statistics discriminator

## Verdict summary

| Claim | Status after sweep |
|---|---|
| C1 factorization | **Not found as stated, in 21 papers.** But it is *corollary-distance* from Åberg's conditional-FT framework (see below) and must be presented as a sharpening of it, not an independent structure. |
| C2 biased past | **Quantitatively anticipated** by Brun–Wilde 1008.0433, Eqs. (4)–(5): the reweighted-and-renormalized prior, with the explicit "changes probabilities set before the P-CTC existed" reading. Residual: only the TPM/thermal framing and the zero-discard emphasis. |
| C3 efficacy | **Anticipated three ways**: Åberg 1601.01302 Eq. (H21) (acceptance-ratio Jarzynski); Malakar–Silva 2504.11664 Eqs. (22),(24) (γ_t ≡ N_R/N_F verbatim in operator form, non-normality noted); Sagawa–Ueda 0907.4914 / Junier et al. 0901.3886 (efficacy lineage). Present as a known-type result applied to P-CTCs. |
| C4 no-go | Operational core (P-CTC ≡ postselected teleportation ≡ heralded device on accepted runs) is **standard** — definitional in Lloyd et al. 1005.2219/1007.2615, explicit in Brun–Wilde 1504.05911. Residuals: the Halmos construction for arbitrary C, TPM completeness incl. the biased marginal, and the discard-count/denominator framing. |
| C5 discriminator | **Cleanest surviving novelty.** 1711.08334 gives only the entropic D-vs-P separation (already cited as complementary); no work-statistics separation found anywhere. Mixture-ambiguity caveat stands. |

## The three closest approaches to C1 (decreasing danger)
1. **Åberg, arXiv:1601.01302** (PRX 8, 011019). Conditional FTs, Eqs. (26)/(28)/(29)/(H16):
   herald-conditioned relations whose Crooks deviation is exactly ln Z(Q^f) − ln Z(Q^i) — a
   genuine two-potential factorization, but with **dynamics-independent boundary potentials**
   Z_βH(Q) = Tr(e^{−βH/2}Q e^{−βH/2}) acting on **sub-normalized** statistics. The
   sub-normalized version of our Theorem 1 is a special case of his global relation (29)
   with sharp energy projectors plus a herald operator; D_ij = u_j − v_i then follows by
   rearranging normalization constants. His Eq. (H21) is the acceptance-ratio Jarzynski.
2. **Junier–Mossa–Manosas–Ritort, arXiv:0901.3886** (classical EFR; also Maragakis et al.
   therein). Eq. (3): subset-conditioned Crooks deviation = βG_S1(λ1) − βG_S0(λ0) — a
   two-potential form with nonzero diagonal deviations, but with **equilibrium branch free
   energies**; collapses to our single-potential forgery class in the ΔF = 0 setting.
   Eq. (2) is the classical acceptance-ratio Jarzynski (2009).
3. **Erdamar et al., arXiv:2309.12393** (Murch group, superconducting qubit). Eq. (6) is
   *literally our conditional object*: P_ij = |⟨i|G(τ)|j⟩|²/⟨j|G†G|j⟩ — per-input-renormalized
   post-selected TPM conditionals of a non-normal operator, measured on hardware, with
   ⟨e^{−βW}⟩ from them. Never forms the reverse process, pairwise ratios, or any factorization.

## External leads checked (beyond the folder)
- **Zeng & Yong, J. Phys. Commun. 1, 031001 (2017)** ("Crooks FT in PT-symmetric QM"):
  works in the **unbroken PT phase where dynamics is unitary** (their stated crucial
  assumption); standard Crooks recovered unmodified; no post-selected deviation, no
  factorization, no efficacy ratio. **Cleared.**
- **Albash–Lidar–Marvian–Zanardi, PRE 88, 032146 (2013)** (arXiv:1212.6589): efficacy
  γ = Tr[ℰ*(ρ_q)] (their Eq. 30) is a **single collective factor** for CPTP maps —
  "unitality replaces microreversibility"; detailed FT (Eq. 29) relates forward PDF to an
  *unnormalized* reverse PDF; **no per-pair factorized deviation, no per-state potentials**.
  **Cleared** for C1; it is the source of Malakar–Silva's efficacy (relevant to C3 lineage).

## Other paper-by-paper results (from the folder scan)
- hep-th/0310281 (Horowitz–Maldacena): background only, no overlap.
- 0908.3023 (BLSS): source of the C5 mixture caveat, no overlap with C1–C4.
- 0911.2666 (Esposito–Van den Broeck): generic FT boundary terms are single-family
  state-function differences; weakly adjacent to C1's form only.
- 1205.4176 (Seifert review): background for the forgery class.
- 1307.5370 (Rastegin–Życzkowski): JE correction for nonunital CPTP maps; adjacent to C3,
  not of N_R/N_F form.
- 1308.4209 (Lloyd–Preskill): postselection unattributable to resource-limited interior
  observers via distinguishability/complexity — thematic C4 kin, different mechanism.
- 1401.4494 (Murashita–Funo–Ueda): absolute irreversibility, already correctly located in note §5.
- 1705.06513 (Murashita et al.): outcome-conditioned FTs with constant outcome-dependent
  RHS; backward post-selection structurally necessary; adjacent to C2/C3.
- 2307.05713 (Manzano et al.): martingale/stopping-time absolute irreversibility; framework
  overlap with §5 only.
- 1711.08334 (Bartkiewicz–Grudka–Horodecki–Łodyga–Wychowaniec): the entropic D-vs-P
  separation. **Correction for the note's reference list: this, not "Czachor et al.," is
  PRA 99, 022304 (2019).**

## Required edits to the note (applied in v4 of the note)
1. §2 priority note rewritten: C1 presented as a sharpening at corollary distance from
   Åberg (29)/(H16); C2 credited to Brun–Wilde Eqs. (4)–(5); C3 presented as known-type
   (Åberg H21; Malakar–Silva γ_t = N_R/N_F; Sagawa–Ueda lineage).
2. §3 remark: cite Brun–Wilde 1504.05911 (simulation equivalence standard) and Lloyd–Preskill
   1308.4209 (complexity-suppressed detectability) as kin of the no-go.
3. References: fix Czachor → Bartkiewicz et al.; add Junier et al., Erdamar et al.,
   Malakar–Silva, Albash et al., Maragakis et al.
4. Reframe the note's headline: the *composition* — all four P-CTC signatures derived from
   one joint, the taxonomy (single- vs two-potential), and the D-vs-P work-statistics
   discriminator — is the contribution; no single formula except possibly C1's operator-norm
   potentials and C5 is individually new.

## Residual risk
This sweep covers 22 papers, chosen well but finitely. Also checked: **Kwon–Kim, PRX 9,
031029** (arXiv:1810.03150) — Petz-recovery-map channel FTs; their pure-state correction
Υ has two-potential form log⟨ψ|e^{βH}|ψ⟩ − log⟨ϕ′|e^{−βH}|ϕ′⟩, but again **Hamiltonian-
dependent boundary potentials, not dynamics-operator norms**; no acceptance renormalization.
**Cleared.** The consistent pattern across the literature: two-potential deviations exist
wherever boundary conditions are exotic, but the potentials are always functions of H and
the boundary operators — never ln‖C|i⟩‖² of the effective dynamics. That specific content
is what C1 owns.

**Citation sweeps completed 2026-09-20.** All 123 Semantic Scholar citations of Åberg
PRX 8, 011019 and all 101 of Albash et al. PRE 88, 032146 screened (~15 deep-dived);
Malakar–Silva and Erdamar citation trees walked; INSPIRE full-text phrase searches
("postselected Crooks", "heralded fluctuation theorem", etc.): zero hits. **C1 not
anticipated anywhere.** New closest approaches, all ADJACENT, now cited in the
manuscript: Prech–Potts PRL 133, 140401 (general measurement/feedback FT — the
umbrella over sub-normalized Theorem 1, alongside Åberg); arXiv:2605.10099 (trapped-ion
non-Hermitian JE, 2026 — uses Theorem 1(a)'s conditionals verbatim, no reverse
protocol; the competing lineage is publishing NOW); Ferri-Cortés et al. PRR 7, 013077
(per-record vs per-pair contrast); Buscemi–Scarani PRE 103, 052111 (reverse-process
foundations — load-bearing for the convention caveat); Maeda et al. PRA 108, L050203
(the only other dynamics-dependent potentials: conditional energies, not operator
norms); Purves–Short PRE 104, 014111 (post-selection energetics). Residual risk now:
Google Scholar's larger index vs Semantic Scholar's, and papers published after this
sweep. The §8(4) debts are discharged in scope.

**Sweep B (non-Hermitian tree), same day — URGENCY FINDING.** Full trees of
Malakar–Silva and Erdamar walked; arXiv phrase sweeps. **arXiv:2607.24961
(Cavina–Quintela Rodríguez–Farina, 27 Jul 2026) ANTICIPATES-IN-PART:** they build the
time-reversed Kraus/adjoint reverse protocol for post-selected no-jump dynamics, derive
a renormalized Crooks relation with a single CONSTANT offset Ξ = ln N − ln Ñ (their
Eq. 9), a constant Jarzynski functional (their Eq. 13 ≡ our C3, now anticipated a
FOURTH way — concede (d) entirely), and tie the quantum part to non-normality. What
they lack, and where C1 survives: renormalization is global with reverse prior ρ(t),
so there is no pairwise D_ij, no operator-norm potentials, no per-transition
D_ii ≠ 0 ⇔ non-normality. Our factorization is what their constant unfolds into under
per-input/thermal-prior renormalization — the manuscript now says exactly this. They
cite neither Åberg nor Malakar–Silva nor Erdamar: the communities are disconnected and
the present paper connects them. Also: 2605.10099's symmetry-enforced JE is the
efficacy=1 special case of (d) — noted in the manuscript.
**Urgency verdict: circulation is time-sensitive on the order of weeks-to-months.**


**OpenAlex walk + Google Scholar phrase checks (2026-09-21).** 299 unique citing
works of six anchor papers screened via OpenAlex (Åberg 127, Albash 120, Kwon–Kim 64,
Funo–Murashita–Ueda proxy 43, Erdamar 8, Prech–Potts 4; the Murashita PRE-90 record
has an OpenAlex indexing gap, covered by the earlier Semantic Scholar sweep).
Whole-index exact-phrase searches (~250M works): zero hits for every C1-shaped
phrase. Google Scholar queried directly in a browser: zero hits for "postselected/
post-selected Crooks", "heralded fluctuation theorem", "post(-)selected fluctuation
theorem", "non-normal + Crooks + postselection". Two new ADJACENT works found and
now cited (Potts–Samuelsson PRL 121, 210603 — precursor of the Prech–Potts umbrella;
Salazar arXiv:2608.10118 — global dynamical-asymmetry bound for apparent second-law
violations). Published versions of 2504.11664 (PRB 113, 064304) and 2605.10099
(Front. Phys.) confirmed as already-dispositioned. **C1 survives; nothing
anticipates; urgency verdict unchanged.** Residual index risk is now limited to
very recent postings and non-indexed venues.


**Final sweep — orthogonal vocabulary (2026-09-24).** Targeted the blind spot of all
citation-tree/phrase sweeps: the same mathematics in other communities' language.
Searched heralded quantum optics, weak-value/arrow-of-time, pseudo-Hermitian,
quantum filtering/Belavkin, and structural phrase combinations (OpenAlex full text).
FINDING: the Dressel/Jordan measurement arrow-of-time lineage holds C1's two
ingredients SEPARATELY — Manikandan–Elouard–Jordan PRA 99, 022117 (boundary terms =
log record success probabilities, no post-selected TPM pairs) and Manikandan–Jordan
Quantum Stud. 6, 241 (post-selected forward/reverse ratios, unnormalized) — never
combined into per-pair renormalized Crooks deviations with the non-normality witness.
Both now positioned; PRA 99, 022117 and Crooks PRA 77, 034101 (canonical reverse-map
construction) added to the manuscript (metadata Crossref-verified). All other
categories cleared; Aug–Sep 2026 postings swept (a post-selected quantum-battery
paper 2609.25699 confirms the field is crowding but contains no FT). Cavina et al.
and Malakar–Silva have ZERO citers to date — nobody has yet unfolded the constant
into the ledger. **Search phase closed: three independent index families converge on
the same anchor set, and structural queries return our own anchors as top hits.
Residual risk is temporal only (future postings), handled by weekly monitoring of the anchor papers' citation feeds and the arXiv query set, plus the standard preprint-v2 mechanism.**
