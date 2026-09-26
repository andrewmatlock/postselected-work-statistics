# Adversarial Review — Findings and Dispositions (2026-09-20)
An adversarial-referee pass was run against the manuscript before submission,
instructed to produce the rejection letter now rather than after. Seventeen findings
resulted (2 rated fatal-as-framed, 9 major, 6 minor). This document records each
finding and its disposition, as part of the audit trail the paper references.
A parallel pass verified all 36 bibliography entries against arXiv/Crossref/journal
pages, and a third pass verified in-text characterizations of cited works.

| # | Finding (abridged) | Severity | Disposition |
|---|---|---|---|
| 1 | Ideal experiments cannot fail — "confirmed" overstates; only noise is at stake | FATAL (framing) | Reframed throughout: "demonstrate/characterize"; abstract rewritten; noise-fragility promoted as the genuinely contingent content |
| 2 | "Verified thermal preparation" false; biased marginal computed by reweighting, not observed | FATAL (framing) | Reweighting construction stated explicitly; "verified" removed; equivalence argument (classical mixture diagonal in measured basis) given in text |
| 3 | Central theorem had no written proof ("direct computation" + numerics) | FATAL | Appendix A written: recorded-circuit derivation, conditionals, cancellation, efficacy, edge cases; Halmos block and reverse-dilation argument added |
| 4 | Theorem 1(c) silently resolves the nonlinear-theory mixture ambiguity (BLSS) the paper flags for D-CTCs | MAJOR | Appendix A displays the per-branch alternative and argues for the global rule from Lloyd's definition; text notes (a),(b) are prescription-independent, (c),(d) global-rule |
| 5 | v2 registration calibrated on primary device; gates could scarcely fail; repo unnamed | MAJOR | fez arm labeled fully out-of-sample in text; abstract no longer implies risky gates; data-availability statement added with DOI placeholder (publication of repo pending, user action) |
| 6 | 182σ shot-noise-only figure indefensible; deep run predates registrations but abstract implied otherwise | MAJOR | 182σ demoted; deep run labeled retrospective; argument replaced by effect size: 0.70 acceptance span vs job-record readout errors 0.6–1.0% (≤0.01 shift) |
| 7 | Control failed its gate; pre-registration discipline applied asymmetrically | MAJOR | Control re-run under registered v3 (hash b63f6a5c…): failed 2/3 (noisier layout). Both attempts now reported as failures in abstract and text; control reframed as measured noise floor (|D_ii| ≤ 0.303 vs signal 1.02) |
| 8 | "Coherent errors cannot break factorization" wrong for separately-run forward/reverse circuits | MAJOR | Weakened to common-error statement; forward–reverse coherent mismatch acknowledged in the residual; logical-circuit basis of predictions stated (same D_th gates both devices) |
| 9 | Theorem 2 scope oversold (single-use TPM vs in-principle unattributability) | MAJOR | Statement and remark scoped to single-use TPM under the global prescription; BLSS-style multi-use/mixture discriminators acknowledged; anthropic rhetoric cut |
| 10 | "Two potentials ⇒ post-selection" stated as theorem; only converse proved | MAJOR | Taxonomy's third arrow marked as inference over tested ensembles; conclusion's "can establish" softened to "can display" |
| 11 | CTC banner + corollary-distance novelty invites desk rejection | MAJOR | Title and abstract lead with post-selected work statistics; CTCs presented as application |
| 12 | Cell error bars assume independence; VE undefined in paper; r robustness unshown | MAJOR | Appendix B: covariance quantified and bounded, VE defined; leave-one-out ranges and Spearman added to text; parameter count stated |
| 13 | Unstated hypotheses (degeneracy, C\|i⟩=0, basis, ΔF) | MINOR/MAJOR | Hypotheses added to Theorem 1 statement; edge cases in Appendix A |
| 14 | Number inconsistencies (13.8σ vs 17.8σ; attenuation coincidence; "four lineages"+5 cites; 1.21 vs 1.02) | MINOR | All labeled/reconciled; attenuation quoted per-run to 4 decimals (0.6931 marrakesh, 0.8414 fez — device-dependent, resolving the coincidence concern); "several lineages"; ideal-vs-measured labeled |
| 15 | No repo URL/DOI, no data availability, no appendices, bloated abstract | MAJOR | Appendices added; data-availability statement added; abstract halved; DOI insertion pending repo publication |
| 16 | Herald readout error unbudgeted | MINOR/MAJOR | Measured-qubit assignment errors pulled from job-linked calibration records (deep: 0.6–1.0%; v2 arms: 0.3–0.5%) and quoted with the bound |
| 17 | AI disclosure + unwritten proofs inverts burden of proof | MINOR | Neutralized by written appendices; disclosure retained |

## Bibliography verification (parallel pass)
All 36 entries checked against arXiv abs pages, APS/IOP pages, and Crossref:
28 verified, **8 wrong — all corrected**: fabricated co-authors (erdamar2024: three
nonexistent names replaced by Ha, Chen, Muldoon per PRR record), wrong first authors
(malakarsilva2025 → Manali Malakar, + published ref PRB 113, 064304 (2026);
zengyong2017 → Meng Zeng; maeda2023 → Kenji Maeda, Tharon Holdsworth;
bartkiewicz2019 → Małgorzata Bartkiewicz, Jacek Wychowaniec; cavina2026 → Frank
Ernesto Quintela Rodríguez), wrong titles (ferricortes2025, purvesshort2021).
Both legacy eprint-ID suspicions resolved in favor of the current bib
(1401.4494 and 0908.3023 correct). funo2017 gained its eprint (1412.5891).
A third pass re-verified the corrected fields letter-by-letter (all exact, incl.
Crossref given/family splits) and checked all 13 substantive in-text
characterizations of cited works against the sources: ALL VERIFIED, no corrections
required. Two soft flags were fixed: the composition/adaptive-use attribution to
BLSS was narrowed to what that paper supports (mixture inputs), and the
Horowitz-Maldacena citation was added at the final-state-projection mention.

## Residual risks that remain open (stated, not hidden)
1. The repository is not yet public; the pre-registration hashes are self-attested
   until a DOI/public host exists.
2. Literature completeness: CLOSED as a search problem 2026-09-24. Four sweep
   families (citation trees x2 indexes, whole-index phrases, direct Google
   Scholar, orthogonal-vocabulary communities) converge on one anchor set; the
   nearest structural kin (arrow-of-time lineage) found and cited. Residual risk
   is temporal only (future postings) and is handled by weekly monitoring of the anchor citation feeds plus the arXiv v2 mechanism.
3. No independent human domain expert has reviewed the manuscript. Adversarial
   AI passes reduce risk; they do not substitute for one. (2026-09-26: the
   author, informed of this risk, elected to proceed to submission without
   prior expert outreach; the risk transfers to arXiv moderation and eventual
   journal peer review, which is a legitimate — if less protected — path.)
4. Closed 2026-09-21: the Deutsch–Lloyd separation is now proven analytically
   for an explicit circuit (cycle invariant ln(5/3), exact, unique fixed points)
   under BOTH mixed-input conventions (exp9_deutsch_analytic.py; paper Sec. VI).
   Genericity across circuits remains numerical.
5. Closed 2026-09-21: the drift-immune paired design (both arms in one job,
   ratio gates, registration c789e667…) passed ALL registered gates
   (R_D 4.55 vs >=3, R_A 5.58 vs >=1.5; control arm cleanly single-potential,
   VE 0.986). The two failed absolute-gated attempts remain reported in the
   paper as the measurement of why absolute gates were the wrong instrument.
