# Work statistics of post-selected quantum dynamics
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22966846.svg)](https://doi.org/10.5281/zenodo.22966846)

Reproducibility package for:

> **"Work statistics of post-selected quantum dynamics: exact two-potential
> structure, an attribution no-go theorem, and pre-registered tests on quantum
> hardware"** — Andrew Matlock (2026). *Preprint; arXiv link to be added. Archived snapshot: [doi:10.5281/zenodo.22966846](https://doi.org/10.5281/zenodo.22966846)*

Post-selection reshapes accepted-run two-point-measurement work statistics in
an exact way: the pairwise Crooks deviations factorize as D_ij = u_j − v_i
into potentials that are operator norms of the effective dynamics, with
nonzero zero-work deviations iff the effective operator is non-normal. A
constructive no-go theorem shows these statistics cannot attribute the
post-selection to any particular selector. Both results, a forgery-class
taxonomy, and a work-statistics discriminator between the Deutsch and Lloyd
closed-timelike-curve models are demonstrated on IBM Heron-r2 hardware under
pre-registered, hash-committed criteria.

## Layout
- `paper/` — LaTeX source, bibliography, and the script that generates every
  figure directly from the raw data in `data/`
- `experiments/` — self-contained numerical and symbolic verifications
  (numpy/scipy/sympy; each prints the claim it supports; `exp9` proves the
  Deutsch–Lloyd separation in exact arithmetic)
- `hardware/` — runners and frozen analysis scripts for every hardware
  experiment (Qiskit; one job per registered run)
- `data/` — raw measurement counts for every run reported in the paper,
  including the runs that returned inconclusive or failed their gates, plus
  the complete recovered job archive
- `notes/` — all five pre-registrations with their SHA-256 files, the code
  and methods audit, the adversarial-review dispositions, and the prior-art
  report

## Pre-registrations
Every hardware claim was gated by criteria committed before data collection,
with the analysis script frozen. The SHA-256 of each registration is stored
alongside it and is bound to a git commit predating its run:

| Registration | Gates | Outcome |
|---|---|---|
| `PREREGISTRATION.json` (v1) | idealized-noise thresholds | Inconclusive (reported) |
| `PREREGISTRATION_v2.json` | audit-corrected, device-calibrated | Confirmed, both devices |
| `PREREGISTRATION_forgery.json` | control, absolute gates | Failed 2/3 (reported) |
| `PREREGISTRATION_control_v2.json` | control, class-gap gates | Failed 1/3 (reported) |
| `PREREGISTRATION_paired.json` | paired ratio gates, one job | Confirmed, all gates |

*Timestamp caveat, stated plainly:* this public repository was created after
those runs, so the pre-data commit timestamps for the experiments above are
attested by the authors' working repository rather than by a third party.
Registrations for all future runs will be pushed here before execution.

## Reproducing
```bash
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt
python experiments/exp4_quantum_factorization.py     # etc.
python hardware/analyze_paired.py data/paired_counts.json
cd paper && python make_figures.py                   # regenerates all figures
```
Hardware reruns require an IBM Quantum account; every runner submits a single
job and prints its job ID.

## Provenance
This work was produced collaboratively by the human author and an AI
assistant (Anthropic's Claude): conjectures, proofs, code, adversarial
review, and prose. The human author takes responsibility for all claims. The
audit trail of that process — including the findings of adversarial review
and every statistical correction between registrations — is in `notes/`.

## License
Code: MIT. Text, figures, and data: CC BY 4.0. See `LICENSE`.
