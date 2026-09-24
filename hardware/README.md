# Hardware protocols

Because the attribution no-go guarantees a heralded circuit reproduces P-CTC statistics
exactly, every prediction here is testable on present-day quantum processors.

| Script | Purpose |
|---|---|
| `run_prereg.py` | Runner for both registrations: shallow 4-CX, 16 cells, 20k shots/input, all 8 circuits as ONE job. v1 executed 2026-09-20 (job dao3sfo2fm4c73f4tkmg, INCONCLUSIVE by 0.004). |
| `analyze_v2.py` | **The v2 analysis** (audit-corrected statistics): proper biased-past errors, uniformity χ², attenuation coefficient; factorization gated on VE + r(D,D_th) + predicted potentials. Criteria in `notes/PREREGISTRATION_v2.json`. |
| `ctc_experiment.py` | Deep-dilation protocol (Haar operator + Halmos dilation). Establishes the noise threshold: the fingerprint dies at ~15% depolarizing. |
| `ctc_shallow.py` | Noise-robust redesign — choose a shallow circuit (4 CX) and read the operator off its heralded sub-block. ~25x better signal-to-noise. |
| `run_highstat.py` | Batched runner (superseded: burns quota on per-job overhead) |
| `salvage_jobs.py` | Retrieve completed job results by ID (free; consumes no QPU). Used 2026-09-20 to recover the full July dataset after account access was restored → `data/ctc_shallow_salvaged.json`. |
| `analyze_counts.py`, `analyze_weighted.py`, `analyze_shallow.py` | Analysis with shot-noise error bars, cell filtering, partial-data tolerance |

**Standing rule (learned the hard way).** Raw counts JSONs are copied into `data/` in the
same session they are measured, and committed immediately. The prior session's 12/16-cell
hardware dataset was lost when IBM deleted the inactive instance — job records and all.

**Restoring hardware access (manual, ~2 min).** Log in at quantum.cloud.ibm.com; create
(or confirm) an Open Plan instance; if the saved key still fails, generate a fresh API key
and run `QiskitRuntimeService.save_account(channel="ibm_quantum_platform", token=..., overwrite=True)`.
Then: `python run_prereg.py hardware`.

**Budget note.** IBM Open Plan: 10 minutes per 28-day window. Per-job overhead, not shot
count, dominates — submit all eight circuits in a single `sampler.run([...])` call.

**Measured to date** (ibm_marrakesh, Heron r2): biased past confirmed at 114.6 sigma;
factorization accounts for 98% of deviation-matrix variance on the shallow circuit
(12/16 cells, partial run); fingerprint destroyed by ~15% incoherent noise on the deep circuit.
