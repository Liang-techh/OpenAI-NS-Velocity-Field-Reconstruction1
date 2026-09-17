# Kokuno Agent 3 — compact radial-stress reconstruction

## Scope of this increment

This increment implements exactly one new correction primitive: the source moment-complement radial inverse, applied to a radial projection of Agent 3's already-measured **real phase-mean NS defect**. It does not yet convert the reconstructed stress into a public velocity increment and does not run a finite correction cycle.

The implementation is `openai_ns_reconstruction.kokuno_radial_stress`.

## Public source / provenance

Source: KokunoYumeto, *Forced Navier–Stokes blowup: reconstruction and validation reader*, corrected 208-page edition, version `2026.09.09-consolidated`, DOI `10.5281/zenodo.22678406`.

The corrected reader's radial-stress audit records, for `e=2` in the angular residual and `e=1` in the axial residual,

`M_e F = integral_0^infinity r^e F(r) dr`,

`P_e F = F - b_e M_e F`, with `integral r^e b_e dr = 1`,

`T_e F(r) = -r^(-e) integral_0^r s^e (P_e F)(s) ds`,

and therefore

`(d/dr + e/r) T_e F = -F + b_e M_e F`.

The same audit states the support reason explicitly: below the annulus the primitive is zero, and above the annulus it is the weighted moment of `P_e F`, which is zero. Without the equal-moment subtraction an arbitrary nonzero moment leaves an `r^(-e)` tail. The source support space is smooth and annular, bounded away from the axis, and the source moving profile is scaled by its positive `q(z,t)`.

The finite-cycle chapter uses the same physical interface with `F_2` equal to the auxiliary-averaged angular defect and `F_1` equal to the auxiliary-averaged axial defect. The source reconstruction cancels the supported averaged source while retaining the positive bump/moment term. None of the source exponent gains are treated here as measured numerical contraction.

## Repository-local adapter and its limitation

Agent 2's current executable oscillatory field is deliberately autonomous: it uses a Cartesian affine phase and a compact Cartesian box envelope rather than Kokuno's full moving annulus/auxiliary coordinates. Therefore the source radial inverse cannot honestly be applied as though the current field already lived in the source annular hierarchy.

This increment makes the adaptation explicit and fixed before reading residual values:

1. evaluate `PhaseMeanDefectContract` on the actual public base candidate plus Agent-2 correction;
2. sample deterministic physical cylindrical rings;
3. take uniform physical-angle means of the measured `e_theta` and `e_z` defect components;
4. multiply them by a fixed C4 annular window `(1-s^2)^5_+` on `r in [0.12, 0.72]`;
5. apply the source moment-complement inverse separately to `e=2/theta` and `e=1/z`.

The annular gate is an **autonomous engineering adapter**, not a fitted defect shape. Its raw-vs-gated RMS capture ratio is reported so the gated source cannot be relabeled as the whole mean defect. The current increment does not consume Agent 1 PR #221's native `q/X/eta` bridge because Agent 1 still does not expose a complete 3-D Kokuno leading velocity. `moving_q_scale_used=false` is recorded in the report.

## Executable numerical contract

`RadialMeanDefectProfile` keeps both the raw ring-average values and the gated compact source, source channel, exponent, physical annulus, time and angular quadrature. The scientific report path constructs it only from the real `PhaseMeanDefectContract`; no hand-entered residual tensor is accepted as success evidence.

`CompactRadialStressInverse` builds the weighted moment, normalized interior bump and stress. Independent centered finite differences of the public `stress(r)` evaluator check

`(d/dr + e/r) sigma = -F + b_e M_e F`.

The report also records:

- weighted source moment `M_e`;
- weighted moment of the complement;
- stress RMS/max;
- inner/outer/outside support values;
- the tail that the same sampled source would have produced **without** moment subtraction, `|M_e|/r_outer^e`;
- independent FD identity RMS/max;
- raw-to-gated defect capture.

## Real artifact path

The command

`python -m openai_ns_reconstruction.kokuno_radial_stress --radial-count 65 --angular-count 16`

loads `artifacts/bipolar_joint_capped/candidate.json`, uses the existing Agent-2 `KokunoCompleteCurlCorrection`, and computes the phase-mean defect through the independent CR006 residual operator. It then reconstructs theta/e=2 and axial/e=1 compact stresses. The generated report is `artifacts/kokuno_agent3/radial_stress_report.json`.

This is a correction **primitive** only. A small radial-identity error proves that the discrete stress implements the intended inverse on the gated measured defect. It does not prove that adding a later velocity correction will lower the full NS residual.

## Next finite-cycle seam

The next Agent-3 increment may use this stress receipt to construct one bounded public velocity correction, then recompute the real NS defect and compare held-in and disjoint held-out residuals before/after. That step must record correction norm, divergence, support, nontriviality and observed contraction. It must reject a correction if held-out residual grows materially. The reader's source exponent gains remain provenance/diagnostic information only; they are not numerical gain factors for this repository.

## Truth boundary

- `public_reconstruction_source = true`
- `consumes_real_phase_mean_defect = true`
- `surrogate_defect_used = false`
- `radial_stress_inverse_applied_to_measured_defect_projection = true`
- `mean_velocity_correction_applied = false`
- `finite_correction_cycle_run = false`
- `residual_reduction_claimed = false`
- `formal_pde_gate_assessed = false`
- `pde_validated = false`
- `paper_exact = false`
- `openai_field_identified = false`
- `blowup_proved = false`
