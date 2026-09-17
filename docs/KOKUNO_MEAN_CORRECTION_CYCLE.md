# Kokuno Agent 3 — first bounded mean-correction cycle screen

## Scope

This increment takes the real phase-mean defect contract from PR #215 and the compact theta/e=2 radial-stress reconstruction from PR #223, constructs one bounded public mean-velocity correction family, selects one scalar amplitude on held-in data only, and recomputes the real defect on disjoint held-out samples. It is intentionally one finite step, not an infinite iteration and not a PDE/blow-up claim.

Implementation: `openai_ns_reconstruction.kokuno_mean_correction_cycle`.

## Source provenance and a necessary boundary

Source: KokunoYumeto, *Forced Navier–Stokes blowup: reconstruction and validation reader*, corrected 208-page edition, version `2026.09.09-consolidated`, DOI `10.5281/zenodo.22678406`.

The reader's compact radial inverse reconstructs `sigma_e` from the auxiliary-averaged tangential source by weighted moment subtraction. In the finite correction cycle it then sends the physical radial stress through a **signed covariance-amplitude map**. The source does **not** state `delta u = sigma` for a mean velocity.

That distinction matters here. Agent 2 owns the oscillatory complete-curl/amplitude lane and PR #230 is already screening that lane. This Agent-3 increment therefore does not duplicate it. Instead, to make one falsifiable mean-correction experiment, it uses the already reconstructed theta/e=2 stress only as a fixed radial shape for an autonomous axisymmetric pure-swirl mean velocity:

`delta u_mean = a * S(r) * chi_z(z) * e_theta`,

where `S` is the normalized compact stress profile and `chi_z=(1-(z/H)^2)^5_+` is a fixed C4 axial gate. This is a repository-local engineering lift and is explicitly `source_exact_velocity_lift=false`.

Because the field is axisymmetric and purely azimuthal, it is analytically divergence-free away from the axis; its radial support is `[0.12,0.72]`, so the axis is excluded. An independent Cartesian finite-difference divergence audit is still reported.

## Typed real-defect cycle

The shape input is not hand-entered. It is rebuilt from the current serialized candidate and Agent-2 correction through `PhaseMeanDefectContract`, cylindrical ring averaging, the fixed annular gate, and `CompactRadialStressInverse`.

The scalar amplitude is selected from the fixed nonzero grid

`[-.04,-.02,-.01,-.005,.005,.01,.02,.04]`.

Zero is the baseline and cannot be relabeled as a successful nontrivial correction. Selection uses seed `9172813`, 64 cylindrical support-interior points, and only `t=.5`. The selected amplitude is frozen before validation.

Validation uses disjoint seed `9172819`, 128 new points, and `t=.375,.5,.625`. For each time the report records the nonlinear phase-mean momentum increment relative to the **original** base operator before/after, theta component norms, raw phase-mean zero-pressure/zero-force operator norms, correction norm, and composite divergence. A third seed `9172829` is used for an independent Cartesian divergence audit of the correction itself.

The report decision is deliberately stricter than an optimization receipt: the nonzero amplitude must improve held-in mean-defect RMS, improve aggregate disjoint held-out mean-defect RMS, keep the held-out raw phase-mean operator from worsening by more than 0.1%, and keep the independent correction divergence max below `1e-5`. If it fails, the attempt remains recorded and the grid is not widened in this run.

## Correction-cycle gain

The reader's stage audit reports exponent-class gains `0.39999` (wave residual), `0.17` (tangential mean), and `0.89996` (defects) for its own hierarchy. They are stored only as provenance/diagnostic values. They are **not** interpreted as numerical contraction factors and are not used as pass thresholds for this repository.

## Reproduce

`python -m openai_ns_reconstruction.kokuno_mean_correction_cycle --radial-count 65 --angular-count 16`

Output: `artifacts/kokuno_agent3/mean_correction_cycle_report.json`.

## Truth boundary

- consumes a real serialized candidate and real measured phase-mean defect;
- no surrogate defect and no free forcing cancellation;
- one public nonzero mean-velocity correction is actually evaluated;
- held-in and held-out samples are disjoint;
- stress-to-velocity lift is autonomous, not source exact;
- no full pressure/force fit is performed;
- the normalized full-NS `1e-3` gate is not assessed;
- `pde_validated=false`, `paper_exact=false`, `openai_field_identified=false`, `blowup_proved=false`.
