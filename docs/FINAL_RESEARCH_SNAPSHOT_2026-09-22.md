# Final research snapshot — 2026-09-22

This snapshot records the state of the Navier–Stokes reconstruction workspace when the scheduled NS research agents were paused on 2026-09-22. It is a consolidation record, not a new fit, proof, or validation result.

## Executive status

- Repository: `Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1`
- Main head immediately before this snapshot: `e6d7b1dea94bc36d61e56a230b0d1ac6c7a8be20`
- Scheduled NS reconstruction agents: paused by the project owner on 2026-09-22.
- Original full momentum target: `<= 1e-3`
- Current scientific status: **target not met**
- `pde_validated=false`
- `source_correspondence_verified=false`
- `paper_exact=false`
- `blowup_proved=false`

The repository contains several distinct research lines. They should not be collapsed into one “best field” without preserving their different objectives, validation samples, and provenance.

## 1. Constrained numerical candidate line

### ST006 — historical compatibility baseline

ST006 remains the backward-compatible Python API baseline. Its historical report records approximately:

- sampled momentum maximum: `0.1082289305`
- spatial volume L2: `0.1075843288`

It was an important early step but is not the latest numerical result and did not pass the original momentum target.

### ST061 — residual-oriented controls

ST061 retained two useful controls rather than hiding the tradeoff:

| Seed | Candidate | Sampled full-vector max | Spatial volume L2 |
|---|---|---:|---:|
| 9206291 | ST061-D | 0.0267880201 | 0.0333346590 |
| 9206291 | ST061-P | 0.0258078222 | 0.0334743264 |
| 9206292 | ST061-D | 0.0272044692 | 0.0342696874 |
| 9206292 | ST061-P | 0.0262079522 | 0.0344404411 |

Within those paired samples, D is the lower-L2 alternative and P is the lower-sampled-peak alternative. Neither passes the original `1e-3` momentum gate.

### ST063 — latest documented geometry line

ST063-G2R is the latest documented geometry experiment. It lengthens the moderate-strength axial rotation structure without simply graphically stretching the field.

| Seed | Candidate | Sampled full-vector max | Spatial volume L2 |
|---|---|---:|---:|
| 9216391 | ST061-P parent | 0.0262792551 | 0.0340356007 |
| 9216391 | ST063-G2R | 0.0243853533 | 0.0342574982 |
| 9216392 | ST061-P parent | 0.0256646731 | 0.0335885450 |
| 9216392 | ST063-G2R | 0.0225035591 | 0.0342060870 |

On these paired samples, G2R reduces sampled maxima by 7.21% and 12.32%, while volume L2 increases by 0.65% and 1.84%. This is a geometry/peak-residual tradeoff, not a universal improvement.

In the declared geometry window, the axial-to-single-transverse RMS aspect increased by about 3.5%–7.3% at the three principal checked times. The hoped-for 10%–20% improvement and a continuous high-threshold strong core were not achieved.

## 2. ST052 morphology/diagnostic line on main

The latest main history contains stable-identity-bound ST052 diagnostic work, including:

- stable streamline trajectory diagnostics;
- core contraction and speed-trend diagnostics;
- stable semantic identity/provenance binding for candidate-side measurements.

These additions improve the auditability of morphology and temporal behavior. They do **not** by themselves promote a new PDE-validated candidate. The corresponding main merge messages deliberately preserve that truth boundary.

## 3. Kokuno-derived reconstruction line

The Kokuno-derived work is preserved as a long, stacked PR/branch lineage rather than force-merged into `main`.

At the pause frontier, the newest integration PR is #1246, with adjacent active frontier PRs #1242–#1245. The integration record states that the latest genuinely self-contained project `u_lead + u_osc` candidate remains the A2 #1117 / A1 #1107 lineage through `xi=11`.

The current integration states remain:

- `leading_ready=false`
- `oscillatory_ready=true`
- `correction_ready=false`
- `velocity_export_ready=false`
- `pde_validated=false`

The remaining closure chain includes, at minimum:

1. resolve the legal exterior/terminal bridge geometry;
2. execute the selected exact runtime rebind and real RF30→RF39 correction path;
3. materialize RF44 post-update correction, Cartesian delta-u, nonlinear remainder and a finite correction cycle;
4. finish a global leading field with matched pressure/grad-p and deterministic Python/MATLAB export;
5. form a matching self-contained global oscillatory composite;
6. bind only preregistered non-residual-defined forcing;
7. independently validate the complete NS defect on held-out/canonical resolutions before any PDE promotion.

The fixed validation target remains normalized momentum max/L2 `<=1e-3`, with the separately governed divergence target and canonical multiresolution checks. Residual-defined/free forcing and post-hoc threshold relaxation remain forbidden.

## 4. Why the open PRs were not mass-merged

There are many stacked experimental PRs. A large fraction target other research branches rather than `main`, and later PRs depend on precise earlier heads. Blindly merging all of them would:

- destroy the experimental ancestry that makes the results auditable;
- mix mutually alternative or duplicate routes;
- risk promoting intermediate audit/routing artifacts as final candidate code;
- make it harder to distinguish software CI success from scientific validation.

Therefore the stop-state policy is: preserve all pushed branches/PRs as GitHub-hosted research artifacts, keep the curated validated/documented state on `main`, and do not reinterpret an open experimental branch as an accepted scientific result.

## 5. Artifact inventory

The repository catalog records three previously delivered complete bundles:

| Bundle | Size | SHA-256 | GitHub state at pause |
|---|---:|---|---|
| `NS_ST061_method_volume_progress.zip` | 35,286,603 bytes | `d11f078d26dd78dc5e1136850f9bb1716afe4e1e0cc3acfff9d6324a906eba37` | cataloged, not a GitHub release asset |
| `NS_ST063_axial_core_visual_progress.zip` | 64,668,298 bytes | `588009a4e42b8adad887b3198229080c7ad2efe57121a03c952245124b6d5196` | cataloged, not a GitHub release asset |
| `NS_ST063_MATLAB_Comparison.zip` | 1,850,963 bytes | `4f78f0930a7e7d1ff1921eccbafd8346da40850eb7bed46d0e1b2009f587343e` | cataloged, not a GitHub release asset |

Their hashes and provenance are retained in `docs/research_catalog.json`. The binary bytes were not available in the active workspace during this consolidation pass, so this snapshot does not falsely claim that those three ZIP files were uploaded as release assets.

## 6. What is fully preserved on GitHub

At this pause point GitHub preserves:

- the curated `main` source tree and documentation;
- ST006 compatibility artifacts and APIs;
- ST061/ST063 study records and exact catalog metadata;
- MATLAB/Python visualization source already present in the repository;
- ST052 diagnostic/provenance work merged to `main`;
- Kokuno and constrained-agent branch/PR lineages, including the latest pause frontier;
- commit, PR, CI and provenance histories needed to reconstruct how each research line evolved.

## 7. Scientific truth boundary

This project has produced materially better numerical candidates and substantially richer geometric/reconstruction machinery than the original ST006 baseline, but no current result should be described as an exact OpenAI field, a paper-exact reconstruction, a proof of Navier–Stokes blow-up, or a candidate that passed the original full `1e-3` momentum target.

For the curated numerical comparison, see `docs/RESEARCH_STATUS.md`. For candidate availability and hashes, see `docs/research_catalog.json`. For the runnable/current entry points, see `docs/CURRENT_CHECKPOINT.md`.
