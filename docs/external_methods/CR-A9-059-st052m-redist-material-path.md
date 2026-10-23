# CR-A9-059 — ST052-M frozen redistribution material-path replay

## Purpose

Replay the exact Agent-7 PR #528 ST052-M `kappa=.05` inner/mid swirl redistribution on Agent 9's unchanged 48-path / 24-pair material-trajectory contract. This is one routing increment toward a callable and visualization-ready `velocity(x,y,z,t) -> [u,v,w]`: determine whether the Eulerian winding-control transfer survives as real-time trajectory improvement before any new candidate identity, pressure reconstruction, restricted-force fit, or final export/render work.

## Frozen sources

- repository base: `main@f0193d66c9d92948b4820ebcb70263673995b324`
- ST052-M parent: PR #508 exact head `b3b8bfdbe1077f9ec967d158602951997d81e17d`
- Agent-7 transform: PR #528 exact head `779ffca71066e2864496d37de55a7aafc45d6f57`, task `CR003-ST052M-FROZEN-REDISTRIBUTION-078`
- transform identity: `alpha=2.520520814687742`, `kappa=.05`, inner window `(.30,1.05)`, outer window `(.95,1.85)`
- Agent-7 exact-head normalization receipt: `1.0032534663681094`
- frozen Agent-9 material-path engine: PR #435 exact head `6699a698c9fb2f0a0986e7aefc5fe0ee431cd5ad`

Agent 7's exact-head receipt is Eulerian expression-capacity evidence only. It records clean transfer, inner/mid angular-rate gain, outer redistribution, independent Cartesian-FD divergence, and local response rank; it explicitly records `material_paths_integrated=false` and `pde_validated=false`. This increment independently applies the frozen trajectory protocol rather than inheriting that conclusion.

## External method classification

`scipy/scipy@eff78058ed0d4cb8e1f0b5c5585d62d712e948d2`, `scipy.integrate.solve_ivp` with DOP853, BSD-3-Clause.

Classification: **direct migration / public API only through the frozen Agent-9 path engine**. No SciPy implementation is copied into the repository.

The Agent-7 transform itself is **direct internal method reuse / independent material-path replay**. Its radial profile, gain, normalization rule and source identity are frozen; this task performs no parameter scan.

## Frozen trajectory protocol

- physical time: `.25 -> .75`
- seed radii: `.6/.9/1.2`
- seed axial positions: `z=+-.3`
- eight azimuths per radius/side
- 48 trajectories / 24 axial pairs
- 33 output times
- DOP853: `rtol=1e-9`, `atol=1e-11`, `max_step=.01`

The same exact ST052-M field is the parent control. Results are descriptive routing evidence, not a visual or PDE acceptance rule.

## Constraint and truth boundary

Canonical physics remain unchanged: `nu=.01`, physical domain `R^3`, evaluation box `[-2,2]^3`, smooth zero extension outside `r<2, |z|<2`, time `[.25,.75]`, the preregistered restricted two-parameter forcing family, `E(.25)=1+-0.001`, independent optimization/validation sampling, divergence max/L2 `<=1e-5`, and momentum max/L2 `<=1e-3`.

This increment does not refit pressure or forcing, recompute a held-out full PDE residual, transfer ST052-M's PDE receipt to the transformed child, use `u->0`, define `f=R(u,p)`, relax a threshold, derive a numeric target from an OpenAI image/video, or claim hidden OpenAI time/camera/seed/velocity recovery.

Hard false remains `production_candidate_selected`, `production_redistribution_gain_selected`, `visualization_ready`, `visual_correspondence_verified`, `pde_validated`, `source_correspondence_verified`, `paper_exact`, `openai_field_identified`, and `blowup_proved`.
