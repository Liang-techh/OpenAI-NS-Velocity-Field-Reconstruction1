# CR-A9-049 — frozen winding-deficit radial attribution

## Purpose

Constrained Agent 9 performs one minimal evidence-reuse increment after CR-A9-048 and Agent 7 PR #437. The question is not whether another swirl basis can be invented, but **where the already-measured candidate-side winding deficit lives on the unchanged frozen material-path population**.

This increment consumes two immutable GitHub Actions receipts rather than integrating new trajectories:

- CR-A9-040 / PR #349, ST006 material paths, exact head `3911718884f021f0c0e3ba41ed694b9688764084`, workflow `35322722094`, artifact `10536949231`, artifact digest `sha256:e73165cbc4ae75e2d9b437ed1b8606a4af9c3e7439ffb50d65c315e33378c88d`, report SHA-256 `2595c8dafd0716b3b07cb9dc761b112e2178bf44200b9c50eff3c000c6166772`.
- CR-A9-048 / PR #435, ST048-S temporal-Piola `kappa=.05` material paths, exact head `6699a698c9fb2f0a0986e7aefc5fe0ee431cd5ad`, workflow `35366730180`, artifact `10556164980`, artifact digest `sha256:15488084eea04a0c4428bc07fbb87dad8dc6d11e719817715ba0c16314f644a0`, report SHA-256 `9463a67cb556ccd22ff24f776dca0af51cc3b7d3363e8f8b700a8fba434f0f6c`.

The analyzer fails closed unless both receipts retain exactly the same autonomous path contract: `t=.25 -> .75`, radii `.6/.9/1.2`, `z=+-.3`, 8 azimuths, 48 paths, 33 output times and the original DOP853 controls.

## External method record

No new external numerical method is introduced in this increment. The two upstream receipts were generated with:

- repository: `scipy/scipy`
- screened commit: `eff78058ed0d4cb8e1f0b5c5585d62d712e948d2`
- API: `scipy.integrate.solve_ivp / DOP853`
- license: BSD-3-Clause
- classification: **direct migration / public API only**
- copied implementation: none

CR-A9-049 itself only groups immutable receipt rows and computes arithmetic summaries. Its migration classification is therefore **direct internal evidence reuse**; no solver code, derivative operator or new dependency is copied.

## Why this does not duplicate Agent 7 PR #437

PR #437 tested a specific C-infinity radial ring multiplier that peaks near `r=sqrt(2)` and correctly stopped that direction after it failed the preregistered fixed-probe efficiency rule. It explicitly routed the next swirl-shape decision to Agent-9 trajectory evidence about which radial band actually needs winding.

CR-A9-049 supplies only that missing attribution. It does not screen another basis or widen #437's gain grid.

## Exact immutable-receipt result

On each radius there are 16 paths. Mean absolute turns are:

| seed radius | ST006 | ST048-S temporal Piola, kappa=0 | same field, kappa=.05 | kappa=.05 vs ST006 |
|---:|---:|---:|---:|---:|
| 0.6 | 0.04361853969 | 0.03087754686 | 0.03148588021 | -27.8154% |
| 0.9 | 0.02595649705 | 0.02266311677 | 0.02307026630 | -11.1195% |
| 1.2 | 0.003880821960 | 0.005279600732 | 0.005381022015 | +38.6568% |

The global `kappa=.05` degree raises mean winding by approximately `+1.97%`, `+1.80%`, and `+1.92%` on the three radii respectively: its effect is nearly uniform over this frozen seed population. The remaining net mean-winding deficit relative to ST006 is not uniform. The `r=.6` paths contribute about **89.75%** of the net deficit, `r=.9` about **21.35%**, while `r=1.2` contributes **-11.10%** because that outer group is already a winding surplus and partially offsets the inner deficit.

The same child remains inward on all 16 paths at every radius. Relative to its own `kappa=0` parent, contraction magnitude weakens by about `2.36% / 2.55% / 2.92%` at `r=.6/.9/1.2` respectively.

## Routing consequence

Under this **candidate-side** frozen contract, the useful next swirl-shape experiment should not spend additional energy on an outer ring around `r~1.2` or beyond. If Agent 7 performs another bounded shape experiment, it should first test an inner/mid redistribution direction that raises azimuthal travel around the `.6` seeds and secondarily `.9`, while explicitly measuring normalization cost and inward/axial tradeoffs.

This does not mean `r=.6` is an OpenAI-derived physical target. ST006 is only the retained repository candidate baseline, and these seed radii are autonomous diagnostics. The result is routing evidence for representation design, not a visual acceptance threshold.

## Constraint / truth boundary

Canonical values remain unchanged: `nu=.01`, physical `R^3`, registered box `[-2,2]^3`, smooth compact support `r<2, |z|<2`, `t in [.25,.75]`, restricted two-parameter forcing, `E(.25)=1+-0.001`, independent optimization/validation data, derivative ladder `.02/.01/.005`, divergence max/L2 `1e-5`, momentum max/L2 `1e-3`.

This increment performs no new trajectory integration, changes no velocity coefficient, fits no pressure/forcing, evaluates no PDE residual, defines no visual pass threshold and uses no OpenAI hidden time/camera/seed/numerical velocity. It keeps `production_radial_profile_selected=false`, `visualization_ready=false`, `visual_correspondence_verified=false`, `pde_validated=false`, `source_correspondence_verified=false`, `paper_exact=false`, `openai_field_identified=false`, and `blowup_proved=false`.
