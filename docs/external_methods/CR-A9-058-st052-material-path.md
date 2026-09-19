# CR-A9-058 — ST052-M frozen material-path replay

## Scope

One minimal Agent-9 routing increment: replay the exact frozen `ST052-M` field from root PR #508 and its immutable `ST051-B` parent on Agent 9's unchanged 48-path / 24-pair material-trajectory contract. The goal is to determine whether the sampled-minimax residual-improving child preserves or damages the candidate-side inward-spiral trajectory geometry before further visualization/export work.

This increment does not fit coefficients, add a basis, change pressure or forcing, rerun the project PDE acceptance gate, or define a visual pass threshold.

## Frozen internal provenance

- base branch point: `main@f0193d66c9d92948b4820ebcb70263673995b324`
- candidate source: PR #508, exact head `b3b8bfdbe1077f9ec967d158602951997d81e17d`
- candidate id: `ST052-M`
- parent id: `ST051-B`
- parent raw SHA-256: `0071106ef10a5d77b620b942becc00b58c8a6765fd5fb7bc4dcd195ac65c970d`
- PR #508 reported frozen candidate raw SHA-256: `e078e753fab38ebfa0284d28ba64d26cb8538b700a3afc7705668849c43e12da`
- frozen Agent-9 material-path engine: PR #435, exact head `6699a698c9fb2f0a0986e7aefc5fe0ee431cd5ad`

The #508 recipe and callable reconstruction are consumed directly. Their existing PDE receipt is recorded only as upstream context and is not promoted by this increment: the reported seed-9175291 momentum sampled max is `0.03681687198184801`, volume-L2 is `0.05179943624413842`, and both remain above the unchanged `1e-3` gate.

## Frozen path protocol

- physical time: `0.25 -> 0.75`
- seed radii: `0.6, 0.9, 1.2`
- seed z: `-0.3, +0.3`
- eight azimuths per radius/z pair
- 48 paths / 24 axial pairs
- 33 output times
- DOP853, `rtol=1e-9`, `atol=1e-11`, `max_step=.01`

No seed, interval, solver tolerance, or interpretation threshold is selected after seeing the ST052 result.

## External method classification

`scipy/scipy@eff78058ed0d4cb8e1f0b5c5585d62d712e948d2`, `scipy.integrate.solve_ivp` with DOP853, BSD-3-Clause.

Classification: **direct migration / public API only through the frozen Agent-9 path engine**. No SciPy integrator implementation is copied into this repository.

The ST052 callable/recipe and the Agent-9 path measurement are **direct internal candidate replay / independent frozen-path measurement**. PR #508 already contains a different 36-seed trajectory diagnostic; this increment does not claim that result as an Agent-9 replay and instead uses the unchanged cross-candidate 48-path contract.

## Constraint and truth boundary

Canonical project physics remain unchanged: `nu=.01`, physical domain `R^3`, evaluation box `[-2,2]^3`, compact support connected to `r<2, |z|<2`, time `[.25,.75]`, restricted two-parameter forcing, nontrivial `E(.25)=1+-0.001`, independent optimization/validation sampling, and original `1e-5` divergence / `1e-3` momentum gates.

Hard false remains: `production_candidate_selected`, `visualization_ready`, `visual_correspondence_verified`, `pde_validated`, `source_correspondence_verified`, `paper_exact`, `openai_field_identified`, and `blowup_proved`. No `u->0`, free residual-defined `f=R(u,p)`, threshold relaxation, hidden OpenAI numerical field/time/camera/seed, image-derived numeric target, parent PDE receipt transfer, or visual acceptance score is introduced.

## Direct contribution to final velocity

The replay supplies one apples-to-apples kinematic check on a candidate that made substantial sampled momentum-max progress. A favorable result would justify carrying ST052-M farther into governed export/render comparison; an unfavorable result would identify a PDE-vs-trajectory tradeoff before spending more visualization or representation budget. Either outcome remains descriptive candidate-side evidence, not OpenAI correspondence or PDE acceptance.
