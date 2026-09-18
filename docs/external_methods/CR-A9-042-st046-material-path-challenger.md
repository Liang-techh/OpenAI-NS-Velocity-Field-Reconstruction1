# CR-A9-042 — ST046-A frozen material-path challenger replay

## Scope

This increment does **not** add a new visual metric. It reuses the frozen CR-A9-040 material-path contract to evaluate the first currently open constrained challenger that is sufficiently reconstructible to call as a numerical velocity field: `ST046-A` from PR #366 at exact head `04d2fa64798d9639bddadb289cc9c23f116eec56`.

The comparison is deliberately narrow: does ST046-A improve the two previously identified ST006 material-path weaknesses—small angular winding and spatially nonuniform axial stretching—under exactly the same registered seeds and time window? No camera, OpenAI streamline seeds, hidden frame time, image fit, or visual score is introduced.

## Frozen comparison contract

The replay preserves CR-A9-040 without retuning after seeing ST046-A:

- physical domain: `R^3`;
- registered evaluation box: `[-2,2]^3`;
- registered time interval: `0.25 -> 0.75`;
- initial radii: `0.6, 0.9, 1.2`;
- initial axial positions: `z=-0.3,+0.3`;
- 8 uniform azimuths per radius;
- 48 material paths / 24 paired axial material lines;
- 33 output times;
- SciPy `solve_ivp`, `DOP853`, `rtol=1e-9`, `atol=1e-11`, `max_step=0.01`.

Those solver settings are numerical integration choices, not visual or PDE acceptance thresholds. The ST006 reference is bound to CR-A9-040 head `3911718884f021f0c0e3ba41ed694b9688764084` and report SHA-256 `2595c8dafd0716b3b07cb9dc761b112e2178bf44200b9c50eff3c000c6166772`.

## Challenger provenance and classification

Candidate source is internal repository work, not an external method: PR #366, exact source head `04d2fa64798d9639bddadb289cc9c23f116eec56`, recipe `ST046-A`, original raw-candidate SHA-256 `94f5eeb0d94568587c9f3d88e69876a8051632830ded796f7c8b0db8e6707619`. The Agent 9 workflow checks out that exact commit read-only and invokes its published reconstruction recipe. PR #366 explicitly states that the readable recipe recreates the mathematical field/reference values while regenerated JSON metadata is not byte-identical to the archived original raw candidate; CR-A9-042 preserves that distinction.

For routing, ST046-A is classified **direct callable replay for candidate-side visual screening only**. It is not promoted to the default candidate and remains `pde_validated=false`. PR #366 reports material residual improvement relative to ST045 parents, but its original `1e-3` momentum gates still fail; this replay does not rerun or reinterpret those scientific gates.

## External method screening

External numerical method: `scipy/scipy@eff78058ed0d4cb8e1f0b5c5585d62d712e948d2`, `scipy.integrate.solve_ivp`, using the public `DOP853` API. Classification: **direct migration / public API only**. License: BSD-3-Clause. Migration scope is only time integration of `dX/dt = velocity(X,t)` under the frozen CR-A9-040 contract. No SciPy source code, RK tableau, adaptive controller, or dense-output implementation is copied into this repository.

Differences from the external library are therefore none at the algorithm-implementation level: this project only supplies the right-hand side, fixed initial seeds, fixed output times, and fixed tolerances. The physics field itself comes exclusively from the pinned ST046-A repository recipe.

## Truth boundary

This receipt is descriptive comparison evidence only. It must remain false for `visualization_ready`, `visual_correspondence_verified`, `pde_validated`, `paper_exact`, `openai_field_identified`, and `blowup_proved`. It defines no visual pass threshold and cannot be used to excuse a failed PDE gate.

Forbidden interpretations include: selecting ST046-A because it is visually closer and then calling it PDE-valid; recovering OpenAI hidden seed locations or frame times; rotating/registering the candidate to improve the score; changing seeds or integration window after seeing the result; treating the original raw SHA as the byte hash of the regenerated recipe output; or replacing the preregistered forcing with `f=R(u,p)`.

## Direct contribution to the final `[u,v,w]`

This is the first candidate-specific cross-check that sends a newer constrained field through the same material-path evidence contract already used for retained ST006. It therefore turns the previously qualitative routing question—whether a lower-residual challenger also has stronger inward spiraling and more coherent axial stretching—into a reproducible like-for-like measurement. The output is a deterministic JSON receipt and CI artifact that can be compared with the ST006 evidence card without inventing a new metric after the challenger is observed.
