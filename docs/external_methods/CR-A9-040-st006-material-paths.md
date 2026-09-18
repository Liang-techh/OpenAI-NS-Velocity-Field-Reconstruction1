# CR-A9-040 — retained ST006 material-path observables

## Task and direct delivery purpose

Constrained Agent 9 uses the already published ST006 callable field to produce one deterministic, candidate-side material-path receipt for the public qualitative observables **inward spiraling** and **axial stretching**. This is not another generic pathline implementation: open PR #45 already owns that reusable primitive. The present increment binds a fixed seed/time contract directly to `research_baseline.load_best().at_points` so the current production baseline can be measured now.

## External method screened

- source repository: `scipy/scipy`
- screened commit: `eff78058ed0d4cb8e1f0b5c5585d62d712e948d2`
- public API: `scipy.integrate.solve_ivp`
- upstream file: `scipy/integrate/_ivp/ivp.py`
- license: SciPy BSD 3-clause terms
- classification: **direct migration / public API only**
- copied upstream implementation: none
- dependency delta: none (`scipy>=1.10,<2` is already a project dependency)

The upstream public API solves `dy/dt=f(t,y)` over a declared interval and supports fixed output times through `t_eval`. It exposes `DOP853`, an explicit high-order Runge–Kutta method suitable for high-accuracy nonstiff integration. This increment calls only that public API; no RK tableau, error controller, interpolation implementation, or solver source is copied.

## Frozen autonomous measurement contract

Before inspecting ST006 pathline output, the following candidate-side contract is fixed:

- time interval: exactly `[0.25,0.75]`;
- 33 equally spaced reported times;
- seed radii: `0.6, 0.9, 1.2`;
- seed axial coordinates: `z=-0.3,+0.3`;
- 8 endpoint-excluded uniform azimuths per radius;
- 48 total seeds and 24 matched `(-z,+z)` material-line pairs;
- solver: `DOP853`, `rtol=1e-9`, `atol=1e-11`, `max_step=.01`;
- registered evaluation box: `[-2,2]^3`.

Seed locations and solver settings are repository-autonomous diagnostic choices. They are not OpenAI streamline seeds, hidden frame times, source constants, visual acceptance thresholds, or PDE tolerances.

For each path the receipt records cylindrical-radius change, unwrapped angular travel/turns, absolute-z change, speed change and sampled path length. For every matched `z=-.3/+ .3` pair it records material-line axial-separation change and ratio. The two public-observable evidence entries therefore remain purely candidate-side:

- `inward_spiraling`: radius evolution together with angular travel;
- `axial_stretching`: matched material-line axial-separation evolution, supplemented by individual `|z|` evolution.

No numerical OpenAI target or pass/fail morphology threshold is defined.

## Constraint governance

The canonical project contract is unchanged: `nu=.01`; physical domain `R^3`; evaluation box `[-2,2]^3`; smooth zero extension outside `r<2, |z|<2`; time `[.25,.75]`; preregistered restricted two-parameter forcing only; `E(.25)=1±.001`; separated optimization/validation sampling; 4096 held-out validation points; derivative ladder `.02/.01/.005`; divergence max/L2 `1e-5`; momentum residual max/L2 `1e-3`; failed results retained; any threshold change requires a new experiment version.

The path integration evaluates no pressure, forcing, spatial/time PDE derivative, or residual. Solver tolerances are numerical integration settings only and cannot replace or modify the registered PDE gates.

## Truth boundary

The receipt hard-keeps all scientific/identity promotion states false: `visualization_ready`, `visual_correspondence_verified`, `pde_validated`, `paper_exact`, `openai_field_identified`, and `blowup_proved`.

It also records that no OpenAI seed locations, camera registration, image fit, hidden-time alignment, or hidden parameters are used. Even if ST006 material paths move inward while rotating and paired particles separate axially, that is candidate-side qualitative evidence only; it is not a numerical recovery of the OpenAI field and cannot override ST006's failed held-out NS residual/divergence gates.

## Remaining limitation

This diagnostic samples one fixed autonomous seed ensemble inside the retained ST006 field. It is not a uniform statement over all material trajectories and it does not establish quantitative OpenAI correspondence. The same frozen seed/time contract can later be rerun on a promoted winning candidate for a directly comparable visualization-side receipt, while PDE validation remains an independent gate.
