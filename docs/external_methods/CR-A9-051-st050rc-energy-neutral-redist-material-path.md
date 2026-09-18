# CR-A9-051 — ST050R-C energy-neutral redistribution material paths

## Minimal increment

Replay exactly one already-screened visualization-side degree from Agent 7 PR #458 on Agent 9's unchanged frozen material-path contract. The parent is `ST050R-C`; the child multiplies cylindrical swirl by the preregistered energy-neutral inner/mid redistribution at `kappa_redist=0.025`, then applies the same one positive common normalization used by the source screen.

This round does not add another swirl basis, optimize a coefficient, fit an OpenAI image, modify pressure/forcing, or evaluate/transfer a held-out PDE residual.

## Frozen repository provenance

- repository: `Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1`
- branch base: `main@f0193d66c9d92948b4820ebcb70263673995b324`
- ST050R-C parent: PR #448 exact head `dc1d5476ea9979566b774293add76196083288c5`
- Agent-7 transfer screen: PR #458 exact head `2ed651cb750f347fd009d3c511645627135192a4`, task `CR003-ST050RC-ENERGY-NEUTRAL-SWIRL-TRANSFER-071`
- frozen Agent-9 material-path engine: PR #435 exact head `6699a698c9fb2f0a0986e7aefc5fe0ee431cd5ad`

Agent 7's frozen radial windows are inner `(0.30,1.05)` and outer `(0.95,1.85)`. Its transfer screen computed the balance coefficient only from frozen ST050R-C at `t=0.25`, preregistered gains `[0,.015,.025,.035,.05]`, and found the first clean Eulerian crossing at `.025`. This Agent-9 round does not widen or retune that grid.

## External method screening / migration

### SciPy DOP853

- source repository: `scipy/scipy`
- pinned source commit inherited from the frozen Agent-9 path engine: `eff78058ed0d4cb8e1f0b5c5585d62d712e948d2`
- API: `scipy.integrate.solve_ivp(..., method="DOP853")`
- license: BSD-3-Clause
- classification: **direct migration / public API only**
- copied implementation: none
- migrated scope: numerical integration of the already-defined material ODE only
- difference from upstream: none at the API level; the seed/time/observable contract is repository-defined

### Agent-7 ST050R-C redistribution

- source repository: this repository
- source PR/head: #458 / `2ed651cb750f347fd009d3c511645627135192a4`
- classification: **direct internal method reuse**
- migrated scope: exact ST050R-C replay, inner/outer compact bumps, first-order energy-neutral balance, `.025` swirl redistribution, common reference-energy normalization
- difference: Agent 7 evaluates target-free Eulerian angular/vorticity proxies; CR-A9-051 evaluates actual cumulative material trajectories with the frozen Agent-9 contract
- external-license migration: none; repository-local source

No OpenAI image/video numerical values, hidden frame time, camera pose, seed locations, or velocity data are consumed in this round.

## Frozen trajectory contract

Unchanged from Agent 9 PR #435 and subsequent replays:

- `t=0.25 -> 0.75`
- initial radii `.6/.9/1.2`
- initial `z=+/-0.3`
- 8 azimuths
- 48 paths / 24 axial pairs
- 33 output samples
- DOP853, `rtol=1e-9`, `atol=1e-11`, `max_step=.01`
- registered evaluation box `[-2,2]^3`

The implementation fails closed if the exact source checkout, path contract, parent lineage, redistribution windows, or preregistered gain grid drift.

## Constraint governance / truth boundary

Canonical CR001 values remain unchanged: `nu=.01`, physical domain `R^3`, evaluation box `[-2,2]^3`, smooth zero extension outside `r<2, |z|<2`, time `[.25,.75]`, restricted two-parameter forcing only, `E(.25)=1 +/- .001`, separate optimization/validation samples, FD ladder `.02/.01/.005`, divergence max/L2 `1e-5`, momentum max/L2 `1e-3`, and retain-failure/no-post-hoc-threshold policy.

Hard false in this increment: canonical velocity changed, production candidate selected, production redistribution gain selected, pressure/forcing changed or transferred, held-out PDE residual evaluated or transferred, hidden OpenAI data used, OpenAI numeric target inferred, visual pass threshold defined, `visualization_ready`, `visual_correspondence_verified`, `pde_validated`, `source_correspondence_verified`, `paper_exact`, `openai_field_identified`, `blowup_proved`.

The ST006 comparison, if reported, is descriptive repository-baseline context only. ST006 is not OpenAI truth and is not a visual acceptance target.

## Direct contribution to final `[u,v,w]`

ST050R-C currently offers a lower held-out momentum-residual challenger than ST048-S on its frozen random samples while retaining useful morphology. PR #458 shows the energy-neutral redistribution transfers cleanly at the Eulerian level but explicitly leaves material paths to Agent 9. This replay decides whether that lower-residual backbone also gains the real inner/mid cumulative winding needed for the visualization route before anyone spends work on a compatible pressure/restricted-forcing refit or 3-D render promotion.
