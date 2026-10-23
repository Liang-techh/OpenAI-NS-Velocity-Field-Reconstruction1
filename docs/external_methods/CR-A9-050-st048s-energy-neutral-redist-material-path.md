# CR-A9-050 — ST048-S energy-neutral swirl redistribution material-path replay

## Scope

This is one visualization-side screening increment. It does **not** select a production field and does not alter the canonical PDE problem. Agent 7 PR #447 identified `kappa_redist=.025` as the first clean Eulerian proxy crossing for an inner/mid-weighted, first-order energy-neutral swirl redistribution on frozen ST048-S plus the existing temporal Piola schedule. CR-A9-050 asks the missing direct question: does that same fixed coordinate improve **actual cumulative material-path winding** on Agent 9's already frozen 48-path contract, while preserving inward contraction and axial material separation better than a global swirl gain?

## Frozen physical / governance contract

The repository contract remains unchanged: `nu=.01`, physical domain `R^3`, registered evaluation box `[-2,2]^3`, compact support connection `r<2, |z|<2`, time interval `[.25,.75]`, reference energy `E(.25)=1`, restricted preregistered two-parameter forcing, separate optimization/validation samples, and the original divergence `1e-5` / momentum `1e-3` acceptance thresholds. This increment neither fits nor transfers pressure/forcing and does not evaluate held-out PDE residuals.

Material paths retain the preregistered Agent-9 diagnostic contract without tuning after seeing this candidate: radii `.6/.9/1.2`, `z=±.3`, 8 azimuths, 48 paths / 24 axial pairs, `.25→.75`, 33 output samples, SciPy DOP853 with `rtol=1e-9`, `atol=1e-11`, `max_step=.01`.

## Migrated methods

**Agent-9 frozen path source:** PR #435 exact head `6699a698c9fb2f0a0986e7aefc5fe0ee431cd5ad`. Classification: **direct internal evidence reuse**. Only its frozen ST048-S loader, temporal-Piola transform, energy quadrature and material-path engine are reused; the earlier global-kappa result is not treated as acceptance truth.

**Agent-7 redistribution source:** PR #447 exact head `527e20e7eaa506d144c725d2c90c468da6aa5d34`, task `CR003-ST048S-ENERGY-NEUTRAL-SWIRL-REDISTRIBUTION-065`. Classification: **direct internal method reuse / independent reimplementation**.

Migrated scope only: C-infinity inner radial bump on `(.30,1.05)`, outer bump on `(.95,1.85)`, `h(r)=g_inner(r)-alpha*g_outer(r)` with `alpha` recomputed from the frozen candidate at `t=.25` so the first derivative of swirl kinetic energy vanishes, and the first clean Agent-7 crossing `kappa_redist=.025`. One positive common normalization restores `E(.25)=1`. The independent profile formula is cross-checked against the pinned Agent-7 source in CI; Agent-7's Eulerian proxy pass is not imported as truth.

**External source:** `scipy/scipy@eff78058ed0d4cb8e1f0b5c5585d62d712e948d2`, `scipy.integrate.solve_ivp / DOP853`, BSD-3-Clause. Classification: **direct migration / public API only** through the pinned path engine; no SciPy implementation is copied.

## Candidate provenance

Frozen ST048-S source: PR #390 exact head `97695a86f85ce68fb4ae70c41fc81c904d655183`, raw candidate SHA-256 `6d9ce8407e29aca20d32599332ae3ec654f45678e17783179267665428ec8e09`. Parent PDE evidence is not inherited by the transformed child.

## Truth boundary

The report is descriptive candidate-side kinematics. It hard-codes false states for production selection, visualization readiness, visual correspondence, PDE validation, source correspondence, paper exactness, OpenAI-field identification and blow-up proof. It also records that no hidden OpenAI time/camera/seed/numerical velocity, no public-image numeric target, no free `f=R(u,p)`, no threshold change and no visual acceptance threshold are used.

## Direct contribution to `[u,v,w]`

A positive result identifies a lower-cost azimuthal redistribution direction that allocates swirl capacity where the frozen trajectory deficit actually lies (`r≈.6-.9`) instead of spending normalization budget at the already-surplus outer diagnostic band. A negative result rules out Agent-7's first proxy crossing before pressure/forcing reconstruction effort is spent on it. Either outcome directly narrows the field family to materialize for MATLAB/Python visualization.

## Remaining limitation

Even a material-path improvement is not PDE evidence and is not proof of correspondence to OpenAI's hidden field. Any promoted child still needs compatible pressure/restricted forcing reconstruction and fresh independent held-out PDE validation under the unchanged gates.
