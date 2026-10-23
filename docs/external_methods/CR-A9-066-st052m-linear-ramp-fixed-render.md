# CR-A9-066 — ST052-M linear-ramp fixed render

## Scope

One target-free visualization replay only. The velocity is Agent-7 PR #587 exact head `0b93819095f6c8576a7571bdc2d2fbef4154944d`: ST052-M plus the frozen `kappa=.05` redistribution and the exact #559 taper/shoulder pair activated by the single autonomous profile `g(t)=2*(t-.25)`. No coefficient, window, basis, pressure, forcing, scientific threshold, seed, camera, or image target is fitted here.

The comparison reuses Agent-9 CR-A9-064 exact head `ee80de11775bd45157d114000378b6ee3422ff46` for the fixed 3-D protocol and evaluates `t=.25/.375/.50/.625/.75` on the same `33^3` box grid. The ramp is required by representation identity to equal the redistributed control at `t=.25` and the static #559 child at `t=.75`; the `1e-11` numerical equality guard is an implementation identity check, not a visual/PDE acceptance gate.

## Sources and classification

- `Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1` PR #587, exact head `0b93819095f6c8576a7571bdc2d2fbef4154944d`: **direct internal reuse**, temporal child under test.
- Same repository PR #594, exact head `b8d45d5d04a03e57739b0b3847c0994f4ade8642`: **supporting internal evidence only**, resolution-stable threshold-free time-evolution receipt; it is not used as a numerical oracle.
- Same repository PR #583, exact head `ee80de11775bd45157d114000378b6ee3422ff46`: **direct internal reuse**, frozen render/grid/streamline/vorticity implementation.
- SciPy `RegularGridInterpolator`, pinned by CR-A9-064 to `scipy/scipy@b12c772edbc1fe0d3db9481cdfcc2e311569cb25`, BSD-3-Clause: **direct migration / public API only**.
- Matplotlib 3.10.6, CR-A9-064 source commit `5cd38c3edcdf0792d0e6aded280a9b7a7de6146f`, Matplotlib License Agreement: **direct migration / public API only**.

No third-party implementation is copied into this increment.

## Frozen observable protocol

- times `.25/.375/.50/.625/.75`;
- Cartesian `33^3` velocity grids over `[-2,2]^3`;
- 48 fixed central seeds inherited from CR-A9-064 (`r=.6/.9/1.2`, `z=+-.3`, eight azimuths);
- deterministic bidirectional normalized-velocity RK4 streamlines and fixed camera inherited unchanged;
- evaluator-only Cartesian finite-difference vorticity;
- top `1.5% |omega|` point cloud for image display only;
- threshold-free enstrophy axial/radial RMS, aspect, fixed tip-band radial RMS, and central-band radial RMS retained for descriptive routing;
- static #559/control metrics are evaluated on the same grids to report the ramp's fraction of the static morphology effect.

There is deliberately **no visual pass threshold** and no OpenAI image-derived numerical target.

## Truth boundary

This replay cannot select a production child or promote `visualization_ready`, visual correspondence, pressure/forcing compatibility, held-out NS momentum validity, source correspondence, paper exactness, OpenAI-field identity, or blow-up. The CR001 domain/viscosity/support/time/forcing/nontriviality/train-validation split and `1e-5` divergence / `1e-3` momentum gates are unchanged. No `u->0`, no residual-defined free `f=R(u,p)`, no threshold relaxation, and no parent PDE receipt transfer are permitted.
