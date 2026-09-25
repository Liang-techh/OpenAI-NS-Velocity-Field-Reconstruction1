# Localized solenoidal poloidal repair

The largest current momentum residuals lie in the axial transition
collar and are predominantly radial. A single prior poloidal mode was
ineffective. This experiment adds eight independent compact
streamfunctions with radial Legendre degrees `0,1`, axial degrees
`0,1`, and even/odd parity across `z=0`. Each mode uses
`psi=15 r^2 (1-(r/R)^2)^5 P_i(2(r/R)^2-1) * 1024 s^5(1-s)^5
P_j(2s-1)` within the collar, optionally multiplied by `sign(z)`.
The velocity is `u_r=-psi_z/r`, `u_z=psi_r/r`, `u_theta=0`, so it is
exactly divergence-free where smooth and remains compactly supported.
The implementation evaluates these formulas directly. A finite-
difference divergence sanity check at interior points decreases by
about fourfold when the step halves, consistent with second-order
truncation.

An unconstrained eight-mode least-squares fit made Gauss-8 L2 look
much better (`129.66` to `93.18`), but worsened Gauss-10 (`314.44`
to `326.34`) and its maximum (`269,257` to `335,294`). This is
**not a robust improvement**. A second optimization constrained each
Gauss-6/8/10 L2 and maximum not to exceed the incoming candidate.
It found a nonzero direction with L2 reductions of about `1.33%`,
`1.10%`, and `0.91%` while maxima stayed at their previous levels.
At full strength that direction fails three edge cone points near
`(r,z)=(0.012,0.0037)` at later times.

At **10%** of that direction, direct full momentum recomputation at
`tau=0.0084` gives:

| Grid | Incoming L2 | 10% L2 | Incoming maximum | 10% maximum |
| --- | ---: | ---: | ---: | ---: |
| Gauss-6 | 159.900 | 159.607 | 94,068.012 | 94,067.912 |
| Gauss-8 | 129.662 | 129.445 | 177,827.212 | 177,827.161 |
| Gauss-10 | 314.442 | 314.095 | 269,256.836 | 269,256.834 |

The local cone is again **105/105 passing** across the existing five
sampled times, with maximum sampled ratio `0.794619`. This is a
small but simultaneous finite-sample improvement in both metrics,
with a nonzero exact-solenoidal velocity correction. The maxima
change only marginally, and the spatial quadrature orders remain
strongly inconsistent. Nothing here establishes the required
`1e-3` bounds, a continuum cone, time-uniform momentum control,
pressure compatibility, or the paper's supported nonaxisymmetric
wave construction. The result remains an experimental candidate.

At box center `(0.01025,0.003625,0.0084)`, the refreshed source has
cone ratio `0.3779` and a nonnegative two-ray covariance fit. The
existing simple exact-curl wave still has no lenient envelope/damping
frequency overlap: required carrier multipliers are `411` and `165`,
while viscous damping allows only `3` and `1` over the tested time
halfwidth. The wave construction remains a separate obstacle.

Reproduce the core screens:

```powershell
python experiments/root_st073/check_local_poloidal_divergence.py
python experiments/root_st073/local_poloidal_basis_screen.py
python experiments/root_st073/local_poloidal_edge_screen.py
python experiments/root_st073/radial_pressure_time_validate.py --local-poloidal-10pct
python experiments/root_st073/local_poloidal_10pct_volume.py
python experiments/root_st073/radial_peak_cone.py --field local-poloidal-10pct --radius 0.01025 --z 0.003625 --tau 0.0084 --output-name local_poloidal_10pct_source.json
python experiments/root_st073/curl_wave_scale_audit.py --source-name local_poloidal_10pct_source.json --output-name local_poloidal_10pct_wave_audit.json --radial-halfwidth 0.00275 --axial-halfwidth 0.000075 --time-step 0.00015
```

Machine-readable reports are under
`experiments/root_st073/compact_potential/`. The 10% candidate is
loadable with `local_poloidal_basis_screen.load_robust_candidate(.1)`.
