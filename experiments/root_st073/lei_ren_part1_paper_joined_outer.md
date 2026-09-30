# Corrected inner to shared outer heat join

Run from the repository root with:

```powershell
python experiments/root_st073/lei_ren_part1_paper_joined_outer.py
```

The prototype constructs `CorrectedInnerField(build_candidate())`, then
dispatches at its terminal `Rh` to the exact shared profile stored at
`inner.provider.reshape.switches.comparison.bundle['profile']`. The public
coordinate is `y = log(R / inner.r)`. The physical wrapper is a
`JoinedLocalCoreExitField` subclass of `LocalCoreExitField`, so its scope
label identifies the outer/heat continuation.

The terminal values use `inner.evaluate_x(exp(1), Z)`, which retains the
five-bump correction. The reference mass is `4 Z Rh` with Z derivative
`4 Rh`. At `Z = .3`, the receipt preserves the actual terminal differences:

- `mass_offset / Rh = +3.0296144408125e-176`;
- `mass_offset_Z / Rh = -1.8764211050504e-186`;
- `pressure_offset = -4.0566932235042e-98`;
- `pressure_offset_Z = -1.4523644574512e-4`.

Those values are carried forward into the outer axial average, radial
transport, pressure, and pressure Z derivative. The pre-offset join relative
errors are about `6.96e-111` for `F`, `1.37e-16` for `FZ`, `3.39e-262` for
`Uz`, `2.77e-177` for `Ur`, and `1.39e-110` for `P`. After carrying the
mass offsets, the mass and mass-Z join errors are about `3.39e-262` and
`1.56e-292`; pressure offsets are reproduced by construction.

The five receipt rows are `Rh_join`, `Rref`, `heat_tail_start`,
`heat_connection`, and `exact_heat`. Each heat row uses the shared source
profile's actual axial average. The full five cumulative moments are marked
unavailable beyond the early exterior; no synthetic moment values are
returned. Outer Z averages use a fourth-order `h = 1e-6` stencil through the
profile API, whose Z argument is float-backed, so the tail coefficient is a
sampled diagnostic rather than a continuous derivative certificate. The
outer angular correction Z jet is also unavailable in that API; schedule
log-Z jets are used at the sampled outer points.

At the exact heat point, the carried tail transport coefficient is nonzero.
The receipt stores it and its ratio to `Rh` in full signed-log and arbitrary
MP form. Its interpretation is `V/R ~ C/R`, giving a radial transport energy
term proportional to `C^2/R` and hence logarithmic divergence when `C != 0`.
This is recorded as an obstruction; the prototype makes no finite-energy,
Navier--Stokes, stress-cone, or scale-recursion claim.

The `LocalCoreExitField` physical roundtrip at the exact heat point has zero
reported relative velocity and radius error at the working precision and a
Z relative error below `8e-446`. This checks the callable chart only.
