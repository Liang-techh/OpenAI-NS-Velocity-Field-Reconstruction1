# Coherent waiting and coupled angular targets â€” 2026-09-30

The actual continuous angular schedule now determines its own zero-Z H=1 waiting length before pressure/core reconstruction. This is an opt-in construction path, not a late moment overwrite or pressure gauge shift.

## Implemented

- `build_joined_field(continuous_pressure=True, coherent_waiting=True)` rebuilds the schedule, continuous pressure anchor, core and dependent inner connection in order.
- `ContinuousWaitingMatch` retains epsilon and the positive collar integral separately, avoiding subtraction of an almost-unit K0.
- `angular_difference_from_axis` propagates the Z-dependent difference directly. Independent direct subtraction at y=.5 and 2 agrees within 7.26e-122 at 120 digits; Z=0 returns exactly zero. These checks cover those preheat points, not the full extreme-scale domain.
- `CoupledAngularTargets` retains actual inner offsets, unbumped preheat differences and heat corrections separately. It does not query corrected cumulative moments recursively.
- `solve_coupled_angular` handles signed, zero and tiny targets with analytic tangents, small-branch and positivity diagnostics. Four independent coefficient fixtures pass.

## Actual source result at Z=.3

The initial target extractor reported r=-1077.580145807... . A direct replay exposed a missing tail-amplitude conversion: the preheat baseline is normalized by Cinf*Rtail^-a*(1-epsilon), while heat atoms use Cinf*Rtail^-a. This factor is now retained. That initial oversized value must not be interpreted as a certified physical field defect. With the coherent waiting identity, the updated path gives retained r=2.31629269314e-837 and rZ=1.47949417366e-836. Its |r|/mu^29 is 4.56070569695e-45. The waiting length is approximately 471.1046619852684520963565252967655; changes at digits beyond ordinary floating precision materially affect the amplified terminal target.

Nominal terminal pressure is 1.18757636672e-98, versus 2.15130102859e-69 for the preceding continuous-pressure core. This is another factor about 1.81e29 improvement, not exact pressure closure.

Reproduce:

```powershell
python experiments/root_st073/lei_ren_part1_paper_coupled_angular_rebuild.py --coherent-waiting
```

Receipt: `experiments/root_st073/lei_ren_part1_paper_coupled_angular_coherent_waiting.json`. The non-coherent replay remains separately available without the flag.

## Limits and next actions

1. The retained angular target uses the solved zero-axis H=1 identity. After the amplitude conversion, direct forward subtraction gives -5.67419145461e-108, far above the retained target. It is an unresolved finite-precision/quadrature floor and is also reported; cancellation at finite precision is not an enclosure for the tiny retained target. Enclose waiting/collar inputs before claiming exact terminal moments.
2. The actual coupled solve gives d1=-8.66472474236e-838 and d2=6.40241372031e-837, within the selected smallness bound and with positive bump multiplier. Its branch `accepted` flag is not row closure: the pressure target-relative residual is enormous. Read the coupled coefficient receipt with its pressure-target scale warning. The canonical heat pressure target is exponentially below the angular target; scalar floating arithmetic cannot certify that cancellation merely from a small branch or a nominal zero residual.
3. Reconstruct the core using complete analytic preheat pressure *components*, including the tiny post-Rv atoms. The present core still uses the prefix datum, so its existing-core pressure target disagrees with the canonical restoration target. Do not install the new bumps as a globally closed field before resolving that compatibility.
4. Regenerate dependent incoming/axial corrections and all five terminal functional moments after installing a compatible coefficient representation.
5. Resolve the actual nonzero 1/r radial tail for finite energy, then establish the admissible cone, lower-order recursion, two-family oscillatory cancellation and flat remainder.

No global finite-energy, stress-cone, scale-recursion or full residual certificate is claimed.
