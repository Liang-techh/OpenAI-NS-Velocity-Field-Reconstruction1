# Second axial derivatives through the common source exit

The reconstruction now has a separate second-order axial ring and comparison,
bridge and analytic continuation adapters. They propagate actual derivatives
f, f_Z and f_ZZ over the retained P9/W2 pressure/width ring. First-order data
are not promoted by assuming a missing second derivative is zero.

## Source and derivative contracts

- `axial_second_jet.py` carries explicit actual derivatives, not factorial-normalized
  Taylor coefficients, through arithmetic and exp/log/sqrt chain rules.
- `pressure_width_second_axial_comparison.py` converts the common source Taylor
  rows and inherits the Section 9.23 comparison equations. Four Taylor
  coefficients per radial row are required: a first-Z driver carrying its own
  second derivative uses the third source derivative.
- `pressure_width_second_axial_bridge.py` propagates the same finite RK equations,
  exposes moments_ZZ and recovers Ur_Z from the radial moment identity. Ur_ZZ
  stays explicitly unavailable because it would require third moment derivatives.
- `pressure_width_second_axial_continuation.py` propagates those slots through
  the same exact exponential-polynomial primitives used by the first-Z branch.
  It also exposes second derivatives of both raw quadratic energies and P0.
- `pressure_width_second_axial_switches.py` propagates complete three-slot
  source states through the R100--R110 switches. It rejects first-only source
  states and retains fields, moments, raw energies and pressure datum in the
  same second-order ring.

The original first-Z adapters and original caches remain usable. A separate
cache is generated from an actual degree-18, Z-depth-3 core; a depth-2 R100
snapshot cannot supply the required higher derivative row.

## Completed checks

The resolved integration fixture replays shifted scalar finite-core data and
the finite scalar/first-Z RK paths. It checks second derivatives of comparison
F, Uz, P, I_theta, I_z, D and E, exit fields and all five moments, continued
fields, all five continued moments, both raw quadratic integrals, and R110
switch fields/moments/raw energies. It also checks Ur_Z against axial stencils.
Maximum scaled discrepancy is approximately `1.513e-32` using an outer stencil
step `1e-8`. This verifies differentiation of
the declared finite construction, not the error of the exact source ODE.

A fresh actual same-preheat source with third-depth Taylor data completed at
Z=.3. The actual inner comparison endpoint receipt stores value, first-Z and
second-Z atoms of all five moments, fields and stresses, with an explicit
pressure datum from that source. It does not infer an axial interval bound.

Actual second-Z collar and analytic continuation to R100 also completed.
The post_collar_R100 receipt preserves value/first/second derivatives of F,
Uz, P, all five moments and both raw quadratic integrals, plus P0_ZZ and Ur_Z.
Its independently generated second-Z cache is available for the actual switches.

Actual second-Z switches to R110 completed as well. The R100 inlet value,
first-Z and second-Z fields and all five moments agree with the continued
source to scaled error below `1e-200`. The R110 source is saved in a separate
local cache and in the R110 JSON section, with all three axial slots for the
fields, five moments and raw energies. P0_ZZ and Ur_Z are explicit.

## Commands and evidence

```powershell
python experiments/root_st073/lei_ren_part1_paper_axial_second_jet_fixture.py
python experiments/root_st073/lei_ren_part1_paper_pressure_width_second_axial_fixture.py
python experiments/root_st073/lei_ren_part1_paper_pressure_width_second_axial_check.py
python experiments/root_st073/lei_ren_part1_paper_pressure_width_second_axial_check.py --resume --post-collar
python experiments/root_st073/lei_ren_part1_paper_pressure_width_second_axial_check.py --resume --switches
```

These continuation commands reuse only the new internally generated Z=.3 caches
and extends it to the collar and R100. It never extrapolates to another Z.
Its completion must be read from the process result and post_collar_R100 JSON
section; existence of the earlier inner-endpoint receipt is not completion.

## Remaining work

The actual source now reaches R110 with second-Z data, but the centered defect
source labels still need the second-Z flat-kernel composition. Only then can
the five-bump inverse and physical recovery
consume the complete second-Z source. Uniform interval/Taylor bounds,
finite core/jet/RK/quadrature remainders, and endpoint Z extension remain open.
The new derivatives alone do not certify cone margins, five-moment functional
closure, finite energy, exact heat matching, or temporal n-recursion.
