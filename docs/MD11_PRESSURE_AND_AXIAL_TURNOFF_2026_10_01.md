# C120-F17: Md>1 pressure datum, fresh core, and axial turnoff

This companion branch resolves the explicit Md=0.5 mismatch by constructing
a new Md=1.1 outer pressure source. It does not relabel old core/matching
receipts or assert that all paper smallness hypotheses are known.

## New fourteen-stage pressure source

`experiments/root_st073/lei_ren_part1_paper_Md11_pressure_datum.py` recomputes
the nominal waiting root and generates a new schedule with SHA
`15cf393dd5e73823acba88721764491a64eb5df2021e92d6e0a1b9d0d8252ebf`.
Its `.json` contains all fourteen true preheat pressure stages, physical
`Pstar_squared=exp(28)`, complete fixed-beta masses and a positive flatten
mass enclosure. The API `pressure_jets(ctx,z,order,datum)` preserves the
analytic q-power representation and a uniform Cauchy envelope for flatten.

Md=1.1 and logPstar=14 satisfy the explicit Md>1 and logPstar>Td inequalities.
The unchanged initial swirl transition mass is transferred from its existing
directed high-order receipt because Eq. (4.6) on y in [0,1] does not involve
Md, mu, delta, or waiting length. Other masses are regenerated from the new
schedule with analytic negative-slope integrals, directed Taylor-remainder
integrals, or complete positive mass bounds. No finite quadrature fit is
used as the analytic pressure datum.

The companion check verifies fourteen nonnegative stages, exact stage sums,
the independent exact axial-turnoff pressure primitive, thirteen exact axis
Taylor coefficients, new schedule fingerprint and dependency hashes.
The exact primitive comparison uses the stored Decimal endpoint: its
exp(Md) rounding remains outside the parameter-error scope. The difference
from the unrounded endpoint is recorded explicitly. The waiting root is
nominal; unknown sufficiently-small/sufficiently-large paper constants and
full terminal moment closure remain uncertified.

## Analytic core admission and fresh coefficient generation

`lei_ren_part1_paper_Md11_core_tail_admission.py` computes the new complex
pressure, derivative, g and Psi bounds from its new fourteen-stage measure.
Each is strictly below the existing universal majorant. The normalized
pressure slack is approximately 2.21849e-13; the Psi norm slack is about
128.34. Axis data, tube, linear and commuting resolvents, Lambda, j, delta
and logC are unchanged.

In the full Eq. (8.50) size and Lipschitz majorants, external pressure enters
only through the Psi center bound. All powers and coefficients in those
majorants are nonnegative. Thus dominance retains the contraction/self-map
gates for the new fixed datum. This transfers upper estimates, not an old
fixed-point solution or its finite coefficients. The new degree124 mixedC3
tail is 7.467286253e-13 < 1e-12 on scaled R in [0,4.1]. These remain
conditional bounds for stored construction data, not a full NS residual.

`lei_ren_part1_paper_Md11_interval_core.py` generates a separate resumable
core state using all128 new pressure jets and the unchanged exact amplitude
and axis machinery. State stem:
`lei_ren_part1_paper_Md11_interval_core_Z049_Z051`.
It has completed all124 radial orders, with 1174.313 seconds of coefficient
compute time. It does **not** implement temporal n-dependent coefficient
recursion. A partial state is never promoted.

`lei_ren_part1_paper_Md11_comparison_jets.py` loads only a completed degree124
new-source state and the new analytic admission. It inherits comparison
equations, while replacing all source loading. The completed controlled inlet runner
`lei_ren_part1_paper_Md11_comparison_inlet.py` checks that all128 pressure jets
in the state match the new datum (all128 differ from the old saved pressure
coefficient intervals), then encloses the core and first switched comparison
using16 directed cells.

`Md11_transition_pipeline.py` has regenerated comparison cells, actual exit,
R100 continuation, R110 switch, functional five-moment defects, and the
uniform implicit C1 five-bump inverse. Its contraction upper bound is
0.002712694947. `Md11_leading_field.py` connects the new source through the
long reshape, restoration, repaired reference annulus, outer slope and axial
cutoff/buffer to Rd. All seven diagnostic packets preserve the same P0;
the terminal Uz/Uz_y are exactly zero and Mz is retained through the buffer.

## Axial turnoff equations

`lei_ren_part1_paper_interval_outer_axial_turnoff.py` takes an explicit
angular-slope endpoint packet and pressure-schedule Md. On
`phase=log(y)/Md`, `y=log(R/Rref)`, phase in [0,1], it implements

\[
U^z=4Z[1-\sigma(\mathrm{phase})],\qquad
U^\theta=U^\theta(eR_{\rm ref})e^{-(y-1)/2}.
\]

Positive directed rectangles enclose the B and B-squared moment integrals
after the substitution y=exp(Md*phase). Angular, pressure and swirl-energy
increments are exact exponential or linear primitives. All five physical
C1 moments and P0 are inherited. At phase1, Uz and Uz_y are exactly zero;
Mz and the mixed moment remain accumulated through the eleven-unit buffer
to Rd. Ur and Ur_R are recovered using Mz and the nonzero cutoff derivative.

The moderate synthetic fixture checks84 quantities using an independent ODE
integration and independent scalar stress formulas. Exact inlet C1 moment
preservation, flat terminal fields, retained axial mass, nested rectangle
refinement and pressure-schedule mismatch rejection also pass. The numerical
fixture is not a production cone certificate. The module is now connected
to the regenerated new-source repaired annulus. Merely passing a new Md
string alongside the old pressure is rejected.

## Important obstruction: Md1.1 is not an admissible candidate

The production axial midpoint at phase0.5 is **ruled out over the enclosed
local family**, not merely an interval non-certificate. The necessary
negative stress/shear direction fails. `Md11_axial_cone_obstruction.py`
normalizes away the enormous radial placement and proves the inertial
direction per R lies in a positive interval approximately [10.3722,64.3827]
on Z in [0.49,0.51]. The additional shear-square term per R is positive.
The obstruction therefore cannot be removed by increasing Rref, adding
radial coefficients, or more endpoint samples. This does not contradict the
paper, which requires Md sufficiently large, not only greater than one.

`outer_parameter_probe.py` rebuilds a fourteen-stage trial pressure and
nominal waiting root for each Md in 1.1,1.3,2,3. It chooses
logPstar=max(14,exp(Md)+11), respecting the explicit Td inequality with
margin. Md1.1 and1.3 fail the sampled midpoint direction. Md2 and3 have
unresolved midpoint cone bounds on the full local family; their other two
sampled phases pass. These are asymptotic inertial-cone screens with exact
radial homogeneity, not full finite-placement or whole-stage certificates.
The trial sources are saved and reproducible; no additional expensive core
was generated for them.

## Next actions

- [x] Complete the new degree124 state and controlled inlet receipt.
- [x] Regenerate comparison, exit, R110, long reshape and restoration.
- [x] Regenerate functional defects and the uniform implicit C1 inverse.
- [x] Attach the new repaired reference/slope endpoint to axial turnoff.
- [x] Identify a radius-independent admissibility obstruction for Md1.1.
- First resolve the Md2/Md3 midpoint using common axial factors and pressure
  represented by the remaining masses from the **same** preheat datum. Use
  the existing high-order first-transition pressure mass at the inlet rather
  than a coarse rectangle when appropriate. Never add a fitted pressure tail.
- Certify full phase and axial cells before selecting a new parameter source
  and rebuilding its analytic majorants/required degree/fresh finite core.
- Construct the -1/2-mu turn and reserved strong-cone power interval.
- Continue pulse, flatten, angular repair, waiting, collar and exact heat field.
- Establish higher axial smoothness, original parameter errors, full energy,
  flat remainder, true temporal recursion, oscillatory correction and Cartesian
  residual as distinct later requirements.

Reproduce by running the pressure producer/check, tail admission, resumable
core generator, and then the controlled inlet runner. Run the turnoff check
separately for its moderate synthetic fixture. Existing receipts stay frozen.
