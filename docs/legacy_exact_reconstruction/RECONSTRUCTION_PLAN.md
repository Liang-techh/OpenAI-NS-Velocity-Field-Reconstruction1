# Reconstruction plan and truth-status ledger

Updated 2026-09-13. The target remains an executable counterpart of the **entire published
velocity construction**, not merely the blow-up exponent. Every construction layer needs
source-mapped choices and independent verification. The goal has not been downgraded.

## Status convention

- **paper-exact:** direct published formula/algorithm; floating evaluation is still numerical.
- **formal-structure:** the structure is encoded but some profiles, parameters or recursive choices are missing.
- **diagnostic-only:** a numerical checker or experiment, not a proof construction.
- **toy:** illustrative data that must never be called the OpenAI counterexample.
- **pending:** no runtime implementation for the claimed layer.

`references/provenance_manifest.json` preserves the stable source/layer IDs introduced by
the parallel provenance work. The installed CLI additionally exposes stages 0–8 separately.
Tests check source-pin consistency. Neither metadata nor a caller's `paper_exact` flag proves
correctness. A green test run never automatically upgrades a stage.

## Stage 0 — similarity geometry: implemented formula

`coordinates.py` implements Eq. (4.1):

```text
A=1/2+h, D=1/2-h, tau=1-t
z=q^D eta, tau=q(1-eta^2), X=r^2/(2q)
q-z^2 q^(2h)=tau
```

The root solve is dimensionless and relative-scaled. `solve_q_from_tau` and
`similarity_coordinates_from_tau` avoid `1-tau` rounding; `d=tau/q` avoids subtracting
nearly equal numbers. Tests cover manufactured roots down to q=1e-200.
This is not a guarantee for every representable input or arbitrary precision.

## Stage 1 — leading velocity/profile: partial

`profiles.py` and `velocity.py` implement Eqs. (4.3)–(4.7):

```text
u_theta=q^(-A) E, u_z=q^(-A) U, r u_r=V0
AU=(1/X) integral_0^X U(x,eta) dx
V0=X/L [2 eta U-2D eta AU-(1-eta^2) d_eta AU]
L=1-2h eta^2, Pi_X=E^2/(2X)
```

The radial averages now use cached Gauss–Legendre nodes or user-supplied exact averages.
Regular-axis profile data can provide `F` with `E=sqrt(2X)F`; invalid nonzero axis swirl
is rejected instead of silently erased.

The actual schedule now has a staged coefficientwise leading-profile chain in
`axis_coefficient_wide_natural_remainder.py`, `axis_coefficient_wide_first_picard.py`, and
`axis_coefficient_wide_first_picard_remainder.py`. The wide mixed-scale adapters represent the
complete `naturalRemainder(x0)`, the first Picard state `x1`, the full coefficientwise
`naturalRemainder(x1)`, and the second Picard state `x2`; ordinary,
pressure-linear, and pressure-square channels remain explicit rather than being collapsed into
binary64 values. The generic sparse mixed-scale algebra and formal triangular solver are now
implemented. The solver computes each requested formal coefficient directly from lower radial rows
and uses no fixed iteration cutoff. These are local formal coefficient recursions and do not provide
weighted-space membership, a converged fixed point, global `AxisSpace` certification, or a
paper-exact `phi/u` profile.

The one-graph solver `jet_prefix` now feeds `axis_coefficient_profile_prefix.py` through
`formal_axis_profile_prefix`: finite natural-radial `Y` polynomial prefixes evaluate `phi`, `u`,
radial derivatives, eta-derivative jets, radial averages, and the scaled formal pressure primitive
`P = primitive((a phi)^2)` while retaining every split signed-log amplitude/Lambda channel. The
pressure primitive uses the full Bell eta family for `a^2` and remains separate from the axial
remainder pressure forcing. `axis_coefficient_radial_tail.py` supplies a conditional omitted-row
majorant using the pinned `AxisSpace` weight; it does not establish global membership or absorb
coefficient roundoff. This prefix does not provide a converged or paper-exact profile; the
reviewed wide prefix now carries the corresponding sparse physical `Pi` map.
`NaturalProfileAssembly` now exposes the smooth ratio `F = a phi` separately
from the physical swirl `E = sqrt(2 X) F`; its focused normalization target
passed `4` tests in `0.15 s`.

The reviewed `wide_natural_profile_prefix` now retains the actual Decimal
`Lambda` and symbolic amplitude while exposing the maps `F`, `E`, `U`,
`d_eta U`, `Ubar`, and `Pi`. The reviewed conditional `F/P` norm-budget
extension reuses the same actual reference/remainder chain, and its optional
six-field `conditional_truncation_bounds` interface applies exact Fraction
scaling and upward square-root rounding. The focused command
`python -m pytest -q tests\test_natural_axis_wide.py tests\test_axis_coefficient_profile_prefix.py -W error`
returned `5 passed in 0.66 s`; the final wide-only check returned `4 passed in
0.64 s`. These remain conditional formal prefixes: no global weighted
membership, coefficient/phase-quadrature roundoff budget, or physical
certification is claimed.

The formal finite-Picard bridge is now available through
`formal_axis_picard_family_state(solver)`. Its
`jet_prefix(iterations, max_n, m, eta)` starts at the actual reference `x0`
when `iterations=0` and applies the same routed sparse map
`reference + R(input)/(2 Lambda)` for each requested finite step. The wide
physical prefix also retains `d_average_U_deta` and the regular radial
velocity numerator `V0`. Its conditional eight-field tail mapping uses exact
Fraction scaling, including an upward-verified `sqrt(2 X)` factor. The
source-pinned manual finite-filtration proof has been accepted: under the
theorem-side hypotheses, formal row `n` stabilizes after update `n` and
identifies the compatible fixed-point coefficient. This exact-arithmetic
conclusion concerns formal recurrence dependence and does not bound errors in
the evaluated scalar jets or coefficient roundoff; the jet values are finite,
but their numerical errors lack certified bounds. Focused runtime acceptance
passed the assigned bridge targets. Even with the conditional exact-identification
proof, evaluated scalar jet enclosures and
coefficient-roundoff control remain the next Stage 1 blocker.

The exact rational input seam is now available through
`RationalAxisCoefficientData(actual_data)`. It treats selected `h`, `j`, and `sigma`
as exact binary rationals, recomputes `A` and `D` from `h`, and evaluates the fixed
scalar fields, `chi`, and angular reference coefficients with exact Fraction Taylor
algebra. Nearest Decimal96 values and directed scalar/angular enclosures are exposed;
the exact domain is `[-11/10, 11/10]`, so binary64 `1.1` is outside this interface
while the physical `|eta| <= 1` chart is unchanged. `zStar`, pressure/Bell-phase
inputs, and accumulated coefficient errors remain outside this seam. The measured
baseline command `python -m pytest -q tests\test_axis_coefficient_rational_data.py
tests\test_axis_coefficient_formal_solver.py tests\test_axis_coefficient_picard_family.py
tests\test_axis_coefficient_profile_prefix.py tests\test_natural_axis_wide.py -W error`
returned `16 passed in 25.22 s`; the remaining cleanup only removes a discarded
reference-half evaluation and preserves the value formulas.

The rational input provider also supplies exact normalized amplitude-power Bell jets through
`amplitude_power_bell_fraction(power, order, eta, Lambda)`. It uses the selected rational
`normalizedGradient` full derivatives and Bell recurrence without forming amplitude, phase, or
`C`, with nearest Decimal96 and directed enclosure wrappers. The affected-source command
`python -m pytest -q tests\test_axis_coefficient_wide_natural_source.py tests\test_axis_coefficient_rational_bell.py -W error` returned `12 passed in 0.30 s`, using independent Fraction expectations for the first two Bell values. The broader formal/Picard/profile/wide integration command
`python -m pytest -q tests\test_axis_coefficient_formal_solver.py tests\test_axis_coefficient_picard_family.py tests\test_axis_coefficient_profile_prefix.py tests\test_natural_axis_wide.py -W error` returned `13 passed in 24.37 s`; phase and total-product construction and accumulated coefficient-error control remain outside this boundary.

The new `axis_phase_log_enclosure.py` primitive maps an externally supplied phase interval through the
exact Fraction affine relation `log a = Lambda * phase - log_C`, then presents directed Decimal bounds
with `paper_exact=False`. The focused command `python -m pytest -q tests\test_axis_phase_log_enclosure.py -W error`
returned `4 passed in 0.14 s`. This does not replace runtime phase quadrature or certify `C`, scalar
parameter selection, or global reconstruction. Fixed-order Cauchy bisection requires order growth to
meet arbitrary tolerance; the related GitHub `hub281`/`task282` exchange remains unintegrated.

`axis_phase_integral.py` now implements the linked exact-rational default phase path. It expands the
rational kernel on centered cells, requires the exact denominator gate `theta < 1/2`, and adaptively
subdivides and increases Taylor order under finite caps, failing closed when a cap is reached. The
actual `ActualScheduleAmplitudeLogState` defaults to phase tolerance `10^-12`; its
`phase_enclosure` receives exact rational `h,j,sigma,eta` inputs and uses an exact-input cache of
size `16`. `log_amplitude` returns the nearest Decimal96 midpoint of the validated log interval,
and `default_log_amplitude_enclosure` exposes the interval after the selected positive `Lambda`
amplification and symbolic `C` exponent transport. The explicit
`log_amplitude_enclosure` path requests a log tolerance and passes
`absolute_log_tolerance/Fraction(Lambda)` to the phase integrator. `phase_samples` and
`legacy_log_amplitude_diagnostic` retain the former 4001-point trapezoid only as a named diagnostic.
The historical pre-default commands returned `5 passed in 0.17 s` for
`tests\test_axis_phase_integral.py`, `3 passed in 0.17 s` for
`tests\test_axis_amplitude_phase_enclosure.py`, and `6 passed in 0.16 s` for the prior amplitude
compatibility target. The completed bounded command
`python -m pytest -q tests\test_axis_phase_integral.py -W error` returned `7 passed in 0.20 s`.
The combined default regression
`python -m pytest -q tests\test_axis_amplitude_validated_default.py tests\test_axis_coefficient_amplitude.py tests\test_axis_phase_log_enclosure.py -W error`
returned `13 passed in 43.71 s` under the hard `120 s` cap. On the actual selected fixture at
`eta=-1/50`, the prior direct phase measurement at tolerance `10^-12` returned in `10.664 s` with
`225` cells, maximum order `64`, estimate `0.05073499503328993`, and error `7.357907e-13`.
The current default validated binary64-input diagnostic hit its `45 s` cap before completing, so no
precise runtime is claimed. These are selected-fixture observations, not a general benchmark. The
chosen rational phase interval and its selected scalar `Lambda` amplification remain
conditional; they do not certify SchedulePressure quadrature, theorem parameter selection, `C`
selection, coefficient roundoff, global `AxisSpace`, or the full reconstruction. A prior combined
default attempt was interrupted after approximately `9` minutes without a pytest summary and remains
historical; the bounded command above resolves the default-regression performance blocker. Full-suite,
demo, and strict-audit acceptance remain separate. The metadata chain now reaches the assigned
remainder, second-Picard, formal-solver, and wide-profile surfaces. The next Stage 1 boundary is
validated Bell-factor/logarithm enclosures, general coefficient arithmetic, SchedulePressure and
parameter proofs, and compatible global `AxisSpace` certification.

The amplitude log-scale metadata evidence now totals `35` distinct tests from the prior
`24` foundation/integration/mixed/aggregate tests, `5` nonlinear metadata tests in `0.80 s`,
`3` remainder/second-Picard tests in `0.38 s`, and `3` formal/wide-profile tests in `0.22 s`.
The chain reaches all assigned remainder, second-Picard, formal, and wide-profile metadata
surfaces, with shared actual-source and cross-branch identity checks. A separate batch of the
existing default test file returned `4 passed in 20.19 s` and included one new nonzero-eta
(`eta=-0.02`) metadata test, bringing the new metadata total to `36` distinct tests. These
focused results do not certify exact Bell-factor/logarithm enclosures, general coefficient
arithmetic, SchedulePressure or parameter selection, global weighted-space bounds, or the
paper-exact reconstruction.

## Current local checkpoint — 2026-09-14

The full suite immediately preceding the metadata chain returned `1936 passed, 1 failed` in
`1723.45 s`; the sole fixture-argument issue was fixed and its affected test subsequently
passed in `0.20 s`. The full suite has not been rerun after that fix, so current full acceptance
remains pending. The demo exited `0` with all checks; the direct strict-audit module exited `2`
as expected and this was confirmed in its log, while an earlier console-wrapper exit `1` is
historical. All relevant source compilation and review checks passed. The local checkpoint is
being saved without an external push.

**New constructed component:** `heat_exterior.py` implements Appendix A.6, Eqs. (A.32)–(A.38):
H and its derivatives, exterior K(r,t), E_heat(X,eta), and centrifugal pressure normalized
at infinity. It includes derivative/ODE checks and independently tested NS balance for
this exterior-only flow. This component is not smooth at r=0 and cannot stand in for the core.

**Still required:** identify compatible global weighted `AxisSpace` majorants and control
coefficient roundoff, preserve adjacent-eta compatibility and approximation-error budgets, and
complete the compatible weighted profile around the reviewed
`wide_natural_profile_prefix` maps. The scaled formal pressure primitive and sparse physical `Pi`
map are present, remain separate from remainder forcing, and are not a paper-exact converged
profile. The reviewed conditional `F/P` norm-budget and six-field truncation interfaces remain
conditional on the actual chain and compatible `AxisSpace` identification; then translate Theorem 4.6
and Appendices A/B/C into a regular inner profile, matching moments, compact moment corrections,
shear modification and admissible cone conditions. Emit all constants/choices and generated
coefficient hashes. Do not splice a Gaussian into the missing core and mark the result complete.

## Stage 2 — all-order background: assembly only

`background.py` implements lambda_n=2nh, streamfunction/potential assembly and cutoff sums.
The generic cutoff in `cutoffs.py` is C-infinity mathematically, replacing the old C1
smoothstep, but **is not the paper-selected recursive cutoff schedule**.

Still required: solve Eqs. (5.2)–(5.6) for each coefficient, preserve compact moments,
and select cutoffs from the paper's estimates.
Acceptance: reproducible coefficients plus increasing derivative-order evidence for the
predicted super-algebraic residual decay. Finite-depth plots alone are insufficient.

## Stage 3 — dyadic charts and transported waves: geometry only

`charts.py` implements Eqs. (6.1)–(6.6): Q=2^-ell, epsilon=Q^h, S*=ell^2,
chart/physical conversion, the fixed integer covering matrix and its powers, and normalized
velocity/pressure/residual scaling. Q is fixed per chart; do not differentiate it as q(z,t).
Integer matrices use Python integers rather than overflowing int64. Floating fast-phase
accuracy is not certified by exact integer matrix arithmetic.

Still required: slow labels, separated auxiliary-torus supports, transported
wavevectors/polarizations, phase dynamics and rounded frequencies. Arbitrary frequency
lists are not acceptable for `paper-exact` status.

## Stage 4 — oscillatory stress realization: pending

Implement Section 7: admissible cone, positive covariance decomposition, primary waves,
amplitude solve, exact divergence-free realization by vector potentials, curl remainder,
physical evaluation and support separation.
Acceptance: independently computed averaged quadratic flux realizes the requested stress
to the stated order, with support and divergence identities checked.

## Stage 5 — compact mean corrections: pending

Implement Section 8 radial inverses, moment-preserving compact corrections and the
five-equation defect solve. Preserve all moments and supports needed by the iteration.
Acceptance: the zero-auxiliary-average defect is cancelled without breaking those constraints.

## Stage 6 — residual-improvement iteration: pending

Implement Section 9, including the recursive correction sequence and summation (9.21):

```text
Delta u_j=curl(A_j)+B_j
A=A0+sum_j chi(a_j q) A_j
B e_theta=B0 e_theta+sum_j chi(a_j q) B_j
u_loc=curl(A)+B e_theta
```

Choose a_j from estimates, not plot fitting. Acceptance: residuals and every tested derivative
vanish to increasing order with stable truncation studies. Full convergence still needs
analytic or formal bounds, not only a numerical experiment.

## Stage 7 — final compact field and force: composition only

`local_field.py` implements u=curl(cA)+cB e_theta, Eq. (10.4). Apply the cutoff before curl.
The direct swirl and its localized product must be axisymmetric to preserve divergence.
Factories encode that restriction; tests include a non-axisymmetric counterexample.

Still required: the completed local field, paper-specific support/cutoffs, force
f=u_t+(u.grad)u-Delta u+grad p at viscosity one, and the Section 10 smooth extension through t=1.
A force defined numerically as the residual has **not** thereby been shown smooth.

## Stage 8 — independent verification: partial diagnostics

Implemented: finite-input and stencil checks, bounded-time differentiation, manufactured
solutions with independently specified force, space/time refinement, exterior heat balance,
leading singular-path slope checks, deterministic CSV/JSON output, hashes and an audit CLI.

Still required: interval or symbolic/AD backends, kinetic energy/compact support checks,
full corrected-field convergence, trajectories/visualizations of the actual field, and
source-mapped Lean theorem cross-checks. Upstream is pinned, but no Lean build is claimed.

## Next independently reviewable work packages

1. Complete the regular inner/outer leading-profile constructor and matching/cone tests.
2. Implement one real Section 5 recursive coefficient and its residual/moment tests before generalizing.
3. Extend fixed charts to transported phases and verified support separation.
4. Implement the Section 7/8 stress and mean correction operators on independently generated defects.
5. Add analytic/interval error certificates, then integrate the Section 9/10 correction sequence.

These are uncompleted work packages, not claims that background jobs have been started.
Avoid parallel edits to the same path; re-read main and preserve concurrent work before committing.
