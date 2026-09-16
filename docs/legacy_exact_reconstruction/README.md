# OpenAI NS Velocity Field Reconstruction

An independent, **partial** executable reconstruction of the velocity construction in
OpenAI's September 2026 *Finite Time Blowup for Navier–Stokes*. This is not an OpenAI repository.

**Agent handoff:** read the [current checkpoint](docs/CURRENT_CHECKPOINT.md) and
[detailed task board](docs/AGENT_TASKS.md) before claiming work. Task completion,
coordinator acceptance, and PR merge status are recorded separately.

**Current result:** tested similarity kinematics; natural-axis/range structure, an admissible
ideal-prefix pressure witness, the constructive outgoing scalar schedule, executable
`finalAngular/clockWeight`, the actual outgoing `SchedulePressure.axisPressure`, an analytic
no-sampling low-|Z| `H^2` margin instantiating the theorem-side `sigma`, the theorem-faithful
`Lambda/C` selection algebra, a fail-closed AxisContraction-style norm-ledger adapter, the
canonical analytic-input propagation that turns a certified common complex radius `rho` and common
value bound `B` into `epsilon=rho/2`, `radiusLoss(1/2)=12`, coefficient norm `<=12B`, normalized
amplitude `M=12`, and a conservative complete naturalResolvent factorial-series bound, plus an
actual-schedule analytic-neighborhood certificate that constructs a positive common `rho`, the
11-field complex bound `B`, and a conservative `realPartSup(axisPhase)` upper bound without sampled
complex maxima. The actual schedule is refined componentwise: the same one-field Cauchy bound is
applied separately as `||element_k|| <= 12 B_k`, the maximum is retained only for the common
`CoefficientFamily.bound`, and the resolvent uses the chi-specific `||chi|| <= 12 B_chi`. On the
current theorem-admissible regression instance this removes the earlier cross-field resolvent
overflow; the complete factorial resolvent majorant is representable in binary64. The original
binary64 remainder path still fails closed during conservative `remainderBound/remainderLip`
propagation, but `axis_remainder_wide_bounds.py` mirrors the same positive theorem-shaped
Controlled/remainder algebra in 96-digit `Decimal`, rounded toward `+infinity`, and keeps the actual
regression ledger finite even when a conservative final majorant exceeds `sys.float_info.max`.
`natural_scale_selection_wide.py` and `stage1_scale_chain_wide.py` carry those same wide bounds
through the pinned `Lambda=max(1+B+L,1+(14000/9)B)` choice without binary64 down-conversion and
represent `C=exp(Lambda*phaseSup)` by its certified symbolic exponent rather than mistaking a huge
finite threshold for floating-point infinity. `axis_fixed_point_picard.py` takes that actual-
schedule wide chain one theorem step farther: it verifies the pinned scalar contraction gates
`s*remainderBound <= 1` and `s*remainderLip <= 1/2`, `s=1/(2 Lambda)`, with upward-rounded Decimal
arithmetic and exposes the exact Picard update shape plus a geometric tail budget. The actual
schedule now also feeds `axis_reference_pair.py`, which materializes the pinned `referencePair`
zeroth parameter-jet radial coefficient functions on the full coefficient window: the angular
coefficients follow the exact `AxisReference` recurrence/factorial formula, and the axial reference
uses `u0[1]=-(1/2) inverseL*zStar` with `zStar` rebuilt from the actual schedule pressure. The
coefficientwise Stage 1 chain now includes the complete mixed-scale `naturalRemainder(x0)`, its
first Picard state `x1`, the full coefficientwise `naturalRemainder(x1)`, and the coefficientwise
second Picard state `x2`. These adapters preserve ordinary, pressure-linear, and pressure-square
scale channels without collapsing them into binary64 totals. A generic sparse mixed-scale algebra
and formal triangular coefficient solver are now implemented. The solver computes each requested
formal coefficient directly from lower radial rows and uses no fixed iteration cutoff. This generic
formal prefix is an x1-anchored triangular coefficient representation; it is separate from the
coefficientwise first-Picard `x1`/`x2` adapters and is not itself a finite Picard iterate. Its one-graph
`jet_prefix` path now feeds `axis_coefficient_profile_prefix.py`, which evaluates finite natural-radial
polynomial prefixes for `phi`, `u`, radial derivatives, eta-derivative jets, and the radial-average
map while preserving split signed-log amplitude/Lambda channels. The same graph now exposes the
scaled formal pressure primitive `P=primitive((a phi)^2)` as a sparse pressure map and signed-log
view; this is separate from the axial remainder pressure forcing and does not apply an outer
`1/Lambda` factor. The conditional `axis_coefficient_radial_tail.py` majorant is tied to the pinned
`AxisSpace` weight and does not infer a global norm or absorb coefficient roundoff.
`natural_axis_wide.py` now exposes the reviewed `wide_natural_profile_prefix`, retaining the actual
Decimal Lambda and symbolic amplitude while exposing sparse `F`, `E`, `U`, `d_eta U`, `Ubar`, and
`Pi` maps. The reviewed conditional `F/P` norm-budget extension derives its inputs from the same
actual reference/remainder chain. `WideNaturalProfilePrefix.conditional_truncation_bounds` also
provides the six-field conditional mapping with exact Fraction scaling and upward square-root
rounding; its latest focused acceptance passed 4 tests in 0.64 s. These interfaces remain conditional
on compatible global `AxisSpace` identification and exclude coefficient and axis-pressure
quadrature roundoff. This remains a formal execution boundary: no weighted-space fixed point,
converged `phi/u` profile, or paper-exact leading field is claimed. `NaturalProfileAssembly` continues
to expose the smooth ratio `F=a phi` separately from the physical swirl `E=sqrt(2X) F`.

The formal Stage 1 bridge now also exposes
`formal_axis_picard_family_state(solver).jet_prefix(iterations, max_n, m, eta)`:
`iterations=0` is the actual reference `x0`, and each further finite step applies
the same sparse map `reference + R(input)/(2 Lambda)` through externally routed
full derivative families. The wide physical prefix additionally retains
`d_average_U_deta` and the regular radial-velocity numerator `V0`; its
`conditional_truncation_bounds` values use exact Fraction scaling and an
upward-verified square-root factor. A source-pinned manual proof of
coefficientwise finite-filtration stabilization has been accepted: under the
theorem-side hypotheses, row `n` stabilizes after update `n` and identifies the
compatible fixed-point coefficient. This exact-arithmetic argument concerns
formal recurrence dependence only and does not bound errors in the evaluated
scalar jets or coefficient roundoff; the jet values are finite, but their
numerical errors lack certified bounds. These additions remain conditional on compatible global
`AxisSpace` identification and do not claim a converged fixed point.

The exact rational input seam is now available through
`RationalAxisCoefficientData(actual_data)`: selected `h`, `j`, and `sigma` inputs are
treated as exact binary rationals, `A` and `D` are recomputed from `h`, and the nine
fixed scalar fields plus `chi` and the angular reference coefficients use exact
normalized Taylor/Fraction algebra. `jet_decimal` and
`angular_reference_jet_decimal` use nearest Decimal96 conversion, while the scalar
and angular enclosure APIs use directed Decimal96 bounds. The exact domain is
`[-11/10, 11/10]`, so binary64 literal `1.1` is outside this rational interface;
the physical `|eta| <= 1` chart is unaffected. `zStar`, pressure/Bell-phase inputs,
and accumulated coefficient-error control remain outside this seam. The measured
baseline command `python -m pytest -q tests\test_axis_coefficient_rational_data.py
tests\test_axis_coefficient_formal_solver.py tests\test_axis_coefficient_picard_family.py
tests\test_axis_coefficient_profile_prefix.py tests\test_natural_axis_wide.py -W error`
returned `16 passed in 25.22 s`; the remaining cleanup only removes a discarded
reference-half evaluation and preserves the value formulas.

The rational input provider also supplies exact normalized amplitude-power Bell jets through
`amplitude_power_bell_fraction(power, order, eta, Lambda)`, using the selected rational
`normalizedGradient` jets and the full Bell recurrence without forming amplitude, phase, or `C`.
Nearest Decimal96 and directed enclosure wrappers preserve the same exact source value. Focused
acceptance for this Bell extension returned `12 passed in 0.30 s` with the affected-source command
`python -m pytest -q tests\test_axis_coefficient_wide_natural_source.py tests\test_axis_coefficient_rational_bell.py -W error`, using independent Fraction expectations for the first two Bell values. The broader formal/Picard/profile/wide integration command
`python -m pytest -q tests\test_axis_coefficient_formal_solver.py tests\test_axis_coefficient_picard_family.py tests\test_axis_coefficient_profile_prefix.py tests\test_natural_axis_wide.py -W error` returned `13 passed in 24.37 s`; phase and total-product construction and accumulated coefficient-error control remain outside the certified boundary.

The new `axis_phase_log_enclosure.py` primitive maps an externally supplied phase interval to exact
Fraction affine bounds for `log a = Lambda * phase - log_C`, with directed Decimal presentation and
`paper_exact=False`. Its focused command `python -m pytest -q tests\test_axis_phase_log_enclosure.py -W error`
returned `4 passed in 0.14 s`. It does not replace runtime phase quadrature or certify `C`, scalar
parameter selection, or global reconstruction; fixed-order Cauchy bisection requires order growth for
arbitrary tolerance. The related GitHub `hub281`/`task282` exchange remains unintegrated.

`axis_phase_integral.py` now supplies the linked exact-rational default phase path: it forms centered
`P/Q` polynomials for `g=-L*H/(H^2+sigma^2)`, accepts a cell only when the exact Cauchy denominator
gate `theta < 1/2` holds, and adaptively bisects cells and raises Taylor order until the rational
error budget is met. `ActualScheduleAmplitudeLogState` defaults to phase tolerance `10^-12`;
`phase_enclosure` uses the selected exact Fraction `h,j,sigma,eta` inputs and an exact-input LRU
cache of size `16`. `log_amplitude` returns the nearest Decimal96 midpoint of that validated
log-amplitude interval, while `default_log_amplitude_enclosure` exposes the interval after the
selected positive `Lambda` amplification and symbolic `C` exponent transport. The explicit
`log_amplitude_enclosure(eta, absolute_log_tolerance=...)` path requests a log-space tolerance and
passes the exact quotient `absolute_log_tolerance/Fraction(Lambda)` to the phase integrator.
`phase_samples` and `legacy_log_amplitude_diagnostic` retain the former 4001-point trapezoid only
as a named diagnostic; it no longer supplies the default phase value. The historical pre-default
commands returned `5 passed in 0.17 s` for
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
precise runtime is claimed. These are selected-fixture observations, not a general benchmark.
The interval certifies the chosen rational phase kernel and its selected scalar `Lambda` transport;
it does not certify SchedulePressure quadrature, theorem parameter selection, `C` selection,
coefficient roundoff, global `AxisSpace`, or the full reconstruction. A prior combined default attempt
was interrupted after approximately `9` minutes without a pytest summary and remains historical; the
bounded command above resolves the default-regression performance blocker. Full-suite, demo, and
strict-audit acceptance remain separate. The current metadata chain reaches the assigned remainder,
second-Picard, formal-solver, and wide-profile surfaces. The next Stage 1 boundary is validated
Bell-factor/logarithm enclosures, general coefficient arithmetic, SchedulePressure and parameter
proofs, and compatible global `AxisSpace` certification.

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

Section 5 includes pressure/Omega recurrence rows, the Eq. (5.7) singular Picard primitive, a
theorem-shaped Lemma 5.1 / Eq. (5.8) Picard-term and complete-tail truncation certificate, Lemma 5.2
compact five-moment repair, a finite `eta`-jet repair adapter that propagates supplied moment and
`p=e_*f` derivatives through exact Leibniz/quotient recurrences and the same fixed repair matrices,
and a function-level repaired-profile adapter that materializes repaired `U_n(X,eta)`,
`E_n(X,eta)`, and `d_eta U_n` from those supplied jets while hard-coding `paper_exact=False` and
never silently inheriting unrepaired pressure. `background_positive_axis.py` fixes the pinned
PositiveAxisSystem positive-order coefficient algebra: from typed
`BaseJet(phi, axial, beta)` / `SourceJet(angular, axial, pressureProduct, omegaQuotient)` inputs it
constructs the actual `A0`, sparse `A1`, and `f_n` formulas at `lambda_n=2nh` and sends the structured
forcing through the existing singular inverse for a genuine theorem-shaped first `G f_n` step.
`background_lower_history_source.py` now removes the need to hand-supply those four SourceJet
scalars once genuine strict lower-order second jets exist: it assembles the `i=1,...,n-1` lower
convolutions, evaluates the preceding `Z_(b-D) Z_b F_(n-1)` term analytically, and uses the landed
regular `Omega_(n-1)/X` path. Those lower-history jets remain upstream dependencies and are not
paper-derived until they come from the genuinely solved/repaired hierarchy. Eq. (5.15) forward
reconstruction from supplied repaired analytic coefficient data, a finite-prefix
SlowBorel/DiagonalScale cutoff-scale constructor, exact finite-prefix
plateau/transition/zero-tail support certification, a theorem-shaped fixed-prefix
truncation/tail-order arithmetic gate, and an exact finite SlowExpansionResidual
recurrence/truncation bridge are also present. The cutoff scheduler preserves the same theorem
inequality when a valid scale exceeds binary64 by using arbitrary-precision integer power-of-two
witnesses and log-space checks; `reciprocal_support_edge` fails closed instead of underflowing such
a support edge to fake zero. Nonzero retained recurrence defects remain explicit, and the first
omitted slow order is factored only after independently established recurrence cancellation.

Sections 6–7 include fixed dyadic geometry, the slow-label/base-jet bridge, executable squared slow
partitions, a constructive auxiliary-slot centers/`r0` witness, the physical-support-to-
`SlotColoring.Adj` bridge, pointwise and uniform phase-normal adapters, the pinned BasePhaseGeometry
frame/damping implication bounds with the corrected family-level
`phaseConstant(M)=normalConstant(frequencyBound(M))`, the exact signed two-slot reference stress-cone
algebra, a fail-closed actual-vs-reference 2x2 covariance perturbation certificate, and a separate
Eq. (7.28) pulse-error-budget implication that turns independently certified Lemma-7.4 ratio/frame
and normalized first-moment inputs into a `C/sqrt(S_*)` normalized-column envelope. The
coefficient-level Section 7 Formula (30) curl algebra records `B=|n|^-2(n cross a)`, the
tangential-projection identity, `inverseCarrier=i/K`, and the displayed `(i/K)curl(B)` derivative
remainder. `cylindrical_curl_jet.py` now narrows the former arbitrary `coefficient_curl` interface:
it constructs the genuine `r>0` cylindrical curl from supplied component derivative jets, including
the `B_theta/r` moving-frame connection term, and the tangent adapter fixes `B` itself from the
pinned `|n|^-2(n cross a)` formula. The derivative jets are still not derived from the actual
localized paper coefficient. A formal-structure primary-amplitude ODE layer records the pinned
MovingFrameODE/GrowingMode/PrimaryODE pointwise algebra: moving-frame `a,b,c`, modal errors, `j^2`
viscosity damping, the exact 2x2 modal operator, projected forcing transform, and `x=p+q`,
`y=h(p-q)` with `h'=rate*h`. Its regression independently compares the physical `rhsX/rhsY` route
with differentiated modal-basis reconstruction. `primary_amplitude_solution.py` also executes the
finite-interval modal IVP for supplied time-dependent datum/forcing and separately computes the
Volterra integral defect; a constant-coefficient specialization has a closed-form oracle, so the
diagnostic is not the production integrator comparing against itself. These adapters are **not**
the actual paper pulse/wave: Proposition 5.5/phase-frame/pulse data, rigorous numerical error
control, analytic localized derivative jets, support/zero-germ/common-glue hypotheses, and
downstream weighted estimates remain unresolved.

Section 10 includes spatial/time localization, the pinned closed-past `zeroBefore` branch, endpoint
Taylor–Borel right-extension infrastructure, and a `close_left_open_past` adapter for upstream
residual evaluators that exist only for `t<T`: it leaves `t<T` untouched, fills `t=T` from the
degree-zero endpoint jet, rejects `t>T`, and preserves the formalization's closed-past branch at the
join without asserting continuity. The full-spacetime-to-normal endpoint-jet adapter, SpatialBorel
scale arithmetic with both per-degree `2^-j` bounds and the exact omitted series budget
`sum_{j>N}2^-j=2^-N`, a value-level CandidateFromLimits traced-residual/Borel glue bridge, the exact
support-cylinder/fixed-time energy implication, a fail-closed late-origin certificate showing final
localization preserves any upstream-certified origin curl on `3/4 <= t < 1`, and an endpoint-limit
majorant that converts an independently proved integrable bound
`sup_x ||partial_t D^n R|| <= C(T-t)^(-alpha)`, `alpha<1`, into a locally uniform endpoint-tail
modulus are present. The endpoint-localization bridge gates that majorant to the paper's endpoint
regime: it requires `T=1` and a start time `>=3/4`, checks `chi=1`, `chi'=0`, and `chi^2-chi=0` on the
late plateau, and only then transfers the upstream residual bound unchanged; transition-collar and
non-`T=1` uses fail closed. `section10_endpoint_ladder.py` now packages degrees `0..N` into one
contiguous finite family: every degree must appear exactly once, all members share one spatial
window, every member passes the official `T=1`, `t>=3/4` gate, and a common finite-order Cauchy/tail
modulus is obtained by maximizing the exact antiderivative budgets. This is a prerequisite interface,
not evidence that the actual residual satisfies the bounds or that the infinite all-order family is
closed. A separate endpoint-force support implication checks exact vanishing outside the official
cylinder branch-by-branch: past residual for `t<T`, degree-zero endpoint trace at `t=T`, and exactly
the locally finite active future normal coefficients for `t>T`. It uses no support tolerance and
does not turn manufactured or finitely queried data into a global support proof. Reproducible
numerical diagnostics remain explicitly diagnostic rather than proof.

**Not yet delivered:** the materialized coefficient-space fixed-point `phi/u/average/pressure`
fields and complete regular leading profile; the actual-schedule wide `Lambda/C`, scalar contraction
gates, zeroth-parameter-jet `referencePair` radial coefficients, and coefficientwise
`naturalRemainder(x0)`, `x1`, `naturalRemainder(x1)`, and `x2` adapters are now connected, and a
generic mixed-scale algebra/triangular solver now computes requested formal coefficients directly
from lower radial rows without a fixed iteration cutoff. Weighted-space membership and fixed-point
convergence remain unavailable. Finite `jet_prefix`/`formal_axis_profile_prefix` polynomial
evaluation, radial-average/derivative maps, and the scaled formal pressure primitive
`P=primitive((a phi)^2)` are available; the axial remainder pressure forcing remains a separate
quantity. The conditional radial-tail majorant does not establish global membership or cover
coefficient roundoff. The reviewed `wide_natural_profile_prefix` now retains Decimal Lambda and
symbolic amplitude in `F`, `E`, `U`, `d_eta U`, `Ubar`, and `Pi`; the reviewed conditional `F/P`
norm-budget extension uses the same actual chain. Its optional six-field truncation mapping passed
the latest 4-test exact-scaling, square-root, axis-tail, and solver-mismatch acceptance check. The
reviewed wide prefix supplies the physical assembly maps; a converged paper-exact profile remains
unavailable, and
possibly tighter theorem-side conservative majorants are also still needed; the profile-derived
converged all-order background solve, genuine strict-lower-history second jets feeding the landed
SourceJet bridge, actual positive-order Picard constants, hierarchy-derived moment/`p=e_*f` jets
with interval nonvanishing control, true eta-dependent repaired hierarchy/support closure, uniform
`C[j,m]` bounds, independently proved residual-cancellation identities, and completed infinite
recursive cutoff/all-jets-flat argument; the paper-exact Proposition 5.5 base-field provider and
actual LocalBase/vector hypotheses on every active box; the actual Eq. (7.27) pulse-integrated
covariance columns and genuine Lemma-7.4/Gaussian-concentration inputs needed to instantiate the
landed Eq. (7.28) budget and perturbation certificates, followed by actual paper-data instantiation
and rigorous certification of the landed finite-interval primary-amplitude solution and analytic
derivation of the localized cylindrical coefficient derivative jets/support/glue hypotheses needed
for a genuine supported-curl wave; mean corrections and the convergent correction sequence; or the
final compact field with actual closed-past residual full-spacetime derivative bounds/limits at every
order on the official late plateau, closure of the infinite endpoint-majorant family beyond the
landed finite ladder, genuine analytic template majorants, a proved endpoint trace/normal-jet match,
proof that the actual residual/jets vanish outside the support cylinder, a force proved smooth and
compactly supported through the singular time, the upstream origin blow-up premise, and the paper's
uniform bounded-energy conclusion along the actual blow-up limit.

A successful demo or a green test suite is **not** a reconstruction of the full counterexample.
The CLI reports `full_reconstruction: false` / `paper_exact_velocity_available: false` and refuses
`--require-paper-exact` requests.

## Run from a clean checkout

Requires Python 3.10+; NumPy and SciPy are installed as dependencies.

```bash
python -m pip install -e '.[dev]'
python -m pytest -q -W error
ns-reconstruct status
ns-reconstruct audit
ns-reconstruct demo --output artifacts
```

`python -m openai_ns_reconstruction` is equivalent to the console command.
The `status` and `audit` commands derive their stage truth from the same fail-closed runtime source,
while retaining their respective output schemas. The demo produces labelled diagnostic/toy
artifacts only; those outputs are not samples of the final counterexample.

To demand a complete paper reconstruction:

```bash
ns-reconstruct audit --require-paper-exact
```

This currently exits **2**, intentionally. It is a conservative completion gate, not a formal
proof checker. `demo --require-paper-exact` also refuses before generating toy data.

## Implemented components

| Component | Implementation and boundary |
|---|---|
| Similarity coordinates, Eq. (4.1) | Relative-scale root solve, finite-input checks, direct `tau` API |
| Leading velocity, Eqs. (4.3)–(4.7) | Caller-supplied profiles, regular-axis handling, pressure evaluation; natural-axis/range formulas, outgoing schedule, actual `axisPressure`, analytic low-|Z| margin, theorem-faithful `Lambda/C` selection, theorem-shaped conservative `remainderBound/remainderLip` propagation, and canonical analytic-input norm/resolvent propagation are present. The actual schedule has an analytic no-grid common-neighborhood certificate producing `rho`, an 11-field complex bound, and a conservative `realPartSup(axisPhase)` upper bound. A theorem-faithful componentwise adapter applies `radiusLoss(1/2)=12` per field, retains the fieldwise maximum only for the family ledger, and sends the chi-specific norm to AxisResolvent. The current regression resolvent majorant is representable. The original binary64 remainder path overflows, while the landed 96-digit Decimal adapter carries the same positive conservative `remainderBound/remainderLip` ledger to finite values beyond float range; the landed wide scale selector then carries those values through `Lambda` and represents `C=exp(Lambda*phaseSup)` symbolically by its exponent without binary64 down-conversion. The landed Picard gate verifies `s*remainderBound<=1`, `s*remainderLip<=1/2`, `s=1/(2Lambda)`, and exposes the theorem-shaped update/tail budget. The landed `axis_reference_pair` path now materializes the pinned angular/axial reference radial coefficients at parameter-jet order zero from the actual schedule, on the full coefficient window. The coefficientwise Stage 1 chain now represents complete `naturalRemainder(x0)`, the first Picard state `x1`, full `naturalRemainder(x1)`, and `x2` with explicit mixed-scale channels. A generic sparse mixed-scale algebra and formal triangular solver now computes requested formal coefficients directly from lower radial rows without a fixed iteration cutoff. The one-graph `jet_prefix`/`formal_axis_profile_prefix` path evaluates finite `Y`-radial polynomial prefixes, radial derivatives, eta-derivative jets, averages, and the scaled formal pressure primitive `P=primitive((a phi)^2)` with split signed-log channels; `wide_natural_profile_prefix` now retains Decimal Lambda/symbolic amplitude in `F`, `E`, `U`, `d_eta U`, `Ubar`, and `Pi`, with reviewed conditional `F/P` norm budgets and a six-field truncation interface. **These are formal coefficient recursions; the weighted-space fixed point, adjacent-eta compatibility/error budget, global `AxisSpace` certificate, interval/Lean proof, and paper-exact leading field remain unavailable. The axial remainder's pressure forcing is not the final profile pressure.** |
| Profile averages | Cached 32-point Gauss–Legendre rule; optional exact-average callbacks |
| Stage 1 formal profile prefix/tail | The generic `jet_prefix`/`formal_axis_profile_prefix` is an x1-anchored triangular formal prefix, separate from the finite first-Picard `x1`/`x2` adapters. It evaluates finite `Y`-radial `phi/u`, radial and eta-derivative jets, radial averages, and scaled formal `P=primitive((a phi)^2)`. Reviewed `wide_natural_profile_prefix` maps retain Decimal Lambda/symbolic amplitude in `F/E/U/d_eta U/Ubar/Pi`, and reviewed conditional `F/P` norm budgets feed the six-field truncation interface. The focused command `python -m pytest -q tests\test_natural_axis_wide.py tests\test_axis_coefficient_profile_prefix.py -W error` returned `5 passed in 0.66 s`; the latest exact-scaling, square-root, axis-tail, and solver-mismatch check `python -m pytest -q tests\test_natural_axis_wide.py -W error` returned `4 passed in 0.64 s`. Global `AxisSpace` identification and coefficient/axis-pressure roundoff control remain open. |
| Exterior heat swirl, Appendix A.6 | Adaptive evaluation of H and derivatives, swirl and centrifugal pressure; **r>0 only** |
| Section 5 background rows | Eq. (5.2) radial flux, Eq. (5.27) streamfunction/vector potential, Eq. (5.5) pressure row, regular Eq. (5.6) `Omega_k/X`, Eq. (5.7) singular inverse/Picard primitive, and a pinned PositiveAxisSystem adapter that constructs the real `A0`, sparse `A1`, and `f_n` from typed BaseJet/SourceJet inputs and evaluates a theorem-shaped first `G f_n` step. The landed lower-history adapter builds the order-zero BaseJet and four SourceJet entries from genuine strict-lower `ProfileSecondJet` history, including analytic preceding double-Z diffusion and the regular `Omega/X` source, so hand-supplied source scalars are no longer required once that hierarchy exists. Also present: the Lemma 5.1 / Eq. (5.8) term majorant and fail-closed complete-tail truncation certificate, Lemma 5.2 compact five-moment repair, finite `eta`-jet and function-level repaired-profile adapters with `paper_exact=False`, Eq. (5.15) forward `F_n/V_n/Pi_n` reconstruction from supplied repaired analytic data, finite-prefix recursive cutoff scheduling from supplied analytic `C[j,m]`, exact finite-prefix cutoff support/truncation certification, exact fixed-prefix tail-order arithmetic for `2^-J q^(h(J+1)-m)` / physical exponent `h(J+1)+b-2M`, and an exact finite recurrence/truncation identity that retains nonzero `R_n` defects and conditionally exposes the first omitted slow-order factor. A valid cutoff scale beyond binary64 is preserved as an arbitrary-precision integer witness with the original log-space inequality; an unrepresentable reciprocal edge fails closed instead of becoming zero. **The strict-lower second jets still need to come from the real solved/repaired hierarchy; actual analytic constants needed to instantiate Picard convergence, genuine hierarchy-derived moment/patch-factor jets and nonvanishing control, the converged coefficient hierarchy, true eta-dependent support/stress closure, true uniform `C[j,m]`, completed infinite schedule/theorem-level local-finiteness argument, independently proved order-by-order residual cancellation identities, and Proposition 5.3 all-jets-flat decay are not yet constructed.** |
| Dyadic geometry / phase, Sections 6–7 | Fixed chart scaling and exact covering matrices; active-shell/slow-label bridge; squared partitions; constructive 2250-color/rational-center/common-`r0` witness; physical slow-support-to-`SlotColoring.Adj` bridge with cross-band common-point handling; Section 7.1 phase/tangent-frame algebra; pointwise rounded-normal adapters; a fail-closed UniformLocalBase implication bridge; and pinned BasePhaseGeometry frame/damping consequence bounds, including the enlarged family-level phase constant. **The paper-exact Proposition 5.5 provider and actual LocalBase C1/C2, unit/orthogonality, normal-closeness and slot-normal derivative-closeness certificates remain missing; the landed implications cannot be promoted until those true hypotheses are supplied.** |
| Oscillatory stress / curl realization, Section 7 | Exact signed two-slot reference covariance solve for target `(-m,t)`, strict cone `|a t| < b m`, positive squared amplitudes, the manuscript ratio specialization, exact `epsilon*mask^2` covariance scaling after amplitude localization, a fail-closed 2x2 perturbation implication that turns independently certified column errors into determinant/inverse-norm/positive-amplitude margins, an Eq. (7.28) error-budget implication `|e_sigma| <= E_ratio + u_*sqrt(1+c_*^2) M1`, coefficient-level Formula (30) curl algebra with exact tangency gating, a cylindrical derivative-jet adapter that computes the real `r>0` cylindrical curl including `B_theta/r` instead of accepting an arbitrary curl vector, the pointwise MovingFrameODE/GrowingMode/PrimaryODE modal algebra with independent physical-vs-modal reconstruction regression, and a finite-interval PrimaryODE execution adapter with an independent Volterra-defect diagnostic and closed-form regression oracle. **The actual Eq. (7.27) pulse functions/integrals, genuine uniform Lemma-7.4 ratio/frame plus weighted first-moment estimates, paper-derived ODE datum/forcing, rigorous ODE/initial-value/weighted solution certification, analytic differentiation of the genuine localized coefficient into the cylindrical derivative jet, and the real smoothness/support/zero-germ/common-glue hypotheses are still missing; caller-supplied algebraic/IVP/derivative data do not establish a paper-exact divergence-free wave.** |
| Section 10 localization | Official spatial/time support/plateau geometry represented by explicit C-infinity bumps, analytic spatial cutoff gradient, time-switch derivative, activated-field adapters and independent residual-identity cross-check; pinned closed-past `zeroBefore/pastVelocity/pastPressure` branch with diagnostic residual evaluation; endpoint Taylor–Borel infrastructure, plus `close_left_open_past`/`glue_to_left_open_past` for a `t<T` evaluator and a separately supplied degree-zero endpoint candidate without claiming continuity; full-spacetime endpoint-jet adapter; exact-rational SpatialBorel scale scheduling from supplied analytic derivative bounds; per-degree `2^-j` derivative-tail certificates and exact conditional all-omitted-degrees budget `sum_{j>N}2^-j=2^-N`; plus a value-level traced-residual/Borel glue bridge for supplied endpoint tensors. The support certificate records the exact radius-`1/4`, height-`1/2` cylinder of volume `pi/32`, hence `E(t)<=pi M^2/64` from an independently certified fixed-time `|u|<=M`. A late-origin certificate preserves an already-certified incoming origin curl on `3/4<=t<1`. A separate endpoint-limit majorant proves the elementary implication from an independently certified integrable time-derivative bound `C(T-t)^(-alpha)`, `alpha<1`, to a uniform Cauchy/endpoint-tail modulus. The official endpoint-localization transfer requires `T=1`, starts no earlier than `3/4`, checks `chi=1`, `chi'=0`, `chi^2-chi=0`, and then carries that independently proved majorant unchanged; it rejects transition-collar inputs. The finite `section10_endpoint_ladder` requires contiguous degrees `0..N`, one common spatial window, and the same official gate for every degree, then returns a common finite-order modulus; missing/duplicate degrees and mixed windows fail closed. A pointwise endpoint-force support certificate requires exact zero on the past/endpoint/locally-active-future branches outside the official cylinder and never infers support from tolerances. **Transition-collar point values are not claimed equal to Mathlib's noncomputable bump; the endpoint candidate is not yet proved to be the actual left limit; actual closed-past residual derivative bounds/limits and zero-fiber facts for all true endpoint jets, closure of the infinite all-order ladder, true analytic template majorants that instantiate the tail certificate, all normal-jet matching, globally compact smooth force gluing through `t=1`, the upstream origin blow-up premise, and a uniform bounded-energy estimate along the actual blow-up limit are still missing.** |
| Verification | Independent manufactured-solution and refinement checks; bounded-time finite differences |
| Audit/provenance | Source-pinned ledger, explicit blockers, synchronized `status`/`audit` truth surface, fail-closed completion gate |

### Near-singularity evaluation

Use `tau=1-t` directly when working close to the singular time. Forming `t=1-tau` in binary64
can round to 1 and irreversibly lose small positive tau.

```python
from openai_ns_reconstruction import similarity_coordinates_from_tau
from openai_ns_reconstruction.heat_exterior import HeatExterior

s = similarity_coordinates_from_tau(r=1e-50, z=0.0, tau=1e-100, h=0.005)
print(s.q, s.X)

exterior = HeatExterior(h=0.005, c_inf=1.0)
print(exterior.swirl_from_tau(r=1.0, tau=0.2))
print(exterior.pressure_from_tau(r=1.0, tau=0.2))
```

`HeatExterior` implements a parameterized published exterior component, **not** an admissible
smooth-axis `LeadingProfile` and not the full field. The numerical code accepts `0<h<1/2` where
the coordinate formulas make sense; the paper's complete construction imposes the stricter
`0<h<1/100` and additional parameter conditions.

## Verification boundaries

Finite differences, quadrature, and sampled convergence are diagnostics, not interval bounds or
Lean certificates. Setting `f=R(u,p)` and checking `R-f` with the same stencil is tautological;
our regression suite additionally uses independently specified manufactured forces. A generic
`B(x,y,z,t)e_theta` need not be divergence-free: the direct swirl must be axisymmetric, and
localization must preserve that property. Prefer the `LocalField.from_axisymmetric` and
`LocalizedField.from_axisymmetric` factories.

The regular inner core and materialized coefficient-space fixed point, and weighted-space
membership/convergence remain unresolved. The generic mixed-scale algebra is implemented, and the
formal triangular solver now
computes requested formal coefficients directly from lower radial rows without a fixed iteration
cutoff. `jet_prefix`/`formal_axis_profile_prefix` evaluates finite `Y`-radial polynomial prefixes,
radial and eta-derivative jets, radial averages, and the scaled formal pressure primitive
`P = primitive((a phi)^2)` while preserving split signed-log channels. The reviewed wide profile
now exposes sparse `F`, `E`, `U`, `d_eta U`, `Ubar`, and physical `Pi` maps, together with a
conditional `F/P` norm budget and six-field truncation interface. Its latest focused acceptance
check passed `4` tests in `0.64 s`; these interfaces remain conditional on compatible global
`AxisSpace` identification and do not absorb coefficient or axis-pressure quadrature roundoff.
These coefficientwise stages are formal local recursions; they do not establish a complete
compatible weighted coefficient state or a converged paper-exact profile. The axial remainder
pressure forcing is not the final profile pressure. Real
leading/lower-order instantiation of the Section 5
BaseJet plus the landed strict-lower-history SourceJet adapter from genuinely solved/repaired second
jets, converged hierarchy and actual Picard constants, genuine hierarchy-derived moment/patch-factor
jets and interval nonvanishing control, true uniform coefficient-template bounds and completed
infinite recursive cutoff schedule/residual cancellation, paper-exact base-field/slow-box and vector
hypotheses, actual pulse-integrated covariance columns and genuine Lemma-7.4/concentration inputs
needed to instantiate the landed Eq. (7.28) analytic budget and finite-dimensional perturbation
layer, actual-paper-data certification of the landed finite-interval primary-amplitude solution,
analytic derivation of the genuine localized cylindrical coefficient derivative jets plus
support/zero-germ/common-glue closure, mean correction, infinite-order summation, completed compact
field, actual closed-past residual full spacetime jets and independently proved all-window endpoint
derivative bounds/locally uniform limits on the official late plateau, closure of the infinite
all-order family beyond the landed finite endpoint ladder, proof that the degree-zero endpoint
candidate used by `close_left_open_past` is the actual left limit, genuine analytic compact-template
derivative majorants, all normal-jet matching, proof that the actual residual and endpoint
coefficients have the required zero fibers outside the official support cylinder, smooth compact-
force extension, upstream blow-up-path premise, and uniform bounded-energy conclusion remain
explicit blockers. The endpoint-limit majorant, finite endpoint ladder and official `T=1`,
`t>=3/4` transfer are conditional sufficient-condition adapters and do not infer derivative bounds
from samples. The endpoint point-support adapter is also conditional: exact zeros at a queried point
do not prove the actual field's global compact support. The late-origin localization identity
preserves a proved incoming origin divergence if one is supplied; it does not infer blow-up from
finite samples. The fixed-time compact-support energy implication likewise does not replace the
missing uniform energy estimate. The exact SpatialBorel `2^-N` omitted-series budget is conditional
on genuine all-order template bounds and is not itself a smooth-force proof. See
[the reconstruction plan](docs/RECONSTRUCTION_PLAN.md),
[the measured validation report](docs/VALIDATION_2026-09-10.md), and
[the machine-readable manifest](references/provenance_manifest.json).

## Sources

- Paper: https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf
- Official Lean source: https://github.com/openai/NavierStokesAndEuler
- Pinned upstream commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Announcement: https://openai.com/index/navier-stokes-solution/

Pinning the source commit does not mean its Lean code has been built or its theorem-to-runtime
mapping verified here. The paper was inspected in rendered pages; its binary hash has not been
computed. These limitations are recorded rather than filled with guessed provenance.
