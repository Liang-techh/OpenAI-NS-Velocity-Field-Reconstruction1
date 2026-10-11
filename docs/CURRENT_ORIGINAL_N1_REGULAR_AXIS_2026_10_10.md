# Genuine n=1 regular operator and coupled axis solution — 2026-10-10

The current original hierarchy has progressed from known forcing to its first positive-order axis solution. Temporal hierarchy n=1 is fixed; a separate radial Taylor degree runs from one through three. F1, Uz1, K1 and P1 have original-source finite axis coefficients through rho^3, with ordinary Z derivative orders 2, 1, 0 respectively. Regular V1/R is recovered through rho^2. The exact six-coordinate regular PDE operator is directly callable on rho in [0,4.1], Z in [-1,1], including the axis and beyond R_a. The full goal remains **ACTIVE / INCOMPLETE**: the finite axis solution has no certified radial remainder or full-interval field, no order-one moment repair, and no completed temporal recursion.

Source commit: [b2d79f50](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/b2d79f50766aad299255bccbafeed249d82defbe). Source quartet: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_regular_system` with `.py`, `.json.gz`, `_check.py`, `_check.json`.

## What is now solved

All positive-order axis traces vanish. If R coefficients are written F1=aR+bR^2+..., Uz1=cR+d2R^2+..., P1=eR+gR^2+..., the equations force

```text
a  = N1theta(0,Z)/4
c  = N1z(0,Z)/2
e  = N1p(0,Z)
g  = F0(0,Z)*a + (partial_R N1p)(0,Z)/2
```

The code computes in rho=Lambda R and converts ordinary R derivatives by k! Lambda^k. The first rho coefficients thus include epsilon_core=1/Lambda. The complete paper equations determine the next coefficients recursively; they are not copied from the leading fixed point's radial rows. In particular, P1 degree two and three include the original 2F0 F1 coupling, and its pressure gradient contributes an F0^2 sector to Uz1 degree three. That sector is retained separately instead of being rounded into the much larger non-swirl part.

For the positive-order recovery, q1=1+delta and

```text
K1_k = -k/(k+1) * Uz1_k
(V1/R)_k = [(2-q1/(k+1))*Z*Uz1_k
             -(1-Z^2)*partial_Z Uz1_k/(k+1)]/(1-delta Z^2)
```

Here k is radial degree, not temporal hierarchy order. V1 itself begins at degree R^2. No division by an interval containing R=0 occurs. Degree-three V1/R needs an additional Z derivative beyond the current total-four leading source jet and is deliberately absent.

The formal scale A_anchor=F0(Z_anchor) is held constant in jet algebra. The original ordinary mixed derivative rows already include derivatives of the anchored amplitude. Output terms are coefficient times A_anchor^p, with log scale p*logF0; the extreme exponential is never materialized. They are enclosures of finite coefficients forced by the admitted original leading function, not selected parameter values or a substituted Bessel model. JSON cannot hydrate the live source/packet owners.

## Full regular operator for the next solve

Paper (14.6)-(14.8) is transformed exactly to x=sqrt(rho). The state is W=(F1,Uz1,K1,P1,partial_x F1,partial_x Uz1), and

```text
partial_x W + D W/x = B0 W + B1 partial_Z W + g
D = diag(0,0,2,0,3,1)
```

The delivered B0, B1 and g retain original leading functions, delta, all epsilon conversions and the pressure feedback. The axis source is regular. B1 maps the first four coordinates into the last two and annihilates the last two; therefore B1(x) D0 B1(s)=0 and the same identity holds when the right B1 is differentiated in Z. Separated derivative blocks remain active.

The regular inverse to implement for the full solve is

```text
(G_x g)_i = x * integral_0^1 t^D_i * g_i(t*x,Z) dt
W = G_x(B0 W + B1 partial_Z W + g)
```

This checkpoint supplies the actual operator and finite starting coefficients. It does not yet perform or certify that whole-interval quadrature/Volterra solve. The x=0 singular term is a regular-system definition; the API never evaluates W/x at the axis.

## Focused checks and scope

- 54 exact rational polynomial rows solve the direct paper equations independently and check all available radial/Z coefficients, including varying F0 amplitude and pressure feedback.
- 96 directed coefficient rows replay the canonical original source at 800 digits, over the whole axis and at Z=0 and .371.
- 57 independent explicit slope, second-coefficient and third-order pressure-feedback comparisons use interval overlap as a consistency diagnostic. They are separate from the directed replay's containment checks.
- 48 uncollected paper-equation contractions compare the full regular operator on the whole real core, axis, inner point and beyond R_a. These interval-overlap comparisons are consistency diagnostics; they do not alone prove equality or tightness. Formal amplitude test values check algebra only and are never selected as original field amplitudes.
- Owner-issued packet identity, source/context/family guards, out-of-domain input, mutation, and unavailable derivative controls pass. The source quartet and inherited inputs are hash-bound. The public accepted constructor/axis/operator API is exercised separately.

The common holomorphic extension, full n=1 inner solution, radial Taylor remainder, five-moment repair, temporal recursion, regional stress, flat remainder and full corrected NS gates remain false. Do not evaluate this degree-three polynomial on the entire core or assert convergence from its finite coefficients. The existing original pressure 4C/first-interface proof remains closed and is not a new blocker.

One existing read-only worker is reused: **GPT-5.6 Luna / max**. No new worker or Astra child is spawned.

## Live use

Starting from the accepted `CurrentOriginalCoreHierarchySource` owner `hierarchy_source`:

```python
from lei_ren_part1_paper_compliant_current_original_n1_regular_system import CurrentOriginalN1RegularSystem

n1 = CurrentOriginalN1RegularSystem(hierarchy_source)
axis_packet = n1.axis('.371')
operator_packet = n1.system('2', '.371')
axis_report = n1.report(axis_packet)
assert axis_report['current_original_n1_regular_system_and_coupled_axis_jets_installed']
assert not axis_report['n1_full_inner_interval_solution_certified']
```

Use axis([-1,1]) to enclose the whole real axis and system([0,'4.1'],[-1,1]) for the real operator rectangle. Source coordinates must stay inside the exact supported domain.

## Detailed next actions

- [x] **CURRENT-RP-N1-REGULAR-AXIS-INITIALIZATION** Zero axis data, actual first slopes, n=1 average/K/radial recovery, original scales and source ownership are installed. Finite coefficients go through degree three; this task completion concerns initialization only.
- [x] **CURRENT-RP-N1-REGULAR-OPERATOR** Deliver all six regular state rows, known source, D, original n=1 weights, x conversion and derivative block support on the same real inner interval.
- [x] **CURRENT-RP-N1-COUPLED-FINITE-AXIS-SOLUTION** Solve the finite axis equations through rho^3 with the original pressure coupling and its degree-three Uz backreaction. Preserve available Z orders; do not fabricate additional jets.
- [ ] **CURRENT-RP-CORE-COMMON-HOLOMORPHIC-DOMAIN** Bind the admitted core's analytic tube and each required radial derivative majorant to the canonical original family and R_in=4.1/Lambda. Supply a named common Z neighborhood, bounds for the regular matrix/source and their needed derivatives, and a proof that it extends beyond R_a. Close Assumption 14.1 only with those same-source bounds; real rectangle data alone is insufficient.
- [ ] **CURRENT-RP-N1-REGULAR-VOLTERRA-SOLVE** Implement G_x with its exact diagonal kernel and zero axis condition. Propagate all six coordinates and B1 partial_Z W, using analytic jet/Cauchy bounds or an equivalent controlled parameter method. Keep amplitude sectors correlated and retain the original coupled pressure. Produce a typed full-interval F1/Uz1/K1/P1/V1 source before setting coefficient-solved flags.
- [ ] **CURRENT-RP-N1-VOLTERRA-TAIL** Separate quadrature/discretization uncertainty, source uncertainty and omitted iteration terms. Use the non-adjacent derivative-block structure and the paper's factorial simplex/Cauchy estimate to certify the omitted series tail on a smaller common tube. Do not use a finite-grid residual alone or leading radial coefficients as the tail certificate. Report any parameter-dependent loss explicitly.
- [ ] **CURRENT-RP-N1-WHOLE-INTERVAL-EQUATIONS** Check both axis traces and the whole common interval, the regular average identity, pressure recovery and all three hierarchy equations. Supply radial/Z derivative enclosures required for independent residuals. Keep a distinct source/equation/truncation gate so successful assembly does not stand in for a solved field.
- [ ] **CURRENT-RP-N1-FIVE-MOMENT-DEFECTS** From the full n=1 field, compute all five actual cumulative defects as Z functions and bind their source/hash provenance. Retain original pressure and the n-dependent recovery; pointwise fits are insufficient.
- [ ] **CURRENT-RP-N1-FIVE-MOMENT-REPAIR** Use the same five bump controls and a controlled inverse/conditioning bound. Cut streamfunction/vector potential before curl. Verify repaired moments, joins, high derivatives and pressure compatibility as Z-function identities/enclosures. Publish the repaired field and checks.
- [ ] **CURRENT-RP-N2-GENUINE-FORCING** Only after n=1 solution/repair, build the n=2 known sums, shifted viscosity and pressure source from genuine lower hierarchy orders. Use e2=4delta and the correct recovery q2=1-delta+e2 on the same R_in. Repeat regular solve and per-order repair; radial degree two is not hierarchy n=2.

Continue regional cone, actual core widths/material winding, physical norms/finite energy, all-chart/axis support, flatness and oscillatory correction from [the multitime task list](CURRENT_ORIGINAL_RP_MULTITIME_REMAINDER_2026_10_10.md). The long-term objective remains active.
