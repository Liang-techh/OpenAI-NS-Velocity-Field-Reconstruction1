# Actual micro2 to macro0 mixed source interface — 2026-10-14

Status: **ACCEPTED finite source interface**. All seven scoped groups pass
in 137.890 seconds, with 246 actual own mixed row pairs and 240 independent
physical transport/history rows. Independent source-only restoration passes with
exit code 0 and six two-sided source queries in a fresh Python process. The full project goal is
**ACTIVE / INCOMPLETE**.

## What this checkpoint constructs

The original second-micro and frozen-macro charts now have an explicit
coordinate-changing source interface. The micro seam is native `s=2` and
the macro seam is native `m=0`. Their physical radius is the same:

```
R_micro = Ra exp(hb s)
R_macro = Ra exp(2 hb + length m)
length = Y - 2 hb
D_y|micro = hb^-1 D_s
D_y|macro = length^-1 D_m
```

The new macro derivative factory uses `D_m phase = N length`; the accepted
micro factory uses `D_s phase = N hb`. Other integral/angle binders hold the
issued phase fixed. Differentiation is applied to the periodic lifted
composites. A raw fractional-part function is not differentiated across a
wrap.

The actual saved macro `U_y` is differentiated to physical `y4` with six
ordinary Z rows. The original mean ODE `m_y=U-m` supplies genuine higher mean
rows. Leading Q uses `1-delta`, with five Z rows. `N1p`, its first two
physical derivatives, pressure derivatives through `y3`, and signed `m15`
derivatives have four ordinary Z rows. Source-only functions are defined;
original-scale numerical enclosures are not supplied by this interface.

The archived graph has **674777 nodes**. The exposed seam neighborhoods are
second_micro `[15/8,2]` and frozen_macro `[0,1/8]`. These are new endpoint
source scopes. The previous acceptance near native 1 is not silently
treated as an acceptance at native 2.

The second-micro endpoint packet is rebuilt from the original complete
second-micro chart, with its own fresh derivative factory and 18 genuine
high-derivative requests. The previous near-one packet is not copied and
assigned a new domain. A separate check derives both leading Q recoveries
and replays the N/reduced-S operators before those source identities are
used as inputs to the own seam comparison.

## Why the background physical jets match

The checker replays the original macro comparison history, comparison
direction, actual field, actual six-history and normalized generic-field
recipes on the actual saved source references. Macro phi/U initial values
are exactly micro phi/U evaluated at 2. Every actual macro history starts
from the corresponding complete micro history at 2.

The comparison micro smoothing uses `alpha=1-sigma(s-1)`. Its value and first
four native derivatives vanish at 2. Leibniz therefore gives zero comparison
profile derivatives through physical order 5. The macro comparison profiles
are those same endpoint functions and are constant in macro radius. Both
comparison history families obey the same physical ODE:

```
J_y = density - rate J
rates: H=2, M=1, K=2, A=1, B=2, C=1
densities: 2 phi, U, 2 phi U, U^2, phi^2, phi^2
```

In particular the frozen B equilibrium is `phi^2/2`; its forcing is still
`phi^2`. Nonzero inlet values are retained. The actual micro control
`1+(hb-1)sigma(s)` has the same value and required jets as macro `hb` at
the seam. Radius jets agree after the separate Jacobian conversions.

The actual fields then solve the same physical phi/U equations and six
history ODEs with the same inlet. This gives actual background physical
jets through order 5 as functions of Z. Five orders are needed because loop
roots use `E_y` and `V_y`. The common root formulas, positive denominators,
unique monotone inverse, common global logarithmic phase, accepted periodic
lift, and exact total-native calculus then give the leading U jets through
order 4. The original mean history and mean ODE give the mean jets; the
unchanged leading recovery gives Q.

This is a finite source identity argument. It is not a certification of a
global mixed norm or a numerical interval enclosure.

## Pressure and own-state memory

For both `m13` and `m15`, frozen-macro incoming signed integrals are the
complete second-micro endpoint integrals. A macro local integral is zero at
native 0, while its incoming contribution remains nonzero. Both finite
global pressure boundary terms remain in the original formula. No pressure
offset is appended to improve a residual.

The before-Rc support theorem applies to the actual original chart packets.
Their own velocities vanish outside their supports, while incoming own mass
and pressure remain arbitrary nonzero functions. The original recovery is

```
mean1 = M1/R, K1 = mean1
Q1 = (-(1+delta) Z mean1 - (1-Z^2) mean1_Z)/(1-delta Z^2)
V1 = R Q1
```

Consequently mean1/K1/Q1 physical derivatives alternate in sign; M1/V1
derivatives vanish there. This does not zero the leading velocity. The
checker compares the saved own first derivatives and the original complete
mass, pressure and five cumulative histories. High own pressure rows use
the common leading N source as an input, not assigned own output values.

## Reproduction and limits

Inspect or independently restore the accepted source with:

```
python scripts/run_original_micro_macro_mixed_source_join.py
python scripts/run_original_micro_macro_mixed_source_join.py --restore-source
```

The focused checker has seven groups: actual source contract; actual background
recipe replay; independent physical transport recurrence; smooth source
seam; actual leading recovery/common N; actual macro pressure/FTC; actual own
mixed seam. The independent transport fixture derives 240 physical
derivative/history rows from the
original constructors with separate symbolic Jacobians and nonzero inlets.
Both the Duhamel branch and the frozen comparison branch with its B/2
equilibrium are checked explicitly.
The own seam contract contains 246 state/pressure/cumulative mixed rows.

Numerical evaluation/bounds, a common admitted global N, global n1, energy,
axis, analytic preheat pressure, exact heat exterior, admissible stress,
genuine temporal n>=2 recursion, smooth summation, oscillatory cancellation
and corrected Cartesian NS residual gates remain open.

## Detailed next-agent queue

Keep `docs/PROJECT_GOAL.md` authoritative. Claim one task by adding an owner,
branch and starting commit; change `[ ]` to `[x]` only with a commit, check
command/result and scope limitation. A completed source identity does not
complete its numerical or physical layer. Preserve every nonzero incoming
history, the common core/datum/frequency, and actual source bytes.

1. [x] DONE — owner: root, branch `codex/st073-transition-next`, starting commit
   `9b5c0c54`; delivered in the commit containing this checkpoint. All seven
   groups pass in137.890s; 246 actual own row pairs, 240 independent physical
   transport/history rows and674777 archive nodes are recorded in the pinned
   receipt. `python scripts/run_original_micro_macro_mixed_source_join.py
   --restore-source` passes in a fresh process with exit0 and six two-sided
   source queries. Global n1 and numerical gates remain false.
2. [x] DONE — complete second_micro `[1,2]` and frozen_macro `[0,1]` source
   accessors are delivered by [CURRENT_ORIGINAL_TEN_CHART_HIGH_SOURCES_2026_10_14.md](CURRENT_ORIGINAL_TEN_CHART_HIGH_SOURCES_2026_10_14.md),
   along with all ten full original domains. Seven groups pass; original
   phase/Jacobian sources are retained. Quantitative bounds remain open.
3. [x] DONE — genuine first_switch high U/mean/Q and N/pressure sources are
   delivered in the same ten-chart checkpoint with its original moving
   controls, physical measure and complete predecessor histories.
4. [x] DONE — macro1 -> first_switch0 at radius 100 now matches all 246 own
   mixed row pairs. Actual coordinate maps and distinct width nodes consume
   the shared K1 source theorem. Signed m15 pressure memory is retained.
5. [x] DONE — first_switch -> second_switch is closed with 246 actual own
   mixed rows, original flat controls and nonzero histories. See
   [CURRENT_ORIGINAL_SWITCH_MIXED_SOURCE_JOINS_2026_10_14.md](CURRENT_ORIGINAL_SWITCH_MIXED_SOURCE_JOINS_2026_10_14.md).
6. [x] DONE — second_switch -> post_power is closed with 246 actual own
   mixed rows and the actual -2/5 power slope, retaining the separate
   post-length Jacobian. The shared five-group receipt and independent
   four-endpoint archive recovery are recorded in the same checkpoint.
7. [ ] Close post_power -> long_reshape. Complete long_reshape high sources
   are already delivered; retain genuine controls and separate Jacobians.
8. [ ] Close long_reshape -> reference with pressure and all five histories
   as Z functions. Complete reference high sources are already delivered.
9. [ ] Close reference -> restoration and retain the unchanged original
   terminal-repair branch. Complete restoration high sources are delivered.
10. [ ] Close restoration -> postrestore. Complete postrestore high sources
    are delivered; identify exact Rm endpoint histories without replacing
    them by interval centers or forcing target values.
11. [ ] Extend the original Rm patch to genuine y4/C3 state and y3/C3
    pressure rows; prove both patch interfaces with its actual angular
    correction, support and nonzero moment data.
12. [ ] Complete remaining preRc outer chart high sources and source joins.
    Reuse accepted O3/Rc/compact/pulse/outer receipts where they apply rather
    than recalculate identical universal operators.
13. [ ] Assemble a single checked global n1 source provider. Enumerate every
    original interval and interface, reject uncovered queries, and require
    actual row pairs for every derivative contract.
14. [ ] Recover numerical evaluators for the actual high sources from the
    unique analytic core and compatible datum. Do not claim numerical
    restoration from source-only DAG recovery.
15. [ ] Admit one common N with the original Z domain and complete control
    denominators. Use directed intervals and separate truncation/quadrature
    errors; reject callback-selected branches and enclosure-center values.
16. [ ] Validate the global numerical five-moment terminal identities as Z
    functions over a complete domain cover, with an independent operator and
    explicit numerical error budgets.
17. [ ] Recover compatible analytic preheat pressure and exact heat exterior
    consumers, retaining radial tail and finite-energy information. Do not
    append a pressure tail to an inconsistent core.
18. [ ] Assemble physical-coordinate axisymmetric velocity and pressure
    evaluators. Preserve divergence by streamfunction/vector potential;
    check axis regularity, support, units and nonzero amplitude.
19. [ ] Compute admissible stress and flat remainder separately. Report
    cone margins by actual region, maxima and physical-volume L2 norms;
    establish the required flat/high-order decay.
20. [ ] Implement the actual n>=2 coefficient recovery equations with their
    n-dependent constants and the common core definition interval. Distinguish
    temporal coefficients from the radial derivatives in this checkpoint.
21. [ ] Build independent moment repairs for each coefficient order and
    establish finite-order remainder estimates before smooth summation.
22. [ ] Implement the original smooth summation and mean/oscillatory pulse
    correction, then test realized averaged quadratic stress cancellation.
23. [ ] Measure multitime radial/axial widths, velocity/vorticity growth,
    finite energy and material winding using the computed field. Separate
    fitted dynamics from imposed coordinate exponents and instantaneous
    streamline geometry.
24. [ ] Validate the final corrected field with an independent Cartesian
    residual and fixed forcing, varying numerical errors independently.
    Require Linf and physical-volume L2 below 1e-3 only at this final stage.
