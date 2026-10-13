# Complete original n1 pulse pressure C1 functions

The actual original pulse now supplies both complete pressure forcing
moment functions, m13 and m15, with genuine first axial derivatives on
Rp<=R<=Rv. All main, inactive gap and compact end integrals use the checked
own physical velocity and cumulative mean. This milestone, dated
2026-10-13, completes the pulse C1 function adapter needed for the full
five-defect composition. The full original goal remains
**ACTIVE / INCOMPLETE**; the true core/lower-collar and inner-cutoff
function interfaces still precede actual linear n1 repair.

## Checked interface

The source quartet is
`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_functions`
(`.py`, `.json.gz`, `_check.py`, `_check.json`). Its local gates are
`current_original_n1_Rp_to_Rv_pulse_pressure_C1_functions_installed` and
`current_original_n1_full_pulse_C1_function_adapter_installed`.

```powershell
python scripts/run_original_n1_pulse_pressure_functions.py
```

`CurrentOriginalN1PulsePressureFunctions()` restores the checked own-field
owner and verifies its complete source/receipt chain. The common graph has
35,558 nodes. `.functions` contains m13/m15 C1 pairs; `.pieces` contains
the corresponding three chart functions; `.bindings` exposes the ten
signed products per chart, complete integrals, summed integrals and finite
endpoints. `.endpoint_certificate` records the true incoming mean, exact
selected-row identities and support theorems. `.assert_graph()` preserves
the accepted prefixes, bodies and source hashes.

These are actual source-bound signed function definitions. Their receipt
does not claim uniform numerical pressure accuracy at the original
microscopic parameters or install the complete five defects.

## Physical measures and forcing signs

The pressure API's q is the own adapter's auxiliary
`Q=Ur/sqrt(R/2)`. No factor converting Q into literal Ur/R belongs in this
API. Physical U and Q already contain Pstar once; the pressure builder
does not apply another physical scale.

The ten products are
`I0_q_Z`, `I0_U_q_Z`, `I0_U_squared`, `I0_U_U_Z`,
`I1_q`, `I1_q_Z`, `I1_q_squared`, `I1_U_q_Z`,
`I1_U_squared`, and `I1_U_U_Z`. All expose genuine value/first-Z rows,
using Q_ZZ where a q_Z product is differentiated. Each is integrated over
the entire original chart:

| Chart | I0 physical measure | I1 physical measure |
| --- | --- | --- |
| Main and gap, xi=mu log(R/Rp) | (R/mu) dxi | (R^2/mu) dxi |
| End, s=log(R/Rv) | R ds | R^2 ds |

The inactive gap retains Q, Q_Z and Q^2 contributions from the genuine
nonzero cumulative mean. Zero local U does not reset the pressure history.
Absolute radii remain exact log/exp bindings.

For `y=log R`, `h=1-delta`, `d=1-Z^2`, `L=1-delta Z^2`, the independent
direct original forcing is

```
S_direct=(Q+h*Z*Q_Z/2+Q_y)/L
        +Q*(Q/2+Q_y)
        +U*(-2Z*Q+d*Q_Z-2Z*Q_y)/L
N1p=-S_direct/2+(Q_yy+Q_y)/R
m13=integral N1p dR
m15=-integral R*N1p dR
```

The installed reduced-S operator returns the paper m15 with this minus
already included. Sum its signed pieces directly. A provider returning
the literal weighted integral needs the conversion once.
The cumulative swirl/base pressure is a different field and never defines
these n1 forcing moments.

## Exact endpoints and selected moment identities

At Rp, the flat startup gives U=U_y=U_Z=0 while the native own mean remains
nonzero. Therefore Q_y=-Q and Q+Q_y=0; the native RQ and RQ_y terms are
retained. At Rv, both beta supports are empty at s=0 and every backward
mean integral has identical zero endpoints. Thus U,m,Q,Q_y vanish as
whole axial functions. These support identities simplify the global
endpoints without selecting zeros from numerical bounds.

Internal endpoints telescope by the same selected coefficient functions.
Let I_i be the original scaled incoming rows, Pi_i the scaled main rows,
and C=u+ap*v from the original divided inverse. Exact symbolic algebra
gives

```
A[C]=-mu*(I1+ap*Pi1)
(A+mu*D)[C]=-mu*(I2+ap*Pi2)
```

The finite-mu second row is exactly A+mu D. In the normalized beta measure,
`A1*D0-A0*D1` equals the installed positive absolute determinant
`2*exp(-2+5mu)*J1*J2*exprel(2mu)`. These are pointwise function identities,
so the internal `[Q+Q_y]` and `-[R Q_y]` terms cancel between matching
charts. The final functions use the summed signed integrals and the two
global endpoints; per-chart functions remain available for inspection.

## Independent evidence

The focused checker integrates the original radial differential N1p
directly for a nonconstant own-mean family satisfying `m_y=U-m`. Its
nonunit fixture has Pstar=3/2, mu=3/20, delta=1/5 and physical R endpoints
37/10 and 74/10. All three chart Jacobians and all ten value/first-Z
products are checked. Both finite viscous endpoints are nonzero in this
fixture, so dropped terms cannot pass silently.

Both complete moments and their first-Z derivatives equal the independent
direct integrals symbolically; each chart agrees separately, and summed
chart functions equal the global result. The admitted axial interval is
[-1,1], where `L=1-Z^2/5>=4/5`. Rational cancellation before integration
excludes irrelevant symbolic branches at the domain-external poles
Z=+/-sqrt(5); exact function equality is still required throughout the
admitted interval.

Actual-source checks replay the original pulse selector and preserve its
accepted value/Z handles. They rebuild every product from the genuine own
U/Q, inspect both derivative integral rows and complete domains, retain
nonzero gap terms and verify the explicit source endpoint certificate.
The checked restore CLI and exact staged source hashes complete this
milestone's handoff.

## Next executable tasks

- [x] Bind all ten pulse pressure products and genuine first-Z rows.
- [x] Install complete main/gap/end I0/I1 integrals with physical measures.
- [x] Retain signed finite endpoints and prove internal telescoping and
  original selected moment identities.
- [x] Install complete m13/m15 C1 functions and pass direct independent
  physical forcing, derivative and checked-source tests.
- [ ] `N1-CORE-FIXED-POINT-FUNCTIONS`: bind the admitted nonlinear core
  fixed point to the actual twenty-term map, its analytic limit and
  enclosing evaluation oracle. Keep the anchored gauge primitive and
  `Uz=4Z+j+epsilon_core*Psi` without extra Pstar.
- [ ] `N1-LOWER-COLLAR-FUNCTIONS`: consume the same true core and its
  nonlocal histories from Rin to r_minus; recover actual pressure C1
  fields and finite lower endpoint derivatives.
- [ ] `N1-INNER-CUTOFF-FUNCTIONS`: bind the original solenoidal cutoff,
  its own cumulative M/V/P and signed pressure functions. Tiny caps
  remain bounds and must not become selected function values.
- [ ] `N1-PULSE-PRESSURE-BOUND-TRANSFER`: map the new complete bodies to
  the existing same-source C5/signed pressure receipts. Carry their
  nonzero numerical errors; no further historical grid refinement is
  needed before function composition.
- [ ] `N1-FULL-FIVE-DEFECT-FUNCTIONS`: compose inner/collar, compact/
  active/quiet, pulse and exact exterior pieces once in a common graph.
- [ ] `N1-ACTUAL-LINEAR-REPAIR`: select true alpha/beta functions with
  the installed positive-order operator and fresh disjoint supports.
- [ ] Recover repaired own M/V/P and prove five terminal zero identities
  as axial functions. Then verify matching, finite energy, global stress
  and flat remainder before true n>=2 temporal recursion.
- [ ] Complete oscillatory quadratic cancellation and independent full
  corrected Cartesian residual validation with fixed/restricted forcing.

The detailed later-stage queue remains in
[the own-field handoff](CURRENT_ORIGINAL_N1_OWN_PULSE_FIELDS_2026_10_13.md).
