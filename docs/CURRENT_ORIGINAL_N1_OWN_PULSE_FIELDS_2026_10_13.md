# Original n1 pulse C3 functions and own physical fields

The selected original pulse coefficients now have genuine ordinary axial
derivatives through order three. They feed physical main, gap and end pulse
fields with their own complete cumulative means and radial recovery. The
checked interfaces restore saved ancestors without rerunning pressure grids.
This milestone, dated 2026-10-13, supplies the fields needed for both pulse
pressure moments. The full original goal remains **ACTIVE / INCOMPLETE**.

## Checked interfaces

The two source quartets are
`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_pulse_C3_functions`
and
`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_pulse_own_fields`
(`.py`, `.json.gz`, `_check.py`, `_check.json`). Their local gates are
`current_original_n1_actual_future_and_pulse_C3_functions_installed` and
`current_original_n1_main_gap_end_own_velocity_mean_Q_functions_installed`.

```powershell
python scripts/run_original_n1_pulse_own_fields.py
```

`CurrentOriginalN1PulseOwnFields()` verifies the saved receipts, source
hashes, source family, accepted graph prefixes and actual function bodies.
The common graph has 34,295 nodes. Its `.jets` owner exposes genuine C3
native, heat, history, angular, future and selected pulse functions. Its
`.fields` dictionary contains `main`, `gap` and `end` chart interfaces:

| Interface | Meaning | Ordinary axial rows |
| --- | --- | --- |
| `Utheta` | Physical swirl velocity | 0 through 3 |
| `U`, `U_y` | Physical axial velocity and log-radius derivative | 0 through 3 |
| `m`, `m_y` | Own physical axial mean per R and its log-radius derivative | 0 through 3 |
| `normalized_linear_histories` | Both original partial moment histories | 0 through 3 |
| `Q`, `Q_y` | Auxiliary radial quantity and log-radius derivative | 0 through 2 |
| `Ur` | Physical radial velocity | 0 through 2 |
| `absolute_swirl_base_pressure` | Original pressure datum plus swirl history | 0 through 3 |

`Q=Ur/sqrt(R/2)` is exactly the auxiliary quantity called `q` by the
existing n1 pressure APIs. Feed `Q` directly to `nonradial_integrands`,
`endpoint_functions` and `pressure_functions`. The separate physical ratio
`Ur/R=Q/sqrt(2R)` must not introduce another factor into those APIs.
`absolute_swirl_base_pressure` is the leading swirl/base pressure, not
either n1 pressure forcing moment `m13` or `m15`.

## Original histories and physical units

Let `y=log(R/Rp)`, `bp=.5+mu`, `lambda_i=.5-i*mu`, and
`Etheta=Pstar*u_native*exp(-bp*y)`. The native frame is already in the
accepted common normalization; physical Pstar multiplies it once. With the
actual selected functions `ap`, `C1`, `C2` and
`endScale=exp(log_end_scale)>0`, the three charts retain:

```
main: xi=mu*y, 0<=xi<=11
B=ap*gp(xi)
mi=native_mi*exp(-lambda_i*y)
   +ap*integral_0^y exp(-lambda_i*(y-v))*gp(mu*v) dv

gap: 11<=xi<=13-4*mu
B=0
mi=-sum_j Cj*row_ij*exp(log_end_scale+lambda_i*(13-xi)/mu)
row_1=A, row_2=A+mu*D

end: s=log(R/Rv), -4<=s<=0, y=13/mu+s
B=endScale*sum_j Cj*beta_j(s)
mi=-endScale*sum_j Cj*integral_s^0 exp(lambda_i*(v-s))*beta_j(v) dv
```

The finite-mu row `A+mu*D` is the actual second row, not its zero-mu limit.
The gap has zero local axial velocity while its own mean and radial
velocity generally remain nonzero. Main forward history and end backward
history are two defining representations joined by the selected original
linear moment equations. End support is empty at s=0, giving exact terminal
zeros for the axial pulse and these two histories.

The physical fields are

```
U=Etheta*B
U_y=Etheta*(B_y-bp*B)
m=Etheta*m1
m_y=U-m
Q=[2Z*U-(1-delta)*Z*m-(1-Z^2)*m_Z]/(1-delta*Z^2)
Q_y=the same linear recovery applied to U_y and m_y
Ur=sqrt(R/2)*Q
```

Because `bp+lambda_1=1`, the physical mean satisfies `m_y=U-m`.
Its third axial derivative supplies the second derivative of Q. A missing
fourth mean derivative is not padded with zero; Q has only its genuine C2
interface. GP and beta radial derivative nodes carry the accepted flat
source AST bindings and explicit exact continued formulas.

Pressure preserves

```
P=Pstar^2*[native_P0+native_Mp
         +native_u^2/2*(1-exp(-(1+2mu)*y))/(1+2mu)]
```

through each pulse boundary. It is not reset when axial support ends.

## Gamma derivatives and focused evidence

The C3 owner preserves every accepted value and first-Z handle. It imports
all eight ordinary C3 native jets, rather than reading factorial-divided
Taylor rows as ordinary derivatives. With `a=delta/2`, `S=1/Rtail>0`, and

```
M_k=1/Gamma(1+a)*integral_0^infinity
    exp(-v)*v^(a+k)*(1+2*(1-Z^2)*S*exp(-t)*v)^(-a-k) dv
```

the continued higher heat-defect derivatives are

```
Dhat_ZZ=-4*exp(-t)*M1-16*(a+1)*Z^2*S*exp(-2t)*M2
Dhat_ZZZ=-48*(a+1)*Z*S*exp(-2t)*M2
          -64*(a+1)*(a+2)*Z^3*S^2*exp(-3t)*M3
```

These remain finite at Z=+/-1. Positive Gamma moments through `v^(a+3)`
dominate the differentiated v integrals. The theta exterior decays with
strict positive rate a; pressure and energy tails decay faster. Both
improper integrations retain their entire domains and exact nonzero S.
The DAG differentiator rejects unseeded sources and unsupported moving
axial endpoints.

Focused checks compare Gamma derivatives to independent differentiation
of the original expectation and compare mixed C3 calculus to symbolic
derivatives. The own-field fixture uses `mu=.07`, `delta=.3`, `Pstar=2.3`
and `logRp=1.2`, with nonconstant true C3 inlet/coefficient functions.
An independent two-row matrix solve and full radial ODE verify both mean
histories and ordinary axial rows through order three. An independent
Taylor quotient verifies Q and Q_y through order two. Maximum normalized
differences are approximately `1.54e-10` for Q/Q_y/Ur, `1.29e-10` for the
physical mean, and `1.44e-13` for joins; terminal ODE mean error is below
`6.11e-14`. The checks also retain nonzero gap history, exact empty-support
end zeros, the pressure datum, and one nonunit physical scale.

These are finite mathematical fixture errors. The actual original
microscopic graph parameters remain unchanged. Uniform original numerical
accuracy, full n1 pressure/five-moment repair, stress and temporal recursion
are separate open gates.

## Completed work

- [x] `N1-PULSE-TRUE-C3-FUNCTIONS`: actual Gamma/future/angular/pulse
  ordinary Z2/Z3 functions with accepted value/Z handles preserved.
- [x] `N1-PULSE-OWN-FIELDS`: original main, inactive gap and compact end
  fields consume the actual selected coefficient functions.
- [x] Carry both native incoming histories through all three charts;
  keep finite mu and the exact nonzero scaled end coefficients.
- [x] Recover the own physical mean, Q, Q_y and Ur with correct derivative
  orders and one physical Pstar multiplication.
- [x] Retain the original cumulative swirl/base pressure separately from
  the n1 pressure forcing moments.
- [x] Pass independent complete ODE, two-row matrix, nonunit scale,
  C3/C2 derivative, gap, join and terminal-support checks.

## Next executable tasks

Complete tasks in dependency order. A task is done only after its source,
result archive, focused checker and checked receipt are committed together.
Update the checkbox and link its concrete result; never mark a generated
graph as accepted merely because it was saved.

- [ ] `N1-PULSE-PRESSURE-INTEGRANDS`: call the existing
  `physical_nonradial_pressure_fields.nonradial_integrands` on each
  chart's genuine physical U C3 and auxiliary Q C2. Preserve all ten
  signed integrands and their true first-Z rows. Use Q_ZZ for derivatives
  of q_Z products; do not request or fabricate Q_ZZZ.
- [ ] `N1-PULSE-PRESSURE-MEASURES`: install the full main/gap/end I0 and
  I1 integral bodies with `dR=(R/mu)*dxi` on main/gap and `dR=R*ds` on
  end. I1 uses one additional R. Keep huge radii as exact log bindings
  and retain all signed cross-products, including nonzero gap terms.
- [ ] `N1-PULSE-PRESSURE-ENDPOINTS`: construct the actual finite
  endpoints using the unchanged auxiliary Q convention. Telescope
  matching internal endpoints. Retain the native Rp data and signed
  viscous terms, and derive terminal zeros from empty support.
- [ ] `N1-PULSE-PRESSURE-FUNCTIONS`: reuse
  `preRp_moment_functions.pressure_functions` for both complete m13/m15
  C1 functions. Add a direct physical-integral fixture and independent
  axial derivative comparison. Keep base swirl pressure distinct.
- [ ] `N1-PULSE-PRESSURE-BOUND-TRANSFER`: prove that the new integral
  bodies are the same original pulse functions covered by the existing
  C5 and signed pressure receipts. Transfer bounds as bounds, retain
  nonzero error budgets and avoid another historical grid refinement.
- [ ] `N1-CORE-FIXED-POINT-FUNCTIONS`: bind the admitted nonlinear core
  fixed point to its real map and enclosing evaluation oracle. Preserve
  all twenty nonlinear terms and the anchored gauge primitive. The
  physical core Uz is `4Z+j+epsilon_core*Psi`, without extra Pstar.
- [ ] `N1-LOWER-COLLAR-FUNCTIONS`: consume that same core in the true
  Rin-to-r_minus collar and its own nonlocal means. Recover genuine
  first-Z pressure fields and finite lower endpoint derivatives.
- [ ] `N1-INNER-CUTOFF-FUNCTIONS`: apply the original solenoidal cutoff
  to the streamfunction/potential; bind the true cutoff velocity and
  own cumulative M/V/P. Its tiny bounds are not selected field values.
- [ ] `N1-FULL-FIVE-DEFECT-FUNCTIONS`: compose each complete signed
  Rin-to-infinity defect once from inner/collar, compact/active/quiet,
  pulse and exact heat pieces in the same accepted source graph.
- [ ] `N1-ACTUAL-LINEAR-REPAIR`: feed those C1 defects into the installed
  positive-order linear operator; select fresh alpha/beta functions on
  its disjoint reserved supports and retain signed inverse actions.
- [ ] `N1-REPAIRED-OWN-HISTORIES`: recover the corrected family's own
  M/V/P over the whole radius range and prove five terminal zero
  identities as axial functions, including their required derivatives.
- [ ] `BACKGROUND-HEAT-STRESS`: verify compatible analytic pressure,
  exact heat matching, finite-energy tail and global admissible stress
  with the separately identified flat remainder.
- [ ] `TEMPORAL-N2-AND-HIGHER`: implement actual n-dependent recovery,
  independent moment repair, shared inner domain and smooth summation.
- [ ] `OSCILLATORY-QUADRATIC-CANCELLATION`: implement both pulse
  families and mean corrections; verify their averaged flux cancels
  the admissible stress before claiming a corrected NS field.
- [ ] `CORRECTED-CARTESIAN-VALIDATION`: check the complete nonzero
  three-dimensional time-dependent field independently, with forcing
  fixed or restricted in advance and both full residual norms <1e-3.

The final five-moment identities, global stress/remainder, true temporal
recursion and corrected Navier-Stokes validation remain open.
