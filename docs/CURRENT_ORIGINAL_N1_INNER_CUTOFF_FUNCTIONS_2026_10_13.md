# True original n1 radial jets, flat cutoff and five inner endpoints

The actual full regular n1 solution now supplies ordinary radial functions
for its two velocities and own mean/K. The original flat cutoff is applied
to F1/Uz1 only. Its own primitive, mean, K, Q/V and pressure are recovered
through the complete core-to-first history. All five signed cumulative
moment functions and their common-inner Rin endpoints have actual bodies.

This milestone supplies exact functions. It does not execute the original
numerical n1 point oracle, compose the full exterior defects or close the
five moments. The full objective stays **ACTIVE / INCOMPLETE**.

## Restore and source interface

```powershell
python scripts/run_original_n1_inner_cutoff_functions.py
```

The two checked source quartets under `experiments/root_st073` are
`lei_ren_part1_paper_compliant_current_original_n1_inner_radial_functions`
and `lei_ren_part1_paper_compliant_current_original_n1_inner_cutoff_functions`
(each `.py`, `.json.gz`, `_check.py`, `_check.json`). Restore saved ancestors;
do not run the old nonlinear hierarchy or historical pressure-grid suites.

`CurrentOriginalN1InnerRadialFunctions` exports `.rows` for F1/Uz1/mean1/K1
through rho order two and ordinary Z order three, `.axis` for their exact
regular limits, `.recovered` for Q1/V1 through rho and Z order two, and
`.pressure_rho` for the first pressure derivative. It retains 52,487 graph
nodes. Second pressure derivatives and odd-component rho-axis jets are not
provided. Higher mixed source derivatives denote analytic functions; they
do not extend the old numerical core oracle's supported derivative orders.

`CurrentOriginalN1InnerCutoffFunctions` exports `.fields`, `.integrands`,
`.cumulative`, `.endpoint`, `.cutoff_rows` and `.history_bodies` in the
common graph of 53,165 nodes. `.endpoint` contains m11/m12/m13/m14/m15 C3
functions at first phase sc/4, the physical Rin. These are signed functions,
not selected enclosure centers or the final full five defects.

Local gates are
`current_original_n1_true_inner_even_radial_functions_installed` and
`current_original_n1_true_inner_cutoff_and_five_endpoint_functions_installed`.

## Complete equation and regular axis derivatives

Write x=sqrt(rho), W=(F1,Uz1,K1,P1,F1_x,Uz1_x), and
H=B0 W+B1 W_Z+g for the accepted actual raw six-row operator. For rho>0,

```
F1_rho  = F1_x/(2x)
Uz1_rho = Uz1_x/(2x)
F1_rhorho  = H5/(4rho)-F1_x/(rho*x)
Uz1_rhorho = H6/(4rho)-Uz1_x/(2rho*x)
```

Every B0/B1 term remains. In particular H includes derivatives of the
complete analytic n1 function rather than derivatives of its error bounds.
At the axis, regularity supplies exact limits: F1_rho=g5(0)/8 and
Uz1_rho=g6(0)/4. The rho coefficients of H5/H6, including all even slopes,
their Z derivatives and both odd diagonal contributions, give ordinary
F1_rhorho=H5_rho(0)/12 and Uz1_rhorho=H6_rho(0)/8. No evaluation of a
singular 0/0 formula is used at the axis.

For the own mean, mean1_rho=(Uz1-mean1)/rho and
mean1_rhorho=(Uz1_rho-2mean1_rho)/rho. Their axis values are
Uz1_rho/2 and Uz1_rhorho/3. K1=mean1-Uz1. Pressure uses the regular formula
P1_rho=epsilon_core*(2F0 F1+N1p), retaining the unchanged leading source.

## Original cutoff and own histories

```
phase = sc/4, rho(s)=4 exp(hb*s), R=epsilon_core*rho
q=s/phase
chi=1-sigma(4q-2)
chi=1 for q<=1/2; chi=0 for q>=3/4
```

The original source sigma and its derivatives remain exact. Ordinary
sigma'' is twice its second Taylor coefficient. The rho derivatives are

```
chi_rho  = -4 sigma'/(phase*hb*rho)
chi_rhorho = 4 sigma'/(phase*hb*rho^2)
             -16 sigma''/(phase^2*hb^2*rho^2)
```

F1_cut=chi F1 and Uz1_cut=chi Uz1, with the complete rho product rules.
The core is unchanged, including its exact even axis continuation.
The first-chart mean and pressure are reconstructed as

```
mean1_cut(s) = [4 mean1_core(4)
                + integral_0^s hb*rho(v)*Uz1_cut(v) dv]/rho(s)
M1_cut = R*mean1_cut
K1_cut = mean1_cut-Uz1_cut
Q1_cut = [2Z Uz1_cut-(1+delta)Z mean1_cut
           -(1-Z^2)*mean1_cut_Z]/(1-delta Z^2)
V1_cut = R*Q1_cut
P1_cut(s) = P1_core(4)
             + integral_0^s epsilon_core*hb*rho(v)
               *[2F0(v)F1_cut(v)+N1p(v)] dv
```

The temporal exponent c1=-1+delta gives c1+2=1+delta. The n0 formula has
a different mean weight and must not be substituted here. The primitive
M1 stays constant after the axial cutoff vanishes; mean/K/Q/V still have
their resulting tail. Pressure continues to integrate the leading N1p
after both cut velocities are zero. Removing these tails requires repair.

## Five signed physical moment functions

The densities per drho are

```
m11: epsilon_core*Uz1_cut
m12: 2epsilon_core^2*rho*F1_cut
m13: epsilon_core*(2F0 F1_cut+N1p)
m14: 2epsilon_core^2*rho*(Uz0 F1_cut+Uz1_cut F0)
m15: 2epsilon_core*Uz0 Uz1_cut
       -epsilon_core^2*rho*(2F0 F1_cut+N1p)
```

Core integrals run from zero to rho. The first collar adds the complete
core endpoint and its own integral with drho=hb*rho ds. The m11 primitive
and m13 pressure use their exact FTC forms. The minus sign in m15 is
applied once. It includes the leading pressure source even after velocity
support ends. Its pressure term cannot be dropped as a velocity-only term.

These functions stop at Rin. The full m13 additionally needs
integral_Rin^infinity N1p dR. The full m15 additionally needs
-integral_Rin^infinity R N1p dR. Existing lower/pulse/exterior function
adapters must supply these signed tails once, with their original datum.
The actual partition is Rin-to-r_minus, r_minus-to-Rp (active, compact,
quiet), and Rp-to-Rv (main, gap, end). The accepted leading radial-support
theorem makes N1p zero beyond Rv. Gamma/future functions are upstream pulse
inputs, not an additional pressure segment. Do not count them twice.
After the canonical zero extension of the cut velocities, m11/m12/m14 are
already their complete velocity-support moments; the old numeric
three-lower-target task is replaced by these true endpoint functions.

## Evidence and limits

An independent cubic manufactured solution checks all raw B0/B1 feedback,
ordinary radial and axial derivatives, exact first/second axis limits,
mean/K/Q/V units and pressure first derivative with nonunit delta/epsilon.
It covers 122 exact mixed identities. A separate nonconstant cutoff fixture
integrates its own incoming mean and pressure and all five signed moments;
208 C3/radial identities use nonunit hb/epsilon/delta. Setting cut velocity
to zero still leaves the leading m13/m15 source. The original flat-source
fixture checks keep/interior/zero branches, both rho derivatives and the
physical geometric-mean zero radius; its finite-fixture difference is
below 1.88e-84. None of these errors claim original n1 point accuracy.

Actual-source/archive guards and checked CLI restoration bind the complete
original functions and exact nonzero hb/sc. Uniform numerical n1 accuracy,
full exterior moment composition, repair and the global NS construction
remain unproved.

## Next executable tasks

The full signed C1 defect functions, selected five C1 coefficients and
physical C1 repair velocity patches are now installed in
[CURRENT_ORIGINAL_N1_FULL_DEFECT_LINEAR_FUNCTIONS_2026_10_13.md](CURRENT_ORIGINAL_N1_FULL_DEFECT_LINEAR_FUNCTIONS_2026_10_13.md).
The algebraic value/Z1 response identities passed; quantitative coefficient
bounds, higher derivatives and repaired own terminal fields remain open.

- [x] Recover true even n1 rho1/rho2 functions and regular axis limits.
- [x] Apply the original flat cutoff to F1/Uz1 with complete rho jets.
- [x] Recover own mean/M/K/Q/V and pressure, preserving incoming histories.
- [x] Export all five signed cumulative C3 functions and Rin endpoints.
- [ ] `N1-NUMERICAL-FUNCTION-EVALUATOR`: attach a same-source directed
  evaluator for the full regular functions. Prove the exact physical
  conjugation, parameter/P0/amplitude and coordinate identities before
  reusing normalized enclosures. Retain every complete tail and Cauchy
  allowance; never define a field by enclosure centers.
- [x] Replace the three velocity-support targets with true m11/m12/m14
  endpoint functions; do not consume serialized numeric signed sectors.
- [x] `N1-EXTERIOR-PRESSURE-FUNCTIONS`: bind the exact lower collar,
  preRp active/compact/quiet and pulse main/gap/end m13/m15 functions,
  retaining the beyond-Rv zero-support proof and all physical measures.
- [x] `N1-FULL-FIVE-DEFECT-FUNCTIONS`: add these Rin functions to the lower,
  preRp and full-pulse m13/m15 contributions once; beyond Rv is exact zero.
  Verify endpoint partition and original signs; do not negate m15 again.
- [ ] Extend genuine pressure defect derivative functions beyond C1 where
  required for repaired own Q/V and stress. Preserve ordinary derivative
  factors and source Cauchy losses; do not differentiate source envelopes.
- [ ] Transfer existing same-source pressure/C5 bounds to the actual new
  function bodies. Include nonzero source errors and higher axial losses.
- [x] `N1-LINEAR-REPAIR-FUNCTIONS`: solve the actual two axial/three angular C1
  coefficient functions with the reserved linear operator; use its true
  determinant/domain data and the complete signed defect functions.
- [ ] Rebuild repaired own M/V/P, check zero axis data and prove all five
  terminal zero identities as axial functions and derivative functions.
- [ ] Complete global smooth matching, finite energy, exact heat exterior
  and admissible stress/flat remainder before true n>=2 temporal recursion.
- [ ] Add oscillatory quadratic cancellation and independently validate
  the corrected Cartesian NS residual with fixed/restricted forcing.

The full objective remains in [PROJECT_GOAL.md](PROJECT_GOAL.md).
