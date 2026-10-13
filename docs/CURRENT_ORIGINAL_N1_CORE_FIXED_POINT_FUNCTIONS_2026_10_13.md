# Actual original n1 nonlinear core functions

The admitted nonlinear core now has a function interface in the same
graph as the complete original pulse pressure functions. It preserves
the full ten angular and ten axial source terms, the unique analytic
fixed point, the root-anchored amplitude and the compatible native
pressure datum. Its lazy evaluation oracle returns directed ranges
from actual nonlinear radial coefficients plus the complete tail.
It never selects an enclosure center as a function value.

This milestone is dated 2026-10-13. The full original goal remains
**ACTIVE / INCOMPLETE**. The lower collar and inner cutoff, complete
five defects and repaired own histories still precede global matching,
stress, genuine temporal recursion and oscillatory cancellation.

## Interface and restoration

The source quartet is
`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions`
(`.py`, `.json.gz`, `_check.py`, `_check.json`). The local acceptance gate is
`current_original_n1_true_core_fixed_point_functions_and_enclosure_bridge_installed`.

```powershell
python scripts/run_original_n1_core_fixed_point_functions.py
python scripts/run_original_n1_core_fixed_point_functions.py --evaluate --rho 2 --Z .371
```

The first command restores checked definitions and source hashes without
constructing the numerical oracle. The second attaches the checked
actual nonlinear core and reports Phi, Psi and mean_Psi as enclosing
ranges, each with its finite polynomial and full radial tail. There is
no selected numerical point output. Use the saved source/receipt chain;
do not regenerate ancestors or repeat historical pressure grids.

`CurrentOriginalN1CoreFixedPointFunctions()` exposes `.functions`,
`.components`, `.parameters`, `.terms`, `.definition` and `.rho`.
The shared graph contains 35,819 nodes. `.functions` contains genuine
ordinary C3 data for Phi, Psi, mean_Psi, U0, gradient, chi, G, logF0,
F0, P0, U and m; Q has the available ordinary C2 data. Components
denote the unique solution of the admitted analytic equations, backed
by the twenty-term map and Green inverse. They are distinct from their
numerical enclosures.

`component_enclosure(quantity, rho, Z, radial_order=0, axial_order=0,
degree=6, depth=3)` supports the existing admitted domain
`0<=rho<=4.1`, `-1<=Z<=1`, mixed total order at most five and axial
order at most the requested depth. Its cache shares one nonlinear
radial packet among components and derivatives at the same Z.
`.assert_graph()` protects accepted graph prefixes and definitions.

## Exact source, scaling and primitive

The retained core definitions are

```
logP = exp(40)+11
delta = exp(-4logP-30)
logLambda = 4logP+1000
epsilon_core = exp(-logLambda) = 1/Lambda
rho = Lambda*R
U0 = 4Z+j; d = 1-Z^2; L = 1-delta*Z^2
H0 = (1-delta)*Z/2+d*U0; sigma = j/500
gradient = L*H0/(H0^2+sigma^2)
chi = H0^2/(H0^2+sigma^2)
G(Z) = integral_a^Z gradient(v) dv
F0 = exp(-logCstar-Lambda*G)
F = F0*Phi; Utheta = sqrt(2R)*F
Uz = U0+epsilon_core*Psi
m = U0+epsilon_core*integral_0^1 Psi(rho*t,Z) dt
```

Here a is the unique real H0 root in [-j,0]. Certified root disk
centers only guide root certificates. They do not define a, G or F0.
Gamma's separate epsilon schedule, delta/1000, is not epsilon_core.
Physical Uz and m receive no extra Pstar factor.

The source P0 is `Pstar^2*P0_normalized`, using the same fourteen
continuous atoms and analytic flatten function as the native pulse
datum. The physical cumulative pressure is
`P0+epsilon_core*F0^2*integral_0^rho Phi(v,Z)^2 dv`.
The derivative retains both the derivative of F0 squared and the
integral of `2Phi*Phi_Z`. The evaluator does not recover this primitive
by dividing a scaled pressure enclosure by an F0 squared box whose
lower endpoint can be zero.

The actual equations are

```
2(rho*Phi_rr+2Phi_r)+chi*Phi = epsilon_core*Etheta
2(rho*Psi_rr+Psi_r) = B+epsilon_core*Ez
```

All twenty terms compile from the existing actual original operator.
The complete centrifugal pressure and first-Z primitive appear in
the source. The nonlinear recurrence and Green inverse remain the
source of the analytic solution; a leading Bessel/linear model is
not selected as that solution.

## Evidence and limits

The focused checker independently reconstructs the original angular
and axial equations and compares all twenty compiled terms exactly.
A nonunit epsilon polynomial fixture checks radial averaging, the
ordinary mixed derivative factorial conversion and Psi extraction.
Psi removes the affine U0 row structurally before division by epsilon;
subtracting two independent U0 interval boxes would lose correlation.
A separate symbolic fixture checks the pressure primitive and its
complete first-Z derivative.

The actual-source test cold-attaches the nonlinear oracle, verifies
Phi(0,Z)=1 and Psi(0,Z)=mean_Psi(0,Z)=0, and encloses all three
components through Z0, Z1 and Z3 at rho=2, Z=.371. It requires a
strictly positive, nontrivial Phi enclosure. Each output retains the
full analytic radial tail; mean_Psi conservatively uses the raw Psi
tail because averaging cannot increase that bound.

The attachment guard verifies exact graph recipes, directed
containment by the actual source, same admitted j and selected
logCstar, sigma=j/500, the same amplitude/core/rebuild objects and
the unique H0-root certificate. The native P0 callable AST, fourteen
atoms and datum hash must also agree. Interval overlap alone is
never accepted as equality of two source functions.

These checks install the function/enclosure bridge. They do not
claim uniform point accuracy for the original microscopic core,
complete lower-collar/cutoff functions, five-moment closure or a
complete three-dimensional corrected Navier–Stokes solution.

## Next executable tasks

- [x] Bind actual Phi/Psi/mean_Psi, the twenty-term operator and
  the admitted unique analytic solution in the pulse graph.
- [x] Retain the unique-root anchored G/F0, compatible native P0
  and inverse-Lambda core epsilon.
- [x] Recover nonlinear coefficient enclosures with complete tails
  and structural affine-row removal.
- [x] Guard exact source/parameter identity and provide lazy restore
  and explicit numerical-enclosure commands.
- [x] `N1-LOWER-COLLAR-FUNCTIONS`: substitute rho=Lambda*R into
  these same functions on Rin<=R<=r_minus. Recover the prescribed
  collar velocity, cumulative mean and pressure without resetting
  incoming histories. Preserve physical/logarithmic Jacobians and
  genuine lower Uy/UyZ and pressure endpoint terms.
- [ ] `N1-INNER-CUTOFF-FUNCTIONS`: apply the original solenoidal
  cutoff to the streamfunction/recovery structure. Build its own
  M/V/P and both signed n1 pressure functions. Keep tiny cutoff
  caps as error bounds rather than selected function values.
- [ ] `N1-PULSE-PRESSURE-BOUND-TRANSFER`: tie the installed complete
  pulse function bodies to the existing same-source C5 and signed
  pressure receipts, retaining their nonzero error bounds.
- [ ] `N1-FULL-FIVE-DEFECT-FUNCTIONS`: compose core/cutoff/collar,
  compact/active/quiet, pulse and exact exterior once in this graph.
  Check all five signs, units and cumulative radial histories.
- [ ] `N1-ACTUAL-LINEAR-REPAIR`: use the installed positive-order
  operator to select true alpha/beta functions with the prescribed
  fresh disjoint supports. Rebuild the repaired own M/V/P.
- [ ] Prove the five terminal identities as axial functions and
  establish compatible matching, radial energy tail and heat exterior.
- [ ] Build admissible stress and the separate flat remainder across
  all regions before genuine n1/n>=2 temporal coefficient recursion.
- [ ] Implement realizable mean/oscillatory stress cancellation,
  then independently validate full corrected Cartesian NS residuals
  with fixed/restricted forcing and both norms below 1e-3.

The complete goal and deliverables remain in [PROJECT_GOAL.md](PROJECT_GOAL.md).

The lower-collar bridge is now checked; continue from
[CURRENT_ORIGINAL_N1_LOWER_COLLAR_FUNCTIONS_2026_10_13.md](CURRENT_ORIGINAL_N1_LOWER_COLLAR_FUNCTIONS_2026_10_13.md),
including the newly explicit true inner n1 Volterra prerequisite.
