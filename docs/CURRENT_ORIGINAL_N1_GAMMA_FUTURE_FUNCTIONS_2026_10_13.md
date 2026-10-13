# Exact original Gamma future C1 functions, 2026-10-13

The following adapter is now consumed by the checked [complete future and selected pulse coefficient functions](CURRENT_ORIGINAL_N1_COMPLETE_FUTURE_FUNCTIONS_2026_10_13.md). Read that newer checkpoint for current completed/open tasks; the queue below records the earlier dependency stage.

The three complete Gamma heat-defect functions now have actual integral
bodies and first axial derivatives in the original n1 repair graph. The
adapter also imports the genuine native angular inlet and logRp, transports
the nonzero pulse memory, and defines the original scalar waiting length
and exact nonzero S=1/Rtail. This removes another function dependency for
the complete future energy and actual pulse coefficient selection.

The full original goal remains **ACTIVE / INCOMPLETE**. These are exact
function definitions with focused numerical fixtures and retained source
receipts. They are not completed five-moment matching, scale recursion,
an original-parameter point evaluation, or a small full NS residual.

## Installed interface

The source quartet is
`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions`
(`.py`, `.json.gz`, `_check.py`, `_check.json`). Its local acceptance gate is
`current_original_n1_exact_Gamma_three_future_C1_functions_installed`.

```powershell
python scripts/run_original_n1_Gamma_future_functions.py
```

`CurrentOriginalN1GammaFutureFunctions()` loads checked saved ancestors,
constructs the small adapter and verifies its own receipt. It does not run
ancestor producers or pressure grids. Its `.g`, `.pulse.g` and
`.pulse.preRp.operator.graph` are the same graph. Use:

- `.functions`: genuine C1 pairs `Dhat`, `theta`, `pressure`, `energy`;
- `.pieces`: collar and complete infinite-exterior C1 pairs separately;
- `.histories`: imported Xp and actual Xv/Xf/Xrel/steep/Xt C1 histories;
- `.schedule`: actual waiting, logRtail, S, delta, epsilon, Ts, Lrel and
  both energy/pressure preheat W and W-squared atoms;
- `.assert_graph()`: preserves every accepted ancestor prefix and the
  new defining bodies while allowing downstream graph appends.

The checked graph contains 31,917 nodes. Current C1 handles are
`Dhat=(31833,31840)`, `theta=(31895,31896)`,
`pressure=(31905,31906)`, `energy=(31915,31916)`; actual S is node 31794.
Prefer API references to hard-coded indices. If relocating this new graph
later, explicitly register operands for `exact_original_source_heat_phi`
and the operand-free `mathematical_positive_infinity`, plus the optional
Gamma normalization reference on an integral. The existing parent importer
does not know those new operations. Same-graph downstream use needs no import.

## Exact definitions and normalization

Let a=delta/2, d=1-Z^2 and x=2dS exp(-t)v. The stable defining integrals are

```
Dhat = 2d exp(-t)/Gamma(1+a)
       * integral_0^infinity exp(-v) v^(a+1)
         log1prel(x) exprel(-a log1p(x)) dv
Dhat_Z = -4Z exp(-t)/Gamma(1+a)
         * integral_0^infinity exp(-v) v^(a+1) (1+x)^(-a-1) dv
```

Here log1prel(x)=log(1+x)/x and exprel(x)=(exp(x)-1)/x,
both with their exact continuation 1 at zero. Thus
`1-H_delta(2dS exp(-t))=a*S*Dhat`; no subtraction of rounded H from 1,
division by a numerical S cap, or finite Gamma series defines the function.

With the original flat sigma and phi(t)=exp(-4/(3-t)^2) for t<3,
zero beyond 3, define

```
C=1-epsilon*phi
pre=1-epsilon*(1-sigma+sigma*phi)
D=sigma*C*Dhat
K=pre-a*S*D
Theta_hat=a*integral_0^infinity exp((1-a)t)*D dt
Pressure_hat=1/2*integral_0^infinity exp(-(1+delta)t)*D*(2pre-a*S*D) dt
Energy_hat=integral_0^infinity exp(-delta*t)*D*(2pre-a*S*D) dt
```

The stored integrals split exactly at t=3. The exterior reaches infinity;
it is not a finite cutoff or the leading asymptotic tail. On the exterior
sigma=C=pre=1. All weights and the pressure factor 1/2 are retained.
The square derivative is `2*K*D_Z`, not `(2pre-a*S*D)*D_Z`.

For 0<=d<=1, positivity of the Gamma measure gives
`0<=Dhat<=2d(1+a)exp(-t)` and
`abs(Dhat_Z)<=4abs(Z)(1+a)exp(-t)`.
These dominate both value and derivative integrals; the angular tail has
positive decay rate a and the other tails decay faster. This justifies
first-Z differentiation across the complete improper integral. At Z=+/-1
all three values vanish, but their first derivatives do not vanish.

The radius recipe retains the actual native origin:

```
Xv=1/(1-mu)+(Xp-1/(1-mu))*exp(-13*(1-mu)/mu)
logRtail=logRp+13/mu+102-30log(mu)+Ts+waiting
S=exp(-logRtail)>0
delta=exp(-4*(exp(40)+11)-30)
epsilon=delta/1000
```

The actual Md=40 selects that strict delta branch of the accepted minimum.
Different directed contexts can give different interval endpoint tuples
for this same recipe; tuple inequality is not evidence of different data.

## Evidence and limits

Focused tests passed: direct original Gamma expectation versus the stable
DAG body, all three complete tails and first-Z rows, symmetry/endpoint
conditions, and independent ODE transport through flatten and both steep
transitions. Native X/logRp imports match the checked frame bodies and all
ancestor graph prefixes are preserved. The checked restore/CLI passed.

Fixtures use ordinary nonzero parameters, including S=.12, so the
quadratic deficit and full exterior have a measurable contribution.
The direct Gamma value maximum discrepancy was about 1.04e-17.
These fixture values are not original-parameter field samples. Existing
accepted Gamma/future source receipts remain bound as numerical enclosure
evidence; their caps are not promoted to defining function values.
The read-only Luna/max review found no missing factor or chain-rule error.

## Next executable tasks

- [x] `N1-GAMMA-DHAT-FUNCTION`: genuine positive Gamma expectation and true DZ.
- [x] `N1-GAMMA-THREE-FUTURES`: complete collar plus infinite exterior C1 pairs.
- [x] `N1-ACTUAL-ANGULAR-HISTORY`: import native Xp/logRp and retain pulse memory.
- [x] `N1-ACTUAL-WAITING-RTAIL`: genuine scalar waiting root and nonzero inverse radius.
- [x] Retain both epsilon and epsilon-squared preheat atoms and original parameters.
- [x] Pass focused Gamma/ODE/first-Z/native-source and checked-restoration checks.
- [ ] `N1-FUTURE-ANGULAR-SMALL-BRANCH`: construct exact A/B/D/E/F beta integrals,
  use these actual heat defects and Xf(0)-Xf(Z), and solve the admitted small
  angular quadratic branch as genuine C1 functions. Do not choose the
  centers of the Newton coefficient boxes. Preserve original P0.
- [ ] `N1-COMPLETE-FUTURE-ENERGY`: compose flatten, reference power, signed
  angular bump energy, steep transitions, steep power, waiting and preheat
  tail. Subtract exactly `Nrel*Ntail/(1-epsilon)^2*(a*S)*Energy_hat` once.
  Nrel must not occur twice. Keep the infinite energy baseline 1/delta and
  both epsilon atoms. Weight the complete total by mu*exp(-26)/2.
- [ ] Bind that exact future C1 function to
  `.pulse.original_coefficient_functions(future)`. Prove the accepted
  positive discriminant and (.9,1.2) amplitude branch using the same source
  before claiming actual pulse selection.
- [ ] Recover main/gap/end axial U, its own nonlocal mean and q, then both
  full pulse pressure moment C1 functions with real radial measures and
  finite viscous endpoints.
- [ ] Bind the true admitted nonlinear core fixed-point functions and their
  enclosure oracle; preserve all twenty nonlinear terms and the anchored
  gauge primitive. Build the lower collar and inner-cutoff C1 functions.
  Physical core Uz=4Z+j+epsilon_core*Psi, with epsilon_core=1/Lambda;
  there is no extra Pstar multiplying this physical axial velocity.
- [ ] Compose all five complete signed n1 defects once, compute the fresh
  linear alpha/beta functions, and recover repaired own M/V/P.
- [ ] Verify five terminal identities and matched heat/pressure/stress.
- [ ] Implement n>=2 order-dependent temporal recovery only after n1 closes.
- [ ] Complete oscillatory quadratic cancellation, flat remainder and
  independent corrected Cartesian residual validation with fixed/restricted forcing.

Later physical stages and the full long-term objective remain open.
