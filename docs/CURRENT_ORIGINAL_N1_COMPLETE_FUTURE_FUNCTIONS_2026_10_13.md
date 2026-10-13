# Actual angular, complete future and pulse coefficient C1 functions, 2026-10-13

The original small angular correction, complete future swirl energy and
selected pulse coefficient functions now share the checked n1 repair graph.
The full future target is supplied to the genuine native incoming functions;
neither Newton-box centers nor numerical kernel centers define coefficients.
A fresh whole-axis interval calculation proves the original positive pulse
branch. This is a prerequisite completed for the full pulse field and its
pressure functions. The full original goal remains **ACTIVE / INCOMPLETE**.

## Checked interface and evidence

The source quartet is
`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_complete_future_functions`
(`.py`, `.json.gz`, `_check.py`, `_check.json`). Its local gates are:

- `current_original_n1_actual_small_angular_C1_functions_installed`;
- `current_original_n1_complete_future_energy_C1_function_installed`;
- `current_original_n1_actual_selected_pulse_coefficient_C1_functions_installed`.

```powershell
python scripts/run_original_n1_complete_future_functions.py
```

`CurrentOriginalN1CompleteFutureFunctions()` consumes checked saved Gamma,
pulse and native-frame ancestors, verifies their receipts and preserves
their graph prefixes. It does not rerun pressure grids or ancestor producers.
The restored common graph has 32,320 nodes. Its API exposes:

- `.angular`: scaled small-root x/y and actual d1/d2 C1 functions;
- `.functions`: complete future swirl energy and the Section 7.34 weighted
  target, both as genuine value/first-Z pairs;
- `.pieces`: all signed energy pieces and original fixed prefactors;
- `.pulse_functions`: selected ap/C1/C2 value/first-Z functions;
- `.pulse_equations`: the original quadratic, target and linear rows;
- `.source_identity`, `.branch_certificate`, `.pulse_branch_certificate`:
  explicit source joins and whole-Z branch evidence;
- `.assert_graph()`: accepted prefix and function-body preservation.

Current energy handles are `complete=(32184,32185)` and
`weighted=(32188,32189)`. Selected pulse handles are
`ap=(32310,32311)`, `C1=(32314,32315)`, `C2=(32318,32319)`.
Prefer API references over hard-coded indices. Pulse Cj are the original
scaled end coefficients: the source end coefficient is
`cj=exp(log_end_scale)*Cj`, with a formal exact nonzero scale. Keep that
factor in the field; its smallness does not authorize replacing it by zero.

The fresh whole-Z amplitude enclosure is approximately
`1.0100453347053863 <= ap(Z) <= 1.0100510646773290`, strictly inside
(.9,1.2). Its discriminant and derivative denominator are positive, and
the two scaled end coefficients have the original negative/positive signs.
These numbers are enclosures of the actual original branch, not selected
point values. First-Z function definitions are installed; a uniform numerical
derivative accuracy claim is not made here.

Focused checks passed: five direct finite-mu bump integrals, nonconstant
angular solutions against an independent implicit matrix derivative, and
direct integration of a complete piecewise swirl profile with physical
radial measure. The latter tests signed bump energy, steep/waiting segments,
both epsilon atoms and the complete positive Gamma exterior. It separately
checks the nonzero heat term and compares the first-Z energy function to an
independent derivative. Fixtures use ordinary nonzero parameters and do not
substitute those parameters into the actual original construction.
The checked restore/CLI and exact staged source pins are part of this handoff.

## Exact small angular root and source transfer

Let r=1-mu, prate=1+2mu, p=exp(-2r), q=exp(-2prate),
kappa=D/(2B), scale=mu^(30r). Original A/B/D/E/F are integral nodes
in the same normalized beta measure. With the exact Gamma future defects,

```
r_scaled = Xf(0)-Xf(Z) + exp(logtheta-logRtail-logscale)*Theta_hat
s_scaled = exp(logpressure-logRtail-logscale)*Pressure_hat
b1=exp(r)*r_scaled/A
b2=exp(-3prate)*s_scaled/B
p*x+y=b1
x+q*y+kappa*scale*(x^2+q*y^2)=b2
```

Eliminate y=b1-px. With

```
Aq=kappa*scale*(1+q*p^2)
Bq=1-p*q-2*kappa*scale*q*p*b1
Cq=q*b1+kappa*scale*q*b1^2-b2
x=-2*Cq/(Bq+sqrt(Bq^2-4*Aq*Cq))
d1=scale*x, d2=scale*y
```

the exact rationalization avoids division by microscopic Aq. The source
ball certificate proves Bq/discriminant/derivative positivity and that both
roots lie in the already admitted unique small branch.

The source transfer is explicit. Native Xp is Z-independent. Actual pulse
memory is bounded in logarithms and Xv is inside the archived symmetric
candidate interval. With Q=1+Z^2 and k=1-sigma(t/100),

```
Xf(0)-Xf(Z)=2*Xv*exp(-100r)*Z^2/Q
 + integral_0^100 exp(-r(100-t))*2^k*(1-Q^-k) dt
1-Q^-k=Z^2*integral_0^1 k*(1+s*Z^2)^(-1-k) ds
```

This proves the new exact value and first-Z RHSs are the functions enclosed
by the original divided flatten source. Production heat conversion
assignments, the complete Gamma integrands, actual radius/waiting recipe,
strict Md40 delta minimum branch and exact mu recipe are bound separately.
Interval endpoint tuples from different contexts are not used as an identity.

The future contains flatten + power + signed angular change + complete
steep/waiting/preheat tail, minus exactly
`Nrel*Ntail/(1-epsilon)^2*(delta/2)*S*Energy_hat` once. Weight the entire
total by `mu*exp(-26)/2` for the pulse selector. C5 source assignments and
independent physical integrals check these normalizations.

For the fresh pulse proof, the interval inverse reuses the accepted
`linear_inverse(A,D,det,r1,r2)`; the second numerator is
`r2*A0-r1*D0`. Actual native incoming and complete future bounds, original
basis/saddle integrals and proven factor caps feed the same quadratic as
the function graph. A nu enclosure may have lower endpoint zero: it is an
upper bound on an extremely small true value. Exact
`nu_j=mu*Kj*exp(2log_end_scale)>0` remains in the graph. The positive cap
itself is not a lower bound for nu and must not be selected as its value.

## Next executable tasks

- [x] Build the actual small angular root as C1 functions and transfer the
  original whole-Z RHS/ball bounds by explicit function identities.
- [x] Compose the actual complete future energy and weighted target,
  including the entire Gamma tail and the single Nrel heat prefactor.
- [x] Select actual original ap/C1/C2 functions from the genuine native
  inlet and complete future. Freshly prove the positive amplitude bracket,
  derivative denominator and signed end coefficients over [-1,1].
- [x] Pass direct normalized-bump, implicit-derivative, physical-energy,
  source-containment and checked-restoration checks.
- [ ] `N1-PULSE-TRUE-C3-FUNCTIONS`: extend the same defining Gamma/future/
  coefficient bodies to ordinary Z2/Z3 functions where pressure recovery
  requires them. Reuse existing C5 enclosures as bounds of those same
  functions. Do not pad missing derivatives with zero or differentiate caps.
- [ ] `N1-PULSE-OWN-FIELDS`: use
  `Utheta=Ep*exp(-(.5+mu)*y)`,
  `Uz=Utheta*ap(Z)*gp(mu*y)` on the main chart, the two original normalized
  end bumps with nonzero cj, and actual zero-axial gaps. Carry native
  incoming histories through the main, gap and end charts. Use the correct
  incoming common normalization; no additional Pstar division.
- [ ] Build each own partial mean by its defining forward integral or
  selected backward equivalent, retaining nonzero past history. Recover
  q from this own mean with the original delta law. Distinguish zero local
  axial velocity from nonzero cumulative mean/radial velocity.
- [ ] Integrate both full pulse pressure C1 functions with correct radius
  measures and signed finite viscous endpoints. Reuse the selected energy
  identity near the terminal pulse to avoid subtracting unrelated enclosures.
- [ ] Bind the true core nonlinear fixed-point functions and enclosure
  oracle, then the lower collar and inner cutoff. Retain all twenty actual
  nonlinear terms and the anchored gauge primitive.
- [ ] Compose all five complete signed n1 defects once, select fresh linear
  alpha/beta functions, recover repaired own M/V/P and verify terminal zeros.
- [ ] Verify heat/pressure matching and global admissible stress/remainder.
- [ ] Implement actual n>=2 temporal coefficient recovery and smooth summation.
- [ ] Complete oscillatory quadratic cancellation and independent corrected
  Cartesian residual validation with fixed/restricted forcing.

Full pulse pressure, repaired five terminal identities, stress, true temporal
recursion and the final corrected Navier-Stokes field remain open.
