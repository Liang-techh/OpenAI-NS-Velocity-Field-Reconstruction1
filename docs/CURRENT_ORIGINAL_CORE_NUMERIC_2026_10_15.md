# Original normalized core numerical approximations — 2026-10-15

## Delivered

`OriginalCoreNumeric` computes numerical approximations of the **leading
nonlinear core profiles**, accompanied by directed error bounds relative to
the same original admitted solution. It exposes Phi, Psi, mean_Psi, leading
Uz, radial recovery Q, and compatible P0/Pstar^2. These are numerical values
with explicit uncertainties, rather than source-DAG handles or selected
interval midpoints.

The calculation shares the exact three-part source family with the current
n1 source atlas. The existing n1 physical/Cartesian exports are unchanged.
This milestone does not evaluate F0, assemble n0+n1 or implement higher
coefficients, global matching or complete physical velocity.

```python
from openai_ns_reconstruction.original_core_numeric import OriginalCoreNumeric
solver = OriginalCoreNumeric(precision=70, quadrature_nodes=192)
rows = solver.coefficients('0.3', degree=18)
result = solver.evaluate(rows, '2')
result['profiles']['Phi']['approximation']
result['profiles']['Phi']['absolute_error_upper']
```

Original analytic coordinates are rho=Lambda R, R=r^2/(2lambda^2),
Z=z/lambda^(1-delta), lambda^2(1-Z^2)=tau=1-t. The numerical profile domain
is 0<=rho<=4.1, -1<Z<1; the original closed core occupies rho<=4. Normalized
viscosity and blow-up time remain nu=1, T=1.

## Numerical values and uncertainty

At Z=0.3, rho=2, the refined calculation gives:

| Original profile | Numerical approximation (rounded) | Absolute error of saved decimal |
|---|---:|---:|
| Phi | 0.5767248077568733872 | 3.39e-36 |
| Uz | 1.2000000000000000000001747835... | 1.23e-70 |
| Q | -3.279999999999999999999947565... | 3.37e-70 |
| P0/Pstar^2 | -2.7898516370796503211 | 7.53e-27 |

The error column applies to the **full decimal stored in the JSON**, not to
the shortened display in this table. Phi=F/F0; Q is the auxiliary profile
with Ur=sqrt(R/2) Q, not the physical radial velocity. Uz here is the leading
similarity profile, before the temporal lambda weight. P0/Pstar^2 is the
compatible normalized axis pressure, not the complete pressure field.

Psi can be enormous in its normalized units. Keep Uz=U0+epsilon*Psi and
mean_Uz=U0+epsilon*mean_Psi factored, using the same original epsilon. This
epsilon=1/Lambda is a construction parameter, not a temporal coefficient
weight. The output includes an independently bounded U0 value as well.

These error bounds are **profile approximation errors**, not full momentum
residuals or measured blow-up diagnostics. No small-NS-residual conclusion
follows from their size.

## Actual source inputs and two arithmetic branches

1. Recompute the original active j minimum branch exactly. The fixed integer
   constants CA=344, CQ=1081, CS=57749553 and KN=18000 imply
   `j=1/(16000000*KN*CS*CA)=1/5721363714816000000000`.
   Exact rational comparison verifies containment in the original required_j
   source and dominance over both other minimum branches. No j midpoint or
   interval endpoint is selected.
2. Retain logPstar=exp(40)+11, logLambda=4logPstar+1000, epsilon=exp(-logLambda)
   and the active delta=exp(-4logPstar-30)<1e-200 branch. Arbitrary-exponent
   mpmath arithmetic prevents ordinary floating-point underflow.
3. Compute a **nominal early-pressure atom proxy** from the original cutoff
   exp(-1/y^2). A Chebyshev approximation to its primitive is numerical
   quadrature machinery. The true pressure source is not redefined by it.
   Its error is bounded against the already accepted original atom interval;
   the new quadrature does not independently certify that interval.
4. Execute the existing production `advance_scaled_one` radial recurrence,
   with all twenty coupled source formulas. The point branch uses the early
   atom proxy, an asymptotic turnoff approximation and S=0 as explicit
   **numerical proxies**, not exact choices for the original source.
5. Independently execute the same recurrence in directed interval arithmetic.
   This branch retains the six aggregated beta-2 atoms, seven aggregated
   beta-0 atoms and separate variable-beta flatten atom. It retains the
   actual nonzero swirl source through S=epsilon^2 F0^2 and its Cauchy bounds,
   using the original selected-amplitude bound exp(-6logLambda-2000).
   Individual atom correlations are not reconstructed by the aggregate.
6. Bound the infinite original tail using the admitted analytic fixed-point
   correction norm. Subtract the known model first: Phi has a factorial
   Bessel-model tail, whereas Psi's leading model is linear in rho. This
   avoids losing the large pressure/model correlation in a full-ball tail.
7. Subtract the saved numerical decimal from the directed original-solution
   enclosure. Thus even a proxy outside the finite original interval is
   accompanied by a valid distance bound. No original source is declared
   equal to its proxy, including zero swirl or omitted late pressure masses.

The output stores exact binary endpoint tuples for authoritative intervals,
a raw directed error tuple, and a displayed error padded by 1e-30 relatively.
Rounded interval displays are labeled as rounded. Numerical point values and
their approximation status remain explicit throughout.

## Reproduction and evidence

```powershell
python scripts/run_original_core_numeric.py --write-report
```

The final focused run exits 0 in 7.937 s. It evaluates Z=-0.3,0,0.3 and
rho=0,1,2,4 at two levels:

- 50 decimal digits, 96 quadrature panels, radial degree 12;
- 70 decimal digits, 192 quadrature panels, radial degree 18.

All 72 cross-level value/error comparisons are consistent. Twelve independent
Bessel special-function evaluations agree with the angular approximation
within its actual-source error budget. The refined maximum relative bounds
are 6.30e-30 for Phi, 2.70e-27 for Psi/mean_Psi/normalized axis pressure and
1.03e-70 for Uz/Q on these samples. Refinement comparisons use 120 digits to
avoid cancellation in near-saturated error bars. The numerical pressure-atom
proxy is 0.479462706996512896080768591422... with an inherited-source error
bound 8.94e-27. This is not a new validated-quadrature claim.

The checked output is
`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_core_numeric_approximation.json`.
It records both levels, original source pins, implementation hashes, exact
interval endpoints, numerical approximations and explicit limitations. The
calculation does not rebuild any archived high-source/Volterra owners.

## Remaining work

Task 15 is IN PROGRESS. A same-source normalized leading-core approximation
now exists, while complete original pointwise physical recovery remains open.

1. Recover the original anchor/primitive G and tagged logarithmic F0/S point
   representation with usable error propagation. Do not materialize infeasible
   double exponentials or promote the point proxy's S=0 to a source choice.
2. Resolve the remaining original pressure atoms/waiting function where finer
   relative source information is needed. The current error branch covers
   them; their individual point values have not been computed.
3. Numerically solve the actual regular n1 equation and all nonzero histories,
   using the same source data. Combine leading and coefficient fields through
   their common streamfunction, with the correct temporal weights.
4. Complete actual higher inner interfaces, physical z inversion and original
   radial chart selection, original-scale error budgets and common N.
5. Continue global moment/energy/heat/stress layers, genuine n>=2 repair and
   summation, oscillatory cancellation, measured contraction/material winding
   and independent fixed-forcing corrected Cartesian NS validation.

The full long-term goal remains ACTIVE / INCOMPLETE.
