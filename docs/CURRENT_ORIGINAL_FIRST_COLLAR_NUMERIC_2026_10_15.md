Subsequent progress: [actual first-collar Z5 jets, own-mean Q4 and pressure5](CURRENT_ORIGINAL_FIRST_COLLAR_JETS_2026_10_15.md) now evaluate the differentiated data. This earlier value-only report retains its original scope.

# Original leading first-collar numerical fields and own histories — 2026-10-15

## Delivered and usage

`OriginalFirstCollarNumeric` evaluates the original prescribed-shear leading
first micro-flat collar, for an interior exact-decimal Z and phase 0<=s<=1.
The output contains Phi, raw Uz and all six actual own cumulative histories.
It stores inlet values and nonzero finite-collar increments separately:

```python
from openai_ns_reconstruction.original_first_collar_numeric import OriginalFirstCollarNumeric
field = OriginalFirstCollarNumeric(precision=80)
packet = field.evaluate('.3', '.5')
phi = packet['normalized_fields_and_own_histories']['Phi']
phi['inlet']
phi['increment_per_original_hb']
raw_six = packet['original_fields_and_six_histories']
```

Each raw term means

`coefficient * epsilon**epsilon_power * hb**hb_power`.

Every coefficient includes a numerical approximation, serialized-decimal
absolute error and exact directed source endpoints. The increment coefficient
encloses `(actual field at s - actual inlet)/hb` at the **finite original
positive width**; it is not merely a formal zero-width derivative. The nominal
coefficient is the first-width approximation, and its error includes the full
source and finite-width difference. No interval midpoint selects a field.

The source radius is `R=4*epsilon*exp(hb*s)`, rho=4*exp(hb*s), y=hb*s.
The ledger defines hb=cstar*K^-100 and supplies its exact log enclosure;
log(hb) is approximately -1.4305906698640853e+408906090034569676.
Neither hb nor F0 is exponentiated. Pstar^2 is a representable arbitrary-
precision intermediate in the directed pressure branch. This interface does
not extend the existing physical Cartesian velocity API beyond the core yet.

## Actual prescribed shear and coefficient errors

The defining source remains the unchanged same-family
`CompliantActualBridgeIntegrals.prepare/contributions` first chart. Its
comparison direction is known data; it is not the actual field or its own
histories. All hydro, Pstar^2 pressure and F0^2 swirl terms enter the directed
finite-collar increment. The original nonzero implicit swirl and complete
pressure uncertainty remain enclosed.

The first control is chi=(1-sigma(s))+hb*sigma(s). Put
I(s)=integral_0^s(1-sigma(u))du. The original inlet equations give the nominal
increments `4*Phi_rho(4,Z)*I(s)` and `4*Psi_rho(4,Z)*I(s)`, where
Uz=U0+epsilon*Psi, U0=4Z+j. Point coefficients are computed from fresh original
core recurrence rows. They are numerical proxy approximations, with no exact
point-oracle claim. Their complete finite-collar error comes independently
from the accepted actual signed integral sources, not a value-tail bound
misapplied to a radial derivative.

For every prefix u<=s, both I(u) and integral_0^u sigma are nonnegative and
monotone. Absolute sums of the actual signed source rows therefore bound
all earlier increments, before any signed cancellation. The angular
exponential remainder is bounded separately. Remaining hb powers and
unmaterializable F0 products are controlled with the actual log(hb) and
source log before a positive error cap is used. A generic 1e-50000 width
bound alone would give useless errors after multiplication by huge
pressure-normalized norms; it never substitutes for the original width.

## Own histories and their normalization

At rho=4, integrate the fresh coupled core polynomials and directed analytic
model/correction tails. Retain the small axial correction before integration:

- H=(1/8) integral_0^4 rho*Phi;
- Mpsi=(1/4) integral_0^4 Psi;
- Kpsi=(1/8) integral_0^4 rho*Phi*Psi;
- Apsi=(1/4) integral_0^4 Psi^2;
- B=(1/16) integral_0^4 rho*Phi^2;
- C=(1/4) integral_0^4 Phi^2.

Raw histories are exactly M=U0+epsilon*Mpsi,
K=U0*H+epsilon*Kpsi and
A=U0^2+2epsilon*U0*Mpsi+epsilon^2*Apsi; H,B,C are unchanged.
There is no subtraction of rounded raw M/K/A to recover these corrections.

The normalized y-ODEs have rates (2,1,2,1,2,1), and sources
(2Phi,Psi,2Phi*Psi,Psi^2,Phi^2,Phi^2). Thus each history satisfies

`X(s)=exp(-rate*hb*s)*Xa + hb*integral_0^s exp(-rate*hb*(s-u))*source(u)du`.

The reported first-width approximation is s*(source_at_inlet-rate*Xa).
The finite-width error uses the actual prefix bounds and
`|(1-exp(-rate*hb*s))/(rate*hb)-s|<=rate*hb*s^2/2`.
Histories are carried from the actual nonzero core integrals, with no seam
reset and no substitution of comparison moments. C remains the actual
pressure primitive; full physical pressure assembly is still open.

## Executed evidence

`python scripts/run_original_first_collar_numeric.py --write-report` exits 0
in 57.307 s. The saved 1296144-byte report pins 559 original input files.
It contains 16 samples on Z=-0.3,0,0.3,0.17 and s=0,0.25,0.5,1. The exact
inlet, full first exit and an offgrid original inlet are included.

Six exact raw-history recombination/ODE identities, 24 independent Gaussian
inlet moment integrals, 256 directed degree16/degree22 comparisons, 768
original log-product bounds and five invalid-input rejections pass.
Maximum nonaxis Phi/Psi increment relative error in the checked samples is
below 4.854e-8. This is a sampled accuracy statement, not a uniform Z bound.
Fresh saved-report inspection exits 0 in 1.594 s without restoring owners.
The existing read-only Luna/max worker independently reviewed the source
factors, first-width identities, normalized history ODEs and log bounds.

At Z=0.3,s=0.5, `(Phi-Phi_a)/hb` is approximately
-0.229869518383453019, with absolute error below 1.116e-8.
I(0.5) is approximately 0.479366124612204887. At s=1, I(1)=1/2 exactly;
the increment is approximately -0.239764041075505351, with absolute error
below 6.880e-39. This control-dependent change differs from simply
continuing the analytic core with I(s)=s.

A very small normalized Kpsi increment is not resolved against its source
uncertainty; its interval is preserved. Do not infer an exact cancellation
or a nonzero sign from its displayed approximation.

## Remaining work and handoff

Task 15 remains IN PROGRESS. This is original **leading** first-collar
value evaluation; it does not transport temporal n1 through the collar.
The full long-term goal stays ACTIVE / INCOMPLETE.

1. Extend the same normalized inlet and history evaluation to required Z
   jets with differentiated analytic tails and actual signed-source rows.
   Preserve all U0/epsilon correlations when differentiating raw histories.
2. Recover actual Q from M and M_Z, and pressure from the same C, using the
   original amplitude dress and compatible P0. Check the core/first inlet
   and first/second value and derivative joins with independent equations.
3. Feed these actual leading jets and own histories to full B0/B1 n1
   transport. Start from the computed core endpoint and retain the whole
   accumulated integral, including dx=hb*sqrt(4exp(hb*s))*ds/2.
4. Continue the actual second/lower collars, source-owned cutoff and moment
   repair; then extend shared-streamfunction physical assembly and chart
   inversion beyond the original core. Do not reset histories at interfaces.
5. Common global N, functional five-moment numerical closure, heat/energy,
   stress/flat remainder, genuine n>=2 recursion and oscillatory corrected
   Cartesian NS remain open. These first-width powers are collar-width
   expansion machinery, not temporal coefficient recursion.
