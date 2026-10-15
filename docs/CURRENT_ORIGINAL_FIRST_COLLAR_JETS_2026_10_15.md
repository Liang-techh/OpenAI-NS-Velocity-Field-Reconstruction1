# Original first-collar axial jets, own-mean Q and pressure — 2026-10-15

## Delivered and API

`OriginalFirstCollarJets` extends the evaluated actual leading first collar
through native fixed-phase axial Taylor order 5. Its six actual histories
use the same nonzero core inlet, and Q is recovered through order 4 from
that own mean. Pressure uses the compatible original P0 and the same C
primitive with the original amplitude applied once.

```python
from openai_ns_reconstruction.original_first_collar_jets import OriginalFirstCollarJets
field = OriginalFirstCollarJets(precision=90)
packet = field.evaluate('.3', '.5', degree=24, order=5)
packet['normalized_fields_and_histories']['Mpsi']
packet['radial_recovery']['Qpsi_increment_per_hb']
packet['pressure']
```

Inputs are exact decimal strings, |Z|<1 and 0<=s<=1. Requested output
order is 1..5. Inlet order 6 is used internally when output order 5 needs a
Z derivative of comparison moments. Coefficient k means D_Z^k/k!, at fixed
original phase s. Each serialized coefficient has exact directed source
endpoints and an absolute error; multiply both by k! for ordinary derivatives.

This is a finite-original-width numerical approximation with errors, not
an exact point oracle or merely a formal zero-width derivative. Point
coefficients use the first-width approximation. hb, epsilon corrections
and original amplitude factors remain separate throughout recovery.

## Differentiated core inlet and model tail

Fresh original coupled radial recurrence rows provide Phi, normalized
Psi=(Uz-U0)/epsilon, six normalized histories and radial derivatives 0..2.
Integrate the normalized rows directly; never recover them by subtracting
rounded raw fields. Nonlinear differentiated tails still use the accepted
Xh correction and `tail_factor(k)/k!`.

The Bessel-model tail is bounded directly in axial Taylor coefficients.
Write chi=chi0+v, |chi0|<=1, with v0=0, and let b[k,m] be coefficient k of
(sum_{j>=1}|chi_j| z^j)^m. For n=N+1, radial order i and q=rho/2, each
positive tail starts at

`b[k,m]*binomial(n,m)*rho^(n-i)/(2^n*(n-i)!*(n+1)!)`

and has ratio at most

`q*(N+2)/((N+2-i)*(N+3)*(N+2-m)) < 1`.

Sum its geometric majorant over m<=k. This retains the small actual chi
axial derivatives, instead of multiplying every model tail by a large
inverse jet weight. The real chi0 interval is intersected with its exact
[0,1] theorem range; derivative rows are unchanged.

## Finite-width errors and exact source cancellation

The original signed first-collar source supplies uniform prefix bounds.
A weighted jet norm is sum |J_k|*tau^k, with the accepted original
comparison weight tau. The actual hb log must be used before exponentiating
or capping large norms; hcap*norm alone is not adequate for high derivatives.

Higher derivatives of F0-relative jets can contain enormous factors that
cancel analytically. The refinement cancels that gauge before arithmetic.
For g=L*Haxis/(Haxis^2+sigma^2), Gamma=F0_Z/F0=-g/epsilon,

`Dbar = rho/L * [epsilon*(-W+angular_no_Gamma/(2Phi))`

`                 +g*((1-delta)*Z*H/2+(1-Z^2)*K)/(2Phi)]`.

Here H,K are the actual core-inlet normalized comparison moments, while
Haxis appears only in g. Comparison variation from that same inlet is
bounded by its original WidthPolynomial first/second rows and cubic
remainder. The actual Phi/barPhi ratio is handled separately.

Difference norms propagate through the exact source products, derivatives
and variable inverses. For an input a=a0+hb*da, inversion uses

`E_inverse <= norm(inv(a0))^2*E_a / (1-hb*norm(inv(a0))*E_a)`.

The denominator condition is proved with the original log(hb). Product
errors include hb*E_a*E_b. Derivative norms use the input order divided by
tau, retaining one more input coefficient. rho=4exp(hb*s), its variation
and the squared radius in the swirl drive are all included.

The implicit S=epsilon^2*F0^2 and all its original directed source derivatives
remain in the swirl term. No F0 value is selected or exponentiated. The
exact zero-width core-join coefficients 4Phi_rho*I(s), 4Psi_rho*I(s) are then
paired with the finite-width error derived from this actual comparison/
ratio/radius source. No analytic core continuation substitutes for the
actual collar. In particular I(s)=integral(1-sigma) is not replaced by s.

The six own-history Volterra errors use the same weighted norm and preserve
all normalized U0/epsilon splits. Their finite-width kernel and source
variation bounds are applied before Q recovery.

## Q and pressure

Use d=1-Z^2, L=1-delta Z^2, U0=4Z+j. The original radial recovery is

`Qaxis=((1+delta)*Z*U0-4*d)/L`,

`Qpsi=[2Z*Psi-(1-delta)*Z*Mpsi-d*(Mpsi)_Z]/L`.

The same linear map is applied to the finite-width increments and their
Z derivatives. Thus Q=Qaxis+epsilon*(Qpsi_inlet+hb*Qpsi_increment).
Q is not fitted as a separate velocity component.

Pressure is P=P0+R*F0^2*C, with a plus sign and R=4epsilon exp(hb*s).
The output provides P0/Pstar^2 and C dressed by the relative F0^2 Taylor
jet at the base point. Apply the formal F0(base)^2 once, keep exp(hb*s),
and retain the separately factored hb increment. The compatible complete
pressure source is not changed or repaired by an independent tail.
These are native pressure jets; physical p=lambda^-2 P and physical
Cartesian collar velocity assembly remain open.

## Executed evidence

`python scripts/run_original_first_collar_jets.py --write-report` exits 0
in 100.815 s. The 1703013-byte report pins 560 original input files.
It contains nine full Z5/Q4/pressure5 packets on Z=-0.3,0,0.3 and
s=0,0.5,1. It executes:

- five exact own-mean Q/split derivative identities;
- 108 independent Bessel special-function mixed derivative/tail checks;
- 1161 directed degree18/degree24 comparisons;
- 144 independent accepted scalar-enclosure overlaps;
- 36 original compatible-pressure and F0-squared dress coefficient joins;
- 20 distinct original-width log proofs and four invalid-input rejections.

Fresh saved-report inspection exits 0 in 1.300 s. No broad historical
suite or source-owner archive rebuild is run. At Z=0.3,s=0.5, sampled
relative errors are below 1.134e-7 for the Qpsi increment value and
4.854e-8 for the Phi increment Z5 Taylor coefficient. The latter is
approximately -6.8559114030685322e-48, with error below 3.328e-55.
These sampled accuracy statements do not prove uniform relative errors.

The existing read-only Luna/max worker reviewed the differentiated model
tail, source factors, normalization, norm errors and source-cancellation
route. A tiny normalized Kpsi increment remains unresolved; no exact
cancellation or nonzero sign is inferred from it.

## Next work and limits

Task 15 remains IN PROGRESS; full goal ACTIVE / INCOMPLETE. These jets
supply the evaluated actual leading fields/histories required for n1
transport. They do not themselves solve a positive temporal order.

1. Build actual first-collar radial rows from the original chi and core
   comparison equations, preserving hb^-1 sectors in chi_rho and second
   radial derivatives. Do not discard a thin-layer term because hb is small.
2. Feed these leading rows and own Q/pressure to the full six-row B0/B1/g
   n1 equation. Retain original core endpoint values and Z jets and the
   complete accumulated history. In dx=hb*sqrt(rho)*ds/2, cancel width
   powers algebraically against hb^-1 forcing before numerical bounds.
3. The currently accepted n1 first-inner operator domain ends at s=1/2.
   First transport on that domain, then admit/extend the remaining collar;
   leading-field coverage through s=1 does not enlarge that operator domain.
4. Establish quantitative first/second derivative joins, actual second/
   lower collar transport, source-owned cutoff and reserved moment repair.
   Extend shared-streamfunction physical assembly with the evaluated data.
5. Common global N, global five-moment numerical closure, heat/energy,
   stress/flat remainder, genuine n>=2 temporal recursion and oscillatory
   corrected Cartesian NS remain open. Axial/radial derivative orders and
   collar width powers are not temporal coefficient orders.
