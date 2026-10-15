# Original finite-width n1 first-collar transport — 2026-10-15

Subsequent progress: [full original first-half collar Z2, q_Z and q_rho_Z](CURRENT_ORIGINAL_FIRST_COLLAR_AXIAL2_2026_10_15.md)
are now evaluated through a separate native API on 0<=s<=1/2. Physical
mixed-derivative assembly and complete Cartesian C2 joins remain open.

Subsequent progress: [actual native radial derivatives and core source seam](CURRENT_ORIGINAL_FIRST_COLLAR_DERIVATIVES_2026_10_15.md)
are now evaluated on the same admitted first-half collar. Complete physical
Cartesian C2 joins and higher n1 axial numerical data remain open.

Subsequent progress: [physical leading+n1 velocity is now assembled](CURRENT_ORIGINAL_FIRST_COLLAR_VELOCITY_2026_10_15.md)
on this same admitted first-half collar, using native phase and separate
original width/geometry factors.

## Evaluated interface

`OriginalFirstCollarTransport.evaluate(Z,s)` transports the complete-error
original regular n1 core endpoint into the actual first collar, on its
currently admitted domain `|Z|<1, 0<=s<=1/2`. It returns numerical values
and first fixed-phase Z derivatives of all six natural state components,
their directed finite-width errors, the recovered own mean and Q value,
and separately factored collar increments.

```python
from openai_ns_reconstruction.original_first_collar_transport import OriginalFirstCollarTransport
field = OriginalFirstCollarTransport(precision=85)
packet = field.evaluate('.3', '.5')
packet['normalized_original_n1']['p']
packet['factored_first_collar_increments']['ux']
packet['finite_width_remainder_budgets']
```

The natural frame is

`f=epsilon*F1/F0, u=Uz1/epsilon, k=K1/epsilon, p=P1/(epsilon*Pstar^2)`,

with `fx=d_x f`, `ux=d_x u`. Here `x=sqrt(rho)` is the computational
radial variable, not physical Cartesian x. Physical radius is
`R=epsilon*rho`, and `rho=4exp(hb*s)`. The leading normalized axial
correction is `Psi=(Uz_leading-(4Z+j))/epsilon`; it is not temporal n1.

Coefficient k in an axial Taylor jet means `D_Z^k/k!` at fixed phase.
For k=0,1 this equals the value and ordinary first derivative. Q is
recovered as a value only, because its next derivative would require
second Z derivatives of the incoming full n1 numerical solution.

## Exact source and the finite collar response

The implementation builds actual radial rows from the unchanged first
collar identities and same-family original leading inlet. In particular,

`chi=1-(1-hb)*sigma(s)`,

`chi_rho=-(1-hb)*sigma_prime(s)/(hb*rho)`.

The complete `(1-hb)` factor stays inside the inverse-width coefficient;
its difference from 1 is explicitly bounded in the centered remainder.
The second radial rows of Phi, Uz and Q retain their `hb^-1` sectors.
The exact normalized forcing includes the positive pressure head

`Np=-flow/(2Pstar^2)+(2Q_rho+rho*Q_rhorho)/(epsilon*Pstar^2)`

and the negative axial pressure correction

`-4epsilon*Pstar^2*Z*rho*Np/L`, where `L=1-delta*Z^2`.

In `dx=hb*sqrt(rho)*ds/2`, the width cancels the singular source before
numerical evaluation. Freezing the remaining coefficient at the original
core inlet rho=4 gives the finite principal increments

| Component | Principal increment |
| --- | --- |
| f, u, k | 0 |
| p | `-8 Z (Psi_rho/Pstar^2) sigma(s)/L` |
| fx | `epsilon^2 * 32 Z^2 Phi_rho sigma(s)/L^2` |
| ux | `epsilon Pstar^2 * 64 Z^2 (Psi_rho/Pstar^2) sigma(s)/L^2` |

All inlet derivatives in this table belong to the full leading core.
The factor 64 includes both the axial and pressure contributions.
Each exported increment remains separate from the nonzero incoming
core value, so small epsilon-scaled changes survive finite precision.
Zero principal increments for f,u,k do not mean their finite-width
increments are set to zero: their local integral is enclosed in the
remainder around the zero numerical proxy.

## Complete finite-width errors

Write the local exact equation as `W_x=A(W,W_Z)+g0+hb^-1*gminus1`.
Only the local collar integral is used; the original core 0-to-4
contribution is already in `OriginalN1CoreVolterraNumeric.evaluate(Z,'4')`.
It is never integrated a second time.

The accepted common-inner regular all-word solution bound controls
`W_i/x^p_i`, with `p=(2,2,2,2,1,1)`. Since x<3, `9M` bounds each raw
natural component on the required tube. Cauchy bounds at half the
accepted solution gap provide W, W_Z and W_ZZ. Taylor coefficient k is
bounded by `9M/gap^k`; differentiation supplies its factorial normally.
The exact conjugated `linear_rhs` applies the original F0 gauge once.
It retains row-specific epsilon/Pstar support before division by each
increment scale. The accepted all-word bound is not a radial Taylor tail.

For each increment scale `a_i` (1, epsilon^2 or epsilon*Pstar^2), the
remaining regular integral is bounded by

`hb * (3s/2) * (norm(AW_i)+norm(g0_i))/a_i`.

The measure bound x/2<3/2 follows from the archived actual rho<4.1
geometry. The matrix and all finite forcing terms are bounded uniformly
over the phase prefix. The actual leading source uses uniform first-chart
direction boxes and nonnegative, monotone control masses; separate
absolute source norms therefore cover every earlier prefix.

The singular coefficient is centered algebraically before its error is
bounded. A width-change calculation `a(h)=a(0)+h*error` retains radius,
sqrt-radius, actual/comparison ratio, comparison radial derivative and
the `(1-hb)` change. Products include the width-weighted product of
errors; inverses use the proved geometric denominator bound. Its local
integral error is `hb*sigma(s)*E_i`, since sigma is increasing on this
domain. Comparison derivative variations use the pinned joint core
Cauchy domain. Direct integration against the full W norm covers
propagator effects; no frozen Green kernel is silently substituted.

Every large error coefficient is multiplied by the original log(hb)
before taking a representable upper cap. hb and F0 are not materialized
or replaced by a selected bound. Original nonzero S=epsilon^2*F0^2,
pressure data, all linear feedback and source-family bindings remain
in the directed error calculation.

The recovered mean is `mean_u=u+k`, and the n1 recovery is

`Q1_normalized=[2Z*u-(1+delta)*Z*mean_u-(1-Z^2)*D_Z mean_u]/L`.

This uses the n1-dependent weight, rather than the leading mean weight.

## Execution evidence

The accepted focused run exits 0 in **151.777 seconds**. It computes
nine packets at Z=-0.3,0,0.3 and s=0,0.25,0.5, with:

- 12 exact independent raw differential-source and inverse-width measure identities;
- 36 exactly zero phase-zero increment coefficients;
- 78 comparisons across precision 60 and 85, using the complete errors;
- four rejected unsupported requests and 20 distinct original-width log proofs.

The first comparison run exposed cancellation in the comparison calculation
itself. Differences of serialized decimals are now computed at at least
160 digits, above both output precisions; the reported source-error
thresholds are unchanged. The final run passes in full.

At Z=0.3,s=0.5 the normalized pressure increment is approximately
**-2.68132896789416667552**, with absolute error below **2.062e-20**.
The fx increment divided by epsilon^2 is approximately
-0.17263010957436385302 (error below 4.842e-44); the ux increment divided
by epsilon*Pstar^2 is approximately 6.43518952294600002126 (error below
4.948e-20). These are sampled accuracy statements, not uniform relative
bounds. They must not be mistaken for physical pressure or unscaled
velocity increments.

The report is 764025 bytes and pins 1477 original inputs plus the
implementation files. A fresh-process saved-report inspection exits 0.
The existing read-only Luna/max worker found no material source-algebra,
pressure-sign, scale or regular-operator error; the same-domain norm and
uniform-prefix assumptions are stated above.

Driver: `scripts/run_original_first_collar_transport.py`.
Report: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_first_collar_transport.json`.

## Remaining construction

Task 15 remains IN PROGRESS; full goal ACTIVE / INCOMPLETE.

1. Admit and evaluate the remainder of the first collar and the second/
   lower collar. Leading coverage to s=1 does not enlarge this n1 domain.
2. Establish quantitative first/second derivative joins, then apply the
   original cutoff and reserved own five-moment repair to evaluated data.
3. Extend shared-streamfunction physical leading+n1 assembly beyond the
   core, including actual histories and appropriate chart inverses.
4. Complete global five-moment numerical closure, heat exterior and
   finite-energy diagnostics, admissible stress and flat remainder.
5. Construct genuine temporal n>=2 coefficients, their independent moment
   repair and final smooth summation, then oscillatory corrections and
   independent corrected Cartesian NS residual validation.

This delivery evaluates one positive temporal order through a restricted
actual collar. It does not establish global matching, a physical collar
velocity API, genuine n>=2 recursion, or a blow-up theorem.
