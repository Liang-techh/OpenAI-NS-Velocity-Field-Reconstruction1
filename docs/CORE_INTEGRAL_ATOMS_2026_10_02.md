# Actual core inlet atoms and axial derivatives

Six normalized core atoms H,M,K,A,B,C are now integrated coefficientwise
from freshly rebuilt coupled Phi/Uz rows, through axial order6. This
replaces the missing numerical core integrals with controlled source
enclosures. It does not yet replace the old bridge providers or resolve
their signed source integrals.

## Source definitions and units

All radial integrals below are over rho in[0,4], with rho=Lambda*R.

```
H=(1/8)*integral rho*Phi
M=(1/4)*integral Uz
K=(1/8)*integral rho*Phi*Uz
A=(1/4)*integral Uz^2
B=(1/16)*integral rho*Phi^2
C=(1/4)*integral Phi^2
```

These six atoms represent five physical moment conditions because the
axial quadratic primitive has distinct A and B contributions. Raw bridge
V means Uz, not the pressure primitive V_pressure=int Phi^2. At the exit,
V_pressure=4C and P-P0=epsilon*F0^2*V_pressure=R*F0^2*C, with R=4epsilon.
The same mean M is used by radial velocity recovery.

Each coupled radial row already includes the affine U0=4Z+j. Linear
integration retains it in M; product convolutions retain it in K and A.
Axial coefficient k is derivative/k!, so coefficient products have no
binomial factors. Exported ordinary derivative arrays multiply by k! once.

## Finite coefficients and controlled tails

The degree24/depth6 rows are freshly generated from current selected
Cstar, original positive swirl-source enclosures and analytic preheat
pressure. No old finite rows or chosen parameter midpoints are read.
Integrating rho^n uses 4^(n+1)/(n+1), or 4^(n+2)/(n+2) with a rho weight.
Products are integrated over their full finite radial convolution, not
truncated again at degree24.

The infinite Phi tail includes its original factorial Bessel model and
admitted nonlinear Xh correction. The Uz model is affine in rho; only
epsilon times nonlinear correction contributes beyond degree24. At the
exact shared H root, fresh rows preserve H(a)=0 before enclosure. Since
chi has axial valuation2, the model radial tail beyond24 vanishes through
axial order6. Nonlinear tails remain nonzero.

The general Bessel tail now supports axial order6. For real Z, chi0 is
in[0,1]. If M=sum |chi_k| for k1..6, coefficient k of chi^n is bounded by
(k+1)*(n+1)^k*(1+M)^k. The factorial radial terms have a checked contracting
ratio. Xh derivative tails are divided by k! before coefficient convolution.

For finite absolute coefficient bounds P_i,U_i and tail bounds E_i,F_i,
the product tail coefficient k is bounded by

```
sum_i[P_i*F_(k-i)+E_i*U_(k-i)+E_i*F_(k-i)].
```

The normalized radial mass is one for H,M,K,A,C and one-half for B.
Independent positive Phi bounds sharpen H/B/C values only. The result is
a directed source enclosure; no cover endpoint is promoted to a source
value. Independent absolute tail bounds are conservative and may lose
some correlation.

At Z=.5, representative displayed values are

| Atom | Value | Directed interval width, approximately |
| --- | ---: | ---: |
| H | 0.4795280821510107 | 1.08e-44 |
| M | 2.0000000000000000 | 3.77e-183 |
| K | 0.9590561643020214 | 2.15e-44 |
| A | 4.0000000000000000 | 1.51e-182 |
| B | 0.1272028385626529 | 2.58e-44 |
| C | 0.4006067501478635 | 5.15e-44 |

Displayed decimals omit smaller source corrections, including j in M.
Use interval artifacts for calculations, not these rounded examples.

## API and evidence

```
from lei_ren_part1_paper_compliant_core_integral_atoms import CompliantCoreIntegralAtoms
provider=CompliantCoreIntegralAtoms()
atoms=provider.core_atoms('.5',z_order=6,degree=24)
root_atoms=provider.root_atoms(z_order=6,degree=24)
```

The API accepts current source intervals Z in[-1,1]. Sample receipts cover
Z=.3,.5,-.5,0 and the exact shared root. Focused reproduction:

```
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage coreatoms
```

The ordered pipeline has119 modules. The checker passes210 source atom
coefficients,35 shared mean/pressure primitive coefficients,42 independent
exact weighted polynomial integrals,42 independently integrated signed
nonlinear product tails and14 independent differentiated factorial-tail
sums through order6. The independent finite fixture is not the paper source
or an NS test. Same-source file hashes, family and datum remain bound.

## Next dependency and limits

Inject these actual core atoms into the original comparison namespace;
transport all six histories through both smoothing intervals. Then solve
the signed actual bridge and first-switch integrations without resetting
their source histories. Keep actual and comparison moments distinct and
preserve V100/V110/E identities. Existing correlated covers are not those
point values; caps must only enclose final errors. Actual implicit bump
values and full physical point/chart assembly follow this step.

Global radial morphology, measured dynamics/winding, required-domain energy,
admissible stress/flat remainder, genuine n-dependent coefficient recursion,
mean/oscillatory correction and full forced-NS residual remain unfinished.
Original unlocalized whole-space energy remains infinite.
