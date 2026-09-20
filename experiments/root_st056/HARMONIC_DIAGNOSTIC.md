# Higher-order harmonic compatibility audit

This supplements RESULTS.md. It is a post-freeze diagnostic of the same fields, not another fit or altered acceptance norm.

## Identity and combined witness

For any harmonic polynomial H and phi=grad(H), under this experiment's smooth compact u/p and compact solenoidal prescribed force:

```
d_H = Integral grad(H) dot R = -Integral u^T Hess(H) u
```

Time derivative and force vanish by incompressibility/integration by parts; pressure vanishes because Delta(H)=0; viscosity vanishes because Delta(grad(H))=0. These conclusions depend on compact support and the specified solenoidal force; they are not asserted for arbitrary noncompact pressure or forcing.

Use axisymmetric solid harmonics of degrees 2 through 8. The first even members are H2=z^2-s/2 and H4=z^4-3*s*z^2+3*s^2/8, s=r^2. Let G_ij=Integral grad(H_i) dot grad(H_j) over r<2, |z|<2. The exact Gram matrix is obtained by symbolic polynomial integration. The normalized linear combination chosen from d gives

```
true spatial ||R||_L2 >= sqrt(d^T G^{-1} d).
```

This is a projection lower bound, not a redefinition of the residual. Numerical moment integrals do not constitute an interval-certified bound. The normalized degree-8 Gram condition number is about 5.425, and its entries were checked by independent 16-point quadrature. Nested degrees strengthen the same weak-witness estimate; the velocity family was not enlarged.

## Actual initial-time estimates at 96-point quadrature

| Largest test degree | Q2 and J2 | Moment-balanced M1 |
|---|---:|---:|
| 2 | 0.004203256891 | approximately zero |
| 4 | 0.005556315217 | 0.006373575226 |
| 6 | 0.010572731043 | 0.009131573266 |
| 8 | 0.010761766379 | 0.010185209786 |

Thus the simple initial energy-anisotropy balance is not the only compatibility condition. M1 eliminates degree-2 initial defect but leaves substantial higher-degree defects. Its larger full L2 was already reported as a negative result, not a successful repair.

For Q2/J2 at t=.25, degree-8 estimates at quadrature48/72/96 are .010761769417/.010761766355/.010761766379. For M1 they are .010185207776/.010185209818/.010185209786. These are converged numerical estimates, not interval certificates. At t=.75, degree-8 estimates are Q2 .009348608058, J2 .009348179604 and M1 .010976071045. All are above .001 but below the actual full residual norms near .047-.061; this does not explain every part of the plateau.

## Independent cross-check

For each field at t=.25 and .75, construct the normalized witness from the velocity moments, then independently integrate phi dot the original analytic momentum residual. All six witness norms equal one to about 7e-15. The largest disagreement between the velocity-moment and residual-integral routes is 3.98e-8. The raw diagnostics and cross-check records are in the complete user package. This compares two formulas; it is not rigorous error enclosure.

## Consequence for the next construction

The exact Q2 initial velocity fixes all these initial moments. Initial-flat swirl corrections and pressure-only changes cannot remove them. J2 is therefore a useful small-improvement checkpoint, not a complete route to true L2<.001 while that entire initial velocity remains fixed. The numerical size of the obstruction still needs interval enclosure for a certified exclusion statement.

Keep E0=1, support, viscosity, prescribed force and required core geometry, but permit nontrivial initial velocity redistribution. Control several harmonic moments through time and jointly adjust poloidal/swirl dynamics and compatible pressure. A single global amplitude rebalance is insufficient. No full-family impossibility or source-paper theorem is claimed.

Actual final focused tests: 10 passed in 6.12 seconds, including the independent Gram and nested-projection check. Full historical tests, native MATLAB and Lean were not run. No new field was fitted after any of these diagnostics.
