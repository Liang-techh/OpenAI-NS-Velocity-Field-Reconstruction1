# Continuous bump integral enclosures

`ContinuousBasisEnclosure` encloses the exact real continuous bump integrals
for declared decimal parameters. It uses an independent interval context,
outward arithmetic, and an analytic midpoint quadrature remainder. It does
not infer an error bound from agreement between two numerical quadratures.

For \(\phi(x)=\exp[-1/(1-x^2)]\), extended by zero at the endpoints, put
\(a=(1-x^2)^{-1}\ge1\). Direct differentiation gives

\[
\phi'=-2xa^2e^{-a},\qquad
\phi''=(4x^2a^4-8x^2a^3-2a^2)e^{-a}.
\]

Since \(a^me^{-a}\) has maximum \(m^me^{-m}\) for \(a\ge1\), valid global
bounds are

\[
D_0=e^{-1},\quad D_1=8e^{-2},\quad
D_2=1024e^{-4}+216e^{-3}+8e^{-2}.
\]

For \(g=e^{kx}\phi\),

\[
\sup|g''|\le e^{|k|}(D_2+2|k|D_1+k^2D_0).
\]

The corresponding bounds for \(\phi^2\) are \(D_{0,2}=D_0^2\),
\(D_{1,2}=2D_0D_1\), and \(D_{2,2}=2D_1^2+2D_0D_2\).
All endpoint derivatives tend to zero, so the midpoint error theorem applies
on the closed interval. For \(n\) equal panels on \([-1,1]\),

\[
|I-Q_n|\le \frac{\sup|g''|}{3n^2}.
\]

The mesh midpoint \(-1+(2j+1)/n\) is constructed with integer and rational
interval operations. No ordinary floating-point mesh is converted after
the fact. The interval midpoint sum includes arithmetic rounding; adding
the symmetric analytic remainder includes quadrature error. Ratios enclose
the normalized weighted bump moments and the squared-bump energy Gram atom.

The determinant is bounded through its correlated analytic identity,

\[
\det M=-2e^{-2+6\mu}\sinh(\mu)B_1B_2.
\]

For positive \(\mu\), \(\mu\le\sinh\mu\le\mu e^\mu\). This avoids
subtracting nearly equal exponentials or independently enclosed matrix
products. The shared candidate's determinant interval is strictly negative.

At 4096 panels the normalizer interval width is approximately
`1.216e-6`, the first weighted moment width `5.706e-6`, and determinant
width `1.477e-33`. These are conservative certified intervals, not a claim
of 60-digit integral accuracy. Increasing panel count from 1024 to 4096
reduces the analytic midpoint remainder by 16.

Every serialized interval retains exact binary endpoint sign, mantissa and
exponent. Decimal values are display approximations and are not used as the
certificate endpoints. Parameter scope is the declared real `mu` and `ell`;
inherited source parameter uncertainty is excluded.

The live axial runtime exposes the optional `basis_enclosure(...)` method.
Its CLI tests containment of the actual live basis. Existing flags for the
complete axial quadrature and global energy remain false: pulse quadrature,
incoming moments, coefficient sensitivity, the full target energy and exact
terminal closure are still open.

Run:

```powershell
python experiments/root_st073/lei_ren_part1_paper_continuous_basis_enclosure.py
python experiments/root_st073/lei_ren_part1_paper_continuous_axial_runtime.py
```
