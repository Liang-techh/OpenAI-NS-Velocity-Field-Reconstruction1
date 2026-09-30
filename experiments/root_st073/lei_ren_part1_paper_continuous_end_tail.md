# Continuous end-tail diagnostic

This diagnostic reads `lei_ren_part1_paper_continuous_incoming_outer.json` and
its nested `continuous_solve` receipt. It instantiates the latest materialized
`ContinuousAxialCorrection` directly, without rebuilding the joined outer
field or regenerating incoming rows.

The direct tail is

\[
 T_\lambda(s)=\int_s^\ell e^{\lambda t}\beta(t)\,dt
              =P_{-\lambda}(-s),
\]

using the even bump. Therefore

\[
 T_\lambda'(s)=-e^{\lambda s}\beta(s),
 \qquad P_\lambda(s)+T_\lambda(s)=P_\lambda(\ell).
\]

The flat-endpoint example uses 100-digit arithmetic at `s=.1499` and
`lambda=.5`. Subtracting partial from full rounds to zero, while direct
integration retains a positive tail with log value about `-763.2991417`.
Factoring the endpoint maximum keeps adaptive quadrature from stopping
merely because the unscaled integral is below its absolute error scale.

The run checks those identities for positive tails approaching `ell = 0.15`,
using the reflected primitive derivative as the direct MP derivative identity
and a finite-difference check at resolved endpoint gaps,
then checks the last bump endpoint at offset `-0.85` and points immediately
below and above it. The endpoint sample has an outward `1e-180` MP guard so
the decimal `-0.85` is classified on the post-support branch despite the
finite MP representation of `0.15`. For the cumulative row it uses the retained materialized
terminal balance

\[
 N(s)=\rho-T_\lambda(s),
 \qquad
 \rho=b+a p+\sum_j c_j M_j.
\]

At the support endpoint the direct tail is zero, so the cumulative value is
the actual `rho` supplied by the materialized coefficient and incoming atoms.
The diagnostic reports this residual and never masks it with zero. It also
copies the installed normalized `mean_Z` residual from the terminal receipt
to show that the source derivative residual remains retained.

`quadrature_enclosure_certified`, `incoming_uncertainty_enclosed`, and
`finite_energy_certified` remain false. The direct-tail checks establish
functional identities for the evaluated MP quadratures; they do not provide a
global quadrature enclosure, energy proof, or recursion certification.

Run from the repository root with:

```powershell
python experiments/root_st073/lei_ren_part1_paper_continuous_end_tail.py
```

The generated evidence is in
`lei_ren_part1_paper_continuous_end_tail.json`.
