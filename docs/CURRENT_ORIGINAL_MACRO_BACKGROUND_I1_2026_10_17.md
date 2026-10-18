# Signed native macro background I1 and macro/first-switch RQy source

Tasks 14/15 and the long-term goal remain **ACTIVE / INCOMPLETE**. The new
numeric layer encloses the actual current C6 native **pre-finite-N macro
background**, uniformly for m in [0,1] and Z in [37/100,93/250]. It does
not promote these values to full finite-N Q4 I1 values. A separate new
module supplies the conditional full finite-N source RQy seam at R=100.

## Actual native macro arithmetic

Old macro numerical samples cover only Z=-.3,0,.3 (and a distinct old
report has .17). Their modal provider is not identified with the current
C6 source. They are not numerical inputs to this adapter.

The comparison phi/U freeze at the second-micro endpoint. Each of its six
macro histories is its nonreset incoming value times exp(-rate*ell*m),
plus the equilibrium times 1-exp(-rate*ell*m). Since the decay is independent
of Z, the new ordinary derivative intervals are convex hulls of the
incoming and equilibrium rows. Comparison rows through Z4 supply the true
Z3 angular and axial drive bounds. Missing top derivatives are not set to
zero. The actual native fields remain separate from comparison fields.

The adapter uses the current C6 angular equation and complete hydro,
pressure and swirl axial drive. F0_Z/F0=-Lambda*gaxis is retained, and
the external F0^2 cancels its matching internal quotient before boxing the
dressed swirl. Pstar^2 occurs in its prescribed pressure sector; no extra
Pstar factor is put in the I1 measure. The actual macro angular and axial
integral errors are bounded by new directed logarithmic arithmetic, using
the saved loghb. Native M follows M_y=U-M with the original incoming M,
so its differentiated rows are positive-kernel hulls, not reset histories.

The same background U_b and M_b define

```
Q_b=(2*Z*U_b-(1-delta)*Z*M_b-(1-Z^2)*(M_b)_Z)/(1-delta*Z^2).
```

Native U/M rows through Z3 produce Q_b through Z2 and all six ordinary C1
densities. The derivative in this formula is ordinary Z, not physical y.

Write h_b=exp(loghb), Ra=4*exp(-(4*logPstar+1000)),
Y=log(100)-log(Ra), ell=Y-2*h_b, and R=Ra*exp(2*h_b+ell*m).
The saved field named `hb` in scale-log packets is loghb. The physical
positive weight is R^2*ell*dm, exactly once. With Wmacro=10000,

```
K=integral_0^1 R^2*ell*dm/Wmacro=(1-exp(-2*ell))/2.
```

The table gives outward-rounded background I1/Wmacro enclosures. Exact
binary interval tuples in the report are authoritative; these are source
function enclosures, not selected cap values, midpoints or Z samples.

| pre-finite-N macro background integral / Wmacro | ordinary Z0 | ordinary Z1 |
|---|---:|---:|
| q | [-1.455368, -1.443495] | [2.952000, 2.984001] |
| q_Z | [2.952000, 2.984001] | [7.999999, 8.000001] |
| q_squared | [4.167361, 4.236193] | [-17.371273, -17.044800] |
| U_q_Z | [4.368960, 4.440193] | [23.648000, 23.840001] |
| U_squared | [1.095200, 1.107073] | [5.920000, 5.952001] |
| U_U_Z | [2.960000, 2.976001] | [7.999999, 8.000001] |

The signed reduced volume is D_b/Wmacro in [5.823035, 5.976471], with
ordinary Z derivative in [37.371251, 37.996241]. It retains the same six
pressure coefficients, including -I[q_squared]/2. These are background
volume terms. Full pressure also requires finite-N differences and the
retained endpoint sources.

## Next endpoint source

The full Q4 finite-N source macro-upper RQy ports 141114/141115 are paired
with new first-switch lower s=0 ports 184687/184688. The new suffix has
16 nodes, extending the complete
184687-node accepted namespace to
184703 nodes. Its factorized
arithmetic uses the accepted common physical Q_y/radius, phase/N and
nonreset-history source premises. Native macro J=66987 and first-switch
J=9 remain distinct; physical y differentiation divides by each own J.

The identity concerns -RQy(macro,m=1)+RQy(first_switch,s=0), only in a
joined adjacent FTC expression. Both standalone endpoints and their old
caps remain, including first-switch upper ports 148514..148517. No first-
switch volume has been added to the forward cumulative pressure here.

## Verification and remaining work

New production exits 0 in 10.884 seconds. Saved-only
inspection, independent bounded Luna/max mathematical review, exact source
pins and graph suffix reconstruction pass. Separate runtime worker model
metadata is unavailable. This turn replays no accepted producer, owner,
checker, norm, recurrence or quadrature. The prior turn's excluded legacy
review-check exception is preserved under inherited execution audit.
Two unaccepted report attempts exposed an archive-key lookup and raw
interval metadata serialization issue; both were repaired before any new
report was written. These interface repairs did not change the arithmetic.

Report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_macro_background_I1.json` is 19,227,724 bytes, SHA256 `5a77d3d8b9b294516cc0be05ed83667e2710d64ff8427c7e917b8aed3379aede`, with
2008 input and 383
implementation pins. Old 212 forward terms, six signed inner/lower sources,
micro finite-N I1 enclosures and all previous error packets are preserved.
The new macro numeric rows are stored separately; the accepted cumulative
ledger is unchanged. P1 and m15 keep their separate source data.

NEXT obtain source-related current-N macro physical delta-U/delta-M
ordinary rows through Z3, then delta-Q through Z2, and add the resulting
six I1 defect radii once to these background integrals. The older whole-Z
R100 history report uses a different N selection and is not a physical
delta-U/delta-Q provider for this source without an all-N function bridge.
Then extend the first-switch pressure volume while using the new endpoint
identity before separate boxing, and continue the reserved-pulse terminal
functional comparison.

Full finite-N macro I1, complete pressure, functional five-moment closure,
exact heat/energy, admissible stress/flat remainder, genuine temporal
recursion/summation, realizable corrections and independent corrected
Cartesian NS validation remain incomplete. Goal ACTIVE.

Implementation: [macro numeric adapter](../src/openai_ns_reconstruction/original_macro_background_I1_numeric.py),
[macro/switch source seam](../src/openai_ns_reconstruction/original_macro_switch_RQy_seam.py),
[runner](../scripts/run_original_macro_background_I1.py).

```
python scripts/run_original_macro_background_I1.py --inspect-report
```
