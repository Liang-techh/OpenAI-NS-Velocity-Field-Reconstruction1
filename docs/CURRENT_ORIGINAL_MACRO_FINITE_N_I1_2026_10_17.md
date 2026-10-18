# Current-N macro six-I1 signed ordinary C1 enclosures

Tasks 14/15 and the long-term goal remain **ACTIVE / INCOMPLETE**. This
stage encloses the current full finite-N macro six-I1 volume terms and D,
uniformly for m in [0,1], Z in [37/100,93/250], in Wmacro=10000 units.
Complete macro pressure still requires its retained endpoint sources.
These local function enclosures do not establish global five-moment closure.

## Source and physical scaling

The current common N remains graph handle 18, the original huge positive
integer, not the number 18. Its saved logarithm is consumed exactly.
Ordinary Z differentiation holds the actual issued phase fixed because
phase=fractional_part(N*offset) has no Z dependence. No phase is selected,
and no physical-y derivative estimate is manufactured.

The new input adapter distinguishes

```
Qdrive = hydro + Pstar^2*pressure + R*dressed_swirl,
fullDrive = R*Qdrive,
a = h_b*R*D_over_R,
b = -sqrt(2*R)*h_b*Qdrive/(F0*barphi).
```

The dressed swirl already contains its F0^2 factor. The actual positive
a floor and the new directed log upper bound on a+b^2/a prove kappa<1,
so the prescribed cut is exactly one on this tube. Both inverse-F0 and
inverse-comparison-phi ordinary derivative products are retained.
Core, comparison and native phi floors are recorded separately.

The macro native fields are kept distinct from comparison fields and
from the micro inlet. A new fourth derivative extension consumes cached
lower rows and proves new angular/axial integral radii, then transports
native M/A/B/C with their original positive kernels and incoming values.
The generic p2 numerator retains native histories, Pstar^2 P0 pressure,
and the differentiated F0^2 B/C sectors. It is not the comparison drive.

The new Z/angle logarithmic engine propagates the exact Poisson denominator,
positive lambda floor, implicit inverse and T1 composition through ordinary
Z3. It retains all factorials. The zeroth loop bound is also obtained by
this propagation; no older generic C0 theorem is used to tighten it.

The physical finite-N changes obey

```
deltaU = Pstar*B_over_Pstar/N,
deltaM = positive-kernel transport(incoming second-micro deltaM, deltaU),
deltaQ = (2*Z*deltaU-(1-delta)*Z*deltaM-(1-Z^2)*deltaM_Z)/(1-delta*Z^2).
```

There is no incoming deltaU term. Native mean memory is not reset, and R
is not a mean scaling factor. DeltaU/M through Z3 give deltaQ through Z2.
Six product difference identities then give uniform ordinary C1 error
radii. The positive measure R^2*ell*dm, normalized by Wmacro, contributes
K=(1-exp(-2*ell))/2 once. Each normalized I1 defect is strictly below
2^-1024. The adapter adds the symmetric radius once around the accepted
signed background I1; it never selects an error cap or midpoint as a
function value. D is computed from those full I1 rows with the ordinary
Leibniz product rule and no separate repeated D radius.

## Actual enclosures and verification

The table is rounded outward to six decimals; exact binary records in the
report are authoritative. The tiny finite-N radius does not change these
displayed ranges from the background at this precision.

| full finite-N macro integral / Wmacro | ordinary Z0 | ordinary Z1 |
|---|---:|---:|
| q | [-1.455368, -1.443495] | [2.952000, 2.984001] |
| q_Z | [2.952000, 2.984001] | [7.999999, 8.000001] |
| q_squared | [4.167361, 4.236193] | [-17.371273, -17.044800] |
| U_q_Z | [4.368960, 4.440193] | [23.648000, 23.840001] |
| U_squared | [1.095200, 1.107073] | [5.920000, 5.952001] |
| U_U_Z | [2.960000, 2.976001] | [7.999999, 8.000001] |

D/Wmacro is [5.823035, 5.976471]; its ordinary Z derivative
is [37.371251, 37.996241].

New production exits 0 in 3.451 seconds. Independent
bounded Luna/max reviews, source AST guards, saved-report inspection,
directed binary record checks and exact source pins pass. Separate runtime
worker model metadata is unavailable. Only new stage arithmetic executes;
no accepted ancestor builder, checker, producer, owner, norm, recurrence
or quadrature is replayed. The inherited execution audit is preserved.

Report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_macro_finite_N_I1.json` is 19,955,216 bytes, SHA256 `066c6eaf48b49f836379b1604725f1b7f33ee9372bb2e4521652ddeaff1d7c3c`, with
2009 input and 389
implementation pins. Old 212 cumulative terms, six signed inner/lower
sources, previous errors and standalone outer RQy sources remain unchanged.
The new volume enclosure is stored separately and is not added again to
the cumulative ledger. P1 and m15 remain separate data.

NEXT extend the first-switch pressure volume using current-source finite-N
inputs and the accepted macro-upper/first-switch-lower RQy identity before
separate endpoint boxing. Keep the distinct native Jacobians and nonreset
histories. Then continue the independent reserved-pulse terminal functional
comparison. Complete pressure and global five-moment closure remain open.
Exact heat/energy, admissible stress/flat remainder, genuine time-order
recursion and summation, realizable corrections and independent corrected
Cartesian NS validation also remain incomplete. Goal ACTIVE.

Implementation: [native Z4 extension](../src/openai_ns_reconstruction/original_macro_native_Z4_extension.py),
[macro loop inputs](../src/openai_ns_reconstruction/original_macro_loop_Z3_inputs.py),
[Z/angle engine](../src/openai_ns_reconstruction/original_log_Z_phase_bounds.py),
[loop propagation](../src/openai_ns_reconstruction/original_macro_loop_Z3_bounds.py),
[physical error adapter](../src/openai_ns_reconstruction/original_macro_finite_N_I1_numeric.py),
[runner](../scripts/run_original_macro_finite_N_I1.py).

```
python scripts/run_original_macro_finite_N_I1.py --inspect-report
```
