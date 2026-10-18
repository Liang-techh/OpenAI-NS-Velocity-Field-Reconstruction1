# Native long/reference background I1/B/D/S1 and current-N pressure source

The long-term goal and Tasks 14/15 remain **ACTIVE / INCOMPLETE**. The native
constant-U background pressure sector now has ordinary-C1 enclosures through
Rz, uniformly for Z in [37/100,93/250]. The exact current-N source graph is
extended over long_reshape and reference. The new background numerical rows
do not enclose the finite-N oscillatory correction. Complete pressure and
functional global five-moment closure remain open.

## Same data and source scalars

The inlet is the published post-switch native U/M packet, with all five
ordinary Z rows retained. Both new charts use the same actual U function.
For each local physical log-radius distance y,
M(y)=U+(M_in-U)*exp(-y). The reference incoming M is the actual long endpoint
enclosure U+(M_R110-U)*[0,r], never a reset to U; r=2^-1024. Whole-segment
M uses the positive-kernel hull of incoming M and U.

Long length T is the exact singleton from the pinned C6 geometry certificate.
Reference length is Lref=10*(logC+logP)-T-8, with logP=exp(40)+11.
The actual reference graph uses true_core_logC node 35560 in node 81513,
not native alias node 6 directly. The scalar guard binds that source to
common.core.logC, the inherited selected-logC singleton and upstream
logCstar tuple, with the accepted native/common source identity. No range
endpoint or midpoint is selected as a construction parameter.

Wlong=Rsh^2=110^2*exp(2T), Wref=Rz^2=110^2*exp(2*(logref-8)).
Only their logarithms are materialized. Neither huge radii nor exp(-T) nor
exp(-Lref) are constructed numerically.

## Analytic integrals with nonzero tails

For a chart of physical length L, use
Jk=(1/W)*integral R^2*exp(-k*y)dy. Exactly,
J0=(1-exp(-2L))/2, J1=exp(-L)-exp(-2L), J2=L*exp(-2L).
Directed checks -L<log(r) and log(L)-2L<log(r) prove
J0 in [(1-r)/2,1/2], J1/J2 in [0,r]. These are enclosures of nonzero
function values; zero is an interval endpoint, not a chosen tail value.

Q=alpha+beta*exp(-y), where
alpha=((1+delta)*Z*U-d*U_Z)/Lz and
beta=(-h*Z*(M_in-U)-d*(M_in_Z-U_Z))/Lz,
h=1-delta, d=1-Z^2, Lz=1-delta*Z^2.
Whole-segment Q uses decay [0,1]; endpoint Q uses [0,r]. Both remain
separate from I[q]=alpha*J0+beta*J1. Six I1 and D use ordinary Leibniz
products and the coefficient 2*h*Z^2/Lz^2. C2 in
G=Q/Lz+Q^2/2-2Z^2*U^2/Lz^2 cancels before interval boxing, giving
B/W=2*C0*J0+C1*J1. S1=B+D uses the accepted, source-pinned Q4 identity as
an inherited lemma without re-executing it. No independent error radius is
added on top of the moment tail boxes.

## Actual signed background enclosures

Ranges below are rounded outward to six decimals. Exact directed records
in the report are authoritative. These are axial interval function
enclosures, not selected centers or terminal closure certificates.

### long_reshape, divided by Wlong

| quantity | ordinary Z0 | ordinary Z1 |
|---|---:|---:|
| q | [-1.452400, -1.446463] | [2.960000, 2.976001] |
| q_Z | [2.960000, 2.976001] | [7.999999, 8.000001] |
| q_squared | [4.184516, 4.218932] | [-17.289370, -17.126133] |
| U_q_Z | [4.380800, 4.428289] | [23.680000, 23.808001] |
| U_squared | [1.095200, 1.107073] | [5.920000, 5.952001] |
| U_U_Z | [2.960000, 2.976001] | [7.999999, 8.000001] |
| B | [0.666911, 0.726272] | [-17.958663, -17.657717] |
| D | [5.846315, 5.953164] | [37.460346, 37.906909] |
| S1 | [6.513227, 6.679436] | [19.501684, 20.249191] |

### reference, divided by Wref

| quantity | ordinary Z0 | ordinary Z1 |
|---|---:|---:|
| q | [-1.452400, -1.446463] | [2.960000, 2.976001] |
| q_Z | [2.960000, 2.976001] | [7.999999, 8.000001] |
| q_squared | [4.184516, 4.218932] | [-17.289370, -17.126133] |
| U_q_Z | [4.380800, 4.428289] | [23.680000, 23.808001] |
| U_squared | [1.095200, 1.107073] | [5.920000, 5.952001] |
| U_U_Z | [2.960000, 2.976001] | [7.999999, 8.000001] |
| B | [0.666911, 0.726272] | [-17.958663, -17.657717] |
| D | [5.846315, 5.953164] | [37.460346, 37.906909] |
| S1 | [6.513227, 6.679436] | [19.501684, 20.249191] |

The joined new background uses Wlong/Wref=exp(-2Lref) in [0,r], preserving
the signed long contribution through interval multiplication. Its D/Wref is
[5.846315, 5.953164], with ordinary derivative
[37.460346, 37.906909]. It is kept separate from the parent's
macro/first/switch full-finite-N numerical packet.

## Current-N source graph and proof scope

The full 185308-node parent prefix remains unchanged; the new source appends
382 nodes. Each chart contributes its own LOCAL
S1/2 + lower RQy - upper RQy exactly once. Cumulative S1 is not added.
Native and correction memories, common N=18, global fast phase, physical
Pstar and the distinct native Jacobians are retained.

The new source extension retains the post/long and long/reference
RQy endpoint differences. Their modified finite-N Uy/Qy congruence is not yet
proved. No unproved endpoint is contracted or deleted; all raw expressions,
standalone endpoints and prior caps remain available. Independent review
rejected a candidate which combined generic chain-rule fixtures and source
token bindings without proving the actual physical endpoint derivative
congruence. Those candidate proof helpers and the contraction path were
removed. Every new seam/theorem receipt explicitly remains unproved. The
accepted new source scope is the raw signed graph extension only.

The surviving reference upper RQy is retained and unevaluated. The old
212-term ledger, six signed inner/lower terms, reserved/target sources and
all errors are unchanged; P1 is distinct from m15. A source identity never
selects a cap as a function value. This stage does not supply native phi or
A/B/C numerical histories, finite-N loop norms, physical-y numerical
bounds, an exact heat exterior, stress admissibility, temporal recursion,
realizable corrections or independent Cartesian NS residual validation.

Production exits 0 in 7.786 seconds. Saved-report inspection,
new directed-record checks, independent Luna/max review and source pins pass.
Only the two new builders execute; accepted ancestor arithmetic is not
replayed. Separate worker runtime model metadata is unavailable. Review
fixed the whole-segment/endpoint Q distinction, the scalar logC graph alias,
and a scalar-tail fixture that initially masked the endpoint-box interface.
The current-N seam proof candidate was rejected and removed as described
above; passing source checks certify the raw graph, not that missing theorem.

Report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_long_reference_background_I1.json` is 25,144,415 bytes, SHA256 `7fce21dc6e763f42b99ef7fcf7416decaf848a5124cda99741222dba9619eb1b`, with
2012 input and 404 implementation pins.

NEXT implement source-bound native long/reference controls. Use actual
anchored-G C0 bounds only as envelopes, derive B_Z1..4 from the exact rational
G gradient and the positive saved R110 phi jet, and retain full p2/pressure
and history terms. Then establish current-N signed pressure enclosures and
complete any remaining modified-derivative seam proof before contracting.
Continue restoration/post_restore toward Rm and the independent reserved
n1 own_swirl versus negative full-pulse m15 functional comparison. Goal ACTIVE.

```
python scripts/run_original_long_reference_background_I1.py --inspect-report
```
