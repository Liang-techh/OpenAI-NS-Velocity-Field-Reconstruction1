# Current-N second-switch and post-power I1, B/D/S1 and pressure source

Tasks 14/15 and the long-term goal remain **ACTIVE / INCOMPLETE**. This stage
encloses both native charts' full finite-N six-I1, B, D and S1 through R110,
uniformly in Z in [37/100,93/250], with ordinary Z orders 0 and 1. It preserves
the accepted first-switch parent and the old 212-term ledger. This is a partial
pressure/m15 construction, not global five-moment closure or an NS validation.

## Actual functions, memories and measures

The authoritative chain is the saved R100 C6 -> C5 -> Q4 source, not an older
standalone conditional R110 report. The second switch uses
R=100*exp(hbS*(1+t)), 0<=t<=1, physical Jacobian hbS. U is the first-switch
exit function and U_y=0. Its phi logarithmic derivative is
-hbS*R*D_over_R*(1-sigma)/2-2*sigma/5. The post chart uses
y=(log(11/10)-2*hbS)*p, R=100*exp(2*hbS+y), 0<=p<=1, phi=phi_R2*exp(-2y/5)
and the same U. The native source names are post_switch in C6 and post_power
in C5/Q4. Physical U=Pstar*V, applied once.

Each chart consumes the preceding native phi/U/M/A/B/C rows. The mean obeys
M(y)=U+(M_in-U)*exp(-y), retaining its inherited memory. Post native history
rows cover the whole segment via positive-kernel hulls of incoming M with U,
A with U^2, B with phi^2/2, and C with phi^2. Endpoint transport is retained
separately. The B endpoint convolution is (5/6)(exp(-4L/5)-exp(-2L)); the C
convolution is 5(exp(-4L/5)-exp(-L)). Scalar radial widths have zero ordinary
Z derivatives. Comparison invariant histories are not asserted equal to
actual native histories.

With Lz=1-delta*Z^2, h=1-delta, d=1-Z^2,
alpha=((1+delta)*Z*U-d*U_Z)/Lz and
beta=(-h*Z*(M_in-U)-d*(M_in_Z-U_Z))/Lz, Q=alpha+beta*exp(-y).
The pointwise Q enclosure is separate from the weighted integral I[q].
Second pointwise decay lies in [1-r,1], post decay in [10/11,1], r=2^-1024.
The six integrals use the analytic J0,J1,J2 moments and ordinary Leibniz rules.

Second normalized moments are Jk=exp(2hbS)*integral_0^1 exp((2-k)hbS*t)dt,
in [1,1+r], with Wswitch=10000*hbS. Post Wmacro=10000 moments are
J0=(1.21-exp(4hbS))/2 in [0.105-r,0.105],
J1=1.1*exp(2hbS)-exp(4hbS) in [0.1-r,0.1+r], and
J2=exp(4hbS)*(log(11/10)-2hbS) in [log(11/10)-r,log(11/10)+r].
The directed log-width proof hbS<r/100 and fresh derivative bounds certify
these enclosures. No hbS=0 value or midpoint is selected.

## Correlated B and finite-N errors

Write G=Q/Lz+Q^2/2-2Z^2*U^2/Lz^2=C0+C1*exp(-y)+C2*exp(-2y).
Then (R^2*G)_y=R^2*(2*C0+C1*exp(-y)), so
B/W=2*C0*J0+C1*J1. C2 cancels before boxing. This evaluates the actual
boundary functional without subtracting two wide endpoint intervals or
dividing that difference by tiny hbS. A new arbitrary-symbol identity binds
the actual Q4 reduced-pressure source to S=2G+G_y+D_density, using m_y=U-m;
this proves S1=B+D for finite N without assuming U_y=0. D retains the correct
coefficient 2*h*Z^2/Lz^2.

Source-related loop inputs use exact b=t0=0, actual second a and post a=4/5,
complete native p2 with the meridional Lz denominator, and the unchanged
N graph handle 18. Ordinary Z3 phase inversion retains factorials. Physical
deltaU is local Pstar*B_over_Pstar/N; the positive-kernel deltaM bound carries
first -> second -> post memory. The linear Q relation yields deltaQ through
ordinary Z2. Pointwise product differences are integrated with the true
physical mass once. Each normalized six-I1 error is strictly below r; one
symmetric r is added to each signed background row. D is computed from those
full rows with no extra D radius.

For B, deltaG=deltaQ/Lz+delta(Q^2)/2-(2Z^2/Lz^2)*delta(U^2), using ordinary
C1 products. Its two endpoints are bounded by 2*R_upper^2/W times sup|deltaG|.
The actual second 1/hbS factor stays logarithmic. Both B errors are below r;
B receives r once, then S1=B_full+D_full receives no extra radius. No parent
macro/first B value is fabricated. Physical-y numerical bounds are not claimed.

## Actual directed enclosures

Ranges are rounded outward to six decimals. The report's exact binary
directed records are authoritative; these are function enclosures over the
whole axial interval, not selected signed centers.

### second_switch, divided by Wswitch=10000*hbS

| quantity | ordinary Z0 | ordinary Z1 |
|---|---:|---:|
| q | [-2.907776, -2.889951] | [5.912000, 5.960001] |
| q_Z | [5.912000, 5.960001] | [15.999999, 16.000001] |
| q_squared | [8.351734, 8.455162] | [-34.660690, -34.170316] |
| U_q_Z | [8.749760, 8.868481] | [47.328000, 47.648001] |
| U_squared | [2.190400, 2.214145] | [11.840000, 11.904001] |
| U_U_Z | [5.920000, 5.952001] | [15.999999, 16.000001] |
| B | [1.322203, 1.464165] | [-35.966276, -35.266483] |
| D | [11.669323, 11.929714] | [74.831289, 75.903174] |
| S1 | [12.991526, 13.393879] | [38.865013, 40.636690] |

### post_switch, divided by Wmacro=10000

| quantity | ordinary Z0 | ordinary Z1 |
|---|---:|---:|
| q | [-0.305600, -0.303162] | [0.620000, 0.626561] |
| q_Z | [0.620000, 0.626561] | [1.679999, 1.680001] |
| q_squared | [0.875287, 0.889437] | [-3.647167, -3.580089] |
| U_q_Z | [0.917600, 0.932322] | [4.966400, 5.006081] |
| U_squared | [0.229992, 0.232486] | [1.243200, 1.249921] |
| U_U_Z | [0.621600, 0.624961] | [1.679999, 1.680001] |
| B | [0.137727, 0.154842] | [-3.781110, -3.698330] |
| D | [1.223064, 1.254843] | [7.848787, 7.978327] |
| S1 | [1.360791, 1.409684] | [4.067678, 4.279996] |

The joined macro/first/second/post six-I1 packet uses Wmacro. It retains the
parent macro/first rows, adds the whole second contribution once through the
actual hbS scale, and adds the full post rows directly. Parent and post radii
are not repeated. Joined D/Wmacro is
[7.046099, 7.231314], with derivative
[45.220039, 45.974568]. Whole-micro
volume is outside this numeric packet; complete joined B/pressure is still open.

## Source joins, audit and next task

The pressure graph preserves its complete 184898-node accepted prefix and
appends 410 nodes. Each new chart contributes its own local
S1/2 + own lower RQy - own upper RQy exactly once; cumulative S1 is not added.
The accepted parent already contains its first-switch upper endpoint. All
standalone endpoints, uncontracted source expressions and prior caps remain.
Both new seams have actual R100 endpoint/memory and C5/Q4 field/phase bindings.
Arbitrary-symbol endpoint slope/C2 identities, the Q4 rational map and the
C5 analytic periodic normalization prove the same finite-N Qy function at
each interface. These are source identities, not selected numerical values.
After replacing each new lower by its proved common upper source, a new
sum/negative-only exact integer linearization makes the first/second internal
upper coefficients zero. The final joined graph has no surviving dependency
on those four ordinary endpoint handles. Products, integrals and substitutions
remain atoms; the raw contracted-before-linearization view is also retained.
The outer R110 RQy is retained and not numerically evaluated.
Old 212 terms, six signed inner/lower contributions and all errors are unchanged.
P1 and m15 remain distinct data.

New production exits 0 in 9.571 seconds. Saved-report inspection,
directed-record validation and exact source pins pass. Only this stage's five
new builders execute. No accepted ancestor builder/checker/producer/owner/norm/
recurrence/quadrature or old exact-algebra computation was replayed. Fixed
Luna/max workers and root cross-review repaired pointwise-vs-integral Q, the
post mean bridge, constant width jets, positive-kernel coefficients and native
whole-segment hulls before acceptance. Separate worker runtime model metadata
is unavailable; no Astra child was created. The previous tautological same-state
endpoint fixture was rejected and is not accepted as an endpoint proof.

The first unaccepted integration attempt exited 1 during record-schema
validation after the new builders: a copied old source geometry record has
no directed flag. The repaired checker matches that complete geometry to its
SHA-pinned saved source and excludes it from new numerical counts; it still
checks endpoint ordering. All newly computed enclosures require directed
records. The successful production above uses the repaired runner. The initial
bridge-to-R100 premise is inherited from the accepted parent's source pins and
cache; this stage proves the subsequent two seams, not a new bridge reconstruction.

Report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_switch_post_finite_N_I1.json` is 24,654,478 bytes, SHA256 `9d0d957bae129853030bad4cd5c484eb27e8b13b0623601e458fc4e95e93bc2b`, with
2011 input and 401 implementation pins.
The prior brief is preserved byte-for-byte in
`docs/continuation/RESUME_FIRST_SWITCH_2026_10_17.json`; the local brief points
here without loading accumulated conversation history.

NEXT carry the same data through the R110-to-Rsh long_reshape, reference and
restoration toward Rm. Keep the R110 RQy interface until a new same-source
identity proves its internal cancellation; evaluate surviving external RQy
with source-correlated finite-N ordinary C1 control. Continue the independent reserved n1 own_swirl
versus negative full-pulse m15 functional comparison, without duplicating any
old ledger term. Global functional five-moment closure, exact heat/energy,
admissible stress/flat remainder, true time-order recursion/summation,
realizable corrections and independent Cartesian NS validation remain open.

```
python scripts/run_original_switch_post_finite_N_I1.py --inspect-report
```
