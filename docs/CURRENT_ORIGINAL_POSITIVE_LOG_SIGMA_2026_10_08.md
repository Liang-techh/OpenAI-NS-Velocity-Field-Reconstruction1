> Successor: [CURRENT_ORIGINAL_NORMALIZED_MICROSCOPIC_Q_2026_10_08.md](CURRENT_ORIGINAL_NORMALIZED_MICROSCOPIC_Q_2026_10_08.md), source [ca0fb994](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/ca0fb99479d436ef5b91d90efc16e38f028f5036). Typed xi/k geometry and original q/ordinary seam source are executable; actual amplitude/inverse/density and full prefix integral remain open.

# Same original positive logarithmic cutoff and transition support scale

Accepted implementation [f6e3f66d](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/f6e3f66dc03f204f516ed99d4b42654e2bf17f3d); predecessor [CURRENT_ORIGINAL_TRANSITION_RIGHT_HALF_2026_10_08.md](CURRENT_ORIGINAL_TRANSITION_RIGHT_HALF_2026_10_08.md). Root owns math/code/compute. Existing read-only GPT-5.6 Luna/max reviewed original logistic derivative, O3 Delta/a/b, mu/eta/dstar and actual phase contracts. No new child or worker edits/tests/compute.

## Original source refinement

Write s for the original O3 transition offset log(R/Rd); source code calls it t. It is not the physical time in the Navier-Stokes blowup scaling. For0<s<=1/2 the unchanged defining cutoff is

```text
L(s) = 1/(1-s)^2 - 1/s^2
sigma(s) = exp(L)/(1+exp(L))
1-sigma(s) = 1/(1+exp(L))
Delta(s) = kappa-2 = 2*mu*sigma(s)
a(s) = 2+Delta(s), b(s)=0.
```

The new source representation stores L as a formal logarithmic offset and the strictly positive bounded coefficient1/(1+exp(L)). The bounded exponential is used only in that positive denominator, so its directed upper tail does not erase sigma's nonzero numerator factor. Monotone endpoint values enclose the original odds over the whole source interval. No sigma/parameter cap midpoint or endpoint is substituted for a defining field value. Choosing exact domain endpoints2^-128 and1/2 is a source partition, not a change in the original field.

Ordinary y=log R derivatives are retained directly:

```text
L_y = 2/(1-s)^3 + 2/s^3
L_yy = 6/(1-s)^4 - 6/s^4
sigma_y = sigma*(1-sigma)*L_y
sigma_yy = sigma*(1-sigma)*(L_yy+(1-2*sigma)*L_y^2)
```

No second Taylor-factorial conversion occurs. Original sigma/Delta/a Z rows are exact zero because the original s and chosen mu are Z independent; all other E/E_Z/p2/mixed roots stay the same native source functions. Original endpoint0 and its flat slow jets are retained, and source intervals containing0 continue to carry explicit tail/refusal status rather than fabricated positive lower bounds. The original selected eta/dstar, q formula, native radius and conditioned inverse/first-jet implementation are unchanged. Accepted defining source modules and prior numerical archives are not edited.

## Genuine source result and exact cumulative operator

At candidate N1024, both Z[.36,.38]/[-.38,-.36] now have whole-source queries on[0,2^-128] and[2^-128,1/2]. The latter has a strictly positive Delta lower above the original eta upper, so the original cutoff returns the flat q branch and every required q slow derivative exactly zero. Actual original phase/conditioned inverse/first-jet/density evaluation then gives zero own C0/Z integral changes. This closes the former positive-cutoff lower-bound refusal over that complete collar.

Joining the accepted right half gives the original quiet interval[2^-128,1], with true log-radius width w=1-2^-128. The complete explicit-input operator is

```text
m_out   = exp(-w)       * m_in
h_out   = exp(-3*w/2)   * h_in
k_out   = exp(-3*w/2)   * k_in
e_out   = exp(-w)       * e_in
p_out   = p_in.
```

It applies identically to ordinary Z rows. Both original source widths enter once; every own integral increment is zero by the exact defining cutoff theorem. Original nonzero incoming histories are explicit arguments; pressure p/p_Z and separate P0 are not reset. The prefix[0,2^-128] remains unresolved and is not skipped. This quiet-interval operator is not yet an entire-transition numerical history because actual incoming at its left boundary still depends on the prefix and axial/buffer route.

2^-128 is approximately2.94e-39 in this radial coordinate. That number is not an overall project completion percentage or physical core width. It is still enormously larger than the source-derived possible active-q support described below.

## Source-derived microscopic support

On0<s<=1/2, L(s)>=-1/s^2 and1+exp(L)<=2, hence

```text
Delta(s) >= mu*exp(-1/s^2).
D = log(mu)-log(eta) > 0.
s >= D^(-1/2)  implies  Delta(s)>=eta, q=0.
```

Therefore the original nonzero-q support lies within s<=W with W=D^(-1/2). The entire original mu/eta parameter covers are retained; W remains a positive formal scale with no unmaterializable exponential or huge denominator allocated. The current source logW upper is approximately **-4.707705336740406e17**. This explains why ordinary decimal/dyadic endpoint sweeps can miss the real activity while zero lower cutoff tails continue to make their broad source boxes appear mixed.

The support theorem is a bound on the original radial cutoff, not a new parameter choice, a zero-field statement on its active prefix, physical blowup evidence or n-dependent scale recursion. The exact active endpoint s0 is retained. A typed coordinate s=W*xi and same-source normalized Delta/eta is the next route to an actual microscopic integral; multiplying an arbitrary cap by an ordinary width does not select a numerical field value.

## Evidence and appropriate checks

All paths under `experiments/root_st073/` use stem `lei_ren_part1_paper_compliant_current_transition_positive_log_sigma`:

- `.py`: unchanged original logistic source/ordinary derivative representation; genuine native source/root/q/actual phase/inverse/density queries; source-derived support theorem; exact two-cell cumulative operator.
- `.json`: same original family/parameters/N, support log-scale, two lossless source archives, original unresolved prefix and exact quiet-interval outputs; all full-construction gates false.
- `_positive.json.gz` / `_negative.json.gz`: original unrefined source record and refined same-source sigma/Delta/a/q records, actual native phase/first-jets/density, true geometry/masses, accepted right-half source reference and composite incoming operator. Total198,904 compressed bytes.
- `_check.py/.json`:39 independent original logistic C0/y1/y2 derivative comparisons including nonmaterializable left-tail factors and midpoint symmetry;9 independent support log-ratio checks; original mu/eta/positive-a certificate bindings; new source roots/ordinary q reconstruction; actual phase/density C0/Z and true mass;20 new left-collar zero integral rows,20 composite rows and independently collected total-width decay. Earlier source owners/integrals/scalar regressions are not rerun.

Focused check PASS;1261 Git-index hashes match. Producer67.469s; saved-source focused checker1.484s.240-digit native intervals/300-digit ordinary arithmetic retained; independent fixture points do not define native field values.

## Detailed next actions

- [x] **POSITIVE-ORIGINAL-LOG-SIGMA:** exact source logistic value and genuine ordinary y2/Z1 rows on positive left-half intervals, original endpoints retained, parameters/field unchanged and source-derived positive lower bounds.
- [x] **ACTUAL-QUIET-COLLAR-C1:** whole[2^-128,1/2] source evaluation joined with accepted[1/2,1], all own density/integral C0/Z zero, true width and explicit incoming/p/P0 memory.
- [x] **ORIGINAL-ACTIVE-SUPPORT-SCALE:** source-bound W=[log(mu)-log(eta)]^-1/2 and directed support proof; no unsupported claim of exact support endpoint or full numerical transition.
- [x] **PARAMETER-NORMALIZED-MICROSCOPIC-COORDINATE (source/geometry component; see successor):** introduce a typed original endpoint/radius identity for s=W*xi using the same original mu/eta expression and its full directed cover. Preserve original dy/ds=1; integration in xi uses ds/dxi=W once. Native ordinary source derivative rows are already y derivatives and must not receive this Jacobian again. Retain actual N*y phase and the exact original s0 active endpoint. Never allocate a denominator with astronomically many digits or pick a parameter midpoint.
- [x] **CORRELATED-DELTA/ETA-AND-LAZY-Q (queried source component; full prefix integral open):** collect the original common log(mu/eta) terms before evaluating L(s)+log(2mu)-log(eta). Keep the original cutoff/seam/active/flat branches and original q formula. On the same active part0<=Delta<eta retain eta<2eta-Delta<=2eta before the positive square root, rather than losing its upper bound through separate giant-log subtraction. Original mixed/ordinary q derivatives and directed error must remain source-derived; do not define a new q.
- [ ] **EXACT-O3-AMPLITUDE-FACTOR:** retain original Utheta/Pstar logarithmic amplitude and its ordinary Z dependence in the original E/p2/root packet before the conditioned inverse. Preserve the common source/ledger/family and source amplitude identities. Avoid the old separate conservative logR/amplitude bounds turning a positive original p2 into a sign-indefinite inverse source. This is a source factor refinement, not an added fitted field.
- [ ] **ACTUAL-MICROSCOPIC-PREFIX-C1-INTEGRAL:** original normalized coordinate/phase/inverse/primitive/density value-Z evaluation and directed positive original Duhamel weights, with true W and derivative-under-integral justification. Cover the real active/seam/flat support and retain a source-derived remainder/refinement error; no prefix reset. Compose the accepted quiet operator after the actual prefix.
- [ ] **CORRELATED-AXIAL-b:** retain same-source b/amplitude/radius factors and actual signed crossings before Delta=b²/2. a2 exact cancellation stays. Fresh native source boxes must reduce actual branch refusals without forcing positivity or erasing zero-shear edges.
- [ ] **FULL-AXIAL/BUFFER-SOURCE-INTEGRALS:** source-bound period means/slow variation and retained actual boundary phase for enormous widths; genuine original[0,9]/[9,11] buffer with actual incoming. The accepted one-period8/32 result remains a local numerical component, not a full buffer integral.
- [ ] **ACTUAL-NUMERICAL-RC-TARGETS/FIVE-FUNCTION-CONTROLS:** complete original route and compatible analytic pressure/heat targets with same physical radius/unit correlations and controlled errors; original background once/P0 separate; source-dependent independent controls and continuous-Z terminal identities.
- [ ] **WHOLE-AXIAL/HIGH-JETS/GLOBAL-N:** cover original Z0/signed-crossing/entire required axial domain and higher mixed source derivatives before selecting one globally compatible finite N. This two-tile candidateN1024 stage does not admit that gate.
- [ ] **MATCHED-HEAT/STRESS/FLAT-REMAINDER/TRUE-RECURSION/CORRECTED-UVW:** pressure/five moments/tail/joins/exact heat first; original cone and flat remainder; true n-dependent coefficient recovery/independent repairs/smooth solenoidal summation; both oscillatory quadratic stress-cancellation families and independent Cartesian physical field/residual/dynamics. Geometric scaling or local radial cutoff bounds cannot substitute for any of these layers.

Mark a named task DONE only with its actual source-bound result and receipt. The original positive cutoff, quiet operator and support scale are completed components; the microscopic prefix and all dependent physical construction layers remain open.
