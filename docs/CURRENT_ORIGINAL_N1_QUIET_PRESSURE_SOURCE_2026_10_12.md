# Actual source-scaled signed quiet pressure moments

The actual quiet contribution to m13/m15 is now computed with a directed signed coefficient enclosure. Its maximum coefficient relative width is **0.0002803156569**, below the 0.001 target. The absolute original Rc and Rc squared remain exact physical source scales. This is a quantitative contribution on the genuine repaired 2Rc to Rp segment, not a complete pressure-tail defect or full n1 repair. The full reconstruction goal remains **ACTIVE / INCOMPLETE**.

Read [the pressure-tail reduction checkpoint](CURRENT_ORIGINAL_N1_PRESSURE_TAIL_REDUCTION_2026_10_12.md) for the accepted all-17-chart mean interface, exact forcing endpoint Rv and quiet physical pressure identities. The previous generic quiet primitives are reused. The new calculation binds their coefficient to the actual native-frame mean kernel, keeping the correlations before numerical evaluation.

## Exact original source reduction

The accepted analytic native inlet gives common mean kappa*Z. The accepted native physical constant satisfies M=Pstar*kappa and Pstar=exp(logP). With Y=exp(40), the actual original turnoff kernel is

```text
A(Y) = integral[1,Y] exp(y-Y)*sigma(1-log(y)/40) dy
F = A(Y)+exp(1-Y)
kappa = 4*exp(-Tw-12-logP)*F
Rp = Rc*exp(Tw-2)
K = Rp*Pstar*kappa = 4*Rc*exp(-14)*F.
```

The implementation identifies this integral node in the actual accepted native-frame graph, checks its defining integrand/lower endpoint/cutoff and verifies the current native physical M identity. It uses the true compact-repaired quiet source; its incoming mean has not been reset. The graph's `original_logRc` is the absolute physical radius origin. The older quiet transport's relative phase offset must not substitute for it.

Let L=log(Rp/(2Rc)), the same accepted positive quiet length. The physical primitive then gives

```text
m13_quiet = +Rc * [2*exp(-28)*F^2*(1-exp(-L))]
m15_quiet = -Rc^2 * [4*exp(-28)*F^2*L].
```

These are signed coefficients times exact physical radial scales. Both ordinary Z1 rows vanish because the actual kappa and source endpoints are Z independent. The pressure-forcing integrals have the m13 plus/m15 minus signs prescribed by the n1 moment definitions.

## Numerical calculation and error control

The original integral is changed to s=Y-y. On its relevant window,

```text
x(s)=-log(1-s/Y)/40
D(s)=1/x(s)^2-1/(1-x(s))^2
sigma(x)=exp(-D)/(1+exp(-D))
H(s)=s+D(s).
```

An independently located approximate saddle is only a quadrature coordinate. Directed derivative signs bracket the true saddle. Exact derivative identities and a rational lower-bound argument prove H is convex for 0<x<=1/100. No approximate saddle is installed as a field value.

An autonomous numerical choice uses 2048 trapezoid cells over the center plus/minus 12 times sqrt(center), at 140 decimal digits of interval precision. Each cell retains the directed trapezoid error `h^3*sup(abs(f''))/12` for normalized f=exp(H(center)-H(s)). Endpoint monotonicity or the certified saddle bracket bounds f; its derivative formula bounds f''. Positive left/right convex tails and the remaining far tail are retained. The sigma denominator error and original inlet-memory term exp(1-Y) are retained as directed nonzero terms.

The result exports log coefficients, their signs and the original Rc/Rc squared function factors. It materializes no physical radius, huge N, Pstar or inverse-mu radius. The quoted 0.0002803156569 is a **normalized quiet coefficient enclosure width** (about 0.028%), not a full physical field error, full-tail width, NS residual or percentage of project completion.

## Callable artifact and checks

Module: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_quiet_pressure_source.py`.

```python
owner = CurrentOriginalN1QuietPressureSource()
view = owner.evaluate(cells=2048)
```

Scoped gate: `current_original_n1_quiet_source_scaled_signed_pressure_moments_installed`.

The independent checker verifies the H derivatives/convexity proof, actual physical scale cancellation, signs and Rc/Rc squared factors against the physical pressure primitives. The production kernel enclosure contains a separately implemented direct original-sigma quadrature. Both production pressure coefficient enclosures contain independently calculated primitive coefficients. The 0.001 coefficient-width target, positive retained tails and checked-owner constructor pass. Existing ancestor checkers are reused through their accepted receipts.

## Next construction work

- [x] Bind the actual turnoff scalar kernel and native Pstar normalization.
- [x] Cancel Pstar/logP and Rp/Tw exactly before evaluation; retain absolute original Rc.
- [x] Evaluate the full positive kernel with directed quadrature error and tails.
- [x] Compute actual signed quiet m13/m15 coefficients, both below the 0.001 width target.
- [ ] Compose physical U/m/q and pressure moments through all 17 pre-repair charts using the newly accepted mean atlas.
- [ ] Preserve signed correlations in those integrations; coarse independent source caps alone do not give the required full signed defects.
- [ ] Include unchanged Rin to r_minus and the genuine Rc to 2Rc repair band with its own histories.
- [ ] Include actual selected pulse charts Rp to Rv with the original absolute radius geometry and independent P0.
- [ ] Add the signed quiet contribution here, apply lower endpoint terms and the exact zero upper terms at Rv, and close the complete m13/m15 forcing tails.
- [ ] Combine all five defects, bind positive I2/Astar and normalized 2+3 repair matrices, and reconstruct the n1 own M/V/P after repair.
- [ ] Finish axis remainder/joins, genuine n2 and temporal summation, admissible stress, oscillatory cancellation, flat remainder and independent corrected NS diagnostics.

All full signed pressure-tail, n1 repair, temporal-recursion, full corrected NS and absolute-radius materialization gates remain false. The linked preceding checklist remains authoritative for the other source segments.
