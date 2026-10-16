# Original cumulative macro m15 through R=100 - 2026-10-16

`OriginalMacroM15.evaluate(Z,q)` continues the checked seventeen-sector second-micro s=2 moment over the original macro chart, q in[0,1]. It adds an analytically integrated signed background proxy and complete directed source/geometry/selected-N errors. The output has28 cumulative sectors. The previous seventeen sectors are preserved exactly once.

The full goal is **ACTIVE / INCOMPLETE**; tasks14/15 remain IN PROGRESS. The signed leading pressure approximation now reaches physical R=100. Resolved finite-N signed centers, useful total relative accuracy, later active/compact pressure and independent Rp numerical cancellation remain open. This is a spatial continuation, not genuine temporal n=2 recursion.

## Exponential mean and original reduced pressure identity

The source is the accepted original macro pressure row owner with its genuine meanZ3/QZZ and ordinary Uy/Qy Z1 error rows. The nonzero second-micro mean is not reset. For the signed frozen-velocity proxy,

```
U = U0 + epsilon*Psi_c
M = U + epsilon*(M_c-Psi_c)*E
Q = Qaxis + epsilon*(Qc+Qa*E)
E = exp(-Y*q)
Rproxy = 4epsilon*exp(Y*q)
```

All source functions and parameters remain the same accepted family. The proxy does not become the selected finite-N field definition. The actual leading field continues changing according to the original macro source modes.

The implementation uses the original identity

```
m15_increment = half(S1) - [R*Q_y].
half(S1) = [R^2*F(U,Q)] + integral R*dR*G(U,Q,Q_Z,U_Z).
```

The exact half-S1 coefficients come from `pressure_functions` in the accepted original preRp moment source. Both the boundary polynomial F and volume polynomial G have separate epsilon degrees0,1,2 and mean exponential modes0,1,2. This preserves epsilon-weighted terms before adding them to an O(1) background.

Writing x=Y*q, the normalized positive volume kernels are

```
W0 = -expm1(-2x)/2
W1 = exp(-x)*(-expm1(-x))
W2 = x*exp(-2x).
```

W2 is the exact resonant mean-squared mode. Because Rproxy^2*E^2 is constant, the F2 boundary term cancels identically. The remaining half-S1 coefficient is integrated as

```
(2F0+G0)*W0 + (F1+G1)*W1 + G2*W2.
```

This avoids numerical subtraction of almost equal lower/upper boundary values at very small q. It does not drop the quadratic incoming-mean volume term. Physical factors are exp(2logRproxy-degree*logLambda), separately for each epsilon degree.

## Actual source and geometry errors keep hb

The actual geometry remains

```
L = Y - 2hb
Rtrue/Rproxy = exp(2hb*(1-q))
R*dR = Rtrue^2*L*dq.
```

The discrepancy between this geometry and the Y-proxy is charged explicitly. The original hb remains a formal log atom on each source/geometry error sector. A numerical cap is used only on the remaining hb in quadratic error products and positive exponential bounds.

The uniform original modal source norm is reconstructed before its one hb factor, including all first/second/macro history and hydro/pressure/swirl parts. The difference between actual and proxy mean contains the second-micro incoming variation, the evolving-Psi Volterra source, and the difference between exp(-L*q) and exp(-Y*q). The large mean relaxation itself remains a signed exponential mode; it is not bounded as a microscopic variation.

For example, a product difference divided by hb is bounded from

```
a*error_b + b*error_a + hcap*error_a*error_b.
```

Both ordinary C0/C1 rows, coefficient derivatives, QZZ and all quadratic cross terms are retained. Radius and R^2*L measure errors are likewise bounded before hb. Three signed epsilon sectors and three corresponding hb-factored source/geometry error sectors are added. The physical macro volume has no extra hb factor.

## Both viscous and finite-N endpoints

The frozen proxy has Q_y=-epsilon*Qa*E. Its Rproxy*Q_y is generally nonzero but constant, so its increment cancels analytically. This statement applies only to the frozen proxy. The actual changing macro Uy source and true mean/geometry retain separate upper/lower leading viscous error terms with their physical R*epsilon*hb factors.

The lower source is independently obtained from the original second-micro pressure owner at s=2. Its exact interval endpoints are transferred into the macro Taylor context before algebra; no bounds are replaced by midpoints. The lower actual finite-N pressure derivative uses bridge_second, while the macro upper one uses bridge_macro. Both endpoints remain nonflat and nonzero.

Three further error sectors retain the finite-N half-S1 boundary, volume and RQ_y contributions, including ordinary Z1. The volume uses complete prefix source covers and the physical bound

```
integral R*dR <= Rupper^2 * min(Y*q, 1/2).
```

The selected integer and all first/second/macro histories are unchanged. Enclosure sectors do not select a signed correction value. The genuine M_ZZZ/Q_ZZ error can remain too large for a useful total relative accuracy.

## Reproduction and focused evidence

Run `python scripts/run_original_macro_m15.py --inspect-report` to inspect the saved source-pinned report. Use `--write-report` only after relevant source changes. The focused driver independently executes the accepted S1 source, compares physical pressure-density integrals with the analytical macro proxy for value/Z1, checks all three kernels and the exact viscous proxy identity, and checks hb-factored product/radius/length errors. It then evaluates actual original packets through R=100, verifies the unchanged inlet, serialized error bounds, exact empty-prefix rows and both physical endpoint budgets.

The final focused run exited0 in160.262 seconds; saved-report inspection exited0. It produced19 original packets,19 exact complete-inlet replays,418 serialized increment-error conditions,66 exact empty rows,98 hb-unit/R100 conditions,35 history/resonance conditions,64 retained lower-endpoint budget conditions and five invalid-query rejections. Independent checks passed: two actual S1 source-density identities, six original physical S1 value/Z1 integrals,24 additional macro pressure-density value/Z1 comparisons,12 modal integrals, two exact viscous-proxy identities, two small-prefix/resonance checks,64 ordinary hb-quadratic product conditions and12 true radius/length error conditions.

The4,831,445-byte report SHA256 is `076900a16aaea965e229d652b311dd1bcbd2fc73645061fe56aa58575acaf54c`; it pins1756 scientific inputs and45 implementations. Existing read-only worker `/root/n1_atlas_map`: **GPT-5.6 Luna / max**; no new agent. Its review confirmed the kernels, half-S1 coefficients, leading/finite-N endpoint retention and the stated proxy/error scope, with no material underbound found.

## Next construction

1. Continue into the original first-switch chart from R=100. Preserve **all28** macro sectors as the next inlet, exactly once.
2. Use `OriginalFirstSwitchNumeric` and its actual prescribed R100/Rm source. The first-switch source controls change; do not reuse the frozen-macro Uy mode outside its chart. Retain its original width log atom, mean history and both pressure endpoints.
3. Use the original switch_first generic derivative bounds and complete first/second/macro/switch correction prefix. Keep source caps separate from selected signed values. Integrate the actual reduced-S1 volume with the one correct native measure.
4. Continue remaining active charts and compact finite-N pressure, then preserve reserved inverse/source correlations for an independent Rp value/Z1 comparison. Interval overlap and a terminal theorem alone do not measure cancellation.
5. Complete all five numerical histories/global mixed joins/common N, original flatten/preheat/exact heat/energy, admissible stress/flat remainder, genuine temporal n>=2/smooth sum, oscillatory correction and independent corrected Cartesian NS.
