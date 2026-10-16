# First-micro m15 from the actual rminus inlet - 2026-10-16

`OriginalFirstMicroM15.evaluate(Z,coordinate,selected_sc_multiple=False)` computes the cumulative pressure moment through the missing first-micro interval [sc/2,1/4]. The six signed core/cutoff/lower terms at q=2 are preserved once. The leading first-micro increment is numerically integrated with full directed source errors, and the actual selected-N correction is included as three factored absolute-error sectors. Its unknown sign is not replaced by a selected interval center.

The full goal is **ACTIVE / INCOMPLETE**; tasks 14/15 remain IN PROGRESS. This closes the first-slice leading numerical integral with complete factored correction errors. It does not resolve the finite-N fast phase, establish useful relative cancellation or close the complete active/preRp pressure join.

## Original source and native measure

Write s0=sc/2, epsilon=Lambda^-1, rho=4exp(hb*s), R=epsilon*rho. The numerical integral uses native s directly; no impossible fraction s/(sc/4) is materialized. The own cut velocities are zero throughout this interval, so their local axial/cut-F increments vanish while the known leading pressure source and cumulative histories remain.

The regular leading increment is

```
epsilon^2*Pstar^2*hb * integral_s0^S -rho(s)^2*Np_regular(s) ds.
```

Directed cells retain local sigma(s) ranges and the actual whole-prefix field/comparison errors. Numerical midpoint approximations charge their complete distance to these integrated source enclosures. The interval width S-s0 retains its exact original lower endpoint even when it is much smaller than the point precision. The empty S=s0 branch returns the original inlet plus exact-zero increment terms.

## Singular pressure without subtracting rounded pressure endpoints

After exact inverse-hb cancellation, the original singular density is epsilon^2*A(s,hb)*d-sigma(s). With L=1-delta*Z^2,

```
A = 2rho^2*Z*(1-hb)*(Phi_actual/Phi_comparison)*Psi_comparison_rho/L,
A0 = 32Z*Psi_core_rho(4)/L.
```

The original width-change estimates give |A-A0|_norm<=hb*E. The increment is represented by two signed endpoint sectors

```
+epsilon^2*A0*sigma(S),
-epsilon^2*A0*sigma(s0),
```

and the complete error epsilon^2*hb*E*sigma(S), using 0<=sigma(S)-sigma(s0)<=sigma(S). The entire 1-hb factor remains in this proof. Pstar^2 cancels exactly in this singular contribution.

Each nonzero sigma value stays exp(odds)*kappa, odds=(1-s)^-2-s^-2. The unmaterializable lower sigma(s0) is retained as its own log atom; its pressure history is not reset. When the original flat cap applies, kappa is enclosed using sigma's positive tail bound; otherwise its moderate exponential is evaluated directly. The first-radial U_y/Q_y fields from the preceding milestone are complete; the inverse-hb term occurs in second radial pressure-source rows and is accounted for here.

## Actual finite-N correction from the reduced S1 identity

The independent actual-source formula is m15=S1/2-[RQ_y]. The modified-minus-leading correction is bounded before any huge exponential is materialized. Let DeltaX=X_N-X0 and E_X bound its absolute ordinary derivative rows. Every product includes

```
E_(ab) <= |a0|*E_b + |b0|*E_a + E_a*E_b,
```

including its ordinary Z1 Leibniz terms. The endpoint half-S1 coefficients for DeltaQ, Delta(Q^2), Delta(U^2) are 1/(2L), 1/4, -Z^2/L^2. The volume coefficients for DeltaQ, DeltaQ_Z, Delta(Q^2), Delta(UQ_Z), Delta(U^2), Delta(UU_Z) are

```
-1/(2L), (1-delta)Z/(4L), -1/4,
(1-Z^2)/(2L), (1-delta)Z^2/L^2, Z(1-Z^2)/L^2.
```

The boundary, volume and viscous endpoint units are respectively epsilon^2*rho^2, epsilon^2*rho^2*hb*ds and epsilon*rho. The first original source cell covers the entire [sc/2,1/4] interval; its local B bounds and accumulated positive-kernel mean errors therefore bound every volume prefix. The source-flat lower correction is explicitly checked for U, Q and Q_y, including all exported derivative rows, before lower endpoint errors are removed.

The output has thirteen cumulative sectors: six unchanged inlet terms, four leading increment sectors and three finite-N error-only sectors. Each correction sector has proxy zero, a symmetric C0/C1 coefficient and a separate positive log-bound atom. The common atom covers independently differentiated value and Z1 errors; no derivative is taken of a cap. Correction bounds may be too large to yield a useful relative accuracy. Their physical factors, cross errors, undivided B_phi/B_Z_phi endpoint contributions and original selected integer remain present.

## Evidence and next construction

The focused driver executes the actual accepted reduced-S1 function body using a symbolic scalar algebra, verifies its physical density identity and six independent native-coordinate value/Z1 integrals, executes the existing actual regular/singular Laurent-source and cut-F checks, and checks all quadratic-product error signs. Original packets preserve the six-term inlet, exact empty prefix, nonzero lower sigma atom, finite-N error-only semantics and directed cell refinements. The focused run exited 0 in 87.247 seconds; saved-report inspection also exited 0: 15 packets, 12 unchanged inlet partitions, 168 serialized increment-error checks, 42 exact empty-prefix rows, nine nonzero lower sigma atoms, 36 finite-N error-only sectors, 32 signed product-error cases, six cell-refinement comparisons and four rejected invalid queries. The 1,872,712-byte report SHA256 is `b289c8b4bf8811db05bce0f7568ada485ce4865c4d0ecc0fcdad63dde947546a`; it pins 1750 inputs and 33 implementations.

Use `python scripts/run_original_first_micro_m15.py --inspect-report` to inspect saved evidence; `--write-report` reconstructs it after relevant source changes. Existing read-only worker `/root/n1_atlas_map`: **GPT-5.6 Luna / max**; no new agent. It confirmed the singular split, half-S1 coefficients, measures, quadratic cross errors and checked lower endpoint zero.

Continue the remaining original first/second-micro and other active pressure charts, then compact finite-N pressure and the independently transported Rp value/Z1 cancellation. Retain the selected inverse/defect correlations. Complete all five numerical histories and global mixed joins, common N, analytic preheat/exact heat/energy, admissible stress/flat remainder, temporal n>=2 recursion/smooth sum, oscillatory correction and independent corrected Cartesian NS validation. All of these later gates remain open.
