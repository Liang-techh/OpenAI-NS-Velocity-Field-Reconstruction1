# Actual finite-N signed density and core/collar source differences

Tasks 14/15 remain **IN PROGRESS**. The long-term goal is **ACTIVE / INCOMPLETE**. This stage implements exact factored pressure-density differences using actual saved compact and Rm-patch C5 function ports. It also identifies the unresolved function-level bridge for the six signed core/collar terms and restores the formula for their signed singular-width remainder. No numerical finite-N integral is newly evaluated and no existing cap is removed or narrowed.

## Actual C5 source expressions

`original_finite_N_density_kernel.py` implements ordinary Z derivative jets and the actual half-S1 decomposition, with D=1-Z^2 and L=1-delta*Z^2:

```
F = Q/(2L) + Q^2/4 - Z^2*U^2/L^2,
v = -Q/(2L) + Z*(1-delta)*Q_Z/(4L) - Q^2/4
  + D*U*Q_Z/(2L) + Z^2*(1-delta)*U^2/L^2
  + Z*D*U*U_Z/L^2.
```

It computes Delta v and Delta F from one correlated baseline/correction pair before boxing. Every term remains, including 2Q*DeltaQ+(DeltaQ)^2, U*DeltaQ_Z+Q_Z*DeltaU+DeltaU*DeltaQ_Z, 2U*DeltaU+(DeltaU)^2 and U*DeltaU_Z+U_Z*DeltaU+DeltaU*DeltaU_Z. Ordinary C1 output requires the actual U_ZZ and Q_ZZ rows. These are ordinary derivatives, not Taylor coefficients.

The independent exact symbolic check proves S/2=d_y F+2F+v under the source continuity recovery identity:

```
L*(Q+Q_y) = Z*(1+delta)*U + 2Z*U_y - D*U_Z.
```

Thus the chart pressure difference has the correct physical FTC:

```
DeltaM15(x) - DeltaM15(1)
 = [R^2*DeltaF - R*DeltaQ_y](x)
   - [R^2*DeltaF - R*DeltaQ_y](1)
   + integral[1..x] R^2*Delta v dx/x.
```

Both the F and viscous endpoint terms are required. An integral of Delta v plus only the Q_y boundary would omit the R^2*DeltaF difference. The exported expression is a chart increment; the inherited incoming pressure difference is retained and remains numerically unresolved. It is not reset to zero.

`original_finite_N_signed_source.py` reads pinned saved graphs without constructing any owner. Physical baseline U=Pstar*V, Q=Pstar*generic.Q and U_y=Pstar*V_y apply Pstar once. Baseline Q_y is recovered using the same continuity formula above. The full physical U/Q/Q_y rows are the actual C5 source functions; corrections are their signed differences from that baseline. No magnitude bound, midpoint or approximation center is installed as a source value.

| chart | native domain | radius | full physical U C0/C1/C2 | full physical Q C0/C1/C2 | local nodes |
| --- | --- | --- | --- | --- | --- |
| compact | x in [1,2] | Rc*x | 262269 / 262270 / 262271 | 262364 / 262368 / 262374 | 332 |
| Rm patch | x in [1,e] | Rm*x | 204317 / 204318 / 204319 | 204412 / 204416 / 204422 | 359 |

The Rm baseline C6 graph is verified as the identical prefix of its finite-N C5 graph. The local expression graph references actual function nodes by report, SHA256 and node ID. It adds signed DeltaU, DeltaQ, Delta v, DeltaF, DeltaQ_y, endpoint and exact integral expressions without editing the original graph. These are whole-chart functions at arbitrary native x, not merely endpoint substitutions. The common exact N node18 remains N=2^(2^J), J=1358356628656378313; it is not materialized. No further N division is inserted. The original Rm fast phase and the compact selected-N controls remain in the referenced source functions.

## Six signed core/collar contributions

`original_inner_collar_source.py` pins the saved forward six-term partition and canonical inner pair53161/53162 plus lower pair36592/36593. It preserves all six numerical approximation/error records and exports the unresolved typed difference:

```
Delta_six^(k) = sum[i=0..5] F_i^(k) - (C_inner^(k)+C_lower^(k)), k=0,1.
s=(sc/4)*q, rho=4*exp(hb*s), R=epsilon*rho,
d rho = hb*(sc/4)*rho dq, q in [0,2].
```

Read-only equation comparison found the same intended core equation and inlet recipe, not a proved PDE or parameter mismatch. The missing bridge includes the forward P0/M2 numerical source versus canonical native P0 as a function of Z. The forward collar cell path uses endpoint-prepared enclosure rows; canonical fields carry q-dependent rho, the ratio ODE, Psi history and own mean. Their function-level identity is still unproved. Different numerical representations alone are not evidence of different exact functions.

The actual singular Q_rhorho sector has the signed source formula (h=hb):

```
A(s,h)=2*rho(s)^2*Z*(1-h)*(Phi_actual/Phi_core(rho(s)))
       *Psi_core_rho(rho(s))/L,
A0=32*Z*Psi_core_rho(4)/L,
epsilon^2*integral[0..s] A(t,h)*sigma_prime(t) dt
 = epsilon^2*A0*sigma(s)
   + epsilon^2*integral[0..s] (A(t,h)-A0)*sigma_prime(t) dt.
```

Term4 is the A0*sigma contribution; term5 is the signed remainder. The saved bound E for |A-A0|<=hb*E remains only a norm enclosure. The normalized remainder bound |integral(A-A0)*sigma_prime|/(hb*sigma(s))<=E does not choose a signed value. This algebra completes the singular derivative sector, not the full weighted N1p identity or six-source cancellation. The adapter's passed flag covers source pins, geometry and the declared unresolved interface; its exact_identity and numeric_value_available flags remain false.

## Saved views and checks

Three179-entry discrepancy views are byte-for-data unchanged from the accepted quiet-scalar report. The independent negative-pulse target and quiet pair have already been consumed; do not subtract either again. Original raw213/reduced183 archives, six signed core/lower terms, four signed leading leaves19/20/41/42, exact N, pressure tail throughRv and all finite-N caps remain unchanged. Canonical joint slot182 remains -(inner+lower+active+compact). P1/m13 and m15 remain distinct.

New production exits0 in 5.1 s. Report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_finite_N_signed_source.json`: 3,712,161 bytes, SHA256 `3c309218eb3798ba077f2ae0c1c74be9e5f0a0dce023e38c3e1dee3283d309c7`. It contains 1970 input pins and 270 implementation pins. Saved-only inspection verifies all2240 pins and packet invariants. Python compilation passes. Six density/FTC AST guards and nine core/collar guards pass.

Exact controls verify full-versus-factored Delta v and DeltaF C0/C1, the ordinary chain rule, the physical pressure FTC, rational arithmetic and provenance rejection. A nonzero quadratic counterexample has gap -18/847 even when the selected fast-linear correction vanishes. Seven new exact rational port-shaped expression controls cover both actual chart schemas and reject a cap as an original source function. These synthetic controls do not evaluate the original physical functions. No accepted ancestor owner, constructor, producer, fixture, recurrence or quadrature is replayed.

## Next executable work

Evaluate the new signed actual C5 density and endpoint differences in physical factored scales, preserving their common histories before numerical boxing. Start with compact and Rm patch, retain both endpoint terms and nonzero pressure inlet, and compare against the existing cap only as an enclosure. The new local graph exposes the algebra but still references unevaluated implicit-control, Duhamel-history and original function integrals; a source-bound numerical bridge is required. Do not promote the local expression or an error cap to an evaluated value.

For Delta_six, build the exact q-dependent collar source/history bridge and identify the core P0/M2 source function before claiming cancellation. The signed singular remainder formula can be bound to canonical function handles next; its current source formula does not establish the missing forward equality. The saved lower graph has s35848, hb8, epsilon35565, rho35851, sigma_prime35953, ratio35940/35941, Psi_core_rho35870/35871 and L30382/30381. A can be formed from these ordinary-C1 factors. A0 requires substituting core radial node35715 at rho=4 (constant13); existing nodes36065..36068 are radial order zero and must not be used as Psi_core_rho(4). No signed remainder integral is already materialized in that archive.

Global functional all-five/common-N closure, exact heat/finite-energy exterior, admissible stress and flat remainder, true n>=2 temporal recursion/smooth summation, realizable mean/oscillatory correction and independent corrected Cartesian NS remain incomplete. The final corrected-NS 1e-3 threshold is not applied to this source stage.
