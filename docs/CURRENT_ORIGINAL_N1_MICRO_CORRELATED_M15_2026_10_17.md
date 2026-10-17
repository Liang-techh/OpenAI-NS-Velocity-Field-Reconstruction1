# True first/second-micro pressure in a correlated six-I1 basis

Tasks 14/15 and the long-term goal remain **ACTIVE / INCOMPLETE**. The actual modified Q4 source on first_micro [sc/2,1] and second_micro [1,2] is now implemented in an append-only six-I1 pressure representation, with ordinary Z0/Z1 rows. It eliminates total-y Q_y and U_y from the volume integral before any numerical boxing. Actual outer R*Q_y, including its slow and fast finite-N source functions, remains. No numerical full-micro enclosure or five-moment terminal closure is claimed.

## Exact pressure reduction on the actual two charts

Let h=1-delta, d=1-Z^2, L=1-delta*Z^2, y=log(R/Ra), and define I1[f]=integral R^2*hb*f ds on the two native micro charts. The six integrands are Q, Q_Z, Q^2, U*Q_Z, U^2, U*U_Z. They use the current physical U,Q rows from `current_original_modified_Q4_pressure_C3_functions.json.gz`, rather than an older 17-chart active sum.

```
B = Delta(R^2 Q)/L + Delta(R^2 Q^2)/2 - 2 Z^2 Delta(R^2 U^2)/L^2
D = -I1[Q]/L + h Z I1[Q_Z]/(2L) - I1[Q^2]/2
    + d I1[U Q_Z]/L + 2 h Z^2 I1[U^2]/L^2 + 2 Z d I1[U U_Z]/L^2
S1 = B+D
Delta M15 = S1/2 - Delta(R Q_y)
M15(2) = H + S1/2 - (R Q_y)(2),  H=M15(sc/2)+(R Q_y)(sc/2).
```

The U*Q_Z coefficient is d/L, with no extra Z. This follows from physical mean recovery Q=(2ZU-hZm-dm_Z)/L, m_y=U-m and Q+Q_y=((1+delta)ZU+2ZU_y-dU_Z)/L. Differentiating R^2 gives 2R^2. Six fresh symbolic checks prove mean recovery, the reduced R^2 pressure identity and its ordinary Z derivative, and same-function seam cancellation. These are identities of source functions, not sampled interval overlaps.

There are exactly24 new definite-integral leaves: six per chart, each with value and ordinary Z1. Their integrands contain the physical weight R^2*hb exactly once, their lexical binders are unique, and their endpoints are Z-independent. Product and quotient rules use ordinary derivatives without factorial scaling. The base archive has 183,700 nodes; the local suffix has 624 nodes. Direct archived second-micro S1/m15 ports remain references for comparison, never a second additive copy of the same contribution.

## Consumed same-source seam and retained outer boundary

The first/second interface was already proved in the saved active-radial4 receipt and `micro_viscous_m15_numeric.json`. This stage reads those accepted results without executing their fingerprint producer. It checks that the active graph extends the exact Q4 graph prefix, maps all14 saved identity rows to the current U,m,Q physical-y0/1 and RQy ordinary-Z0/1 ports, and checks the fingerprint source scope keeps Z free while substituting only native coordinate s=1.

Both charts share coordinate66988, radius68644, physical measure8, phase97301, original roots and global offset. Second incoming Duhamel histories and pressure integrals are exactly the first endpoint functions. Consequently the internal reduced B and viscous RQy cancel as functions. Only global endpoints are constructed in the reduced output. Any future error cancellation must also use a single shared approximation/error function at the interface; independently boxed errors cannot simply be discarded.

The physical shifted inlet is bound to the accepted parent H and the actual Q4 initial RQy ports. It is not replaced by its approximate center. H is the sc/2 inlet: the full first_micro integral already includes the accepted [sc/2,3sc/4] flat prefix, so that saved prefix increment is not added a second time. The two-chart outer endpoint still contains actual U_y=Pstar*(background V_y+slow_B_y/N+B_phi), and the full Q_y recovery. The fast B_phi term has no 1/N prefactor. No old finite-N cap is set to zero or removed.

## Exact finite-N volume difference and remaining quantitative gap

The new same-source split uses du=Pstar*delta_V and dm=Pstar*history.m from the actual prescribed-loop chart, including the complete predecessor history. It recovers dq=(2Zdu-hZdm-d*dm_Z)/L and defines U0=U-du,Q0=Q-dq. Six polynomial differences retain all linear and quadratic cross terms, with ordinary Z1 propagated by the same graph product rule. The source profiles retain the fixed common N, its phase and its histories; this is not an N-independent coefficient expansion and supplies no numerical error radius.

In particular ordinary Z1 of I1[Q_Z] and I1[U*Q_Z] needs Q_ZZ, hence the actual dm_ZZZ port. The new graph exposes that port instead of silently using only the available lower derivative rows. Low-order small finite-N estimates do not by themselves control this third mean derivative or the retained outer RQy. Quantitative mixed C5 containment and global frequency admission remain open. The next numerical task is to enclose the correlated six-I1 source and its Q_ZZ-dependent defect at the actual selected N, using signed functions and coherent histories; then handle the actual outer boundary. This stage does not raise N or alter the construction family to hide an unresolved radius. This frequency-correction split is not temporal coefficient recursion in order n.

## Execution and scope

New production exits0 in 6.204 seconds. Report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_micro_correlated_m15.json` is 16,101,496 bytes, SHA256 `7868f1899398a54aaeba1b799b8b1015465718445444966648295e66ae6c97aa`, with 1994 input and 346 implementation pins. Focused compilation, new graph arithmetic smoke, six pressure/source identities, signed defect expansion identities, exact24-leaf/measure/binder checks, accepted14-row seam consumption and saved-only inspection pass. A fixed Luna/max worker independently reviewed the pressure reduction, source ports and derivative arithmetic. Its independent synthetic graph evaluation checked dq orders0..2 and all six raw defect values/Z1 against symbolic differentiation exactly (166 nodes,162 appended). Tool responses do not expose separate runtime model metadata; none is claimed.

All inherited source/error packets remain JSON-value identical except stage/hash/contract/check/time metadata. The original six signed inner/lower terms, numerical flat prefix, all previous finite-N errors and the212-term archived forward history are retained. This new representation is an alternative for the same two-chart source, not an extra summand. No accepted owner, builder, fixture, recurrence, norm producer, old quadrature or endpoint fingerprint is replayed. P1/m13 remains distinct from m15.

Whole micro numerical propagation, full forward Rp, independent reserved-versus-negative-pulse terminal matching as axial functions, exact heat/finite-energy exterior, admissible stress/flat remainder, genuine temporal n>=2 recursion/summation, realizable mean/oscillatory correction and independent corrected Cartesian NS remain incomplete. Goal ACTIVE.
