# Original O2-buffer numerical and physical continuation - 2026-10-15

`OriginalOuterBufferNumeric.evaluate(Z,b)` and `OriginalOuterBufferVelocity.evaluate(Z,b,log_tau="-2",azimuth=".25")` now continue the actual original selected-N construction across the complete O2_buffer chart, b in[0,11]. They consume the accepted O2_axial endpoint, preserve its nonzero five histories, independent P0 and n1 primitive J, and export physical cylindrical/Cartesian values. The full reconstruction remains **ACTIVE / INCOMPLETE**; tasks14/15 remain IN PROGRESS.

## Exact leading buffer fields and small sectors

The physical log measure is dy=db. Put A=1/(1+Z^2), B=4Z/Pstar, D=exp40-1, and use the accepted slope endpoint F1/H1/G1/C1 and actual positive axial endpoint M/W. With F0=F1*exp(-D/2) and O=(H1-F1)*exp(-3D/2), the leading fields are

```
E = A*F0*exp(-b/2), V = 0
m = B*M*exp(-b)
h = E+A*O*exp(-3b/2)
k = A*B*(F0*M+O)*exp(-3b/2)
e = B^2*W*exp(-b)
    -A^2*(G1+F1^2*D)*exp(-D)*exp(-b)/2
    -A^2*F1^2*exp(-D)*b*exp(-b)/2
p = A^2*(C1+F1^2)/2-A^2*F1^2*exp(-D)*exp(-b)/2
pressure = P0/Pstar^2+p
```

These solve the original five own-rate density equations m'=-m, h'=E-3h/2, k'=-3k/2, e'=-E^2/2-e and p'=E^2/2 with the actual axial incoming values. The positive M/W, mixed moment and resonant swirl-square memory are never reset when local leading V becomes zero.

`leading_small_sectors` exports each separate signed contribution through TaylorZ2 before adding it to larger terms. It includes h-E, both k sectors, the three e sectors, constant pressure, negative pressure tail and the positive buffer increment

```
Delta p = A^2*F1^2*exp(-D)*(-expm1(-b))/2.
```

The increment remains strictly positive even at b=1e-140. Preserving it separately avoids its disappearance in the rounded large sum; it does not improve the independent incoming scalar uncertainties. The combined P0/Pstar^2 plus pressure constant is also exported, without asserting exact cancellation or fitting the independent P0 to the numerical endpoint.

With V=0 and m=4Z*M*exp(-b)/Pstar, the original shared-streamfunction identity reduces exactly to Q=Ur/sqrt(R/2)=-4*M*exp(-b), independent of Z. Its leading Z derivative is exactly zero. The modified finite-N Q still includes the incoming mass derivative and local V correction, so its Z derivative is not asserted zero.

## Actual finite-N history, geometry and physical field

The exact original source cells are [0,1], [1,5] and [5,11]. The reader merges overlapping axial covers and traverses them by exact rational comparisons, using each actual width db once. It transports the incoming five finite-N error rows with rates1,3/2,3/2,1,0. Each incoming error already contains1/N; each local all-N density driver divides by the unchanged original selected N once. Pressure memory is exactly1.

Local E_N/V_N are current-radius source values rather than additive histories. Exact b=0 reuses the accepted local enclosures; later points use the same-original all-N source bounds. Signed finite-N phase centers remain unresolved, so fields retain a leading point proxy plus directed error. Neither source cap endpoints nor small candidate integers are substituted for original values.

The global log coordinate is y=log(R/110)=10logCstar+10logPstar+exp40+b. The large exp40 and finite b atoms remain separate. Physical assembly includes log110 once, the original lambda/time factors, Cartesian rotation and nonzero n1 meridional history J. Original logPstar=exp40+11 gives log(Rw/Rbuffer_end)=1>0, so the entire buffer remains before the reserved repair. Local n1 axial/swirl zero is source-admitted; no global high-order join is claimed.

## Focused evidence

The focused driver exits0 in50.368s and saved-report inspection exits0. Evidence includes33 actual physical packets across5 axial positions and ordinary/scaled times;45 independent integrals of the raw five-density ODEs;48 shared-streamfunction Q identity conditions;250 numerical and40 physical axial seam conditions;65 original-N/n1/local-inlet conditions;80 physical-measure/pressure-memory checks;120 P0 conditions;86 retained signed-small-sector conditions;4 exact close-boundary classifications;21 scaled-time conditions; and5 invalid requests rejected. No accepted ancestor scientific input is changed. The existing read-only worker `/root/n1_atlas_map` (GPT-5.6 Luna, max) reviewed the separated sectors, streamfunction identity, original source inheritance and open pressure/recursion limits without finding a material algebra error. No new agent was spawned.

Report: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_outer_buffer_numeric.json`;13,948,589 bytes; SHA256 `3183ac356a587248c2875db665baac551c45dd565a0fed6827ba5295f56ef3f9`;1840 scientific input pins and62 implementation pins.

## Next original construction work and limits

1. Continue the actual buffer endpoint through O3_transition x in[0,1] and original O3_power x in[0,2], reaching Rc. Use the original power_flow, not a replacement decay law. Retain the nonzero selected-N five-history inlet even when a later quiet local driver vanishes.
2. Preserve the original mu=exp(-4logPstar)/1000 and its cone-margin contribution analytically or in separate factored sectors. Ordinary exp(-mu*J) enclosure rounding is much larger than mu; a box containing1 cannot certify the needed sign or justify mu=0.
3. Carry the actual signed terminal defects into the reserved repair. Supply same-family S0/S1 and actual total Q_y/RQ_y for reduced m13/m15 pressure, with physical R/R^2 measures. Current Z rows are not total radial or fast-phase derivatives.
4. Complete full n1 pressure and global high-order joins, physical C2/common N, compatible analytic preheat datum/exact heat exterior, physical energy/radial tail, admissible stress and flat remainder. Independent P0 plus numerical p does not yet certify relative cancellation-sensitive pressure accuracy.
5. Construct independent genuine temporal n>=2 coefficients, smooth sum and two-family oscillatory cancellation, then independently validate fixed-forcing corrected Cartesian NS. Four numerical outer charts do not establish temporal scale recursion or a full corrected NS solution.
