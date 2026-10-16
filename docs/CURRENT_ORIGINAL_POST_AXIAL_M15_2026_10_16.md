# Original buffer/O3 pressure and m15 to Rc - 2026-10-16

`OriginalPostAxialPressureRows` and `OriginalPostAxialM15` now continue the actual original pressure through O2_buffer[0,11], O3_transition[0,1] and O3_power[0,2]. The complete195-term axial inlet becomes199,203 and207 cumulative sectors at the successive endpoints, reachingRc. Every inherited sector remains once, with the same selected integer and source family.

Tasks14/15 remain IN PROGRESS and the full goal is **ACTIVE / INCOMPLETE**. This is cumulative pressure progress through three actual outer intervals. It is not global five-moment closure, an independently resolved Rp seam, admissible-stress completion or temporal coefficient recursion. Signed finite-N centers/source correlations and useful total relative accuracy remain unresolved.

## Exact shared original source

Let R0 be the actual axial-end radius and mstar its positive scalar own mean, evaluated in the preceding accepted source. Let tau=b on buffer,11+x on transition,12+x on power. Then

```
R=R0 exp(tau), m=mstar exp(-tau),
U=0, M=4Z m, Q=-4m, Qy=4m.
```

`original_outer_buffer_numeric.py:_leading` gives these identities. `original_outer_o3_numeric.py:_leading` calls that same background at offsets11+x and12+x. Its actual nonzero mu corrections modify E,h,e,p; they do not change leading V,m or Q. The new adapter retains those actual original angular sources in all chart-specific finite-N histories and derivative bounds. No mu=0 substitution is used.

The pressure adapter restores the saved positive original scalar mean and its exact directed enclosure, with its source proof and implementation pins. Because this scalar is independent ofZ, all saved axial endpoint values must agree. It is not replaced by an arbitrary constant or byzero. U/Q/mean/Uy/Qy derivatives are recovered analytically, including genuine mean Z3 zero and Q_Z=Q_ZZ=Qy_Z=0. Their finite-N derivative errors remain separate.

## Half-S1 integral and common radius anchor

With A=1-delta Z^2, the existing half-S1 identity is

```
delta m15=[R^2 F]+integral R^2 G dy-[R Qy],
F=-2m/A+4m^2, G=2m/A-4m^2.
```

F/G already include the factor1/2 of S1. The linear terms cancel as the same source functions. Both leading RQy endpoints equal4R0*mstar and are individually nonzero. Therefore, on a chart of local length L,

```
delta m15_leading=-4 R0^2 mstar^2 L.
```

The implementation pairs the **global axial radius** with the **global axial mean** on all three charts. Equivalently, a local lower radius must be paired with its local mean mstar*exp(-tau_start). Mixing the local O3 radius with the global mean would introduce an erroneous exp(22) or exp(24) factor. The shared atomic source binding prevents that mismatch. The leading increment has exact Z1 zero.

Each chart adds four sectors: the exact signed leading integral and finite-N endpoint, volume and viscous errors. Empty intervals have four exact zero increments. Leading RQy cancellation does not erase finite-N endpoint errors or authorize treating magnitude bounds as selected signed values.

## Actual C3 histories and physical units

All three native coordinates are physical log-radius offsets: dy=db or dx once. Source covers retain four ordinary Z derivative rows for all five histories, using their original rates. Overlapping axial covers do not duplicate radial cells. Incoming histories are already scaled by1/N; each new driver is divided once. The prefix envelope retains all incoming/source magnitude without subsequent decay.

History chart names are O2_buffer,O3_transition,O3_power. The derivative-loop table uses O2_buffer,**O3_slope_mu**,O3_power respectively. The explicit map follows `current_generic_shear_O3_sources.py:SPEC`: O3_slope_mu binds the actual transition tensor. O3_power may use t/Tw as its coverage coordinate, while its stored y rows are already ordinary physical logR derivatives. No extra Tw, native-coordinate Jacobian, hb or epsilon belongs on those stored rows.

Physical local deltaU is Pstar*B/N; the mean error is Pstar*history_m. Slow By/ByZ carry1/N, fast Bphi/BZphi do not. New-chart and predecessor U/Uy budgets are unioned at lower endpoints and for full-prefix bounds. QZZ uses the genuine mean C3 error. Both R^2F and RQy finite-N endpoints remain. The volume mass is factored as Rupper^2*(1-exp(-2L))/2; no huge physical radius is materialized.

## Observed execution

The final original-source run exited0 in82.265 seconds; saved-report inspection exited0. It uses Z=-.3,0,.3; buffer coordinates0,5.5,11; transition0,1; power0,2.

- 21 moment and21 pressure packets,21 exact inlet replays; cumulative endpoints199,203,207.
- 168 serialized new coefficient/error rows,72 exact empty rows;294 pressure rows and420 C3 history rows.
- Nine exact four-row history seams,315 prefix-error conditions,48 lower finite-N endpoint conditions and21 physical-width conditions.
- 144 saved original buffer/O3 velocity overlaps,21 normalized nonzero RQy identities and21 exact leading Z1 zero conditions.
- 42 actual fast-phase unit checks and five rejected invalid queries.
- Independent finite fixtures execute216 original leading-source rows and18 genuinely nonzero-mu O3 cases. Eighteen physical C0/C1 integrals and an exact half-S1 density identity verify the normalization and common radius/mean anchor.

The48,613,025-byte report `lei_ren_part1_paper_compliant_current_original_post_axial_m15_numeric.json` has SHA256 `48272c642ed8a28248aab35e7feda450f6a10aeff534f77b29f4b56a07f913f4`. It pins1796 inputs and112 implementations.

Normally inspect with `python scripts/run_original_post_axial_m15.py --inspect-report`. Reproduce only when needed with `--write-report`; accepted ancestors need no replay. The durable production log is `tmp/original_post_axial_m15_2026_10_16.log`, final session61298 terminal0. Earlier terminal attempts were repaired:98823 exposed the transition loop-name mismatch;1713 exposed a validation helper that tried to decode an exact-zero budget as a finite log. Exact-zero log comparisons now use minus infinity. The physical source implementation was unchanged by that second repair.

The existing read-only `/root/n1_atlas_map` is **GPT-5.6 Luna / max**, reused without spawning an agent. It checked the source path and angular/mu independence. Root retained mathematical decisions, including the distinction between global and local mean/radius factors and the already included half-S1 normalization.

## Next actual source work

1. Continue the full207-term Rc state through `original_compact_numeric.py`. Its native s=log(R/Rc) ends atlog2; leading U=0, M=4Z*mstar*exp(-14-s), Q=-4mstar*exp(-14-s), Qy=4mstar*exp(-14-s). The exact leading increment remains -4R0^2*mstar^2*s. Preserve actual signed source-function bindings and all four finite-N histories across the three supported compact bumps.
2. Recover compact Qy and Qy_Z errors from the actual physical-y local velocity/mean source. The older public packet's Q values or C3 magnitude caps alone do not provide those rows. Use the original `actual_selected_N_profile_y0_y1_y2_C3_bounds`, signed-density magnitude ranges and incoming correction C3 bounds with the correct physical units. These are enclosures, not selected signed centers. Preserve support/radius factors, actual incoming mean and nonlinear/quadratic errors.
3. Continue `original_quiet_numeric.py` toRp with s=qL, L=-60log(mu)-2-log2. Native q width is not physical y width: use Ldq once. The leading mean is mstar*exp(-14-log2-s), so the global pressure primitive is still -4R0^2*mstar^2*s. Inspect whether the **full** quiet U=0/Qy=-Q identity also permits source-correlated cancellation of finite-N endpoint terms; do not assume it merely from leading values.
4. Build an exact source-component crosswalk against `original_preRp_m15.py` and `original_preRp_signed_inner_m15.py`. Their bound terminal identity is inner+lower+active+compact+own_swirl+quiet+pulse=0. The new207-term inlet already covers core/lower/active/outer throughRc; do not append it wholesale to legacy inner/lower/active components. Audit compact, own_swirl, quiet and pulse functions individually before replacing or retaining them.
5. Compare independently defined reserved response and cumulative value/Z1 using the same source identities. `reserved_n1_pressure` is a P1 packet and is not a replacement for the full m15 inlet or its independent target. An overlap of broad intervals is not a resolved signed seam. Prioritize selected finite-N source correlations and actual cancellation when the full forward integral is available.
6. Complete global all-five/common N/high-order joins, flatten/preheat/exact heat/finite energy, admissible stress/flat remainder, genuine n-dependent temporal recursion/summation, realizable oscillatory stress cancellation and independent corrected Cartesian NS.
