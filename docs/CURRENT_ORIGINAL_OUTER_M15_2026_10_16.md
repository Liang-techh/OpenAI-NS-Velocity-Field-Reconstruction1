# Original pressure and m15 through outer reference and angular slope - 2026-10-16

`OriginalOuterPressureRows` and `OriginalOuterM15` continue the original source through Rh_reference, native offset[-5,0], and O2_slope, native coordinate[0,1]. The full168-term Rh pressure inlet becomes172 terms atRref and176 terms at the end of O2_slope. Each chart contributes its actual leading pressure integral and three finite-N error sectors; all earlier terms remain exactly once.

The full goal remains **ACTIVE / INCOMPLETE**; tasks14/15 remain IN PROGRESS. This is cumulative pressure progress through two actual outer intervals, not a global five-moment closure or genuine temporal recursion result. Signed finite-N centers/useful relative cancellation, the changing O2 axial pressure integral, independent reserved Rp value/Z1, all-five/global joins/common N, exact heat/finite energy/stress/flat remainder, temporal n>=2 and corrected NS remain open.

## Exact physical source and pressure

The pinned true_outer_C6_C5 source and smooth Rh join give V=m=4Z/Pstar on both charts. The original angular fields change in O2_slope, while physical U=M=4Z remains unchanged. Therefore My=0, Uy=0 and leading Qy=0. The leading mean third derivative is genuinely zero by this analytic identity; it is not a missing numerical derivative padded with zero.

With A=1-delta*Z^2 and D=1-Z^2,

```
Q = [4(1+delta)Z^2 - 4D]/A.
```

The new pressure adapter evaluates its actual rational axial derivatives through Q_ZZ and exports the mean through Z3 and Uy/Qy through Z1. Angular effects remain in all original history/source and finite-N loop bounds. They are not a missing term in the leading reduced half-S1 pressure identity.

Both native coordinates have unit physical log-radius Jacobian. The radius is R=110*C^10*Pstar^10*exp(offset): offset=-5 atRh,0 atRref and1 at the end of O2_slope. Unit offsets remain separate atomic log-radius terms; the huge radius is never materialized.

## Actual four-row finite-N transport

The Rh inlet reads all four ordinary derivative rows of all five genuine patch histories. O2_slope consumes the complete four-row Rh_reference endpoint. Neither chart uses an old public Z2 history as a substitute. The original rates1,1.5,1.5,1,0 and the actual native cell widths remain. Coincident radial cells across overlapping axial covers are merged, new drivers are divided by N once and already-scaled incoming histories are not divided again.

Derivative bounds bind exactly to the matching `Rh_reference` and `O2_slope` loop keys. Physical deltaU=Pstar*B/N and deltaM=Pstar*history_m feed the pressure error recovery. Slow By/ByZ carry1/N; fast Bphi/BZphi do not. No extra x, hb or epsilon occurs. C3 mean errors supply the required Q_ZZ error, and the actual slow/fast derivative rows supply Qy_Z.

At a shared endpoint the new chart's own local B/By/Bphi bounds remain active. The adapter unions the predecessor U/Uy error budget with those bounds; it does not identify two different conservative chart caps as values of the same function. The lower endpoint m15 error uses this chart-aware union as well. A read-only review caught the need for that lower endpoint union before the production run.

For volume errors, each chart stores a full-prefix envelope by retaining incoming magnitude and all positive source masses without later decay. Local bounds cover every visited cell, plus the incoming endpoint budget. The background is constant over the prefix, so its exact leading enclosure needs no radial approximation. Both R^2F and RQy endpoint errors and all quadratic/C1 product terms remain. Leading Qy=0 never resets the finite-N error.

## Exact pressure integral

For each chart let L be its physical log-radius distance from its lower endpoint. The original F/G coefficients are constant in radius, so

```
delta m15_leading = Rlower^2 * (2F+G) * expm1(2L)/2.
```

This combines the signed F endpoint identity and G integral before subtracting close radii. It uses the existing half-S1 normalization once. The leading signed RQy endpoints are exactly zero; their finite-N errors remain in the viscous error sector. The exact positive volume mass is Rlower^2*expm1(2L)/2. There is no numerical radial quadrature or new free parameter in this leading integral.

Each chart appends four terms: one signed background integral and the finite-N endpoint, volume and viscous errors. Empty intervals have four exact zero terms. Source coefficient/error enclosures and the shared source-function binding remain separate from physical scale atoms.

## Observed execution

The original-source run exited0 in93.876 seconds; saved-report inspection exited0. It evaluates Z=-.3,0,.3 at Rh_reference offsets-5,-4.5,-2.5,0 and O2_slope coordinates0,.25,.5,1.

- 24 moment packets,24 pressure packets and24 exact cumulative-inlet replays;168->172->176 terms.
- 192 serialized new coefficient/error rows, including48 exact empty rows;336 ordinary pressure/mean rows and480 genuine C3 history rows.
- 360 full-prefix error conditions,72 lower finite-N endpoint conditions and24 positive native kernel/width conditions.
- Six exact four-row history seams,96 analytic Uy/Qy zero rows,192 independently saved lower-pressure overlaps and48 current-chart fast-phase unit checks.
- Independent fixtures execute36 original leading-source coefficients, compare18 physical C0/C1 integrals and check20 C3 Duhamel rows against direct integration with finite N=17. Five invalid original-domain queries are rejected.

The47,559,524-byte report SHA256 is `b8dd2f642c4c55ac183a087ed712e6a95a46aaaa7828cb4c70326b44deda8bc8`. It pins1791 inputs and88 implementations. Inspect with `python scripts/run_original_outer_m15.py --inspect-report`; do not rerun accepted ancestors.

Existing read-only `/root/n1_atlas_map`: **GPT-5.6 Luna / max**, reused without spawning another agent. It checked the exact source/loop chart mapping, leading angular independence in half-S1 and local endpoint error scope. The independent rate-zero fixture also caught that interval zero must be compared with r==0 rather than interpreted as a Python boolean. The recovery checkpoint retains terminal attempts, the disappeared earlier run handle and the successful final handle77431. The final run has persistent output in `tmp/original_outer_m15_2026_10_16.log`.

## Next actual source path

1. Build the O2_axial pressure adapter from the complete176-term O2_slope endpoint and all four finite-N histories. In `original_outer_axial_numeric.py`, physical U=4Z*beta(s), M=4Z*Mmemory(s), with beta(s)=sigma(1-s). Mmemory is the genuine positive own convolution from `OriginalAxialMemory`; beta(1)=0 does not imply Mmemory(1)=0.
2. Recover Q_ZZ and Qy_Z from these actual scalar beta/M sources and the analytic Z dependence. Keep beta_s divided by40*exp(40s) for Uy. Its physical log distance is d=expm1(40s); the exact transport cell width is exp(40left)*expm1(40(end-left)), not end-left. Bind the actual O2_axial slow/fast loop key and preserve chart-aware lower endpoint unions.
3. Integrate the changing-source pressure with the original identity [R^2F]+integralR^2G dy-[RQy]. Separate beta/M linear and quadratic scalar modes before introducing large radius atoms. Use source-aware positive/normalized kernels across the exponentially stretched coordinate and retain nonzero terminal M and both RQy endpoints. A uniform native grid alone does not resolve the concentration scales or establish a useful quadrature error.
4. Continue the pressure state through O2_buffer, O3_transition/O3_power and the true compact/quiet intervals toRp. Preserve all finite incoming histories and signed/source/quadrature errors. Then compare source-correlated cumulative value and Z1 with the independent reserved pressure moment.
5. Complete global all-five/common N, flatten/preheat/exact heat/finite energy, admissible stress/flat remainder, genuine temporal n-dependent recursion/smooth summation, realizable mean/oscillatory stress cancellation and independent corrected Cartesian NS.
