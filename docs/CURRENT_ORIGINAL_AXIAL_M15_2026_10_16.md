# Original changing axial pressure and m15 - 2026-10-16

`OriginalAxialPressureRows`, `OriginalAxialRateMemory` and `OriginalAxialM15` now evaluate the original O2_axial turn-off, including its nonzero own mean and viscous endpoint. The complete176-term O2_slope inlet becomes195 terms, with all previous sectors retained once. This is an actual changing-source pressure integral, not a replay of constant outer fields.

Tasks14/15 remain IN PROGRESS and the full goal remains **ACTIVE / INCOMPLETE**. Signed selected-N centers and useful source-correlated cancellation, continuation toRp and independent reserved pressure, global all-five/common N/joins, preheat/exact heat/finite energy/stress/flat remainder, genuine temporal n>=2 and corrected NS remain open.

## Exact original source and nonzero endpoint

Let s in[0,1], d=expm1(40s), beta=sigma(1-s), and let m be the original positive own mean. It satisfies m_d=beta-m, m(0)=1. The physical fields are

```
U=4Z beta, M=4Z m,
A=1-delta Z^2,
Q=8Z^2 beta/A-4m,
Qy=8Z^2 beta_y/A-4(beta-m).
```

The cutoff derivative is beta_y=-sigma'(1-s)/(40exp(40s)). The adapter supplies U/Q through Z2, M through Z3 and Uy/Qy through Z1. M_ZZZ=0 follows analytically from M=4Z*m; finite-N C3 mean errors remain. At s=1, beta=beta_y=0 but m>0, so Qy=4m>0. The terminal original log(m) is approximately-8426695757850.9881764476. Its small magnitude never authorizes zeroing its function or the signed RQy boundary.

The physical radius is represented by separate atoms ln110,10lnC,10lnPstar,1,d. The physical measure is dd=40exp(40s)ds. Each transport cell uses exp(40left)*expm1(40(end-left)) exactly once; native s widths do not replace physical widths.

All five incoming histories retain four ordinary axial derivative rows. New sources carry1/N once; incoming histories are not divided again. The actual O2_axial loop supplies slow By/N and fast Bphi, without an additional physical Jacobian, hb or epsilon. At the lower seam and for prefix bounds, predecessor and current-chart local budgets are unioned. Endpoint and full-prefix error bounds are not selected signed finite-N values.

## Changing-source integral and normalization

The original half-S1 coefficients give

```
delta m15 = [R^2 F] + integral R^2 G dy - [R Qy].
F = 4Z^2 beta/A^2 -16Z^2 beta*m/A +4m^2 -2m/A,
G = 16Z^2 P beta^2/A^3 +16Z^2 beta*m/A
    -4m^2 -4delta Z^2(1-Z^2) beta/A^3 +2m/A,
P = (delta^2+delta)Z^4 -(2delta+3)Z^2 +3.
```

F/G already include the factor1/2 of S1. Do not halve this identity again. The delta-sized G_beta coefficient is evaluated in its simplified factored form; subtracting its O(1) precursor terms would destroy that source sector at original parameters.

For I_f=integral_0^D exp(2d)f(d)dd, the mean equation gives exact identities

```
I_beta*m=(exp(2D)m(D)^2-1)/2,
I_m=exp(2D)m(D)-1-I_beta.
```

The implementation stores the signed pieces as separate physical scale atoms, preserving the shared source-function bindings. There are19 new sectors:5 F-boundary,8 G-volume,3 upper -RQy components and3 finite-N error sectors. The lower leading Qy is exactly zero; its finite-N error remains. Every empty interval has19 exact zero increments.

Finite-N endpoint/quadratic/viscous errors use both actual endpoints. The volume error uses every visited cell and all incoming history magnitude without later decay. Its positive mass is factored as Rupper^2*(1-exp(-2D))/2; no huge physical radius is materialized.

## Positive source kernels and squared mean

`OriginalAxialRateMemory` extends the unchanged pinned rate-one implementation with rates1 and2. Ordinary cells integrate the exponential kernel times a beta Taylor3 polynomial with a fourth-derivative remainder. Near the terminal concentration, peak normalization, Taylor2/fourth remainders and positive tangent tails resolve the true source. Kernel homogeneous and forced contributions remain separate, so no nearly equal homogeneous subtraction defines I_beta or I_beta2.

For q=D-d, the exact identity

```
exp(-q)m(D-q)=exp(-D)+integral_q^D exp(-r)beta(D-r)dr
```

turns the squared mean into the square of a positive future source CDF. Ordinary cells use exact exponential-polynomial square primitives plus a source Taylor error bound. Peak cells use monotonicity between their certified positive endpoint masses. Both source tails, the homogeneous term and finite endpoint mass/length rounding remain. This avoids nested mean quadrature and a uniform s grid. The observed terminal normalized mean-square interval has relative width approximately1.7e-9; this is a kernel enclosure, not a useful relative error claim for the total signed m15.

## Observed execution

The original-source driver exited0 in114.834 seconds; saved-report inspection exited0. The samples use Z=-.3,0,.3 and s=0,.001,.5,.99,.99999,1.

- 18 moment and18 pressure packets;18 exact176-term inlet replays;195 cumulative terms each.
- 684 serialized new coefficient/error rows,114 exact empty rows;252 pressure rows and360 C3 history rows.
- 270 prefix error conditions,60 lower endpoint conditions,18 physical-width conditions and45 positive-kernel conditions.
- 72 overlaps with the independent saved original axial velocity archive,3 nonzero terminal Qy conditions and3 rejected invalid queries.
-Independent fixtures check60 F/G ordinary coefficients,12 changing-source physical C0/C1 integrals, one positive cubic-source nested CDF-square integral and20 physical-width C3 Duhamel rows with finite N=17.

The38,793,438-byte report `lei_ren_part1_paper_compliant_current_original_axial_m15_numeric.json` has SHA256 `57f5495394ab864ace99905fac07aab5bd315ff420329704121de121312b8d4a` and pins1793 inputs/94 implementations.

`scripts/verify_original_axial_kernels.py` separately reads this saved report, without reconstructing ancestors. Direct140-digit tanh-sinh integration of beta and beta^2 at s=.5 and1 lies inside all four saved rate-two source enclosures. It exited0 in1.922 seconds. The check receipt SHA256 is `42fee348f89b3290e2e7f76a8a15e09847684df46a11e7ca5192239ba5f6af36`. Observed normalized kernel widths are approximately8.5e-18 in the ordinary region and1.73e-11 at the peak.

Reproduce with `python scripts/run_original_axial_m15.py --write-report`; normally inspect with `--inspect-report`, then run `python scripts/verify_original_axial_kernels.py`. Do not rerun accepted ancestors.

Existing read-only `/root/n1_atlas_map`: **GPT-5.6 Luna / max**, reused without spawning another agent. It reviewed source equations, units, F/G algebra, exact mean reductions and positive CDF bounds. Root retained responsibility for normalization: a subsequent worker buffer note proposed an extra1/2, which was rejected because F/G already represent half-S1. The accepted driver independently checks the existing normalization.

The durable production log is `tmp/original_axial_m15_2026_10_16.log`, session39048 terminal0. Initial fixture-only failures were repaired before production: the symbolic quadrature needed an mpmath printer to avoid Python float rational constants, and the synthetic transport object needed its chart name. No accepted ancestor was changed.

## Next actual pressure continuation

1. Build O2_buffer pressure from the complete195-term axial endpoint and all four histories. The original source `original_outer_buffer_numeric.py:_leading` gives U=0, M=4Z*mstar*exp(-b), Q=-4mstar*exp(-b), Qy=4mstar*exp(-b), b in[0,11], with physical dy=db. Keep mstar from the accepted positive original endpoint, not a zero surrogate.
2. The leading RQy endpoints are nonzero and equal. The half-S1 identity simplifies exactly to delta m15=-4*Rlower^2*mstar^2*b, independent ofZ. Do not apply another1/2. Preserve all finite-N endpoint/prefix/quadratic errors even though the leading viscous difference cancels.
3. Bind the actual O2_buffer source/loop key and retain local B/N, slow By/N and fast Bphi units. New-chart inlet caps must be unioned with incoming caps. Carry atomic log radii through the11-unit buffer without collapsing order-one offsets.
4. Continue through original O3_transition/O3_power and true compact/quiet intervals toRp. Compare the source-correlated signed cumulative value and Z1 with the independent reserved moment. A large absolute error budget or count of stored sectors is not terminal closure.
5. Finish global all-five/common N/smooth joins, flatten/preheat/exact heat/finite energy, admissible stress/flat remainder, true n-dependent temporal recursion/summation, mean/oscillatory stress cancellation and independent corrected Cartesian NS.
