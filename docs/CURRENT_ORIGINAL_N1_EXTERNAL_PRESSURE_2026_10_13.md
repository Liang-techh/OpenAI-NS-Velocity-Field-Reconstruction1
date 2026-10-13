# Complete original external n1 pressure values — 2026-10-13

The two original leading-pressure forcing integrals now cover the whole physical interval from the actual Rin to infinity, continuously for every real Z in [-1,1]. The 512 checked inner-collar Z cells have been combined with the 256 accepted r_minus-to-infinity models. There are 512 common interval-affine models, with maximum absolute errors in the original fixed positive pressure units of 1.7214685424e-5 for m13 and 2.3006530659e-5 for m15. The positive lower scale ratio, signed lower coefficients, collar comparison errors, both Taylor remainders and nonzero outer source errors are retained.

This completes the external pressure-value dependency. It does not include the n1 inner cutoff pressure moments and does not close the five terminal defects, original 2+3 repair, full pressure derivative accuracy or true temporal recursion. The full original goal remains **ACTIVE / INCOMPLETE**.

## Source and units

The lower provider is `current_original_n1_unchanged_axial_pressure_models`: the accepted actual Rin-to-r_minus pressure prefix with 512 continuous Z cells. Its source basis is [hb,F0_anchor,Pstar_squared,Lambda], and its reference exponents are [1,0,1,-1] for m13 and [1,0,1,-2] for m15. Thus its positive physical units are hb*Pstar^2*Lambda^(-p)*M_p, with p=1,2 and the archived directed positive coefficient norm M_p. No F0 normalization enters these reference units.

The upper provider is the accepted [r_minus-to-infinity pressure](CURRENT_ORIGINAL_N1_REMAINING_PRESSURE_2026_10_13.md). Its fixed units are Phi*A*Rv^p, where Phi=10.99*sqrt(2*pi)*(exp(5)-exp(3)), A=Pstar*U0*exp(logE-13/(2mu)-13), Rv=Rp*exp(13/mu), Rp=Rc*exp(-2-60*log(mu)), and logE=-1/mu-3/(4mu)^(2/3)+(2/3)*log(4mu)-log(6)/2-2*log(mu).

The new provider computes the exact positive unit-ratio logarithm before taking any enormous physical exponential. With lmu=log(mu), log(mu)=-4*logPstar-log(1000), and p=1,2:

```text
log(lower_unit/upper_unit)
 = loghb + logPstar - p*(logLambda+logRc) - logU0
   + log(M_p) - log(Phi) + 13 + 2*p
   + (60*p+4/3)*lmu + (15/2-13*p)/mu
   + 3*exp(-(2/3)*(log(4)+lmu))
   - (2/3)*log(4) + log(6)/2.
```

The inverse-mu coefficient includes the leading -1/mu term of logE. It is -11/2 for p=1 and -37/2 for p=2. The original source interval proves both ratio logs below -1000. The actual ratios remain strictly positive mathematical functions, enclosed numerically by [0,exp(-1000)]; zero is not selected as a ratio or a signed physical field. Multiplication by each archived signed lower coefficient retains its sign enclosure and all original uncertainties.

Both segments use the same pinned source family. The lower owner defines Rin=sc/4 and r_minus=sc/2; the accepted original upper geometry begins at the same r_minus. The physical radial measure and minus sign for m15 were already included by the lower source and are not applied again. Upper active charts, compact band, quiet interval and whole pulse are included once, with exact meridional support beyond Rv. The isolated older end-pulse result is not added separately.

The outer grid is split to the exact lower grid. Each parent intercept is shifted to the common cell center, its slope is retained, and the lower affine coefficients are multiplied by their positive uncertain ratio and added. Lower Taylor errors and inherited outer Taylor/source errors stay separate. These interval models enclose actual values; their midpoints are not chosen as the functions used by the repair equations. Full pressure derivative approximation accuracy remains false. These errors are pressure source units, not the final physical NS residual units.

## Reproduction and acceptance

Run `python scripts/run_original_n1_external_pressure_values.py` from the repository root. It checks shared archive dependency bytes once, composes existing models and runs only the new focused checker. It does not reconstruct parent field owners, recompute the 68 active-chart packets or rerun their mathematical suites.

The public provider is `CurrentOriginalN1ExternalPressureValues()`, with `evaluate(Z)` and `report(packet)`. It supports intervals as well as points in [-1,1]. Its source quartet is `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_external_pressure_values` with `.py`, `.json.gz`, `_check.py`, `_check.json`. The installed gate is only `current_original_n1_Rin_to_infinity_pressure_values_installed`.

Completed work:

- [x] Identify actual lower source exponents and positive directed coefficient norms.
- [x] Bind the original common Rin/r_minus, selected Cstar, Pstar, Lambda and source U0.
- [x] Derive the correlated source-unit conversion, including logE's inverse-mu term.
- [x] Retain the strictly positive tiny ratio as an enclosure, never a selected zero.
- [x] Compose all 512 lower and 256 upper models into 512 exact common cells.
- [x] Preserve signed source values, both remainders, coefficient uncertainty and source errors.
- [x] Check the scale formula symbolically and with independent finite nonunit numerical scales.
- [x] Independently recompose every model and check full coverage, boundaries, seams and interior queries.
- [x] Confirm the accepted public interface works at default mp.dps=15 and restores precision.

## Next tasks: complete n1 defects and repair

Use the existing cutoff and cumulative owners before inventing another inner reconstruction. `CurrentOriginalN1RadialCutoffExtension(checked_cumulative_owner)` already evaluates the actual cut family. Query `evaluate(Z,"1","first",cells=4)` and `report(packet)`, then read `physical_cutoff_moment_values_at_Rin`. Its m13/m15 are signed factored jets retaining the original source exponents and ordinary Z derivatives. The same output includes the three full velocity-support moments m11/m12/m14. The archived inner-exit sample is evidence at Z=.371 only; it does not by itself prove a uniform full-Z bound.

- [ ] Construct the checked live cutoff/cumulative owners and bind their family and independent P0 datum to the new external pressure provider.
- [ ] Evaluate the cutoff's actual inner-exit moments on every required full axial cell, preserving the original source logs, signed sectors, own M/V/P and cutoff-before-recovery rule.
- [ ] Define m13_full=m13_cut(Rin)+the new unweighted external integral.
- [ ] Define m15_full=m15_cut(Rin)+the new already-negated R-weighted external integral. Do not negate the new m15 result twice.
- [ ] Preserve distinct signed source sectors while converting the cutoff pressure moments to common units. Prove any negligible contribution with a nonzero error rather than selecting zero.
- [ ] Compose the resulting two complete pressure defect functions with the original three velocity moments in one consistent n1 five-defect interface.
- [ ] Retain value enclosures separately from genuine axial derivatives required by the repair system. The current affine slopes alone do not certify the full pressure derivatives.
- [ ] Resolve the full pressure derivative source/mean bounds and any wide q_ZZ/C3 dependency before using them in repair coefficient formulas.
- [ ] Bind the original positive I2/Astar controls, normalized bumps and exact two-plus-three repair equations to these full defects.
- [ ] Implement the first repair pair, then the remaining three coefficients, without replacing functions by arbitrary interval midpoints.
- [ ] Apply repair to the correct F1/Uz1 family and recover its own M/V/P, keeping the independent analytic pressure datum.
- [ ] Check all five terminal identities as functions on the original axial domain; sampled near-zero values are not closure.
- [ ] Check axis regularity, full support and high-order joins with the original annuli, flatten region and exact heat exterior.
- [ ] Install the n1 repaired physical layer only after its functional identities and own-field recovery are established.
- [ ] Implement n2 and later orders with the original n-dependent equations and their independent moment repairs.
- [ ] Establish finite-order remainder bounds and smooth summation; distinguish these from coordinate scaling.
- [ ] Restore admissible stress-cone margins and the flat remainder across all physical regions.
- [ ] Implement both oscillatory pulse families and their mean corrections; check averaged quadratic stress cancellation.
- [ ] Independently validate the complete corrected Cartesian NS residual against a force restricted in advance.
- [ ] Mark each completed task `[x]`, record its source/receipt/commit and limits, and update the five entry documents.

Full goal **ACTIVE / INCOMPLETE**. Whole external pressure values are available; inner-plus-outer n1 defect composition and repair are next.
