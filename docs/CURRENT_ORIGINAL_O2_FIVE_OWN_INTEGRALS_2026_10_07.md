# Original O2 five own-rate integrals and explicit incoming history

> Successor: [ordered original source masses and nonzero-Z coefficient cells](CURRENT_ORIGINAL_O2_ORDERED_SOURCE_CELLS_2026_10_07.md). General defining masses and the full source coefficient bridge are completed; nonmidplane phase integrals and functional/global matching remain open.

Checked source: [622563a2](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/622563a264afbcede1aff30177e55de0b6472883). Predecessor: [signed densities and genuine pressure integral](CURRENT_ORIGINAL_O2_SIGNED_DENSITIES_PRESSURE_INTEGRAL_2026_10_07.md). The focused checker passes; the direct Git-index dependency audit covers 1074 hashes. Existing read-only reviewer: **GPT-5.6 Luna / max**. No new workers or ancestor producers were run.

All five genuine original O2 contributions over the complete y in[0,1] window at exact Z=0 and diagnostic N=7 are now enclosed. The calculation keeps both original primitives, the axial increment B/N, the signed cross terms, and each moment's prescribed Duhamel rate. It also transports the original nonzero histories and exposes an affine interface for five explicitly supplied incoming defects.

This closes the previous midplane-all-five task. It does not establish terminal identities as functions of Z, complete all-chart histories, admit a global N, install five controls or stress recovery, or implement coefficient recursion. A current-window contribution initialized at zero is an integral definition; it is not an assumption that the physical preceding-chart defect vanishes.

## Executable interfaces and archived evidence

Production: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_O2_five_own_integrals.py`. Focused checker: the same stem with `_check.py`. The complete producer evidence ships directly as `.json.gz` (about 41.0 MB); its exact decompressed JSON digest and compressed digest are in `_check.json`. The archive retains every density hull, positive mass, signed cell contribution, and quarter-window cumulative history. Source cells point into the accepted pressure archive, whose hash/decompressed digest and full dependency closure are bound.

```python
import gzip, json
import lei_ren_part1_paper_compliant_current_original_O2_five_own_integrals as five

owner = five.OriginalO2FiveOwnIntegrals()
saved = owner.saved['actual_original_pressure_integral_refinements'][-1]
contribution = owner.integrate_cached_level(saved, bits=24)
```

The constructor reattaches accepted source receipts and restores the checked whole-cell source archive. It does not rerun original nested defining quadratures or the pressure-cell producer. To inspect the shipped result without recomputation:

```python
report = json.loads(gzip.decompress((five.HERE / five.NAME).read_bytes()))
record = report['actual_original_five_own_integral_refinements'][-1]
```

`apply_incoming(ctx, record, incoming, source_family=..., original_P0_datum_sha256=...)` requires all five explicit, finite incoming normalized defect enclosures in the same units below. It rejects a missing row, changed family, changed datum, or nonfinite inlet. The interval context should use at least the recorded source precision (260 decimal digits) when transporting the exact saved endpoint bits. No zero inlet is supplied automatically. The return is C0 at exact Z=0 and cannot be passed off as C1/Z-functional input to `GenericMomentRecovery.advance` or its stress recovery.

## Own rates, units and memory

Let S=Pstar, E=Utheta/S, V=Uz/S, e_vel=E*expm1(A/N), and u=B_over_Pstar/N. The accepted signed density graph still executes all five rows, including the quadratic increments:

```text
delta_m = u
delta_h = e_vel
delta_k = V*e_vel + E*u + e_vel*u
delta_e = 2*V*u + u² - E*e_vel - e_vel²/2
delta_p = E*e_vel + e_vel²/2.
```

| Row | Normalized own history | Decay rate |
|---|---|---:|
| m | Mz/(R*S) | 1 |
| h | Mtheta/(sqrt(2)*R^1.5*S) | 3/2 |
| k | Mtheta_z/(sqrt(2)*R^1.5*S²) | 3/2 |
| e | Mztheta/(R*S²) | 1 |
| p | Mp/S² | 0 |

The coordinate is y=log(R/Rref). The original radius normalization is already in these density definitions; no extra R/Jacobian is applied. Physical Pstar/R powers remain formal, so the table of contributions below is normalized, not a physical full NS residual.

For a cell [left,right] of width dy and rate r, the positive final-endpoint mass is

```text
mass_final = exp(-r*(1-right))*dy*exp_average(-r*dy).
```

At r=0 it is exactly dy. The forward prefix recurrence instead uses the local mass `dy*exp_average(-r*dy)` and prior-history decay `exp(-r*dy)` once. Both forms are computed; their enclosures overlap at the terminal point. Multiplying a signed density hull by this positive interval uses all interval corner products.

Original midplane inlet histories are bound to the same original source recipe and are `(m,h,k,e,p)=(0,5/8,0,-5/12,5/2)`. Original V=0, but B and the changed axial velocity generally are nonzero. The original nonzero h/e/p inlet contributions are included in the saved original histories. P0 remains the separately identified analytic datum and is not added into or substituted for cumulative p.

For incoming defect D_j, the actual affine rule is

```text
own_out_j = original_out_j + current_window_contribution_j + exp(-rate_j)*D_j.
```

Pressure rate0 therefore keeps all incoming pressure memory. The other four rows keep their own exponential factor. The checker uses five nonzero diagnostic inlet intervals to check this interface; those diagnostic values are not claimed as the actual preceding-chart defects, which remain to be supplied by the all-chart computation.

## Directed contributions and independent computation

The accepted whole-source f/a/q/positive-eta and true common-N phase unions are restored as function ranges, never selected as field points. The exact Z0 symmetry is bound to the original p2/V/b/t0 templates and pressure parity. It does not zero p2_Z. Both exact original r0 primitives are evaluated with the conditioned inverse:

```text
A = a*q²*sin(2*psi)/(4*pi*(1+2*q²))
B_over_Pstar = -a*E*q*sin(psi)/(2*pi).
```

The graph then returns entire-cell signed density hulls. Refinements at 64, 256, 2048 and 8192 cells contract for all five rows, with common intersections. The complete producer took about212 seconds locally, including source attachment; the focused checker took about9 seconds. The table rounds endpoints; exact directed MPF bits remain in the receipt.

| Row | 8192-cell lower contribution | Upper contribution | Proved sign | Independent ODE reference |
|---|---:|---:|---|---:|
| m | -3.62606e-5 | -1.44459e-5 | negative | -2.53543086e-5 |
| h | -5.33746e-6 | -2.75584e-7 | negative | -2.80650091e-6 |
| k | -1.63915e-5 | 2.36949e-7 | unresolved | -8.07812734e-6 |
| e | 2.64586e-5 | 3.39911e-5 | positive | 3.02246762e-5 |
| p | -1.88671e-5 | -4.95116e-6 | negative | -1.19091267e-5 |

The p contribution exactly reproduces the accepted pressure integral, including its directed endpoint bits. All 10560 pressure density cell hulls are also unchanged. The checker compares 52800 positive masses with an independent high-precision closed-form exponential integral and replays 80 cumulative prefix rows. Original histories, decay factors, and signs are retained separately from the changed contributions.

An independent defining-J plus five Duhamel ODE calculation uses scalar original phase inversion and `history_densities(changed)-history_densities(original)`. Both original histories and changed contributions are integrated with the prescribed rates; quarter-window and terminal values fall within every directed refinement. Two tolerances differ by less than4.1e-13 on these terminal rows. This floating reference corroborates the calculation; it is not the proof of the native directed integral.

The independent reference uses the eta=0 limit at the same fixed native phase origin and N7. Native positive eta is bounded above by1e-520 and is retained in production. The comparison is between two densities, not an absolute density bound: q0 and B0 generally are nonzero. At fixed phase, with s=q² and nu=1+2s,

```text
|dpsi/ds| <= 1/(nu*sqrt(1+4s)) <= 1
|Delta A| <= 5*eta/(8*pi)
|Delta B| <= exp(.1)/pi*(sqrt(5*eta/4)+sqrt(.75)*5*eta/4) < sqrt(eta).
```

For the candidate N7, conservative five density comparison bounds are `(sqrt(eta), eta, 4sqrt(eta), 6sqrt(eta), 2eta)`. The nonnegative-rate Duhamel weights on a unit window preserve those bounds for contributions. B-dependent rows need the sqrt(eta) estimate; the old pressure-only O(eta) estimate would be unsafe for them. Directed square-root/error computations keep this bound valid when serialized.

## Next executable tasks

- [x] **MIDPLANE-ALL-FIVE-OWN-RATES:** all five actual original O2 C0 full-window contributions from accepted whole-cell sources, original B/cross terms, exact positive own-rate masses, true common candidate phase, refinement and independent five ODE comparison.
- [x] **MIDPLANE-INCOMING-MEMORY:** original nonzero h/e/p inlets and cumulative prefixes; explicit affine incoming-defect interface, rate0 pressure memory, separate same P0; no default physical zero inlet.
- [ ] **NEXT GENERAL-ORDERED-SOURCE-MASSES:** integrate the actual defining M0/M1/M2 on ordered original J/f cell covers with positive exponential weights. Return cumulative node integrals and whole-cell f/H/D/P covers, including nonzero inlet constants. Compare against an independent defining ODE, not repeated nested quadrature. Use the source mass rates1.6,.2,1.2 and powers1,2,2, which differ from the normalized own-history decay rates. At Z0, new H/-D/P terminal intervals must agree with these original histories.
- [ ] **NONZERO-Z-COEFFICIENT-CELLS:** bind continuous f/H/D/P covers to full original p1/p2 coefficient templates and the accepted analytic P0/alpha datum. Include both Pstar sectors, delta, L=1-delta*Z² and all pressure tails, with same source family and formal R=Rref*exp(y). Add genuine fixed nonzero-Z and Z-range cells. Enclosure midpoints or log-majorants cannot become fields. Preserve exact Z0 parity while retaining nonzero p2_Z.
- [ ] **NONZERO-Z-INVERSE/INTEGRALS:** enclose signed rho/s/h and the original inverse on each joint source/phase cell; supply both A/B and graph densities for the same N across the complete window. Handle p2=0 as a regular coordinate branch, never by deleting signed terms. Enclose all five own-rate contributions at nonzero Z before claiming all-chart transport.
- [ ] **Z-FUNCTIONAL-MOMENTS:** integrate actual density Z rows with the common source cells and own rates, including p2_Z at Z0. Keep exact correlated pressure/source objects for terminal identities as functions of Z. A midplane scalar or a sampled axial grid cannot establish five functional terminal conditions.
- [ ] **RADIAL/HIGHER-ORDINARY-JETS:** implement slow-y primitive derivatives, add fast A_phi/B_phi terms to total-y modulation, and supply ordinary mixed jets through recovery/stress orders. The fast terms survive the1/N amplitude and cannot be omitted. Maintain original positivity, denominator certificates, tiny endpoint and phase conventions.
- [ ] **ALL-CHART-ORIGINAL/OWN-HISTORIES:** extend the same defining source and own-rate integration through preceding/following charts, quiet intervals, buffer/flatten and exterior joins. Supply actual incoming five defect functions and preserve pressure memory/P0. Install a real continuous source/integral oracle; local C0 summaries cannot be treated as the required C1 jets.
- [ ] **COMMON-N/INDEPENDENT-REPAIR:** solve simultaneous original whole-domain source/stability/tail/frequency constraints for one admitted N. N7 is diagnostic. Compute independent five repair functions and correlated functional terminal closure; a signed local contribution is an input defect, not proof it is canceled.
- [ ] **MATCHED-BACKGROUND/STRESS:** combine repaired five controls, analytic pressure, annular/preheat/heat-exterior joins, finite energy and higher smoothness. Establish the actual admissible stress cone and flat remainder before complete residual admission.
- [ ] **TRUE-COEFFICIENT-RECURSION:** solve the correct n=1 and n>=2 equations on the common inner interval with independent moment repair; truncate through streamfunction/vector potential, then curl. Verify finite-order remainder and flat summation. Scaling an existing velocity field is not this recursion.
- [ ] **OSCILLATORY/PHYSICAL-UVW:** install the two pulse families and mean corrections, validate averaged quadratic stress cancellation, then physical Cartesian uvw and independent full forced-NS residual. Finally report radial/axial core widths, aspect ratio, vorticity/swirl exponents and material winding.

The full long-term goal remains active. This checkpoint advances genuine five-history computation; it does not claim coefficient/scale recursion or a complete forced Navier–Stokes blowup reconstruction.
