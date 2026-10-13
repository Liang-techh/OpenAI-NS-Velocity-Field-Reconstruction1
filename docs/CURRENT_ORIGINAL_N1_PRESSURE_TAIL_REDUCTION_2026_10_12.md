# All-chart own mean, actual pressure-forcing support and exact quiet primitives

The actual finite-N partial own mean now admits all 17 original charts. The known n1 pressure forcing has an exact finite support endpoint Rv, and the genuine repaired quiet segment 2Rc to Rp has exact signed pressure primitives. These remove two large integration obstacles. The full reconstruction goal remains **ACTIVE / INCOMPLETE**: complete signed pressure defects and the n1 moment repair are still open.

The earlier [physical-field checkpoint](CURRENT_ORIGINAL_N1_PHYSICAL_PRESSURE_FIELDS_2026_10_12.md) supplies physical leading U/m/q and ten local nonradial moment inputs. Its two-micro-cell numerical scope has not been silently promoted to all-chart physical moments. The new all-chart interface extends the finite-N mean/Q increment, not the full signed pressure-tail quadrature.

## Accepted interfaces and evidence

All modules below live in `experiments/root_st073/` and use prefix `lei_ren_part1_paper_compliant_`. Each has its own source, raw report, independent checker and checked receipt, bound to the unchanged accepted input graph.

| Module suffix | New accepted scope | Evidence |
| --- | --- | --- |
| `current_original_n1_all_chart_partial_mean` | Own mean Z3, mean_y Z3 and Q increment Z2 on all 17 original charts | 17 actual first-cell queries, 17 checked-owner full-chart union queries, 11 independent moving-endpoint rows |
| `current_original_n1_pressure_forcing_support` | Known n1 pressure forcing and upper endpoint terms vanish from actual Rv to infinity | 30 mixed axial/radial zero rows; exact coverage by nine future charts; actual pressure operator checked with nonzero angular field |
| `current_original_n1_quiet_pressure_primitives` | Exact signed pressure moments on actual repaired 2Rc to Rp quiet continuation | Actual inlet mean is linear in Z; q_Z cancellation; direct physical known-pressure equation and signed primitives |

The scoped gates are:

```text
current_original_finite_N_all_17_chart_partial_own_mean_Z3_mean_y_Z3_and_Q_increment_Z2_installed
current_original_n1_known_pressure_forcing_zero_from_Rv_to_infinity_certified
current_original_n1_quiet_2Rc_to_Rp_pressure_primitives_installed
```

The actual all-chart replay uses the original Z cell [0,1/2]. Its whole-chart queries cover native radial partitions at that Z cell; they do not certify a global Z cover or signed moment accuracy.

## All 17 original charts

`CurrentOriginalN1AllChartPartialMeanQ(checked_partial_owner)` reuses the signed free-endpoint mean functions and the actual incoming N-scaled C3 histories. `range_cell` accepts exact contiguous original native partition unions; `range_chart` uses the full original chart. Symbolic `sc_half` and `Rh` endpoints remain original source endpoints. Arbitrary replacement endpoints are rejected.

The rate-one own-mean kernel, original physical Jacobian, true predecessor memory and selected epsilon_N are retained. Incoming and local-density absolute bounds give directed mean/mean_y/Q-increment enclosures. Quiet local density does not erase incoming mean. Baseline mean and physical S=Pstar still need composition before full pressure moments are evaluated. Magnitude bounds are not signed centers or relative-width accuracy.

## Actual finite pressure-forcing support

```text
Rv = exp(original_logRp + 13/original_mu).
```

The actual end axial bumps and their backward histories already vanish on an open collar ending at Rv. At flatten and throughout outer_power, outer_angular, steep_entry, steep_power, steep_exit, waiting, heat_collar and heat_exterior, the same source defines U=Ur=0 with their required derivatives. The same-family original radius maps cover the entire ray Rv to infinity, including all joins. Therefore q=Ur/sqrt(R/2), q_R, S and N1p0 vanish there, as do every upper boundary term in the two reduced pressure identities.

Both remaining forcing integrals can stop at Rv exactly. The angular velocity and pressure rows remain nonzero where the original source prescribes them. This is a pressure-forcing support theorem, not a statement that the full leading velocity/pressure, global stress or corrected NS solution vanishes. O3 loop flatness before the repair band is not a support argument.

## Exact quiet pressure contribution

The already accepted repaired-exit/native-frame checker proves that the actual compact n0 repair is consumed before quiet continuation. Its actual analytic inlet strings give common-unit mean `kappa*Z`, with kappa independent of Z. Keep the physical constant formal:

```text
K = Rp*Pstar*kappa
m = K*Z/R, U = 0, m_y = -m
q = -K/R, q_Z = 0, q_y = -q
S = -q^2/2
2*q_R + R*q_RR = 0
N1p0 = K^2/(4*R^2).

Integral[a,b] N1p0 dR   = K^2/4*(1/a - 1/b)
Integral[a,b] R*N1p0 dR = K^2/4*log(b/a)
a = 2Rc, b = Rp.
```

The m13 contribution uses the first expression with a plus sign. The m15 contribution uses the second with a minus sign. At fixed original endpoints both Z derivatives are zero. `CurrentOriginalN1QuietPressurePrimitives().primitive(K,a,b)` returns these signed expressions and the physical fields. Original endpoint function nodes and accepted repaired-exit source artifacts remain bound; huge radii, N and Pstar are not materialized.

Preserve this exact cancellation when enclosing the quiet contribution. Independent magnitude bounds on mean rows can lose q_Z=0 and create spurious large R-weighted terms. U=0 does not imply q=0 on this segment: the incoming own mean drives a nonzero radial outflow.

## Detailed next actions

- [x] Extend the partial own-mean numerical interface to all 17 original charts and exact original partition unions.
- [x] Execute the first original cell and checked full-chart query on every chart, retaining symbolic endpoints and original predecessor memory.
- [x] Bind actual Rv and all nine future radius maps; certify known n1 pressure forcing and its endpoint terms are zero from Rv to infinity.
- [x] Reuse the already accepted actual C3 repaired exit, quiet continuation, Rp native frame and local interval inlet. Their later accepted gates supersede old historical missing-callback flags.
- [x] Extract the actual Z-independent quiet kappa and exact physical K/q relation; publish signed quiet pressure primitives with the m13/m15 signs intact.
- [ ] Compose physical baseline U/m Z3 and finite-N increments on all 17 pre-repair charts. Use the same source/P0/N/global phase and apply S once.
- [ ] Export q through Z2 and honest value/Z1 moment functions for every chart. Do not fill absent derivatives with zero or differentiate magnitude caps.
- [ ] Integrate the ten nonradial moments with signed correlation or rigorous signed intervals. Keep distinct free query, inner mean integration and outer pressure integration variables.
- [ ] Compute the unchanged original leading contribution from Rin to r_minus, using the actual physical lower radius and q/q_R at Rin.
- [ ] Bind physical moment contributions on the actual reserved Rc to 2Rc repair band using `current_original_C3_repair_band`. Retain all own repaired histories; do not continue the pre-band source through it unchanged.
- [ ] Bind exact a=2Rc, b=Rp and the same-source K to directed numerical/formal enclosures. Use the quiet primitives instead of quadrature over this enormous interval.
- [ ] Compose the actual selected pulse fields from Rp to Rv using the accepted Rp interval inlet, selected pulse owner and segmented radius. Retain true axial/radial histories and the original Jacobian once.
- [ ] Sum all segments without overlaps or gaps. Keep original log-radius origins and inverse-mu factors symbolic wherever direct floating-point radii are impossible.
- [ ] Apply the independently verified I0(S)/I1(S) identities and lower endpoint terms. Replace the upper terms by zero only using the accepted Rv support theorem.
- [ ] Produce both complete signed forcing integrals and their Z1 rows. Add the unweighted tail to m13 and subtract the R-weighted tail from m15.
- [ ] Combine full m13/m15 with the available full velocity-support m11/m12/m14. Report actual five signed defects, separately from enclosure widths.
- [ ] Bind positive I2/Astar, disjoint normalized 2+3 bump matrices and the actual original repair interval. Compute alpha/beta from all five defects with the original signs and physical units.
- [ ] Reconstruct the n1 own M/mean/V/P after its repair. Verify all five terminal identities, pressure compatibility and derivative joins without applying the cutoff directly to V/P.
- [ ] Finish the common-inner axis remainder and derivative joins; implement genuine n2 forcing and order-dependent recovery/repair before temporal summation.
- [ ] Continue original admissible stress, oscillatory stress cancellation, flat remainder and independent corrected NS diagnostics. Do not substitute scaled copies of n0 for temporal coefficient recursion.

## Handoff and scope

The next concrete blocker is signed physical pressure-moment accumulation over the remaining actual source segments, followed by the full n1 repair. The quiet interval now has an exact formula and the future tail beyond Rv is exactly zero. Missing repaired-exit or local Rp callback is not the current blocker: those source interfaces are already accepted.

Complete signed tails, n1 repair, full-inner remainder, n2/temporal recursion, global stress, flat remainder and full corrected NS remain open. Keep their gates false. Reuse existing accepted receipts; rerun only checks affected by new changes. A response-stream interruption does not erase committed or locally saved artifacts. Read this checkpoint and Git status when resuming.
