# Current pulse-end tensor and five tensor traces (2026-10-06)

Current continuation: [actual gap/gap-end tensors and two tensor joins](CURRENT_PULSE_GAP_BACKGROUND_TENSOR_2026_10_06.md) completes the next gap construction below. Twelve actual tensor regions, eleven adjacent joins and four internal traces are now available. Next is actual main/exit construction. Global and temporal obligations remain open.

Implementation and scoped pulse-end tensor/five-trace receipt: commit [97ceb07e](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/97ceb07e520c610fd18b7ce2949fe1c4aac27efd).

The current pulse-end has a completed physical background stress tensor, its divergence, all three physical remainder components and the Cartesian momentum decomposition. It attaches to flatten at the original end endpoint, and the four exact internal end-support edges now have actual common tensor traces. The current tensor chain has **10 regions, 9 adjacent joins and 4 internal traces**. This completes F57C5a-pulse-end, F57C5b-end-flatten and F57C5b-pulse-end-support in the declared source-enclosure scope.

The separate current22 velocity/absolute-pressure atlas remains 14 adjacent / 8 internal source traces. Full-domain cone admission, independent temporal remainder bounds, prescribed-domain total energy, resolved point values and n-dependent temporal recursion remain open. The full-unbounded Gamma exterior retains its previously established regional physical NS identity.

## Sources and execution

- Producer: `experiments/root_st073/lei_ren_part1_paper_compliant_current_pulse_end_background_tensor.py` and matching `.json`.
- Checker and receipt: matching `_check.py` and `_check.json`.
- Focused controller stage: `currentpulseendtensor`.
- Owner: `CurrentPulseEndBackgroundTensor(flatten_power=checked_current_flatten_power_tensor)`.
- The focused controller ran the real producer/checker pair on one reused checked graph; Python compilation passed.
- Read-only scoped review: **GPT-5.6 Luna / max**, existing `/root/heat_terminal_history_scan`; no material source-ownership or mathematical blocker was found.

The owner requires the identical pulse object from `history.selected.pulse`, `physical.pulse` and `flatten_power.flatten.pulse`. Its context, current atlas, complete-energy history, closed absolute pressure and checked flatten/power owner are shared. The older native dispatcher cannot supply this owner. Reuse the checked current graph; reconstruct an inherited prerequisite only when its defining source has changed.

## Full histories and physical operators

The selected coefficient functions and current controls are replayed through the original `capture_native_end`, `pulse_coefficients`, `pulse_velocity_rows` and `lift_physical_packet`. All five original primitive histories remain present. Native forward X, the complete future energy, the selected backward energy loss, the nonzero incoming first moments and full absolute pressure are retained.

For the original end offset `s in [-4,0]`, the ordinary radial rows obey

```
m1' = Bhat - (1/2-mu)*m1
m2' = Bhat - (1/2-2mu)*m2
e0' = 2mu*e0 - 1/2
J'  = 2mu*J - Bhat^2
P'  = (1+2mu)*P + C^2/2,  C=1/(1+Z^2)
eactual = e0 - D^2*J
```

Axial jets through order 5 supply the original mixed stress3/divergence2/remainder2 operators. The full radial, axial, convective, nonlinear and viscosity contributions survive. Each view includes 100 stress contributions, 60 divergence contributions, 36 completed-diagonal contributions, 36 three-component remainder contributions, 43 Cartesian tensor contributions, 14 Cartesian divergence contributions, 11 Cartesian remainder contributions and 25 momentum-decomposition contributions: **325 total**. These are contributions to the field expressions; their count is not a number of distinct field components.

Physical scales remain exact signed source logs. In particular,

```
log B(s) = logPstar + log(current flatten U)
           -13/(2mu) -13 -(1/2+mu)*s
log H(s) = -13*(1-mu)/mu -(1-mu)*s
```

The original current radius map and selected `log D` are retained. Exponentially enormous radii and inverse tiny amplitudes are not materialized. Directed bounds enclose the source factors; they do not replace defining source functions.

## Stable current absolute pressure

The current completed flatten/power pressure supplies the end datum. Let `a=delta/2`, `d=a-mu`, `L=Lrel`, `pp=1+2mu`, and `C0=2*exp(-d*(L+100))`. The source theorem reduces the full flatten inlet pressure before interval enclosure:

```
signed P_Rv / Ev0^2 =
    -JP_flatten(0)/4
    -exp(-100*pp)*I_pp(L-4)/8
    -exp(8*d-pp*(L+96))*P_angular(-4)/4

I_c(length) = (1-exp(-c*length))/c
P(s) = P_Rv*exp(pp*s) + C^2*expm1(pp*s)/(2*pp)
```

`JP_flatten(0)` is the same directly integrated remaining flatten pressure kernel, and `P_angular(-4)` is the same current full angular datum. The cancelled `C0^2` factors are proved symbolically before numerical enclosure. The implementation does not divide by a tiny amplitude or cap. The inherited P0/Pin/Pstar-squared forward history remains available as diagnostic data; admission consumes the current complete remaining-pressure source and current Cp=0 closure. The physical absolute pressure keeps the signed P and its positive B-squared factor.

## End-flatten join and exact support edges

At end s=0 and flatten t=0, the arbitrary-full-history original tensor join theorem is instantiated with the current complete energy, selected controls, absolute pressure, radius and amplitude-unit bindings from `atlas.interfaces.proof`. The checked generic theorem alone is not the current admission. The meridional support terms vanish at this endpoint, while complete future energy, native X and full pressure remain.

The end's 325-sector layout is grouped into the flatten chart's 65 canonical physical contributions after function equality is identified. Theta and diagonal Cartesian contributions retain separate slots. Common directed upper bounds are then formed from both sides; interval overlap is not used as a function-equality theorem.

The theta-viscosity time powers appear in different algebraic forms on the two source paths. Their six mixed2 expressions are proved exactly equal before enclosure. Rounded byte equality of those equivalent interval expressions is not required.

At each rational support edge, the original source is evaluated with exact flat beta and exact full/empty backward weights before floating-point enclosure. All incoming/full moments remain. The current weighted FTC and endpoint theorem identify both one-sided actual tensors and physical remainders. A vanishing local difference is not substituted for the full common boundary tensor.

| Exact edge | Actual tensor contributions | Nonzero source enclosures |
| --- | ---: | ---: |
| -63/20 | 325 | 229 |
| -57/20 | 325 | 229 |
| -23/20 | 325 | 229 |
| -17/20 | 325 | 108 |

Fresh compact time/axial/viscosity sectors pass at each edge.

## Scope and completed checks

- Eight whole-domain, inlet, center, gap, attachment and fresh views: **2600 physical contributions, 2186 nonzero enclosures**.
- End-flatten common interface: **65 contributions, 62 nonzero**; a fresh sector also passes.
- Four actual support traces: 325 contributions each, retaining nonzero boundary histories; a fresh sector also passes at each edge.
- Current selected owner, complete five-history binding, full pressure reduction and the three-component decomposition pass.
- Invalid axial/end/time/viscosity domains, foreign exact edge data, foreign pulse owners and unchecked flatten/power owners are rejected.
- Working/index source hash audit: **668** files, including both saved outputs; staged whitespace check passes.

The positive-time source domain is R>0, |Z|<1, tau>0 with constant nu>0 and finite compact logtau sectors. Z=+/-1 boxes are source-limit enclosures, not the physical axis theorem. Whole original radial-coordinate boxes are included, but broad enclosures do not establish cone margins or resolved point values. The physical remainder is retained; its independent temporal flatness is not established by these radial/support traces.

The actual regional inventory is pulse_end, flatten, outer_power, outer_angular, steep_entry, steep_power, steep_exit, waiting, heat_collar and heat_exterior. Remaining pulse/core regions and four angular internal tensor traces still prevent global composition. No actual n-dependent temporal recovery step has been completed by this milestone.

## Next executable tasks

- [x] F57C5a-pulse-end: actual full tensor, selected controls, five complete histories, signed absolute pressure and three-component physical remainder.
- [x] F57C5b-end-flatten: current source-function tensor attachment, 65 common contributions and fresh sector.
- [x] F57C5b-pulse-end-support: actual common tensor and remainder at all four exact support edges, with retained nonzero history.
- [x] **F57C5a-pulse-gap-gapend / F57C5b-gap-end-end.** Reuse this checked owner and authoritative selected pulse. Adapt the original `pulse_gap_similarity_C4` and `pulse_gap_physical_C2` operators to the current full histories, units and pressure. Construct actual gap-end and gap tensors with stress mixed3, divergence/remainder mixed2, completed diagonal and Cartesian decomposition. Preserve native forward X and all inherited nonzero moments. Use the original reciprocal coordinates and ordinary physical logR derivatives; do not interpret phase derivatives as radial derivatives. Obtain source-function joins for gap-end to end s=-4 and the reciprocal gap/gap-end attachment. Record each endpoint separately, whole original coordinate domains and at least one fresh compact sector. Reject old/foreign/unchecked owners. Save a focused producer/checker/receipt and controller stage, commit/push and mark only the constructed scopes done.
- [ ] F57C5a-pulse-main-exit / F57C5b-main-exit-gap: current actual main/exit tensors using original `pulse_main_exit_similarity_C4` and `pulse_main_exit_physical_C2`. Retain selected axial controls, five histories, meridional terms and absolute pressure. Prove the main-exit and exit-gap tensor function joins in common physical units.
- [ ] F57C5a-pulse-entrance/incoming: current actual entrance tensor using original entrance operators and the retained incoming attachment. Preserve boundary functions, complete moments and source exponential factors. Supply incoming-entrance and entrance-main completed-tensor joins.
- [ ] F57C5a-core-retained: actual fixed-point core, bridge and incoming tensors from the same checked nonlinear core and analytic pressure. Retain the first-interface traces and all common current boundary functions. List every required chart before global composition.
- [ ] F57C4e-axis: derive physical axis regularity and limiting tensor/remainder bounds for that core. Exterior zero and axial source-limit boxes are not an axis certificate.
- [ ] F57C5b-angular-internal: actual tensor common traces at the four angular support edges using the checked angular owner, retained A/E/P and existing local-difference theorem. Keep adjacent and internal inventories distinct.
- [ ] F57C6a-global: compose all actual regional tensors, adjacent/internal traces and the physical momentum decomposition over the declared full domain. Report every uncovered sector; keep global acceptance false until it is covered.
- [ ] F57C6b/F57E: independently bound the physical temporal remainder and high-order/flat decay. Spatial support flatness and radial Taylor transport do not supply a temporal remainder proof.
- [ ] F57D/F57E/F57F: resolve nonzero u,v,w,p, establish current cone margins and realizable lift, and integrate kinetic energy on the prescribed physical domain with independently controlled tails.
- [ ] F58a-c: implement actual n=1 and n>=2 recovery equations, a common core definition interval, per-order five-moment repair, finite-order remainder and smooth summation. Keep n-dependent recursion false until those actual equations are supplied.
- [ ] F59/F60/F61: both oscillatory families, mean corrections, averaged quadratic stress cancellation, independent corrected Cartesian NS validation, and measured contraction/slenderness/winding.

Advance the next coupled construction on the current graph. Save source hashes and scoped receipts, update the top handoffs and commit/push each finished scope. Re-run inherited checks only when a changed prerequisite or new failure justifies them. The long-term goal remains active.
