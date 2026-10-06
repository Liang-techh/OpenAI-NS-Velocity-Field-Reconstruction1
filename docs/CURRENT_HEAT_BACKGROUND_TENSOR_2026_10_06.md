# Current collar/exterior tensors and heat attachment (2026-10-06)

Current continuation: [actual flatten/power tensors and two joins](CURRENT_FLATTEN_POWER_BACKGROUND_TENSOR_2026_10_06.md) completes the next flatten/outer-power construction below. Nine actual tensor regions and eight adjacent tensor joins are available. Next is actual pulse-end tensor construction and its end-flatten attachment. Global and temporal obligations remain open.

Implementation and scoped heat tensor/join receipt: commit [d04f7cda](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/d04f7cda9efcbfbdb2158e5247c512ee7f677dc1).

The actual current heat collar and exact Gamma exterior now have completed physical tensors, stress mixed3 and physical divergence/remainder mixed2. Waiting-collar and collar-exterior are source-function completed-tensor joins with common physical bounds. The complete current chain now has seven actual tensor regions and six adjacent tensor joins. The separate velocity/absolute-pressure trace atlas remains 14 adjacent and 8 internal.

The full Gamma exterior additionally has a source-exact regional physical Navier–Stokes identity: its stress tensor, divergence, axial-viscosity remainder and all momentum components vanish as functions. Velocity, full energy and absolute pressure remain nonzero. This covers every original finite exterior offset t>=3 and its infinity limit, without a radial cutoff. It does not certify the remaining background, the axis, global admissibility/energy/flatness or temporal coefficient recursion.

Implementation: `experiments/root_st073/lei_ren_part1_paper_compliant_current_heat_background_tensor.py`, matching `_check.py`, producer `.json` and scoped `_check.json`. Focused controller stage: `currentheattensor`. Consume the admitted [current O7 tensor owner](CURRENT_STEEP_WAITING_BACKGROUND_STRESS_2026_10_06.md), retaining the same physical/history/heat and checked `CurrentFullExteriorStress` object. This uses the latest closed current Cp/Dtheta/full-energy source graph, not the older heat companion whose terminal constants were still retained.

## Callable result and domains

```python
from lei_ren_part1_paper_compliant_current_heat_background_tensor import CurrentHeatBackgroundTensor
tensor = CurrentHeatBackgroundTensor(o7=checked_current_O7_tensor)
collar = tensor.chart("heat_collar", Z=".537", t=".417", log_tau="-2.6", theta=".41", viscosity=".8")
join = tensor.interface("waiting_collar", Z=("-.8", ".8"), log_tau=("-5", "-2"), viscosity=".2")
exterior = tensor.chart("heat_exterior", Z=".731", t="4.17", log_tau="-2.6", theta=".41", viscosity=".8")
whole_exterior_identity = tensor.unbounded_exterior()
```

The collar's original offset t is in [0,3]; the exterior finite-offset API requires t>=3. Both return ordinary logR derivative rows. Physical time enters separately through log_tau. Source Z=[-1,1] covers finite physical R>0, |Z|<1, tau>0 and constant viscosity nu>0. Source Z=+/-1 are infinity limits, not the axis. theta=None covers all angles. Compact positive-time estimates are distinct from a uniform critical-time estimate.

Every numerical moment, stress, pressure and physical contribution is a directed source enclosure with exact positive scale logs. No cap or interval midpoint defines a resolved point field. `unbounded_exterior()` evaluates the zero-function theorem over the entire original exterior; a finite t=3 anchor supplies the component layout only. It does not replace infinity by a finite computational cutoff or introduce infinite source logs into the physical row calculator.

## Full current moments and stable stress

The local heat normalization is converted to the same exact KR units used throughout angular/entry/O7. Let LT=a-mu-(1-mu)/2-k*Ts-k/2, k=1-a and a=delta/2. Then

`C = 1/KR = exp(LT-logone)`.

This identity follows from the checked waiting source: LT+logKR=logone. It never divides an amplitude cap. Native collar K, angular numerator A, complete energy E and reduced positive pressure P are converted as K/KR=C*K, A/KR=C*A, E/KR^2=C^2*E and P/KR^2=C^2*P. At the inlet, C*(1-epsilon) equals the exact normalized waiting shape. All nonzero current angular/energy/pressure histories remain included.

Full collar moments use the original `collar_tails` angular numerator, remaining energy times exp(delta*t), and remaining pressure times exp((1+delta)*t). A 45-identity symbolic replay binds the actual native numerator ASTs to `collar_defect_rows`, their unit baselines and the common normalized mixed4 recurrence. The canonical infinite Gamma integral and the separate epsilon/epsilon-squared atoms remain the defining source.

The collar stress is evaluated using the original stable small-defect AST, then converted by C for theta and C^2 for axial stress. Unit baselines cancel before enclosure. Computing large full baseline quantities and subtracting their interval boxes would lose the small stress; that operation is avoided. The checked general-K baseline/homogeneity theorem identifies this stable expression with the full actual moments. K_Z/K_ZZ and the original axial physical operators remain present.

## Same original pressure and current terminal functions

The latest current closure proves Dtheta=0 and Cp=0 for the same native recipes, prescribed analytic P0/Pin, complete selected energy and canonical Gamma future. Native forward pressure and angular diagnostics remain preserved in the output.

The adapter binds the original forward-pressure density K^2*exp(-p*t)/2, p=1+delta, `P3=forward_pressure(3,Ptail)`, the collar remaining-pressure function and the exterior formula involving P3 and the full local Gamma pressure numerator. Complete FTC additivity gives

`I0 = prefix_0_to_3 + exp(-3p)*N3`,

so `P3+scale*exp(-3p)*N3 = Ptail+scale*I0`. These are the same native infinity constant already proved zero by the checked current pressure closure; their interval overlap is not the proof.

The exact waiting/full-heat pressure endpoint exponent is also recorded: p*(Ts+2+wait)+rB=2*(LT-logone). This identifies the O7 future pressure Ih with the local collar inlet after the common KR conversion. Original signed absolute pressure is restored after the positive reduced remaining moment and its pP-K^2/2 recurrence. No pressure constant or gauge is fitted.

## Physical exterior identity and two joins

The checked current exterior's five-history function bridge supplies the actual zero-stress theorem. The original physical heat identity then uses R=r^2/(2*lambda^2), d=tau/lambda^2, xi=2d/R=4tau/r^2. The physical lambda powers cancel exactly. With constant viscosity nu>0,

`u_theta = A_nu*r^(-1-delta)*H(4*nu*tau/r^2)`, `u_r=u_z=0`,

`p = -A_nu^2*integral_r^infinity q^(-3-2delta)*H(4*nu*tau/q^2)^2 dq`,

where A_nu=nu^(1+delta/2)*A and A=2^((1+delta)/2)*Ev0*theta_base*Rtail^((1+delta)/2) is the actual source amplitude. This is not a new free parameter. The full Gamma ODE supplies the heat identity, the same pressure FTC cancels the centrifugal force, and physical axial derivatives vanish. The exterior regional residual and remainder therefore vanish exactly. The original finite interval operator enclosures are retained as diagnostics; exact zero rows evaluate the source-function theorem after the current five histories have been identified.

Waiting-collar and collar-exterior each consume 60 current primitive mixed4 identities, the exact ordinary-coordinate/radius proof, the original arbitrary-terminal waiting stress/pressure AST join or canonical collar/Gamma stress AST theorem, common KR units and the unchanged physical tensor operators. The completed diagonal remains Ttheta_theta=r*partial_z(Trz). Source equality transfers to stress3, divergence2, diagonal2 and remainder2 before valid endpoint envelopes are merged. The collar/exterior completed-tensor join is source-exact zero. Waiting/collar is generally nonzero.

## Scoped evidence and executable next construction

Eight views cover the entire collar, inlet, Gamma boundary, phi endpoint crossing, a fresh collar point, the exterior inlet, a fresh exterior point and a far exterior interval. They check 520 physical tensor/divergence/remainder rows, with 238 nonzero enclosures. Both new interfaces have 65 common contribution bounds plus fresh compact sectors. A further 65 zero rows cover the full unbounded Gamma exterior by function identity. Forty-five original full-moment normalization identities and 120 current primitive seam identities bind the source functions. Invalid domains and foreign/unchecked current owners are rejected. Global and temporal flags remain false.

Read-only review metadata: GPT-5.6 Luna / max. Scoped review passes with no material blocker. The authoritative complete retained-pressure history bridge proves the same-density Ptail/P3 infinity identity and is consumed through the accepted current source/hash chain. The new local exponent/FTC identities are consistent with that bridge. An additional direct assertion of its receipt flag would be optional provenance hardening, not a new mathematical gap. The reviewer confirms the full Gamma physical zero identities and the unbounded function claim; the finite anchor supplies layout only. The actual focused controller passes with one shared graph, seven regions/six tensor joins and the separate 14/8 primitive trace atlas. Compilation and raw working/Git-index SHA256 checks pass for all 660 dependencies; Git whitespace checks pass. Unrelated experimental edits are preserved.

- [x] F57C5a-postpulse-collar: actual current full-moment collar tensor with stable original small-defect stress and complete Gamma future.
- [x] F57C5b-waiting-collar-join: common source units, actual shape/pressure/energy identity and completed-tensor bounds.
- [x] F57C5a-postpulse-exterior: current full Gamma physical tensor, divergence/remainder and exact regional NS identity on the original unbounded exterior.
- [x] F57C5b-collar-exterior-join: source-exact zero completed-tensor/physical-remainder attachment.
- [x] **F57C5a-postpulse-flatten-outerpower and F57C5b-flatten-power/power-angular joins.** Completed in [the current flatten/power tensor layer](CURRENT_FLATTEN_POWER_BACKGROUND_TENSOR_2026_10_06.md): same complete energy/pressure graph, variable flatten K, ordinary stress3/divergence2/diagonal2/remainder2, actual physical tensors and both source-function joins. Nine regions/eight tensor joins are available; global and temporal scope stays open.
- [ ] F57C5a-pulse-end: current actual end tensor from all five complete histories and selected amplitudes, including meridional cross terms. Connect end-flatten and all four internal end support edges; local differences are not actual tensor values.
- [ ] F57C5a-pulse-main/gap/entrance: actual entrance/main/exit/gap/gap-end tensors, reciprocal coordinates, original Pin/P0 and complete future/2. Connect all adjacent pulse tensor joins to end and the retained incoming chain.
- [ ] F57C5a-core-retained/F57C4e-axis: common fixed-point core and original incoming/bridge tensors, axis regularity, first core/bridge traces and all remaining retained boundaries.
- [ ] F57C5b/F57C6a-global: compose every actual chart tensor/decomposition and every tensor interface, including angular support edges; preserve separate tensor and 22-trace velocity-pressure inventories.
- [ ] F57C6b/F57D/E/F: independent global temporal remainder, admissible cone/lift, prescribed-domain kinetic energy and resolved source-defined u,v,w,p.
- [ ] F58/F59/F60/F61: actual n-dependent temporal equations/per-order repair and summation, both oscillatory families/mean corrections, global corrected NS and measured contraction/slenderness/winding.

Save each coupled construction and scoped source receipt, mark only its actual task complete, commit/push and continue. Exact exterior heat flow, spatial moment recurrences and regional zero tensors do not certify temporal scale recursion. The full goal remains active.
