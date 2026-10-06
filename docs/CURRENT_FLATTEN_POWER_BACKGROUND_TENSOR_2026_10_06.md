# Current flatten/power tensors and two joins (2026-10-06)

Implementation and scoped flatten/power tensor/join receipt: commit [c5327243](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/c5327243b5ebd82338f65ac454b6ebf541384851).

The actual current flatten and outer-power tensors now attach to the checked angular-entry-O7-collar-Gamma tensor chain. Nine actual tensor regions and eight adjacent completed-tensor joins are available. This completes F57C5a-postpulse-flatten-outerpower and F57C5b-flatten-power/power-angular for the declared source-enclosure scope.

These are completed physical tensor, divergence, remainder and momentum-decomposition maps of the original current background. They retain the complete nonzero histories. They are not resolved point values, cone admission, global NS completion or temporal coefficient recursion.

## Sources and execution

- Producer: `experiments/root_st073/lei_ren_part1_paper_compliant_current_flatten_power_background_tensor.py` and matching `.json`.
- Checker/receipt: matching `_check.py` and `_check.json`.
- Focused controller stage: `currentflattenpowertensor`.
- Owner: `CurrentFlattenPowerBackgroundTensor(heat_tensor=checked_current_heat_tensor)`.
- Reuse the checked live graph. Do not reconstruct old pressure/terminal companions or the historical quadrature chain merely to obtain the new regional tensors.

The ownership chain is current flatten/power tensor → checked current heat tensor → checked O7 → checked entry → checked angular → checked current22 atlas → the common selected/complete-energy/closed-pressure physical graph. Both regional providers are exactly `physical.history.flatten` and `physical.history.outer`. The actual angular `s=-4` full E/P source supplies the right datum. Native forward X supplies A; there is no exponentially unstable backward X recovery.

## Units, coordinates and full moments

Let `a=delta/2`, `d=a-mu`, `k=1-a`, `p=1+delta`. Actual K/A use KR; actual E/P use KR². The exact source log

`logKR=1+mu/2-3*a/2+k*Ts+log(1-epsilon)`

is already included in the analytically reduced `angular_source_log_parts`. Do not multiply those factors by KR again. Keep the enormous scales as exact source logs; directed caps enclose those sources but do not define the field. No radius or inverse tiny amplitude is materialized.

Outer power retains the full original interval `phase in [0,1]`. Its forward coordinate is `yp=(Lrel-4)*phase`, its offset from angular `s=-4` is `y=-(Lrel-4)*(1-phase)`, and its offset from Rrel is `v=-4+y`. Every derivative row is an ordinary logR derivative, not a phase derivative.

With the checked angular datum `Eright,Pright` and `Kright=exp(-4*d)`,

```
K = Kright*exp(d*y)
E = exp(delta*y)*Eright + K^2*I_(2mu)(-y)
P = exp(p*y)*Pright + K^2*I_(1+2mu)(-y)/2
I_c(L) = (1-exp(-c*L))/c
```

Both suffix terms are positive remaining integrals. They use the same original squared swirl density and current Cp=0 pressure closure. The physical absolute pressure is `-B²*KR²*P`, with its positive source factor retained separately. The inherited P/Pstar² packet is retained as a forward diagnostic; it is not used as a normalized remaining P.

Flatten retains `t in [0,100]`, `y=t-100`, `v=-Lrel+y`, `f=(1+Z²)/2`, `F=f^sigma(t/100)`, and `G=F/f`. At its right endpoint `K0=exp(-d*Lrel)`; E0/P0 are the full current power-inlet moments.

```
K = K0*exp(d*y)*G
JE = integral_t^100 exp(-2mu*(v-t))*G(v,Z)^2 dv
JP = integral_t^100 exp(-(1+2mu)*(v-t))*G(v,Z)^2 dv/2
E = exp(delta*y)*E0 + K0^2*exp(2d*y)*JE
P = exp(p*y)*P0 + K0^2*exp(2d*y)*JP
```

The original `flatten_remaining_kernels` evaluates the correlated positive suffix, with exact exponential cell weights. It avoids subtracting independently enclosed past energy from the complete future. The same original flatten theta and forward X functions are retained. Full radial rows satisfy

```
A' = K-k*A
E' = delta*E-K^2
P' = p*P-K^2/2
```

Flatten K depends on Z. Its true ordinary radial rate is `d+log(f)*sigma_t`; native sigma_y_derivatives already include the 100-unit coordinate conversion. The general stress operator retains K_Z through A_Z, and K_Z/K_ZZ enter the physical viscosity remainder. The Z-independent outer-power reduction is not applied to flatten. F/f cancels exactly only at t=100, where the source sigma endpoint is flat.

The original defect AST is replayed on arbitrary nonzero right-data functions, then the constant unit baselines are cancelled before KR normalization. Ninety full-moment mixed4 identities prove equivalence. The checked general-K stress homogeneity theorem and physical tensor decomposition are reused on these actual full moment functions.

## Completed tensors and seams

Each chart supplies stress mixed3, physical divergence mixed2, the completed theta-theta tensor coefficient, physical axial-viscosity remainder mixed2, Cartesian tensor components, Cartesian divergence/remainder, and the momentum decomposition `R_B=-div(T_B)+E_B`. There are 65 physical contributions per view.

Two actual tensor function seams are source-identified:

1. `flatten_power`: flatten t=100 equals power phase=0.
2. `power_angular`: power phase=1 equals angular s=-4.

For each, the current complete-history primitive theorem supplies 60 mixed4 identities; the original arbitrary-terminal stress/absolute-pressure AST theorem transfers the full moments; exact common-radius and ordinary-coordinate identities transfer the physical tensor operators. Only after those function identities are established are common directed bounds formed from both side enclosures. Interval overlap and rounded byte equality of enormous logs are not function-equality proofs.

## Completed checks and limitations

- Real focused controller producer/checker run with one reused checked graph passed.
- Eight whole-domain/endpoint/fresh views: 520 physical contributions, 496 nonzero enclosures.
- Two seams: 65 common contributions each, 62 nonzero each; a fresh compact time/axial/viscosity sector also passed for each seam.
- Ninety original full-moment normalization identities; 120 current primitive seam identities; native theta/shape/rate and complete source bindings passed.
- Remaining E equals the same native complete half-energy after normalization; the interval check is a consistency diagnostic under the separately established source theorem.
- Invalid regional domains, nonfinite coordinates/time, nonpositive viscosity and foreign/unchecked owners are rejected.
- Python compilation and staged whitespace/source hash checks passed.

The source domain is R>0, |Z|<1, tau>0, constant nu>0 with finite compact logtau sectors. Z=±1 boxes are source-limit bounds, not an axis certificate. Whole-original-coordinate boxes are covered, but their broad directed intervals do not establish cone margins. The full-unbounded Gamma exterior regional identity from the previous heat milestone remains available; it does not close the rest of the background.

Actual tensor inventory: flatten, outer_power, outer_angular, steep_entry, steep_power, steep_exit, waiting, heat_collar, heat_exterior; eight adjacent joins. The separate velocity/absolute-pressure inventory is still 14 adjacent and 8 internal source traces. Do not report either count as overall project completion or as n-dependent temporal recursion.

## Next executable construction

- [x] F57C5a-postpulse-flatten-outerpower: both actual full tensors, current histories, common KR units, variable flatten K and ordinary derivatives.
- [x] F57C5b-flatten-power/power-angular: two actual tensor joins, 130 common rows plus fresh sectors.
- [ ] **Next: F57C5a-pulse-end / F57C5b-end-flatten.** Reuse this checked owner and `physical.history.pulse`; construct the actual end tensor from full five histories and the selected axial amplitudes. Reuse original end stress/physical operators with current units. Do not replace the tensor by previously available local differences. Supply whole end-domain and fresh compact sectors, actual nonzero rows and a source-function end-flatten join to the new flatten t=0 tensor. Save producer/checker/receipt and add a focused stage.
- [ ] F57C5b-pulse-end-support: identify actual tensor one-sided traces at each of the four end-support edges, retaining the shared incoming/full quadratic histories. Compose existing flat local differences only after their actual common boundary tensor has been identified. Record each edge separately.
- [ ] F57C5a-pulse-main-exit: actual main/exit tensors on the current selected-amplitude graph; correct reciprocal radial coordinates, five moments, pressure, and their completed tensor interface.
- [ ] F57C5a-pulse-gap-gapend: actual gap and gap-end tensors; preserve inherited nonzero histories and current pressure. Join exit-gap and gap-end/end with original radius/ordinary derivative conversion.
- [ ] F57C5a-pulse-entrance/incoming: actual entrance tensor and incoming attachment; preserve source exponential factors and all retained boundary data. Record common completed tensor traces rather than only velocity/pressure continuity.
- [ ] F57C5a-core-retained: actual fixed-point core and bridge/incoming tensors from the same checked nonlinear core and compatible analytic preheat pressure. Retain first-interface traces and all current retained boundary functions.
- [ ] F57C4e-axis: separate physical axis regularity and limiting tensor bounds for the admitted core; Z=±1 and exterior zero are not the axis theorem.
- [ ] F57C5b-angular-internal: lift the four existing angular local-difference companions into actual completed tensor common-boundary joins with nonzero inherited A/E/P. Keep internal and adjacent counts distinct.
- [ ] F57C6a-global: compose every actual regional tensor, adjacent/internal trace and `R_B=-div(T_B)+E_B` on the declared full domain. List any uncovered sector and keep the acceptance flag false until every required sector is covered.
- [ ] F57C6b/F57E: independently bound the physical temporal remainder and its high-order/flat decay; radial Taylor transport and support-distance flatness do not establish temporal flatness.
- [ ] F57D/F57E/F57F: resolved nonzero u,v,w,p, current cone margins and realizable lift, required-domain kinetic energy and independently controlled tails. Fix the physical integration domain before claiming finite total energy.
- [ ] F58a-c: actual n=1 and n>=2 temporal recovery equations, common core interval, per-order five-moment repair, finite-order remainder and smooth summation.
- [ ] F59/F60/F61: both oscillatory families and mean corrections; realizable averaged quadratic stress cancellation; independent corrected Cartesian NS validation and measured contraction/slenderness/winding.

For each completed bounded task, save the current defining source/receipt hashes, actual checks and limitations, update the top handoffs, mark only its precise scope done, and commit/push. Advance construction; do not rerun the whole inherited chain unless a changed prerequisite requires it. Global cone/lift/NS/remainder/energy/points and temporal recursion remain false. The long-term goal remains active.
