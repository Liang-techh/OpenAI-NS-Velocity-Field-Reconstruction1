# Current angular local stress and error differences (2026-10-06)

Current update: common 22 source composition and actual velocity/pressure compact-time trace bounds are complete in [CURRENT_INTERFACE_ATLAS_2026_10_06.md](CURRENT_INTERFACE_ATLAS_2026_10_06.md). This supersedes common affected trace-bound-open text below; actual completed tensor and remaining/global/time gates stay open.

Implementation and scoped local stress/error receipt: commit [d7477882](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/d7477882583af4799795d8a1974c8b1310dcda1a).

The four current angular support edges now have local stress3 and axial-viscosity error2 difference bounds, complementing the four previously checked pulse support differences. The current source inventory remains fourteen adjacent and eight internal traces. Both support families now have local stress/error difference companions. Quantitative common interface bounds and the prescribed complete current physical stress/remainder remain next work; global cone/lift/NS, independent temporal flat remainder, prescribed-domain energy, resolved fields and genuine n-dependent recursion remain open.

Implementation: `experiments/root_st073/lei_ren_part1_paper_compliant_current_angular_support_differences.py`, matching `_check.py`, producer `.json` and scoped `_check.json`. Focused stage: `currentangulardifferences`. Construct `CurrentAngularSupportDifferences(support=checked_current_angular_support)` to reuse the admitted current graph; call `interface(edge,h,Z,log_tau,viscosity)` for a fresh local bound. Producer and checker share one field. Do not rebuild historical deep constructors.

## What is constructed

At each exact source edge `center +/-3/20`, the reference removes only the local beta input while retaining the same nonzero boundary A/E/P histories, pressure datum, current signed C5 coefficient functions and complete future. This comparison does not replace the actual field by zero. The input factor F=1+h has current beta forcing h; its squared difference is `F^2-1=2h+h^2`, retaining both signed and quadratic terms.

The current amplitude normalization is derived directly from current thetaR/S/Q/T and theta_base definitions:

`log(KR)=1+mu/2-3*a/2+(1-a)*Ts+log(1-epsilon)`, with `a=delta/2` and `qR=-wait-Ts-2`.

KR is enormous in this parameter regime. It is kept as an exact positive source logarithm, never materialized and never recovered by dividing capped amplitudes. The identities `B*KR*exp((a-mu)*s)=Ev0*thetaR*exp(-(.5+mu)*s)` and the original B source log recipe are bound before enclosure. Theta stress and viscosity differences carry KR; axial stress and quadratic moment differences carry KR^2.

The normalized differences use `deltaK/KR=exp((a-mu)*s)*h` and `deltaQ/KR^2=exp(2*(a-mu)*s)*(2h+h^2)`. Same-boundary full-moment differences satisfy the original equations:

- `deltaA_y=deltaK-k*deltaA`;
- `deltaE_y=delta*deltaE-deltaQ`;
- `deltaP_y=(1+delta)*deltaP-deltaQ/2`.

Both integral orientations are enclosed by `h*exp(abs(rate)*h)*sup(abs(input difference))`. Normalized beta-tail bounds and current C5 controls give the input jets. All inherited histories cancel only in the difference; their actual values are not reset.

The original collar stress AST implements the linear stress difference. Only its unrelated diagnostic ratio K1/K0 is omitted, since a difference K can be zero and has no cone margin of its own. The stress, inertial/shear parts, current inverse-radius bound, axial derivatives and factor derivatives remain original. An arbitrary nonzero-reference AST superposition proof establishes all forty stress mixed3 rows. The original axial-viscosity operator establishes six mixed2 difference rows with total source order at most4. The error is the difference of `Etheta=-nu*partial_zz(utheta)`; radial and axial errors remain zero in this pure-swirl region.

The original physical stress and viscosity maps consume those coefficient bounds and common log factors. No positive amplitude is rounded to zero. Finite log(tau) and positive finite viscosity sectors can be requested. Source coefficient sectors include Z=[-1,1]; actual finite physical traces use the existing coordinate domains and limiting endpoint interpretation. These are local difference bounds rather than a global physical tensor or temporal-flat remainder certificate.

## Evidence

- 2,712 finite normalized local difference bounds over four edges, h=.01/.000001/0 and whole Z.
- 904 exact source-endpoint zeros, 2,016 decreasing-envelope comparisons and 312 original physical stress3/error2 bound rows.
- Forty original stress mixed3 superposition identities and six original axial-viscosity mixed2 superposition identities retain arbitrary nonzero reference A/E/P/K functions.
- Current full-moment transport, original paper stress units and exact current KR source identities are recomputed; changed steep length or waiting normalization are rejected on isolated clones.
- Fresh Z=.631, log(tau)=-2.7, nu=.8 physical difference view and checked receipt loading pass. Controller preserves one shared field and scoped gates; global/temporal gates remain false.

Exact working/Git-index source and receipt audit passes for all 640 dependency files.

Read-only review by GPT-5.6 Luna / max found no remaining material normalization or physical-scale gap within this local scope. B/Qtheta carry one exact KR log and Qz carries two, including the required B*KR factor in axial viscosity. The numerical S cap remains a directed enclosure of the exact inverse radius; it is not an exact field value. Quantitative/global/NS/temporal gates remain open.

## Next executable work

- [x] F57C4c current angular local input/moment differences: retain current signed coefficients, both beta powers, exact KR/KR^2, full boundary histories and ordinary mixed derivatives.
- [x] F57C4d current angular stress/error differences: original stress3/axial-viscosity error2 operators and physical logarithmic maps over all four support edges.
- [x] F57C3 combine the checked fourteen adjacent/eight internal source receipts under one current selected/future/pressure graph; publish a scoped receipt with exact source and coordinate obligations.
- [ ] F57C4a common adjacent physical bounds: transfer spatial4/fixed-position time1 bounds at all fourteen exact joins, retaining shared radius, current units, lambda powers and declared Z/time/angle sectors.
- [ ] F57C4b common internal physical bounds: compose pulse and angular source/difference companions for eight internal edges; distinguish velocity/pressure trace equality, stress/error difference bounds and quantitative sector admission.
- [ ] F57C4e quantitative interface admission: assemble common bounds over covered sectors and axis limits; record every uncovered sector instead of setting a blanket global gate.
- [ ] F57C5a prescribed current chart stress: instantiate the original stress formulas from the current five cumulative moments and actual derivatives in every required chart. Keep whole histories and the analytic pressure datum; do not promote local difference stress to the actual tensor.
- [ ] F57C5b completed physical tensor joins: assemble Trtheta, Trz and the required diagonal completion, with absolute pressure and the correct source factors; verify current interfaces and gauges.
- [ ] F57C6a current residual decomposition: establish residual=-div(T_B)+E_B in physical coordinates with moving-basis, scale/time and viscosity terms.
- [ ] F57C6b independent global remainder: bound each chart/interface/exterior contribution over the prescribed domain, including actual temporal scale dependence; support-distance flatness does not establish this task.
- [ ] F57D1/D2 resolved coefficients and velocity API: solve the actual defining equations at directed accuracy and expose nonzero u(x,y,z,t),v(x,y,z,t),w(x,y,z,t),pressure with domain/error metadata.
- [ ] F57E1/E2 and F57F: current cone margins/admissible lift/temporal flatness and physical kinetic energy over the required domain.
- [ ] F58a-c: true n-dependent temporal coefficient equations, common core, each order's independent moment/pressure repair and finite-order remainder/smooth summation.
- [ ] F59/F60/F61: both oscillatory families and mean corrections, averaged quadratic cancellation, independent corrected NS and measured contraction/slenderness/winding.

The long-term goal remains active. This milestone completes the two local angular difference tasks rather than the global tensor, NS reconstruction or time-scale recursion.
