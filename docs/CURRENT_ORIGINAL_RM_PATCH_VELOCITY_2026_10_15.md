Latest continuation: [Original selected-N Rh-to-Rref numerical/physical values](CURRENT_ORIGINAL_OUTER_REFERENCE_2026_10_15.md) now carries this same source through the next complete reference annulus. The accepted Rm/Rh source and its limits below remain unchanged.

# Physical original finite-N Rm-to-Rh velocity - 2026-10-15

`OriginalRmPatchVelocity.evaluate(Z,x,log_tau="-2",azimuth=".25")` now assembles physical cylindrical and Cartesian leading+n1 velocity values across the actual original patch, 1<=x=R/Rm<=e. Use x="Rh" for the exact endpoint. It consumes the accepted selected-N local/five-history numerical enclosures without re-solving the leading controls or rebuilding ancestor computations. The full reconstruction remains **ACTIVE / INCOMPLETE**; tasks14/15 remain IN PROGRESS.

Values are signed sums of `coefficient*exp(sum(log_scale_terms))`. Every coefficient and logarithmic atom has a directed enclosure. The original huge physical exponentials are not materialized. Native coordinates (Z,x,azimuth,time) are required; there is no direct arbitrary physical (x,y,z,t) inverse for this patch yet. Finite-N centers remain leading proxies with complete unknown-signed corrections, rather than resolved fast-phase values.

## Physical geometry and temporal powers

The native geometry is y=log(R/110), so logR=log110+y. The factor110 is included exactly once:

```
R = 110*exp(10logCstar+10logPstar-6)*x
lambda = sqrt(tau/(1-Z^2)), tau=1-t
r = lambda*sqrt(2R), z = Z*lambda^(1-delta)
ur0 = lambda^-1*sqrt(R/2)*Q
utheta0 = lambda^(-1-delta)*Pstar*E
uz0 = lambda^(-1-delta)*Uz
psi0 = lambda^(1-delta)*R*(Mz/R)
```

Cartesian values use u=ur*cos(theta)-utheta*sin(theta), v=ur*sin(theta)+utheta*cos(theta), w=uz. Original viscosity1, selected N, parameters and source family remain. `log_tau=None,scaled_time="1e18"` selects S=-delta*log(tau)/2; delta*time cancels before interval boxing.

## Retained n1 primitive and local support admission

The post-cutoff terminal primitive J and its actual Z derivatives remain nonzero input data. The checked original lower cutoff and positive T/L place rminus<R110<Rm. The exact outer recipes give log(Rw/Rh)=logPstar+6>0, Rc=Rw*exp(2), and the reserved supports lie beyond Rc. Checked intermediate own-field source functions have u1=f1=0 and mean1_y=-mean1 on this local pre-repair branch. Consequently (rho*mean1)_y=0; J is continued by its source ODE, rather than declared constant just because the inlet was copied.

```
QJ = (-(1+delta)*Z*J-(1-Z^2)*J_Z)/(1-delta*Z^2)
psi1 = J*lambda^(1+delta)*Lambda^-2
ur1 = QJ*lambda^(-1+2delta)*Lambda^-2/sqrt(2R)
uz1 = utheta1 = 0 locally
```

The n1 streamfunction has no radius atom. Zero local axial/swirl values do not erase the meridional history. Source receipts and all their inputs are pinned; this is a local homogeneous continuation, not admission of the still-open global high-order n1 joins. The exported temporal-order-zero pressure includes the same P0 once; full n1 pressure is not assembled.

## Focused validation

`python scripts/run_original_Rm_patch_velocity.py --write-report` exits0 in47.267s. An independent moderate-parameter example differentiates the streamfunction at fixed physical r/z, including inversion of the axial coordinate. Its144 velocity/coordinate/streamfunction comparisons have maximum scaled error4.40079976818e-69. This example verifies conversion algebra and does not replace original-parameter evidence.

The actual source run exports20 physical packets across5 Z positions, interior native cells, Rh and ordinary/scaled time. It checks24 reference-to-Rm physical sum seams in common factored units (including signed zero crossings),51 N/history/local-support identities,68 radius/zero-local-n1 conditions,21 scaled-time invariants and3 rejected invalid requests. Saved-report inspection also exits0. No broad historical suite or ancestor report is rebuilt.

Report: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rm_patch_velocity.json`; 10,061,126 bytes; SHA256 `2469d38d4f78c1887b429ca3df330e5948eba60d89cff1b8779ecf7ba926f8d9`. It pins1837 scientific inputs and51 implementation files. Existing read-only GPT-5.6 Luna/max review found no local physical algebra error; its admission caveat is explicitly retained above.

## Next work

1. Bind the original Rh exit to the same-source outer reference/reshape charts and preserve all five finite-N own memories, independent P0 and J. Do not substitute historical N257/N1024 source families.
2. Implement native outer values through Rc and carry the actual remaining signed defect enclosure into the later reserved repair. Do not assert modified terminal zeros from the leading identity.
3. Derive the total radial/fast-phase and higher-Z rows needed for physical derivatives and joins; current value assembly does not supply gradients, vorticity or global Cartesian C2.
4. Complete global pressure/n1 joins, common N, exact heat exterior, physical energy and admissible stress/flat remainder with quantitative estimates.
5. Construct genuinely independent temporal n>=2 coefficients, the smooth sum and two-family oscillatory cancellation, then independently validate corrected Cartesian NS with fixed forcing. The present leading+n1 value extension does not count as temporal scale recursion.
