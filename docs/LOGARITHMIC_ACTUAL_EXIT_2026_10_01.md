# F22: same-source comparison and actual exit through R110

Date: 2026-10-01. This advances the fresh Md40 degree-110 core into the Section 9 exit construction. Source/datum identities are unchanged; the real axial family is still Z in [0.49,0.51].

## New candidate exit chain

The switched comparison has directed integral enclosures on 32 closed cells over y in [0.005,0.01]; the earlier branch uses the admitted analytic core. It includes the core tail, supplies C2 axial field/moment jets and C1 base-driver jets, and preserves positive Phi throughout the denominator cells. The minimum source Phi lower bound is about 0.0382933.

The actual exit solves g'=chi A_base and Uz'=exp(g) chi B_base. It is not the frozen comparison. Its initial bridge covers 16 core-comparison and 32 switched-comparison cells, retaining C1 axial field and moment enclosures. The actual normalized swirl at the bridge has a positive lower bound about 0.0449678.

After the comparison freezes, the actual exit is continued over the WHOLE interval to physical R100. A single directed driver range controls the log-shape and axial increments. Radial moment weights are integrated analytically, avoiding a log-radius grid across log Lambda approximately 9.4e17. Comparison moments use exact constant-field increments; actual moments use full-domain actual field/axial ranges. This widens correlations but retains integral and directed-arithmetic enclosure for the declared candidate parameters. Paper parameter admission remains a separate missing requirement.

The R100..R110 terminal source switches use 32 closed cells, then exact constant-power moments. The terminal angular log slope is -0.4 and the axial slope is exactly zero. At R110:

| field | directed enclosure summary |
| --- | --- |
| normalized swirl Phi | [0.043413074060995, 0.4992498124076] |
| Uz | approximately [1.96, 2.04] |

Displayed endpoints are rounded summaries, not the authoritative exact endpoints. Physical F0 remains implicit and strictly positive by its analytic definition. Phi widths describe the axial family and interval propagation; they are not truncation-error estimates.

## Scale cancellation and tiny shear

The inlet transfer is rewritten in epsilon ell, epsilon^2 F0^2 and epsilon P0. Its pressure primitive, both swirl contributions and all retained axial derivatives match the frozen physical formulas. A separate fixture checked 75 coefficient comparisons at Lambda=8,500,1e8 and chi=0,0.25,1. Both base drivers require chi; they are explicitly named unmodulated_driver_A/B.

The true microscopic exit shear is

    gamma = exp(-Lambda*(2Grealbar+1)) > 0,

with scalar Grealbar the directed universal real-primitive bound. It is not materialized. The guard Lambda >= 4logLambda+1000 gives

    0 < gamma <= exp(-4logLambda-1000).

Numerical propagation uses [0, positive upper], preserving the true implicit function. A zero lower endpoint does not select zero shear or discard a term. The physical exit cutoff is gamma+(1-gamma) alpha on the first switch, gamma afterward, and the Section 9.4 angular/axial cutoffs are applied separately near R100.

Pressure is recovered throughout from the same swirl primitive in epsilon-P units. The radial velocity is recovered from the axial cumulative moment, not fitted independently. No pressure tail or freely fitted velocity is added.

## Unresolved paper parameter relation

The supplied paper explicitly sets h_b=epsilon_b=c_* K^-100 and uses chi_b=epsilon_b beyond the bridge. This candidate run instead uses h_b=0.005 and an independently chosen microscopic gamma. They are unequal. The receipt flags paper_hb_cstar_K_inverse100_relation_verified=false and full_Section9_parameter_admission=false make this mismatch explicit. Thus the ODE/enclosure implementation and local candidate cone are useful intermediate results, but this is NOT yet a parameter-faithful Section 9 exit. The next action is to recover the K/c_* hierarchy, couple the width and plateau to one epsilon_b, and rerun the connection. The supplied input tests (4.35) and frozen-profile inequalities also need their own certificates. No cone conclusion for this candidate is automatically transferred to the future coupled-parameter field.

## One candidate interface cone

At actual physical R110, Sz=0 and St/F=-0.8, so kappa=0.8 lies on the weak branch of relaxed cone (3.23). The same-source moment transfer gives D=Itheta/F. Exact cancellation reduces its direction/margin checks to D>2. The local axial-family interval bound passes.

This is a certificate at one physical interface over Z in [0.49,0.51]. It is NOT coverage of the whole exit, long reshape, O.2, all axes, pulse or heat collar. It does NOT construct the admissible divergence-form stress lift. F19 O.2 still needs the newly repaired reference inlet.

## Artifacts

Under `experiments/root_st073/`, prefix `lei_ren_part1_paper_`:

- logarithmic_transfer_check.py/.json
- logarithmic_comparison.py/.json
- logarithmic_exit_bridge.py/.json
- logarithmic_exit_continuation.py/.json
- logarithmic_exit_switch.py/.json
- logarithmic_R110_cone.py/.json

Run these scripts in that order using repository-relative paths or from their directory. Each producer checks current source identities and input hashes. Downstream work should consume logarithmic_exit_switch.json for the new R110 endpoint; old Md1.1 transition receipts are incompatible.

## Next actions and unfinished gates

- [x] Fresh comparison integral enclosure with core tails and C2 fields/C1 drivers.
- [x] Scale-cancelled pressure/swirl inlet transfer and required chi modulation.
- [x] Actual initial exit ODE with nonzero implicit microscopic shear.
- [x] Whole-interval actual continuation to physical R100 using exact weights.
- [x] Section 9.4 switches and constant-power moments to R110.
- [x] R110 relaxed weak-branch cone on the current local axial family.
- [ ] Recover K and the c_* restrictions, enforce h_b=epsilon_b=c_* K^-100 (one shared positive parameter), and regenerate the comparison/exit. Retain logarithmic/symbolic scales rather than choosing an unrelated easy width.
- [ ] Certify the supplied exit inputs (4.35), frozen D/E inequalities and core collar domain for the coupled parameters.
- [ ] Build the long reshape from the new R110 endpoint using symbolic logC, segmented radii and normalized moments. Do not exponentiate symbolic log Rref or reuse old amplitude controls.
- [ ] Derive fresh five-moment defect functions and pressure-compatible angular correction from this same source, including the combined fifth moment.
- [ ] Solve the five-bump uniform functional inverse and realize exact reference inlet identities. Then discharge F19 O.2's remaining inlet condition.
- [ ] Prove cone signs through the full new inner exit and long reshape; extend the axial scope. Numerical endpoint positivity is not a full cone certificate.
- [ ] Restore corrected heat pressure and exact heat exterior; verify continuity, axis regularity and radial energy tails.
- [ ] Construct the admissible stress and flat remainder separately, then the true temporal n-dependent recursion, mean and oscillatory corrections and full Cartesian residual validation.

The new candidate exit remains a leading-profile construction with C1 local enclosures and an implicit physical amplitude. It is not a fully sampled physical time-dependent field or a completed matched background. The full goal remains active.
