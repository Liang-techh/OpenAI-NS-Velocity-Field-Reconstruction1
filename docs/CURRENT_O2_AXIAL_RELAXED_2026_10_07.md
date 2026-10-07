# Full original O2 axial turnoff: relaxed input completed

Checked implementation: [1633f6f1](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/1633f6f12f80373ca63996ad3defcdb35957645f). **LEFT4a-axial is complete on phase[0,1], Z[-1,1]. Together with the checked Rh/reference, slope and buffer results, the original upstream outer chain now has its full relaxed input cone.** The axial interior also satisfies the original strict cone where Z is nonzero. Its two flat phase endpoints and midplane remain relaxed; whole closed strict admission is false. This is a source stress result, not coefficient recursion or a corrected NS solution.

## Actual source and the decisive correlation

CurrentOriginalAxialRelaxedCone in experiments/root_st073/lei_ren_part1_paper_compliant_current_O2_axial_relaxed_cone.py consumes checked source packets and receipts, without reconstructing the ancestor graph. query(phase,Z) accepts finite subsets of the declared closed domains. Signed source parameters, original absolute pressure, all cumulative histories and original joins retain the same accepted family.

The ordinary coordinate is y=exp(40*phase)=log(R/Rref), t=y-1. Write U=exp(-1/5)*exp(-t/2), C=1/(1+Z^2), D=Ds*exp(-t), B=sigma(1-log(1+t)/40) and M'=4B-M, M(0)=4. The shared h/U and k/U histories are (1-D)C and (M-4D)ZC. Since B decreases from 1 to 0,

    (exp(t)*(M-4B))'=-4*exp(t)*B'>=0,
    4B<=M<=4.

Using this source correlation before bounding M gives a full theta lower

    Theta >= (1-delta)*Z^2*C^2*(1+8*(1-Z^2)*B)+cD*D+reserve.

The positive reserve includes the selected delta and complete angular radial shear; the exact a=2 shear is retained. Dropping the M/B correlation loses the useful margin at the actual Md40. The separate weaker lower cX*Z^2+cD*da+reserve controls the absolute pressure memory by the same AM-GM geometry.

## Complete axial stress and quantitative cone budget

The unchanged raw_pre_stress_rows operator supplies all six axial sectors: local axial transport, nonlinear meridional transport, retained linear moment, full energy, original absolute pressure and axial radial shear. The linear moment is identically zero here; every nonzero sector stays in the estimate.

The actual kinetic term is

    AZ=16/(Pstar^2*U^2)*(exp(-t)+K2),
    K2=integral_0^t exp(s-t)*B(s)^2 ds.

Its source comparison gives 0<=K2<=1-exp(-t), hence exp(-t)+K2<=1. The angular energy is AQ=AQs-t/2. A0<1 bounds AQs only; the entire -t/2 term enters the dominant budget using t/(1+t)<=1.

The future-pressure function is transported from the actual Rd datum:

    H=1+(H_Rd-1)*exp(t-(exp(40)-1+11)).

The nonnegative backward weights and checked H_Rd range give H_Rd<=H<=1. The signed product uses H<=1; it does not replace pressure by an absolute unit cap. The same tiny absolute pressure memory is bounded separately with its original source amplitude and tail. This is the original pressure FTC and datum, with no adjustable pressure correction.

The local symbolic Pmemory is bound to the inherited exact production formula by current_O2_modified_taper_cone.py:109-118 and current_O3_transition_direction_operator.py:58-65 (both with the lei_ren_part1_paper_compliant_ prefix). Their source hashes and exact theorem bindings are included in this receipt. The read-only reviewer accepted the final kernel and Rd-pressure comparisons without a mathematical or dimensional blocker.

The exact signed shear is a=2, bs=8Z*B_y/(Pstar*U*C), kappa-2=bs^2/2. With r=Tz/Ttheta,

    Dcone/Ttheta=1-bs*r/2,
    Q/Ttheta^2=2*(1+bs^2/4)*(1-bs*r-bs^2/4).

A directed cover of 1024 contiguous rational cutoff cells, including both flat edges, proves sigma'/(1+8sigma)<9/4 and sigma'/y<1e-4 over the whole original cutoff. The complete source budget gives

    bs*Tz/Ttheta+bs^2/4 < .91,
    Dcone/Ttheta > .54,
    Q/Ttheta^2 > .17.

The dominant full AQ term is 16/(40*(1-delta)^2)*(9/4+8delta). Initial AQ and future pressure have separate 1e-4 caps. Kinetic, local/meridional, radial shear, absolute pressure memory and bs^2/4 each have checked source logarithmic suppression. No point sampling substitutes for the whole-domain bound.

## Evidence and next implementation work

Focused checker: 23 exact source identities, 14 strict source budget inequalities, 1024 directed cells, 30 independent signed cone cases, 20 full original stress-sector cases, 6 source queries and 7 invalid guards. The independent stress fixtures also exercise negative conditional H, preventing an accidental absolute-H substitution; these moderate fixtures are algebra checks, not current physical field evaluations. Read-only source/unit reviewer: GPT-5.6 Luna / max.

- [x] **LEFT4a-reference/slope/axial/buffer:** full original relaxed input on the four closed outer charts, including the source joins. See [CURRENT_O2_REFERENCE_SLOPE_RELAXED_2026_10_06.md](CURRENT_O2_REFERENCE_SLOPE_RELAXED_2026_10_06.md) and CURRENT_O2_RELAXED_BUFFER_2026_10_06.md for the other charts.
- [ ] **LEFT4b:** attach the existing current compliant analytic exit collar to the current core-first physical graph. Use compliant_global_exit_certificate.py/.json and compliant_K1_ledger.py/.json, with current core family a4056322..., C-star family 2ed4dc13..., parameter family 6449e4d0..., source 5aec1119... and datum dd4040ee.... The strict interval is Ra<R<=Ra*exp(h_b*sigma_inverse(gamma/(10*K^10))). Pin the current core-first family 3983d0ddb... and prove equality of the actual six inlet moment/pressure atoms, P0 and amplitude/width before claiming physical attachment. The shared legacy exit uses different source/datum hashes and cannot be relabelled. Rh offset-5 and a shifted flat cutoff are insufficient.
- [ ] **LEFT4c:** extend the actual source shear loop across all defective charts from that collar through reference/slope/axial/buffer to the checked right collar. Use ordinary logR derivatives and own-moment radial recovery; preserve divergence, all source joins and the shared pressure datum.
- [ ] **LEFT4d:** recompute the new source family's five correlated defects, independent moment repair, unique implicit terminal control and all-N error envelopes. Existing fixed-N coefficients and minimum frequencies apply only to their checked source definitions.
- [ ] **LEFT4e / BOUND3-global:** prove support-edge directions and complete strict cone after repair, then combine remaining inner regions into a genuinely global sufficient finite N. The existing scoped N>=68,533,403 remains valid only on the checked changed-source taper/O3/pre-repair/quiet domains.
- [ ] **CONT / ENERGY:** remaining physical interfaces and actual required-domain finite energy with volume factors and infinite tails.
- [ ] **REC:** implement actual n=1 and n>=2 recovery equations, independent moment repairs and final smooth sum. Source scaling and cone receipts do not constitute recursive coefficients.
- [ ] **WAVE / PHYS:** oscillatory and mean stress cancellation, corrected u/v/w/p/f, full Cartesian NS residual, quantitative contraction/aspect ratio/scale recursion and material winding.

current_original_O2_axial_relaxed_input_cone_certified=true; whole_closed_axial_strict_cone_certified=false. Global common N, completed global admissibility/physical composition/energy/actual recursion/global flat remainder/full corrected NS remain false. Historical registry counts are unchanged. The full long-term goal stays active.
