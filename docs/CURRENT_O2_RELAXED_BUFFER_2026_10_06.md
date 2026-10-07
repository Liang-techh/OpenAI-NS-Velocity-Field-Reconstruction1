# Original O2 relaxed buffer input and actual strict-collar obstruction

Checked implementation: [574723cf](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/574723cfd8b3a8867518f4ca641598a19db105d7). **The current original O2 shared buffer offset[-11,-2], all Z[-1,1], including both seams, now has a continuous source-bound relaxed input cone proof.** This is input for the paper's shear modification, not a strict admissible stress certificate. Stress is nonzero and the exact shear remains a=2, bs=0, kappa=2 at every frequency. The previous changed-region strict cone/common-N result remains valid on its stated domains, N>=68,533,403; no global gate or original registry count is promoted.

## What was implemented

CurrentOriginalO2RelaxedBufferCone in experiments/root_st073/lei_ren_part1_paper_compliant_current_O2_relaxed_buffer_cone.py reads the checked original O2/O3 source and tensor packets. query(offset,N=1) accepts finite shared offset subsets of[-11,-2] and positive finite integer N. The actual buffer selector is11+offset in[0,9], and R=Rd*exp(offset). Source profiles and all cumulative modulation errors are exactly zero here because modulation has not started. The original nonzero tensor, incoming moments, full kinetic history, absolute pressure datum/memory and radial shear are retained.

The exact defining buffer source is replayed through the accepted O2 ODEs and pressure transport. Its correlated canonical theta lower is Theta>=cX*Z^2+cD*da*exp(-t)+reserve. Extension to the full buffer uses the worst log radius logRd-11 and the full AQ bound abs(EQa/Ua^2)+11/2, rather than the taper's shorter interval bounds. The same future-pressure H(t)=(1-exp(t))+exp(t)*H(0) stays in its checked positive-to1 range for every t<=0. U(t)>=Ua gives a conservative complete theta and constant absolute-memory bound. The remaining original energy/pressure axial expression is preserved as a signed source function and bounded separately.

The canonical reserve and full theta are strictly positive. With F>0, S=(-2F,0), the exact paper criterion is kappa=2, T dot S=-2F*Ttheta<0. Equation(3.23) at kappa2 therefore holds strictly. D/Ttheta=1 and Q/Ttheta^2=2, but kappa-2=0. The new provider exposes relaxed admission as true and strict admission as false. Frequency changes cannot remove this source obstruction.

## Focused evidence

The checker passes 17 exact source/cone identities, 6 strict whole-buffer baseline inequalities, 15 independent full signed stress cases on earlier offsets, 240 exact local/cumulative/radial mixed4 source jet zeros, 28 exact pre-support tensor modification zeros at N=1 and N=10^12, 5 closed buffer/seam queries and 7 invalid guards. It explicitly rejects strict admission for the nonzero degenerate stress and compares the actual flat modulation edge with the checked taper source. Full-domain admission comes from source identities and directed bounds, not these fixtures. Read-only reviewer/model: GPT-5.6 Luna / max.

## Why shifting the cutoff is insufficient

Paper: Lei-Ren v2 p22 equations(3.21)-(3.23); the strict target requires kappa>2 wherever stress is nonzero, while the relaxed input has a special kappa<=2 branch. Section11 also requires already-admissible collars outside the shear modification; its loop needs kappa>=2+eta where it turns off. The supplied paper/cache and repository cone operator provide the exact references.

The turnoff/buffer seam is O2_axial phase1 to O2_buffer offset0, at y=exp(Md). Replacing sigma(t+2) by sigma(t+11) simply moves the flat source correction edge to t=-11. Sigma and its derivatives vanish there, so a2/bs0/kappa2 still holds with nonzero stress. The upstream source does not already supply an admissible collar:

- Rh_reference has Utheta_y/Utheta=1/10, bs0, kappa=4/5.
- O2_slope has a=4/5+(6/5)*sigma(y)<=2 and bs0.
- O2_axial has a2 and bs zero on Z=0 and both flat endpoints.
- O2_buffer has a2/bs0 throughout.

These source facts are bound in pre_pulse_mixed_C4.py:164-200 and current_O2_background_tensor.py:103-106. Rh_reference offset-5 is a source boundary, not an accepted repaired-core/inner admissible collar. A tensor-only adapter cannot alter the source shear kappa.

## Ordered tasks for the next agent

- [x] **BOUND2b-buffer-relaxed:** full original buffer/edge relaxed direction and exact shear proof, with pressure/energy/radial sectors retained.
- [x] **COMMONN-changed / BOUND3-bridge:** prior strict changed-source chain and transported full tensor cone; detailed receipt and task history in [CURRENT_MODIFIED_TRANSPORT_COMMON_N_2026_10_06.md](CURRENT_MODIFIED_TRANSPORT_COMMON_N_2026_10_06.md).
- [ ] **LEFT4a — relaxed input before buffer:** source-replay the complete Rh_reference and O2_slope stress. For kappa<2 prove D>F*(2-kappa), not merely D>0; preserve the absolute pressure, kinetic history and radial contribution. Add the O2_axial full-Z relaxed branches and exact chart seams. Produce a scoped receipt without a strict/global flag.
- [ ] **LEFT4b — construct an actual left collar:** locate or build a source-bound strictly admissible inner collar, or a genuine zero-stress support boundary compatible with the paper's zero case. Do not treat Rh offset-5 or a shifted flat cutoff as this collar. Bind the chosen endpoint, radius, source values/derivatives and complete tensor before deciding modulation support.
- [ ] **LEFT4c — extend actual shear loop:** cover the defective O2_axial, buffer, slope and reference source intervals up to the accepted left collar. Derivatives must use physical y=logR on both phase and buffer selectors. Keep exact divergence through own-moment radial recovery, shared original pressure datum and all four current source joins. Use the already checked right O3 collar only after proving same-source composition.
- [ ] **LEFT4d — repair changed histories:** recompute modulation integrals, their full correlated five-moment defects, the independent terminal repair map and uniqueness/control bounds for the changed source family. Preserve absolute pressure and positive kinetic increments. Do not reuse old h_N vectors or thresholds after changing the defining source; retain the present receipt as its original-family certificate.
- [ ] **LEFT4e / BOUND3-global:** prove full tensor direction/shear/quadratic inequalities on the new complete support and both admissible outside collars; recompute one sufficient N across all new constraints. Add remaining inner region cones before promoting a global gate.
- [ ] **CONT / ENERGY / REC / WAVE / PHYS:** remaining physical interface composition; finite energy with actual physical Jacobian/tails; distinct n=1/n>=2 recursion and independent moment repair/smooth sum; two oscillatory/mean families and averaged stress cancellation; resolved u/v/w/p/f with corrected Cartesian residual and measured scale recursion/material winding.

current_original_O2_buffer_relaxed_input_cone_certified=true. Strict O2, global admissibility/common N/physical composition, energy, actual coefficient recursion, global flat remainder and full corrected NS remain false. Full long-term goal remains active.
