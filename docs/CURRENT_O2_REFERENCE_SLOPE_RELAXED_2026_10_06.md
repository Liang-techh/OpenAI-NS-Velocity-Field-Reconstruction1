# Original Rh/reference and O2 slope: stronger relaxed input cone

Checked implementation: [6640f02c](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/6640f02c12050cafce4ffbeb54f809b12ad4405e). **LEFT4a-reference/slope is complete: the full original stress satisfies the paper's stronger relaxed cone on Rh_reference offset[-5,0] and O2_slope y[0,1], all Z[-1,1], including the reference/slope seam.** The source has bs=0 and kappa in[4/5,2]; strict admissibility remains false. The previously completed buffer relaxed input and changed-source strict cone/common-N certificates are retained. The axial turnoff relaxed input, actual strict inner collar/shear extension and the global target remain open.

## Actual source and correlated angular deficit

CurrentOriginalReferenceSlopeRelaxedCone in experiments/root_st073/lei_ren_part1_paper_compliant_current_O2_reference_slope_relaxed_cone.py consumes checked current original O2 tensor/source packets and selected original radius/delta data. No ancestor graph or fixed-N modified controls are reconstructed. query(chart,coordinate) accepts the two finite declared source domains. The source pressure datum, full energy, all five cumulative histories and completed four O2 source joins stay bound by current_O2_background_tensor.json.gz/_check.json.

The original native0 stress is replayed from the unchanged raw_pre_stress_rows operator, with Utheta/Pstar=U*C, C=1/(1+Z^2), m=4Z, h=U*X*C, k=4Z*U*X*C and Uz=4Z. The signed axial derivative Uz_y is zero on these two charts. The actual energy and absolute pressure remain in the complete axial stress; they do not enter the relaxed branch when bs=0.

On the reference branch X=5/8. On the slope branch U=exp(y/10-3J/5), A=8y/5-3J/5 and X=exp(-A)*(5/8+H), with the same defining angular mass H=integral exp(A). Direct source replay gives:

- (exp(A)*(1-X))'=3/5*(1-sigma)*exp(A).
- (exp(A)*(X-5/8))'=3/8*sigma*exp(A).

Since0<=sigma<=1 and A<=8y/5 on y[0,1], the actual angular deficit satisfies 1-X>=3/8*exp(-8/5)>0 and X>=5/8. This is a whole-source integral comparison, not an endpoint extrapolation or a fit. Public X enclosures intersect the independent mass enclosure with these proved correlated bounds, preserving the same defining function.

## Full theta units and the stronger paper criterion

Let L=1-delta*Z^2, Ageom=(1-delta/2)*C+(1-delta)*Z^2*C^2 and Bgeom=((2delta*Z^2-1)*C+2*(1-Z^2)*Z^2*C^2)/L. The complete inertial theta coefficient is

Theta=(Ageom*X-C)/L+4*C+4*Bgeom*X
     =(Ageom-C)/L+4*(C+Bgeom)+(-Ageom/L-4*Bgeom)*(1-X).

The exact shared whole-Z geometry makes the M term nonnegative and bounds the deficit coefficient below by cD=(15/8-14delta)/4. The initial shape is at least cX*Z^2-delta/(2*(1-delta)), with cX=(1-delta)/4. Thus the actual canonical inertial theta has a positive whole-source lower cD*(3/8*exp(-8/5))-delta/(2*(1-delta)). Its radial shear is retained separately, without adding a second C factor.

With F=Pstar*U*C/sqrt(2R), the complete stress satisfies Ttheta/F=R*Theta/C-a. The paper's kappa<=2 rule is Ttheta/F>2-kappa, with kappa=a here. Therefore its exact stronger reserve is R*Theta/C-2. The checked current radius R>=Rref*exp(-5) makes log(R*Theta/C)>1000, proving a strictly positive whole-domain reserve. Positivity of Ttheta alone would be insufficient for kappa<2; the stronger condition is explicitly implemented and checked.

The reusable relaxed_cone_margins(c,a,bs,Ttheta/F,Tz/F) implements both paper branches in signed normalized units. It checks D/F+kappa-2>0 for kappa<=2, and direction/quadratic conditions for kappa>2. Boxes crossing kappa2 use conservative sufficient bounds on both branches and cannot promote strict admission for the entire box. This helper is not a global field certificate.

## Evidence and ordered next work

Focused checker: 20 exact source/deficit/geometry/normalization identities, 6 strict continuous inequalities, 142 independent signed paper branch decisions (including rejected positive-theta inputs), 30 complete signed source stress-unit cases, 7 source queries, 3 seam equalities and 7 invalid guards. Source hashes pin the current unmodified full tensor and analytic datum. Read-only reviewer: GPT-5.6 Luna / max.

- [x] **LEFT4a-reference:** full relaxed input on Rh_reference[-5,0], actual reference histories, bs0, kappa4/5, retained radial shear and stronger stress lower.
- [x] **LEFT4a-slope:** full relaxed input on O2_slope[0,1], same positive angular deficit, ordinary logR source and exact reference seam.
- [x] **LEFT4a-buffer:** previous original buffer[-11,-2] relaxed input and exact nonzero kappa2 flat edge; see [CURRENT_O2_RELAXED_BUFFER_2026_10_06.md](CURRENT_O2_RELAXED_BUFFER_2026_10_06.md).
- [ ] **LEFT4a-axial:** replay the complete O2_axial stress with the real B_y, energy, absolute pressure and own-moment radial source. The ordinary coordinate is y=log(R/Rref)=exp(Md*phase), not phase. With a2 and bs potentially nonzero, prove relaxed direction at bs0 and the strict direction/quadratic branch where bs!=0. A useful exact reduction is Q/(2*Ttheta^2)=(1+bs^2/4)*(1-bs*Tz/Ttheta-bs^2/4); derive source-bound full-budget estimates before claiming admission. Include both source seams and the full Z domain.
- [ ] **LEFT4b:** construct or locate a source-bound strictly admissible inner collar or a genuine zero-stress boundary compatible with the paper. Existing Rh/slope/axial/buffer source charts do not provide the already-admissible left collar. A shifted flat cutoff alone cannot fix this.
- [ ] **LEFT4c-d:** extend the actual shear loop across every defective upstream chart, preserve all source joins/divergence/shared pressure, then recompute the new correlated five defects, implicit terminal repair/unique control bounds and finite-N envelopes for the changed source family.
- [ ] **LEFT4e / BOUND3-global:** combine all new support/collar stress conditions, remaining inner regions and constraints into a genuinely global sufficient N and admissible stress certificate. Keep the old checked family receipts valid only for their own source definitions.
- [ ] **CONT / ENERGY / REC / WAVE / PHYS:** remaining physical interfaces, finite energy with actual volume/tails, n=1/n>=2 recovery and independent moment repair/smooth sum, oscillatory/mean stress cancellation, corrected u/v/w/p/f and Cartesian residual/scale recursion/material winding.

current_original_Rh_O2_slope_relaxed_input_cone_certified=true. Original strict Rh/O2, global common N/admissibility/physical composition/energy/actual recursion/global flat remainder/full corrected NS remain false. Original registry counts are unchanged. Full long-term goal remains active.
