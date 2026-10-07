# Modified O2 open taper signed cone, with exact degenerate left edge

Checked implementation: [f60f89f8](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/f60f89f8747166009d6e8daac3479e5842c88e19). **The complete modified source satisfies the strict signed two-vector stress cone on the entire O2 modulation taper offset(-2,0], all Z[-1,1], all phases and every finite integer N>=22.** Combined with the preceding closed O3 result, both modulation tapers and the old seam now have complete directional/alignment estimates. The left support endpoint offset=-2 is explicitly degenerate: a=2, bs=0 and vs-2=0. It is not labeled strict. Earlier O2 buffer, quiet repair, global/common-N/energy/actual-recursion/full-NS gates remain open.

## Source-faithful negative buffer baseline

The actual O2 buffer coordinate is v in[0,11]; the modulation phase coordinate is t=v-11. On t[-2,0], the post-turnoff cutoff is zero but the accumulated incoming histories are retained. The producer replays the original pre_pulse_mixed_C4 axial source assignments, their first-order ODEs and their exact backward semigroup. In canonical C(Z)=(1+Z^2)^(-1) units:

- U=Ua*exp(-t/2), D=da*exp(-t), M=Ma*exp(-t).
- h/U=1-D and k/U=M-4D retain the exact common moment correlation.
- e/U^2=Z^2*AZ+C^2*AQ, where AZ=EZa/Ua^2 and AQ=EQa/Ua^2-t/2. The actual EQa may be negative; both energy coefficients remain signed.
- The complete absolute pressure is transported by Pabs_y=U^2*C^2/2 with the same original axis datum. No pressure correction is added.

The source theorem replays the same O2 buffer_offset=11 parent in slope_mu, all five seam histories, the unchanged absolute pressure, and the checked canonical Ua=exp(logAd) endpoint identity. It derives the full original O2 stress directly from the original raw operator. It does not extend a positive-offset O3 interval certificate to negative offsets without a new source identity.

Write Pabs=-U^2*C^2*Hfuture/2+Pmemory(Z). The actual O3 source assignments bind the same constant memory to Ua^2*exp(-1-mu)*(Pv+C^2/(2*(1+2mu)))*exp(-(1+2mu)*(13/mu+Tw)). The checked Hfuture(0) lies in[1/(1+2mu),1]; Hfuture(t)=1+(Hfuture(0)-1)*exp(t) is a convex combination for t<=0. This bounds the future-integral part, not Pabs/U^2. The original signed Pv and its derivative remain in the memory cap.

## Whole signed inequalities and finite-N transfer

The complete theta source splits into its axial shape, the common positive deficit and nonnegative meridional term. Exact rational geometry gives

Theta0 >= cX*Z^2+cD*da*exp(-t)+reserve,

where cD=(15/8-14delta)/8, cX=(1-delta)/4 and reserve=cD*da-delta/(2*(1-delta))-radial_cap>0. The radial cap uses the true minimum radius Rd*exp(-2). U(t)>=Ua gives the published absolute theta log lower. The source energy derivative is replayed in full, including all derivatives of C; AQ is bounded by |EQa/Ua^2|+1 on precisely this interval. The exact ratio U(t)^2/D(t)=Ua^2/da removes the negative time factor from the energy-direction estimate. The same pressure memory is constant and its checked log bound is rescaled by the new theta floor.

All 14 complete native stress error pieces at N=22 lie below exp(-1010) relative to the O2 baseline lower. Each component sum is at most exp(-1000). The checked weighted algebra and native order-zero inverse-N monotonicity apply to this actual negative source as well: no growing higher local derivatives enter. Thus every N>=22 is covered by inequalities, not frequency samples.

The weighted axial ratio is below exp(-390), D/Ttheta>0.9 and Q/Ttheta^2>1.9, uniformly on the closed direction/alignment interval[-2,0]. The actual primitive excess satisfies vs-2>=mu*chi(t)^2/4, which is strictly positive at every t>-2. There is no positive constant primitive floor over the entire open taper; the cutoff becomes arbitrarily small near -2. Exact closure at -2 remains distinct from strict admission. All 14 stress modification pieces are exactly zero at that flat edge, so the nonzero original stress is retained.

## API and evidence

CurrentModifiedO2TaperCone in experiments/root_st073/lei_ren_part1_paper_compliant_current_O2_modified_taper_cone.py consumes checked packets and pure source programs. query(offset,N) accepts finite subsets of[-2,0] and integer N>=22. A box containing -2 returns complete_modified_signed_cone_certified_for_entire_query_box=false, while retaining the closed direction/alignment estimates. Boxes strictly to its right return true. The scoped module gate is current_modified_O2_open_taper_signed_two_vector_cone_certified; current_modified_O2_whole_cone_certified remains false. No fixed-N control midpoint or ancestor constructor is used.

The checker passes 58 exact O2/source/seam/weighted identities, 28 strict baseline/transfer inequalities, 30 independent negative-buffer signed source fixtures with nonzero functional pressure memory, 64 original signed-cone success/failure fixtures, 14 exact flat-edge stress error zeros, 28 larger-N comparisons, 6 current queries and 7 invalid guards. Numerical original-zero fixtures use only a working-precision subtraction roundoff floor; production bounds are unchanged. Read-only source review: GPT-5.6 Luna / max. Producer/checker .py and .json/_check.json share the prefix above.

## Remaining tasks for GitHub agents

- [x] **BOUND2b-left1:** source-exact O2 negative buffer baseline with full signed energy, absolute pressure and incoming histories.
- [x] **BOUND2b-left2:** full modified left taper direction/alignment and primitive positivity, all phases/Z and N>=22 on(-2,0].
- [x] **BOUND2b-left3a:** exact flat-edge modification zeros and explicit non-strict status at -2.
- [ ] **BOUND2b-left3b/buffer:** reconcile the a=2/bs=0/vs=2 support edge and unchanged earlier O2 buffer[-11,-2] with the paper's applicable cone/closure rule and adjacent source joins. Strict vs-2 positivity cannot be asserted there. Any full[-11,0] energy cap needs +11/2 instead of +1.
- [ ] **BOUND2c1:** derive quiet bs=+2*Uz_y/Utheta from original source assignments and bind this signed convention; the earlier repair absolute majorant's negative b label must not be reused as signed bs.
- [ ] **BOUND2c2:** use local signed theta/axial baseline bounds on each disjoint quiet bump and gap, retaining cumulative/complement histories, actual pressure memory, radial cross terms and the unique finite-N implicit control family.
- [ ] **BOUND2c3:** prove every quiet D/Q condition and entry/internal/terminal join for all N above an explicit threshold. Exact final-source equality does not prove its preceding strip.
- [ ] **BOUND3 / COMMONN:** combine all regional completed cone inequalities with independent repair, history, pressure/radial and diagonal constraints. Modulation N22, repair27,303,666 and quiet primitive68,533,403 are distinct partial thresholds.
- [ ] **CONT4f2b and successors:** remaining accepted physical interfaces and global composition, using the correct lambda convention and preserving obsolete full33 views as obsolete.
- [ ] **ENERGY / REC / WAVE / PHYS:** finite physical energy, actual n-dependent recovery/repair/summation, two oscillatory/mean correction families, averaged stress cancellation, resolved u/v/w/p and full corrected Cartesian NS/dynamics.

Mark only each completed scoped task; publish checked source hashes and remaining gates. Preserve unrelated dirty files and historical sections. Previous report: [CURRENT_MODIFIED_O3_TRANSITION_CONE_2026_10_06.md](CURRENT_MODIFIED_O3_TRANSITION_CONE_2026_10_06.md). Full goal remains active.
