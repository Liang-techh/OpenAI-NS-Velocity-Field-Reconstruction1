# Own-moment radial and absolute-pressure error jets for variable N

Checked implementation: [d446dc0e](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/d446dc0ecc8756be088225c162c18bd9b9dbc431). **BOUND1b1c and BOUND1b2c are implemented in the changed native O2/O3 and quiet repair charts.** The provider supplies the missing derivative/error inputs for subsequent signed tensor estimates. Completed stress cones and actual coefficient recursion remain open.

## APIs and scope

CurrentRecoveredErrorMajorants in experiments/root_st073/lei_ren_part1_paper_compliant_current_O3_recovered_error_majorants.py exposes:

- modulation(offset,N): actual shared logR offset in[-11,1], finite integer N>=1. Negative offsets cover the original O2 buffer, positive offsets the O3 slope. Returns current own cumulative history, radial and absolute-pressure error bounds and signed symmetric history enclosures.
- quiet_repair(N,logx=(0,1)): independent five-bump repair with N>=27,303,666, logx=y in[0,1]. Reuses the checked uniform implicit h-ball and exact terminal equation; never reuses saved N=10^12 coefficient boxes at another N.

Every returned error table has ordinary logR orders0..4 and axial orders0..5. The axial native domain is[-1,1]. The query retains the current mu/delta and original feedback family/implicit source/pressure datum. It reads accepted source packets and pure formulas; no ancestor graph is rebuilt.

Upstream initial scalar caps use the accepted all-prefix frequency_uniform_bounds. The kinetic contribution is bounded separately by K_kinetic/N^2, with its positive buffer lower bound retained after offset0; the swirl contribution has its own1/N cap. Pointwise cutoff derivative enclosures are intersected with the global source caps before bounding local profile derivatives. Before the left support all error rows are exactly zero; after the cutoff the cumulative histories remain transported even though local profile increments vanish.

Quiet q=1+y gives t=1+q=2+y, after the modulation support ends at1/2. The provider binds the actual callable scalar source map and converts the checked fixed A/R0 cumulative caps by f2^power*exp(-rate*y), f2=exp(-1-3mu/2). Local normalized bump jets are bounded using the maximum across disjoint supports, while primitive integrals continue summing their support contributions. Shrinking remaining strips and the same implicit terminal identity yield exact zero mixed error rows after the final bump. No incoming moments or control midpoints are substituted.

## Recovery equations retained

The five normalized history derivative bounds use their actual defining ODEs, including axial kinetic squares, old-swirl/increment cross terms and increment squares. C-linear histories m/h use k! axial caps; C-squared histories k/e/p use(k+1)!.

**Pressure is a cumulative value.** The correction delta_p/Pstar^2 equals its own cumulative Cp history. Only its first logR derivative equals Uold*deltaUtheta+(deltaUtheta)^2/2. The absolute pressure remains original_P0+original_p+delta_p. Replacing the value with that local derivative would incorrectly remove the post-cutoff pressure history; the checker explicitly rejects that behavior.

For the radial error, let L=1-delta*Z^2, C=1/(1+Z^2). Factoring Ad from the own-moment source gives:

`Q/Ad = rV*axial_scalar + rM*own_M_scalar`

`rV=2Z*C/L; rM=[-(1-delta)*Z*C-(1-Z^2)*C']/L`.

The actual0<=delta<1 guard supplies a positive denominator. With H=(1-delta*Z^2)^-1, the axial majorants satisfy H0=1/(1-delta) and Hk=[2delta*k*H(k-1)+delta*k(k-1)*H(k-2)]/(1-delta). Leibniz bounds include every polynomial/C/C'/inverse-L contribution. C derivatives through6 are supplied because radial axial order5 differentiates the own-M term once more.

The native radial source is Pstar*Ad*sqrt(R/2)*Qhat. Its ordinary logR rows include the(DlogR+1/2)^j shift exactly once. The returned radial log prefactor is logPstar+logAd+logR/2-log2/2; the pressure log prefactor is2(logPstar+logAd). Modulation logR=logRd+offset; quiet logR=logRw+1+y. These are source radius recipes and positive factored units, not a new Cartesian/time mapping or materialized giant radius.

## Verification and unresolved gates

The focused checker passes61 exact original source/history/radial-shift/C/inverse-L/quiet-unit identities. It independently executes the original IntervalTaylor increment_rows and shifted_rows on40 signed fixtures, including delta0,.03,.5,.9 and both axial closure endpoints, and compares8400 mixed history/radial/pressure rows. Nine current variable-N queries, six invalid guards, exact left/post-repair zero rows, surviving post-cutoff pressure value, shrinking near-exit recovered bounds and a separate positive kinetic mass are checked.917 working/index dependencies agree. Mathematical/source review used **GPT-5.6 Luna / max**, read-only, without constructors/tests/edits/descendants.

Files use prefix lei_ren_part1_paper_compliant_current_O3_recovered_error_majorants: .py, _check.py, .json, _check.json. Run producer and checker directly. These provide rigorous source error enclosures and implicit-family bounds; they do not resolve arbitrary-N signed field coefficients. No completed signed tensor, global physical interface composition, common cone N, energy, coefficient recursion or corrected NS gate is admitted.

## Next tasks for GitHub agents

- [x] **BOUND1b1c:** upstream all-prefix signed history enclosures and mixed logR4/Z5 error jets, retaining separate kineticN^-2 and swirlN^-1 contributions.
- [x] **BOUND1b2c:** own-M radial and own-Cp cumulative absolute-pressure mixed error jets in the current changed native charts, same pressure datum and single radial factor shift.
- [ ] **BOUND2a:** bind the actual completed O2/O3 tensor recovery assignments to these error jets. Keep all original and modified pressure/radial sectors and source correlations. Publish explicit diagonal, off-diagonal, divergence/residual error formulas rather than bounding only the primitive shear.
- [ ] **BOUND2b:** carry the accepted correlated primitive floor through both directional vector inequalities and alignment conditions on both flat tapers. Preserve local chi/chi_y/phase dependence where the source floor vanishes. Check endpoints as source identities; do not divide by a cutoff that is zero there.
- [ ] **BOUND2c:** bound completed signed tensor errors on the quiet repair band relative to the accepted a-2>=3mu/2 reserve. Include cumulative pressure values and their gradient, own-M radial derivatives, C versus C², all cross/kinetic terms, actual mu/delta and logPstar/logAd factors.
- [ ] **BOUND3 / COMMONN:** combine every scalar/repair/pressure/radial/diagonal/direction/alignment constraint into one sufficient finite integer N. Provide each inequality and its source domain. The repair27,303,666 and quiet primitive68,533,403 thresholds alone are insufficient.
- [ ] **CONT4f2b:** implement the reference_restore/restore_buffer physical interface providers from accepted germs and original radii. Complete all original and current-modified adjacent/internal providers, then compose them in one global physical API. Do not relabel the obsolete full33 numerical view.
- [ ] **ENERGY:** physical kinetic cross terms, volume Jacobians and radial/time tails, including the revised fields; terminal signed energy closure does not erase positive kinetic source.
- [ ] **REC:** actual n=1 and n>=2 coefficient recovery, common inner domain, independent per-order moment repairs, divergence-preserving truncation, controlled finite-order remainder and smooth summation. Local source inverse or changing N is not this recursion.
- [ ] **WAVE / PHYS / DYNAMICS:** oscillatory/mean correction and averaged stress cancellation, resolved physical u/v/w/p, independent corrected Cartesian NS residual and measured contraction/aspect ratio/material winding.

The full persistent goal remains active. The controller's last full runtime stage remains currentmodifiedphysicalvelocity. Preserve unrelated dirty files. Previous report: [CURRENT_UNIFORM_REPAIR_BOUNDS_2026_10_06.md](CURRENT_UNIFORM_REPAIR_BOUNDS_2026_10_06.md).
