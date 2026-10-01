# F24: fixed bump constants select a fresh inlet-tolerance core

Date: 2026-10-01. Fixed Section 10 bump constants have directed/analytic majorants C_A=344, C_Q=1081 and C_S=57749553. With the same-source Kp=17 and KN=18000, the selected hierarchy gives e_star approximately 1.39826803516848e-19 and j=eta_tol/8 approximately 1.74783504396060e-22. The old j=1e-14 is about 57,213,637 times larger. It does not satisfy this selected hierarchy. Tighter choices of constants could change the tolerance; no impossibility claim is made about other choices.

The preheat datum is reusable because its fourteen-atom outer definition has no j dependence. The core is not reusable: U0=4Z+j, sigma, the complex tube, the anchored primitive G and the Cstar definition change. F21/F22/F23 remain immutable old-j branches. Their R110 results have not been transferred to the new core.

## Fixed constants and their proof

All weights come from the directed fixed-bump integral receipt, not finite Gauss samples, old functional defects or old inverse controls. The bumps have radius 1/40 and centers 5/4,3/2,7/4. Use h=(c1,c2,xi1,xi2,xi3) and the normalized row order of (10.8).

- C_A: the maximum column sum gives the induced l1 matrix norm. For a fixed finite preconditioner R, directed q=||I-RL||_1<1 yields ||L^-1||_1<=||R||_1/(1-q). The outward integer majorant is 344. This does not reuse the older inverse's weighted row norm.
- C_Q: symmetric ordered-pair output bounds are fg/2, 40*gg, and (ff+ff_over_x)/2. Disjoint supports give exact zero cross terms. The maximum pair sum is enclosed by the integer 1081 in the paper's sum of C1_Z norms, uniformly for Pstar>=1.
- C_S: with N the enclosed raw bump normalizer, B=1/(r*N) and D=8/(r^2*N) bound beta and beta_x. The selected integer majorant ceil(1000*(1+B+2D)) dominates every displayed velocity, radial derivative, pressure and SQ-error coefficient in (10.6)/(10.12). It is intentionally conservative and independent of physical profile parameters.

For C_S*t<=.01, the corrected denominator is at least .99; |a-.8|<=4(2D+.1B)t and the zeta drift is <=2B*t. The norm in (10.6) is bounded by (2.1B+4D)t. Partial positive unit bump mass bounds the average correction by t in C1_Z. Pointwise |(1-delta)Zq+(1-Z^2)q_Z|<=||q||_1, so inherited average error contributes at most 2epsilon0 to W. The pressure correction is bounded by 3t Pstar^2. Expanding SQ gives error <=10epsilon0+(1+2*a_error_coefficient+12B)t, dominated by C_S.

The proof includes x in [1,e]: after x=2 the local bumps vanish, while the cumulative moment and pressure bounds persist. Its stated hypotheses remain inputs: delta<=.001, epsilon0<=1e-6, Pstar>=1, reference velocities before correction, and the inherited average/pressure bounds. The current delta/epsilon0 meet the elementary limits; full inherited average, pressure and defect admission has not been proved. Read-only mathematical review confirmed these distinctions.

## Fresh analytic core admission

The new tube is generated from the selected j, not obtained by shrinking an old numerical field:

    sigma=j/500; eta=sigma/40; h=eta/8.

Its eta is approximately 8.739175219803e-27. The common complex-axis denominator factors, L and pressure poles are excluded with directed bounds. G is anchored at the unique H0 root in [-j,0]; the symbolic logC uses the fresh scalar Gbar. A separate analytic_core_family_sha256 binds the bump constants, j and analytic core definitions in addition to the unchanged outer-pressure identity.

Both linear resolvent bounds have been regenerated from this tube. All twenty coupled nonlinear terms, including pressure and nonzero swirl, remain in the normalized majorant. Directed scaled size and Lipschitz gates pass with the existing logarithmic Lambda choice. No old contraction norm or old finite coefficient tensor is selected as the new solution.

Mixed derivatives of total order up to three for Phi and Uz on scaled radius [0,4.1] meet the 1e-12 target at radial degree 144; the maximum bound is approximately 6.84210593099188e-14. Raw Psi does not meet that target and is not confused with Uz, whose epsilon prefactor is retained. These are scaled radial derivatives, not physical-R derivative bounds or a full K certificate.

## Finite core now running in resumable batches

The fresh factory has generated radial orders 0 through 9 of 144 from 148 initial axial jets on Z in [0.49,0.51]. The first batch took approximately 39.904 seconds of recurrence compute. Its state and paired receipt are saved. Seed acceptance reproduced every fixed jet and checked j/sigma, the fresh tube, the primitive anchor, all fourteen pressure atoms, pressure units and the positive nonzero upper swirl enclosure.

This is incomplete finite radial core generation, not temporal n-dependent recursion. No fresh exit, five-moment closure, heat exterior, stress lift or time-dependent background has yet been derived from this branch.

## Reproduction

Under experiments/root_st073/, prefix lei_ren_part1_paper_, run in this order:

1. shared_bump_constants.py
2. shared_analytic_tube.py
3. shared_linear_resolvent.py
4. shared_commuting_resolvent.py
5. shared_core_majorant.py
6. shared_core_tail_admission.py
7. shared_source_binding.py
8. shared_interval_core.py --seconds 60 --max-steps 144 (repeat sequentially until complete)
9. shared_core_seed_check.py

Do not start two recurrence writers concurrently. The factory resumes only when all dependency hashes and the new family identity match. Its files are shared_interval_core_Z049_Z051.json and shared_interval_core_Z049_Z051_state.json. The paired receipt's completed_target and completed_radial_order are authoritative. Once new coefficients are complete, rerun seed acceptance, then build a fresh exit loader and propagate the new family into shared exit-parameter definitions and all downstream receipts. Do not instantiate the old LogarithmicCoreExit on the new state.

## Ordered remaining work

- [x] Directed fixed matrix l1 inverse majorant C_A.
- [x] Directed quadratic C1 majorant C_Q.
- [x] Analytic C_S majorant under explicitly stated hypotheses.
- [x] Derive t_star/e_star/eta_tol/j and identify old-j incompatibility.
- [x] Fresh j-dependent common analytic tube and Gbar.
- [x] Fresh resolvent bounds and all-20-term core contraction.
- [x] Fresh degree-144 scaled mixed-C3 tail admission and source binding.
- [x] Generate and accept the 148-jet same-pressure seed; start fresh coefficients.
- [ ] Complete finite radial orders 10..144; retain the state and per-step timing.
- [ ] Validate finite-plus-tail positivity and controlled core moments with a fresh core exit loader. Include analytic_core_family_sha256 in every downstream identity check.
- [ ] Rebind the shared h_b=epsilon_b family to the new majorant/core and selected constants; certify K1 and physical C3 K norms without selecting a convenient unrelated width.
- [ ] Certify supplied exit inputs and frozen Df/Ef inequalities, plus radius compatibility; full Section 9 admission must remain false until all pass.
- [ ] Regenerate comparison/actual exit/R110 checks for the new core. Old-j R110 positivity does not establish these conclusions.
- [ ] Build the long reshape and all five functional defects; solve the uniform five-bump inverse and restore pressure-compatible terminal identities.
- [ ] Recover exact heat exterior and global matching; construct admissible stress and flat remainder.
- [ ] Implement true temporal n-dependent recovery and moment repair, oscillatory stress correction and independent Cartesian residual validation.

The complete reconstruction objective remains active and incomplete.
