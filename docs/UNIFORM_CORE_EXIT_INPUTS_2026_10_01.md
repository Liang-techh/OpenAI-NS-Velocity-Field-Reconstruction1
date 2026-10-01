# F25: full-real-axis analytic core exit inputs

Date: 2026-10-01. The fresh F24 inlet-tolerance core now has whole-real-axis analytic bounds for the three core exit requirements (4.35)/(9.9). They use the new j-dependent analytic fixed point, not the local finite coefficient samples or the old-j R110 branch.

On Z in [-1,1] and scaled radius r=Lambda R in [0,4.1]:

| requirement | directed bound, rounded summary |
| --- | --- |
| normalized swirl Phi positive | Phi >= 0.26538107638889 |
| sum of the two C2 axial/average errors | <= 3.4956700879212e-22 |
| signed Ha*partial_Z log F | <= 2.8883589395426e-380 |

The C2 bound uses the sum of separate suprema of derivatives of order 0,1,2, without factorial discount. It is smaller than epsilon0 approximately 5.55525e-11. The signed logarithmic-gradient bound is smaller than 1/20. Exact interval endpoint tuples, not these rounded summaries, are authoritative.

## Uniform positivity and norm embeddings

For real Z, chi=H0^2/(H0^2+sigma^2) lies in [0,1]. The zero-epsilon leading model in (8.38)/(8.51) is

    Phi0(q)=sum_n (-q)^n/[n!(n+1)!], q=chi*r/2 in [0,2.05].

The odd cubic p3(q)=1-q/2+q^2/12-q^3/144 is a lower bound. Its derivative is negative, and the omitted alternating tail beginning at order four decreases. Thus its value at q=2.05 is a uniform lower bound. The even quadratic is an upper bound and <=1 on q<=6. The derivative series gives -1/2<=B'(q)<=-1/2+q/6<0 for q<3.

The full normalized core differs from the leading model by at most the fresh scaled_map_size_upper in the paper's sum X_h norm. The epsilon*beta term is retained in this correction. No finite coefficient tensor is used in the proof.

Dropping the positive denominator weight (n+1)^2 in (8.30) and differentiating the binomial generating function gives the conservative embedding

    E(i,k)=(i+k)!/[20^i h^k (k+1)^2 (1-r/20)^(i+k+1)].

It controls all mixed derivatives through total order three here. The correction value is bounded by its X_h norm times E(0,0), leaving the displayed positive Phi bound. Uz=4Z+j+epsilon*Psi retains its epsilon prefactor. Radial averaging Mz/R divides each coefficient by n+1, so the same positive embedding controls the average. Each C2 error is bounded by j plus epsilon*||Psi||_h times the sum of E(0,k), k=0,1,2.

## Signed logarithmic-gradient gate

The physical amplitude remains implicit and positive. Its derivative is fixed analytically:

    partial_Z log F = -Lambda*g_axis + Phi_Z/Phi,
    g_axis=L*H0/(H0^2+sigma^2),
    Ha=H0+epsilon*d*Psi.

The main term is -Lambda*L*chi. The model derivative contributes at most M*chi, with M=20*r/(2*Phi_lower), because |H0*chi_Z|<=40*chi and |B'|<=1/2. The potentially large cross term is unscaled: -d*Psi*g_axis. It is not incorrectly multiplied by epsilon or discarded. Since |g_axis|<=sqrt(chi)/sigma, it is bounded by B*sqrt(chi), where B=Psi_bound/sigma.

For A=Lambda*(1-delta_upper)-M>0, complete the square:

    -A*chi+B*sqrt(chi) <= B^2/(4A).

The remaining terms are bounded by H0_bound*DeltaPhi_Z/Phi_lower and epsilon*Psi_bound*(20*r/sigma+DeltaPhi_Z)/Phi_lower. The full Psi bound, including its nonlinear correction, is used. The Lambda hierarchy makes the resulting directed bound extremely small even at the narrow H0-root region. This proves the signed inequality required by (4.35); it does not claim an absolute bound on the logarithmic gradient. Read-only mathematical review confirmed the embedding, model/correction split, moment averaging and square-completion argument.

## Finite computation and prepared downstream code

The degree-144 finite recurrence is still running from the saved F24 seed; its live process has advanced beyond the last closed batch. Do not treat a growing state file as a terminal job or start a second recurrence writer. The paired finite receipt is written at batch completion, and its hash must match the state before consumption.

New companions under experiments/root_st073/, prefix lei_ren_part1_paper_:

- shared_core_uniform_bounds.py/.json: completed analytic bounds described above.
- shared_tolerance_core_exit.py: fresh core/tail/seed loader, requires completed_target=true and matching analytic_core_family_sha256.
- shared_tolerance_comparison_base.py: scale-cancelled comparison helpers bound to the new loader; the unrelated easy-width candidate integration is omitted.
- shared_tolerance_exit_parameters.py/.json: regenerated shared h_b=epsilon_b family, including the selected C_A/C_Q/C_S/e_star/j and new analytic core identity. K1 and full K remain formal.
- shared_tolerance_comparison.py, shared_tolerance_exit_bridge.py, shared_tolerance_exit_continuation.py, shared_tolerance_exit_switch.py, shared_tolerance_R110_cone.py, shared_tolerance_phase_check.py: prepared new-family exit chain. These downstream producers have not been run on an incomplete finite core; their numeric conclusions are pending.

After the recurrence is terminal and complete, run shared_core_seed_check.py, then shared_tolerance_core_exit.py. Validate its finite-plus-tail positivity before running the prepared exit chain in the listed order. Do not transfer the F23 old-j R110 cone to this branch.

## Remaining gates

- [x] Whole-real-axis analytic normalized swirl positivity on the required core interval.
- [x] Full-axis summed C2 axial-velocity and mean-moment closeness.
- [x] Signed Ha*partial_Z log F core input.
- [x] New core identity bound into the prepared exit loader and parameter family.
- [ ] Complete the degree-144 finite recurrence and inspect terminal receipt/state consistency; rerun seed acceptance.
- [ ] Run finite-plus-tail exit loader and new-family comparison/actual exit/R110 producers.
- [ ] Certify frozen-profile Df>0, (Df^2+Ef^2)/Df>=2+4gamma on Ra..110 and Df>=4 on 100..110. These are separate from the core inputs proved here.
- [ ] Construct a fixed Section 9 K1 coefficient ledger for interpolation, product/moment transfer, stress quotients and actual bridge estimates. The paper supplies existence, not a numeric value; unrelated kernel integrals named K1 must not be used.
- [ ] Certify every physical C3 term in K and radius compatibility, then admit the shared parameters. Normalized analytic derivatives are not a physical K certificate.
- [ ] Build long reshape and five functional defect repairs, exact heat exterior, admissible stress and flat remainder.
- [ ] Implement true temporal n-dependent recursion and oscillatory corrections; evaluate the complete Cartesian momentum residual separately.

The global analytic bounds do not build a whole-axis finite velocity evaluator, a global cone, an admissible stress lift, a heat-matched background or temporal recursion. Full Section 9 admission and the full reconstruction objective remain incomplete.
