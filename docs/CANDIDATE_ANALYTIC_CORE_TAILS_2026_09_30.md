# Candidate analytic core positivity and infinite Taylor tails

The Lambda1e48 contraction candidate now supplies actual conditional analytic norm bounds, rather than tail factors with an unknown norm. `lei_ren_part1_paper_candidate_core_tail_budget.py` checks the candidate and its input hashes before applying the bounds.

## Whole-domain normalized angular positivity

The fixed point satisfies ||X-X0||_h<=0.146355076889422. In particular, each component's norm is bounded by its leading-model norm plus that correction. The angular difference Phi-Phi0 has zero axis value; its value series starts at radial degree1.

For scaled_R<=4.1, this gives

    |Phi-Phi0| <= correction_norm * (4.1/20)/[4(1-4.1/20)].

On the real axis chi is in[0,1]. The leading Bessel model has q=chi scaled_R/2<=2.05. Its alternating series gives the cubic lower bound1-q/2+q^2/12-q^3/144 and an upper bound1. Combining these bounds proves

    Phi = F/F0 >= 0.25594623652652049

on the entire rectangle scaled_R in[0,4.1], Z in[-1,1], relative to the accepted stored data. Since F0 is a positive exponential on the real axis, this also proves positivity of the candidate analytic angular profile there. It does not establish the near-root radial derivative sign or a stress-cone condition.

## Infinite Taylor tail budget

The norm bounds are approximately

    ||Phi||_h <= 47312.54852459592,
    ||Psi||_h <= 1.917518615499938e15.

The existing directed analytic tail formula is multiplied by these norm enclosures. The resulting receipt covers every mixed derivative of total order at most3 in(scaled_R,Z). It compares the exact analytic solution with its exact radial Taylor truncation; coefficients from an unverified finite adapter are not automatically identified with that truncation.

A search through the proved majorant finds that radial degree124 is sufficient for all these normalized tail bounds to be<=1e-12 throughout the rectangle. The corresponding coupled recurrence needs at least128 axis Taylor coefficients to retain three axial derivatives. Degree20 remains sufficient for a small value tail but not this global mixed derivative budget.

The1e-12 threshold is an internal normalized core truncation target. It is not the full corrected Navier–Stokes residual target. Physical radial derivatives acquire factors Lambda^i; the axial velocity correction additionally contains epsilon, while angular velocity includes F0 and its derivatives. These conversions and adapter/source errors must be propagated separately.

## Next integration

1. Generate higher directed Taylor jets of the same full physical pressure from true fixed-beta mass intervals and a controlled flatten mixture, with no pressure-parameter truncation.
2. Generate candidate amplitude-factored coupled coefficients through degree124, keeping the same P0 and all F0-squared nonlinear terms. Use resumable computation and exact interval state.
3. Compare finite coefficients with the exact recurrence, then combine coefficient/input errors with the infinite Taylor tails. This is required before reporting a computed field error enclosure.
4. Restore pressure and all five moments from the same candidate coefficients, regenerate its inlet and collar, and continue terminal functional matching.

Original parameter derivation, finite-adapter errors, candidate core-to-collar matching, near-root derivative sign, full stress cone, temporal recursion and corrected residual validation remain open.

Reproduce:

    python experiments/root_st073/lei_ren_part1_paper_candidate_core_tail_budget.py
