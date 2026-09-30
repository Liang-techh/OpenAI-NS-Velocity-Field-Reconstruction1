# Actual initial inner exit: Section 9.25

`ExitBridge` uses the same degree-18 `CorePolynomial`, source-scaled
schedule, dominant axis pressure and Section 9.23 comparison as
`lei_ren_part1_paper_exit_comparison.py`. The candidate has Lambda=1e36,
delta=1e-32 and h_b=.005. Its present domain is only
0 <= log(R/R_a) <= 2 h_b, with R_a=4/Lambda.

The joint ODE advances log(F/F_a), Uz, and five cumulative moments.
Separate axial and swirl quadratic integrals reconstruct `z_theta`
without imposing a terminal target. Pressure is the same axis pressure
plus the accumulated integral of F squared. The velocity equations are
Sections 9.25-9.26 of [the pinned source](https://arxiv.org/html/2609.35406v1).

Moment components have explicit amplitude/radius normalizations.
The axial equation uses I_z and F/barF directly, avoiding an enormous
intermediate I_z/barF. The positive collar multiplier is represented in
arbitrary-exponent arithmetic. Its real-domain amplitude margin is a
numerical policy, not a certified paper constant.

Run `python experiments/root_st073/lei_ren_part1_paper_exit_bridge.py`.
The JSON compares 16, 32 and 64 actual-ODE steps at the exact H0 root
and Z=.3, holding comparison discretization fixed. It records initial
matching, the finite-polynomial slope defect, signed-log moment increments
and normalized refinement errors. At Z=.3, the 32/64 log-amplitude
difference is 3.05e-11; the largest normalized moment difference is
3.85e-11. These quantify the initial exit only.

## Next required construction

1. Propagate actual Z derivatives alongside the five moments, including
   derivatives of the comparison stress driving the ODE. Recover Ur
   from the same Mz and Mz_Z rather than fitting it independently.
2. Evaluate actual Section 3 stress and connection cone with those jets.
   Retain endpoint slope defects caused by finite core truncation.
3. Extend the frozen comparison and actual exit toward R=100..110 and
   the reference matching radius with adequate step refinement.
4. Restore the reference power law and repair actual five-moment defects
   using Section 10. Preserve the corrections' pressure effects.
5. Bound angular-correction coefficient Z derivatives in future pressure.
   The baseline tail envelope does not close the corrected-tail bound.
6. Assemble/localize the complete field, then measure divergence, energy,
   stress/remainder and temporal scale laws on that same field.

Actual exit Z derivatives, its cone, full outer matching, global finite
energy, temporal scale recursion and oscillatory cancellation remain open.
