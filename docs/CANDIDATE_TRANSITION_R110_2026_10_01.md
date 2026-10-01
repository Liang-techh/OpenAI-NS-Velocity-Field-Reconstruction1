# Actual Lambda120 candidate exit through R = 110

The completed 124-order Lambda120 core now drives the actual exit bridge, post-collar continuation, and Section 9.4–9.5 switches through physical R = 110 at Z = 0.3. This replaces the old Lambda1e36 source for this candidate calculation. It does not replace the old stored center inverse/cone receipts or prove whole-axis matching.

Run `experiments/root_st073/lei_ren_part1_paper_candidate_transition_analytic.py` with Python. It uses 32 comparison steps, 16 initial exit steps, and 16 switch steps, with a separate 32-step switch refinement. The same nonzero amplitude, fixed pressure datum, five cumulative moments, axial tangents, and raw quadratic integrals propagate through the chain. Analytic Taylor drivers bypass both old axial finite-difference stencils.

The adapter `lei_ren_part1_paper_candidate_comparison_adapter.py` projects interval input jets to nominal MP values for existing transition classes. Its independent check compares the frozen power formulas against direct comparison evaluation at R = 1 and R = 100: driver relative discrepancy is at most approximately 2.32e-261, moment discrepancy approximately 1.05e-261, and pressure discrepancy zero at working precision.

At the four reported points, angular amplitude remains nonzero. Pointwise code diagnostics show:

| Physical R | Region | Relaxed condition | Stronger admissible condition |
| --- | --- | --- | --- |
| 1 | post-collar exit | passes | passes |
| 100 | exit endpoint | passes | passes |
| 101 | second switch | passes | does not pass |
| 110 | constant-power branch | passes | does not pass |

These are nominal center point tests using the existing Section 3 diagnostics. They do not enclose integration error, infinite core remainders, parameter errors, or the full axial/radial domain. Passing the relaxed condition in a connecting segment is not a certificate of a global admissible stress.

The receipt includes exact MP tuples for the refined R110 endpoint and initial-exit endpoint. These retain actual tiny amplitude and moment values; short printed logarithms cannot reconstruct this candidate amplitude accurately. Both exit and switches now expose the original P0 explicitly, together with P0_Z, for downstream pressure compatibility.

The root comparison precomputes radial products and caches core slopes; exit endpoint evaluation is also cached. This avoids repeated core moment convolution while continuing to use the same algebra and data.

Next: regenerate candidate-specific centered five-row defects, retaining their component atoms and axial tangents; attach controlled integration and radial errors; solve the actual five-bump repair; and extend these results to axial functions. Heat matching, uniform cone margins, temporal recursion and oscillatory correction remain open.
