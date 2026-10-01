# Collar first-width analytic linkage — 2026-09-30

For width-tied epsilon, the zeroth-width switch is chi0(s)=sigma(1-s) on [0,1], and zero for s>=1. The exact primitive is Jchi(s)=s-Jsigma(s) for s<=1 and 1/2 afterwards. At zero width all comparison drivers remain at the fixed core inlet. Thus the first-width collar ODE coefficient is available analytically, without a nonlinear collar integration.

At the accepted endpoint s=2, the normalized state coefficients divided by h_b are:

| State | First-width coefficient / h_b |
| --- | --- |
| g | -D0/4 |
| u | -sqrt(Ra/2) Iz0/2 |
| theta | 4 |
| mz | 2u0 |
| mixed | 4u0 |
| axial | 2u0^2 |
| swirl | 2 |
| p | 2 |

These formulas also determine the first and second Z derivatives by differentiating the same inlet data, with h_b held fixed.

## Integration evidence

The directed switch primitive uses 128 panels and Taylor order12. The existing eight-step RK4 first-width switch quadrature is mathematically exact at s=1 and s=2 by paired nodes and sigma(x)+sigma(1-x)=1. This is a statement about the mathematical first-width quadrature, excluding stored floating-point generation errors.

Interior first-width quadrature error upper bounds are 3.058924740713091e-10 at s=.25, 2.604926237978e-5 at s=.5, and 1.12681095069964e-4 at s=.75. Endpoint symmetry does not certify intermediate collar values or second-width coefficients.

## Actual stored-coefficient discrepancy

The trusted accepted core and R100 source caches were linked at Z=.3 without rerunning their expensive collar chain. All 240 comparisons (eight states, three axial slots, ten pressure powers) are recorded, including directed analytic replacement coefficients. The original caches are unchanged and identified by SHA256.

28 coefficient comparisons exceed the chosen 1e-250 relative generation-roundoff threshold. The largest relative discrepancy is 1.41642315466968342e99 in u, second Z derivative, pressure power9. These are extraordinarily small absolute high-pressure coefficients: a large relative discrepancy is not evidence of a large physical velocity error. Conversely, their small absolute sizes do not establish coefficientwise correctness. The receipt deliberately records `within_stored_generation_roundoff=false`.

The analytic endpoint coefficients bypass the nonlinear log/exp quotient evaluation used by the stored collar integration. A direct inlet diagnostic finds 27 nonzero pressure/Z entries in the finite-MP computation of F/F minus1 or log(F/F), whose mathematical values are identically zero. These artifacts demonstrate a concrete cancellation mechanism, but do not account quantitatively for every cached discrepancy. All entries are preserved in the receipt. An interval implementation preserving independent pressure and width atoms is the next route to enclosing that finite arithmetic. Its contribution to the full field still requires propagation; the current comparison does not enclose the original core-driver error.

## Replay and remaining requirements

Run `experiments/root_st073/lei_ren_part1_paper_collar_first_width_integrals.py` for the switch integrals (no local cache required). Run `experiments/root_st073/lei_ren_part1_paper_collar_first_width_cache_link.py --cache-dir work_paper_cache` for the optional local accepted-cache audit. Only use trusted project pickles. Committed JSON contains every resulting comparison, so readers can inspect the discrepancy without the local caches.

Next: connect directed finite atom arithmetic to the axial collar equations, use the analytic first-width terms as controlled coefficients, and enclose second-width integration and truncation errors. True source-driver errors, full nonlinear ODE errors, global moment-function closure, admissible stress, and temporal coefficient recursion remain open.
