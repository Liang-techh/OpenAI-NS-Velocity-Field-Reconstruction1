# Whole-axis analytic core moment envelopes

The Lambda = 1e120 candidate now supplies ordinary axial jets through order 3 for all five cumulative core moments at scaled radius s = 4, uniformly over real Z in [-1,1]. The prescribed amplitude and its derivatives, quadratic swirl coupling, and accepted analytic pressure datum are retained.

Implementation: `experiments/root_st073/lei_ren_part1_paper_candidate_global_core_moments.py`; generated receipt has the same basename and `.json` suffix. Run with Python from that directory.

The Xh coefficient inequality is summed over radial indices n >= 1 using the binomial generating function. Phi's constant value uses the existing global positivity bounds; positive radial integration weights then give integral envelopes without quadrature. All 20 existing center analytic moment coefficients lie inside these global envelopes. This is a consistency check; the whole-axis justification is the analytic norm bound, not center sampling.

These conservative envelopes are functional input bounds, not a finite whole-axis core atlas, terminal five-moment closure, or a transition certificate. Original construction parameter errors remain outside the receipt's scope.

Next dependency: replace the old transition core data and supply analytic driver derivatives. The old ExitTangents and ExitSwitches use an axial step of 1e-45, far larger than the Lambda120 amplitude variation scale. Merely decreasing the step would not prove a derivative error bound. The preferred path propagates axial Taylor jets through the comparison equations and uses their coefficients directly.
