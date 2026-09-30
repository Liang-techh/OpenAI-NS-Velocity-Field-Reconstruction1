# Continuous axial coefficient tangents

ContinuousAxialAlgebra differentiates the SAME two linear row equations and quadratic energy constraint used by the continuous coefficient solve. For fixed matrix/pulse/Gram atoms and mu:

M c_Z + base_Z + a_Z pulse = 0

2 Kp a a_Z + 2 mu sum(K c c_Z) = target_Z.

The inverse matrix uses the analytic hyperbolic determinant. Let u_Z=-M^-1 base_Z and v=-M^-1 pulse. Then c_Z=u_Z+v a_Z, with a_Z obtained from the differentiated energy equation. The denominator is the positive-branch energy derivative; singular/nonpositive branches raise.

Run `python experiments/root_st073/lei_ren_part1_paper_continuous_axial_tangents.py` for an explicitly declared input direction. It is an algebra diagnostic, not actual Z data: the independent fourth-order difference errors refine 2.19e-17 to 1.37e-18. Derivative row arithmetic replays around 1e-201. A zero rounded energy replay is not a proof of exact physical closure.

Actual nominal inner/source Z data are separately supplied by seeded_input_tangents.py. Integral enclosure, source input uncertainty and derivatives of angular heat bounds remain open; this solver cannot certify them.
