# Continuous axial correction solve

Run `python experiments/root_st073/lei_ren_part1_paper_continuous_axial_solve.py`.

This uses the actual seeded candidate incoming rows and energy target, the shared continuous bump matrix/energy Gram integral, MP full pulse rows, and continuous pulse Kp. The positive quadratic energy branch must remain in the source interval (.9,1.2). It does not replace the incoming data by synthetic targets.

ContinuousEndCorrection provides solved bump point values, first/second derivatives, and cumulative weighted primitives with the same basis. Primitive derivatives have an independent two-step diagnostic in JSON. Old coefficients replayed against the new continuous rows give relative defects about 6.68e-15. Re-solved row arithmetic is about 1e-173; energy algebra replay about 1e-201 at 200-digit working precision. These are algebra diagnostics, not full PDE residuals or integral error bounds.

The energy atom uses 100-digit arithmetic and order80 MP quadrature; observed order80/112 refinement is 2.76e-21 relative. Therefore the tiny algebra replay does not establish 173-digit physical moment closure. Incoming rows/energy target remain inherited, pulse off-window cutoff/startup primitives and Z derivatives remain incomplete, and the global velocity profile does not yet install these coefficients.

Next implement the remaining continuous incoming/pulse primitives and Z provenance, then install pointwise and cumulative definitions together in the joined exterior. Keep numerical terminal tails visible until actual continuous mean/Z-mean cancellation is established. Global finite energy, recursive transitions, full stress/remainder and oscillatory residual closure remain open.
