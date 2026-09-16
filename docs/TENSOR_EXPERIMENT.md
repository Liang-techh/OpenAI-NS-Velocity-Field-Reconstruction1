# Tensor correction experiment, stage 1

Current low-dimensional fits plateau above the target residual. Stage 1 adds
separable polynomial corrections to the existing streamfunction and swirl,
while keeping base geometry fixed. The correction basis uses s=R²/4,
v=Z²/4, w=2(t-0.25), with monomials s^i v^j w^k. The full representation
contains 27 coefficients per field; this first fit activates i+j<=2, k<=1,
24 coefficients in total. Coefficients are bounded in [-1,1]. Pressure and
the same two force coefficients remain fitted within their existing bounds.

Zero corrections must reproduce the saved base field. Spatial derivatives of
the streamfunction correction must be included in velocity and agree with the
legacy numerical curl. The optimizer normalizes initial energy per trial and
retains energy, flow-sign and core-drift losses. These penalties do not imply
hard acceptance: each condition still needs independent checks after fitting.

Configuration: `configs/constraints_tensor.json`. The finite domain, viscosity,
force family and validation tolerances are unchanged. New degrees of freedom
change the ansatz, not the definition of a successful NS candidate.

The existing seed-914027 holdout is a development comparison set: it has not
been fed to the training loss, but its outcomes have informed representation
choices. It must not be the sole final acceptance set. Once a candidate and
configuration are frozen, generate a fresh, independently seeded final audit
and retain the old results. Until then all results remain development evidence.

## Stage 1 actual result

TensorCandidate implemented by Luna/max and integrated with the existing validator. Zero correction, nonzero legacy-curl agreement, divergence, serialization and bounds checks are included. Parent integration run: 10 passed in 0.37 s. First fit used the full 2000-call budget and obtained training loss 0.142697. Independent sampled maximum residual is 2.464107, still above 0.001. The 21-time structure check also fails: core drift 0.071794 > 0.05 and minimum energy 0.099549 < 0.1. Therefore this result is not an accepted improvement and must not replace the previous candidate on acceptance grounds.

The soft structural penalties permit trading constraints against PDE loss; the next optimizer step must improve constraint enforcement while retaining the unchanged acceptance thresholds. All raw results are in `artifacts/constrained/tensor_stage1/`. This is evidence for further constrained fitting, not mathematical success.
