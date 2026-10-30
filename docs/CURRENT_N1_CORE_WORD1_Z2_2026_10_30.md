# Core-exit first-word Z2 Green enclosure

At the canonical core exit `rho=4`, this step encloses all six components of
the first coupled word's ordinary `Z2` derivative for real `Z in [-1,1]`.
It is a same-source, signed-sector **interval enclosure** of the core word,
not a point value, a complete repaired n1 field, a Navier–Stokes residual,
or a five-moment closure. The source intervals are the same saved actual
whole-core enclosures on `rho in [0,4.1]`; the new integration uses only
`rho in [0,4]`. They retain nonlinear-correction and factorial-tail
uncertainty and are not independent per-cell measurements.

## Source, word and integral

The new `n1_core_low_source_matrix_v2.py` forms the actual raw source `g`
ordinary Z0–Z2 and the original dynamic `C0/B1` matrix jets Z0–Z2 from the
same `Phi`, `U` and radial mean intervals. It retains the current-Z positive
`F0(Z)` and outer `Lambda` powers as formal factors; Lambda factors inside
relative-F0 derivative coefficients are enclosed numerically. It does not
use model-only fields or add a second F0 derivative connection.

For each outer radial cell, low-order W0 jets use the signed whole-core
source enclosure times the exact positive Green-kernel mass. For W0_Z3,
the eight already saved source cells are kept individually: the positive
mass from each cell is enclosed by `[0, total-kernel-mass upper]`, with
the negative K kernel handled separately. This is conservative and can be
wide; it does not assert that the cell contributions are independent point
values or use a midpoint. The construction includes the axis without a
division by `rho` or `x`.

The source integrated for W1_Z2 contains all six ordinary-Z product terms:

```
C0_ZZ W0 + 2 C0_Z W0_Z + C0 W0_ZZ
+ B1_ZZ W0_Z + 2 B1_Z W0_ZZ + B1 W0_ZZZ.
```

The actual core Green kernels then integrate these terms over eight outer
cells to `rho=4`. The result retains formal source factors and unmerged
cell contributions. It has 404, 528, 528, 24, 404 and 528 terms for `F1`,
`Uz1`, `K1`, `P1`, `d_x_F1` and `d_x_Uz1`, respectively. Most coefficient
intervals cross zero, so no sign or useful relative-accuracy conclusion is
claimed. In particular `P1` is a state component, not `m15`.

## Checks, failure record and scope

The first frozen v1 invocation failed during the independent synthetic
fixture because an interval was multiplied on the left of a Taylor-jet
object (`NotImplementedError`). No actual-source integration occurred in
that invocation. Its exact source files, entry, failure and hashes are
preserved in `experiments/root_st073/n1_core_word1_z2_v1_failure.json.gz`.
The distinct v2 implementation fixes operand ordering and was entered once.

The v2 actual-data run returned successfully. Its independent non-unit-F0
symbolic fixture passed 119 source, matrix, derivative and axis checks;
an independent polynomial feedback/Green fixture passed six component
checks. The runner also checks the raw source hashes, family, exact F0 and
Lambda log enclosures, each saved g_Z3 source-cell boundary, the matching
outer-cell coordinates and the raw matrix support. Full directed endpoints,
source manifest, entered record, fixture results, W0 cell enclosures and
W1 terms are in `experiments/root_st073/n1_core_word1_z2_v2_receipt.json.gz`.

This is the complete **core segment** through the seam only. It does not
integrate the first-collar feedback path, append any m15 contribution,
or remove the previously saved Cauchy bound, all-word tail, 212 cumulative
m15 terms, six signed inner/lower terms, reserved contribution or inherited
errors. A wide interval or overlap is not evidence of a functional endpoint
identity.

Next, bind this core-exit W1_Z2 enclosure as the no-reset inlet for the actual
first collar. Construct same-source first-collar W0 and raw matrix jets over
the collar cells, integrate the collar feedback with its exact `hb`-weighted
measure and physical scales, then propagate independent tail and history
errors to the signed m15 computation without duplicating any existing term.
