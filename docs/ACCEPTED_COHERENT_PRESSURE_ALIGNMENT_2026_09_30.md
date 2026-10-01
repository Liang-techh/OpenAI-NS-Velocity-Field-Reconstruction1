# Accepted coherent pressure alignment and affine flatten bounds — 2026-09-30

## Accepted source alignment

The pressure-error helper and the accepted continuous-waiting core differ only in the waiting, heat-collar, and exterior-power-tail components. The other eleven components have identical stored masses. A fresh coherent source-adapter construction reproduces every accepted component mass in the existing second-axial core snapshot exactly. No core coefficients, collar continuation, pressure anchor, or accepted snapshot is modified.

The schedule manifest stores exact input strings, hashes, and component masses with exact MP tuples. All schedule inputs except waiting length agree numerically. Before the waiting stage, the pressure formulas therefore use the same stored parameters and switches. Their true-integral enclosures transfer. Waiting is enclosed again with its accepted length, and the collar/exterior integrals are recomputed from the accepted normalization data.

The resulting complete accepted pressure budget has weighted normalized uniform C2 upper 3.079264485429091e-12 on |Z|<=0.8. Its derivative bounds through order24 and pressure perturbation bounds through the degree18 core/five core moments are rebuilt from this accepted source. The uniform core-exit bounds remain Uz <=4.394371482111663e-36 and Ur <=4.516471968020449e-53. These are pressure-driven finite-core errors, not total core approximation error.

## Sharper flatten coefficient enclosure

The former endpoint Riemann coefficient estimate gave a relative uniform pressure-value error upper about0.8592. The affine method integrates exp(2ell) exactly on every normalized u cell using expm1, then bounds the remaining factor 2^(beta-3)*(beta)_k/k! by its monotone endpoint beta intervals. It retains the exact stored length100 and slope -1/2-mu, and verifies that the relevant switches are saturated before using that affine slope.

With K96, 512 panels, and the stored192-node finite adapter, the accepted-source error upper bounds, divided by a directed lower bound on the true Z0 mass, are:

| Quantity | Uniform relative-to-mass upper |
| --- | ---: |
| Pressure value | 3.93004537174433e-16 |
| First Z derivative | 2.518273828566633e-14 |
| Second Z derivative | 5.390000058904187e-12 |

These are uniform bounds on |Z|<=0.8. They are normalized to the Z0 mass, not the local value of a derivative. The first derivative vanishes at zero, so no relative-to-local-derivative statement is intended. Coefficient errors and analytic derivative tails are enclosed separately. This does not certify the later relative-flat remainder conditions.

## Reproduction and canonical input

Use `lei_ren_part1_paper_coherent_pressure_error_transfer.py` to rebuild accepted fixed-beta, complete pressure, high-derivative, and uniform core receipts. It reads the committed alignment manifest, reconstructs its exact schedule inputs, and verifies all finite masses against that manifest; no local pickle is needed for this transfer/replay. Its `accepted_profile()` is the canonical source helper for subsequent error propagation.

The source-alignment audit script separately compares a fresh source to the existing local `second_axial_core_Z03.pkl`. That cache is an optional local audit artifact, not required to replay the committed accepted pressure transfer. It never rebuilds the expensive radial/collar chain.

Use `lei_ren_part1_paper_affine_flatten_coefficients_fixture.py` for independent quadrature, endpoint-beta conventions, analytic tail coverage, and panel contraction. Use `lei_ren_part1_paper_affine_flatten_coefficients_check.py` for the accepted production source and exact-tuple receipt. The sharper flatten receipt is separate; the high-order pressure/core propagation retains the conservative valid positive-mass bound for flatten, rather than claiming unavailable higher derivative coefficient precision.

## Next required work

- [x] Identify source-helper versus accepted-source differences.
- [x] Replay accepted component masses against the accepted core snapshot.
- [x] Rebuild the three changed terminal-stage pressure error bounds.
- [x] Transfer all14 pressure error bounds and order24 derivative inputs to the accepted source.
- [x] Recompute uniform pressure-to-core/five-core-moment propagation on that source.
- [x] Tighten uniform flatten coefficient error with exact affine cell integration.
- [ ] Enclose the infinite core radial-series remainder and analytic continuation.
- [ ] Enclose collar RK continuation and R100/R110 transition integrals on an axial interval.
- [ ] Assemble full five-defect C2 function error bounds and the five-bump inverse tail.
- [ ] Validate exact heat velocity exterior/finite energy and the full admissible stress cone.
- [ ] Implement the true temporal coefficient recursion and oscillatory correction.

Original parameter derivation errors, runtime evaluation errors of arbitrary caller code, complex analytic norms, total field approximation error, and full corrected PDE residual remain outside these receipts.
