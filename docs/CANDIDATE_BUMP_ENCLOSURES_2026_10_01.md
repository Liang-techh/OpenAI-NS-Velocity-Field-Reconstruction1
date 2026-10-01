# Directed paper-bump weights and refined candidate repair

The actual Section 10.2 exponential bump family now has a directed MP interval integral backend. `lei_ren_part1_paper_bump_integral_enclosures.py` encloses weighted single-bump and squared-bump integrals, including the normalization. Interior integrals use Taylor polynomials with cell-wide derivative remainders; positive endpoint strips use explicit exponential tail bounds. Products on disjoint supports are exactly zero.

## Reproducible runs

1. Run `python experiments/root_st073/lei_ren_part1_paper_bump_integral_enclosures_check.py`. The saved run uses 80 digits, Taylor remainder order 12, requested tolerance 1e-20 and 16 initial adaptive cells. It took approximately 75 seconds. Independent finite Gauss orders 128 and 256 fall within every saved weight enclosure; normalized masses contain one. These comparisons are diagnostics, while the interval Taylor and tail calculations provide the integral bounds.
2. Run `python experiments/root_st073/lei_ren_part1_paper_candidate_bump_enclosed_residual.py` to evaluate the original persisted candidate controls against the directed weights. It reuses exact interval endpoints from the weight receipt and validates source hashes.
3. Run `python experiments/root_st073/lei_ren_part1_paper_candidate_enclosed_weight_inverse.py` to re-solve the same source targets using the directed-weight midpoints. Ten updates preserve values and analytic first axial tangents. The script writes a new inverse receipt and its separate directed integral residual receipt; the old inverse is retained.

## Result and implication

The original finite-Gauss controls have a second-row residual in [1.1373475745e-27, 1.1373476687e-27] against the paper-normalized bump family. This interval excludes zero and has width about 9.41331e-35. Other rows except the normalized mass also exclude zero. Thus the previous 6.23766e-210 algebra residual did not represent the accuracy of the true bump integrals.

Re-solving shifts c1 by approximately +4.4469834632e-27 and c2 by the opposite amount; the three angular controls shift by approximately -9.06263e-49, +2.06621e-48 and -1.14700e-48. These are corrections of the numerical moment map, with the original pressure and source targets preserved.

For the new fixed controls, all five value residual intervals and all five first-Z residual intervals contain zero. Every value residual is bounded in absolute value by 5e-35; every first-Z residual by 1.3e-56. The separate midpoint-map polarization residual is approximately 6.23766e-210 and is not used as an integral error bound.

These are directed bounds on finite, nominal source targets and materialized controls at Z = 0.3. Inclusion of zero does not prove an exact infinite inverse, nor does it bound source reconstruction errors, discarded pressure/width series terms or the full axial functional problem. The refined controls have not yet been installed into the physical correction field. The previous 39-point relaxed-cone scan applies to the previous finite-Gauss field, not automatically to this new field.

## Next implementation work

- Install the refined controls with the same paper-normalized bump definition, and provide controlled partial cumulative integrals across each support. Do not combine new full-support weights with old cached partial quadrature without explicitly tracking that inconsistency.
- Recompute repaired fields, pressure and stress using those consistent partial integrals; retain the fixed P0 and its first derivative.
- Bound source/transition errors separately and lift the center repair to axial functions with uniform inverse and cone control.
- Continue outer matching and the Section 11 shear-modification inputs before claiming final admissibility. Temporal n-dependent recursion, exact heat-exterior matching and full corrected NS residual remain open.
