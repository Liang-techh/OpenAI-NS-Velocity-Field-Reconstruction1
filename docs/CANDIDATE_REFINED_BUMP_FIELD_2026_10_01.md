# Installed refined candidate matching field

The directed-weight inverse from C120-M7 is now installed into the same-source matching field. The implementation is `experiments/root_st073/lei_ren_part1_paper_candidate_refined_bump_field.py`. Its callable `RefinedCandidateBumpField().evaluate_x(x, Z='.3')` supplies F, Utheta, Uz, Ur, pressure, the fixed P0, all five cumulative moments, and available analytic first axial derivatives, for 1 <= x <= e at the supplied axial center.

The new partial backend `lei_ren_part1_paper_bump_partial_enclosures.py` encloses cumulative weighted integrals of the same normalized exponential bump family. Empty lower supports are exactly zero. Once a support is complete, it returns the exact saved full-weight interval, preventing a new full/partial normalization discrepancy. Clipped interior cells use the frozen directed Taylor derivative remainder; endpoint strips use positive analytic tail bounds.

## Reproduction and result

Run `python experiments/root_st073/lei_ren_part1_paper_bump_partial_enclosures_check.py` for representative clipped-support checks. Run `python experiments/root_st073/lei_ren_part1_paper_candidate_refined_bump_field.py` for the actual candidate field. The latter uses the refined inverse and validates its dependency hashes.

Eleven actual-candidate samples cover x = 1, the three support centers and both edges of each support, and x = 2. All sampled nominal stresses satisfy the connecting-region relaxed cone. P0 and its first axial derivative are unchanged. At x = 2, the partial cumulative map equals the refined full-support map atom-for-atom, for values and first axial derivatives. The receipt saves exact atoms of the velocities, pressure and five moments, plus the directed partial-weight bounds.

The existing directed finite-center terminal residual bounds remain in `lei_ren_part1_paper_candidate_bump_enclosed_refined_residual.json`: less than 5e-35 for normalized values and 1.3e-56 for first Z derivatives, with all ten intervals containing zero. The new field uses that same full-support map.

## Scope

Field values and stress values are nominal midpoint approximations of the directed integral data. The pointwise bump uses the midpoint normalization, while the saved intervals bound the paper-normalized integrals. The partial moments are midpoint evaluations with explicit uncertainty; they are not an exact smooth primitive certificate. No continuous high-order matching, independent Cartesian divergence, uniform support cone, source error or formal pressure/width remainder is certified.

The source remains the actual persisted Lambda120 center Z = 0.3 and rejects other axial coordinates. No center jet is extrapolated to the whole axis. Kappa remains near 0.8 in this region; final admissibility still requires the later shear-modification construction and its outer-profile inputs. Temporal n-dependent recursion is not implemented by this correction.

## Next work

Build a new-center/local-interval source factory from the exact same-source amplitude and analytic pressure, using the coupled functional core recursion. Recompute comparison and transition data at those centers or controlled intervals, rather than interpolating the existing center. Then lift the directed moment repair and stress analysis to axial functions. Retain the original accepted exterior/preheat datum and the common construction parameters.
