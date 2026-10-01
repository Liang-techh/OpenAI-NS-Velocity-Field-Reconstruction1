# C120-F12: callable C1 five-bump reference-annulus field

The certified implicit five-bump controls now enter an actual interval field adapter, rather than remaining a coefficient receipt. `IntervalRepairedReferenceField.evaluate_x(x)` accepts an exact rational `x=R/Rm` in `[1,2]` and encloses the entire fresh axial family `Z in [0.49,0.51]`. Controls and source defects retain their interval values and first axial derivatives. No midpoint controls are used.

## Field and cumulative data

The adapter uses the paper's three normalized, disjoint compact bumps at `5/4, 3/2, 7/4`, each of radius `1/40`. The partial map computes all five rows `L(x)h+Q(x,h)` with directed cumulative integrals. Completed axial mass supports use the exact normalized mass identity; other completed weights come from the frozen directed receipt, and clipped supports use directed Taylor cells and endpoint bounds. Point bump values and radial derivatives vanish exactly at support endpoints.

Writing `f=sum xi_i beta_i`, `g=c1 beta_1+c2 beta_3`, the fields are `Utheta=Am*(x^0.1+f)` and `Uz=4Z+g`. Physical cumulative moments are the reference moments plus the original five source defects and partial repair contributions. Pressure is always `P0+Mp` with the accepted analytic axis datum. Radial velocity and its radial derivative are recovered from the same axial moment and its first axial derivative. The radial moment identities are retained for divergence recovery; an independent Cartesian divergence calculation has not been completed.

At `x=2`, the actual partial map is the full map used in the uniform implicit inverse. Thus the five terminal reference identities follow for that implicit solution family. The receipt also checks zero containment in all five centered value/derivative rows. Zero containment alone is not the existence proof: that proof remains the strict self-map and contraction in the inverse receipt. This terminal statement is restricted to this reference-annulus map and accepted fixed-parameter source family, not a complete inner-to-heat-exterior closure.

## Evidence

- The partial-map companion compares completed weights with the parent inverse and checks clipped C1 rows against independent quadrature.
- An independent moderate-scale physical fixture integrates the actual corrected field minus the reference field. All ten physical moment value/first-derivative coefficients are enclosed. All six stress components, pressure value/derivative and recovered radial velocity contain the independent scalar results. Synthetic fixture controls exercise quadratic coupling and are not production controls.
- The production field is evaluated at `x=1,5/4,3/2,7/4,2`. All five radial points certify the relaxed `kappa<=2` cone on the full local axial family. These are point-radial certificates; the entire radial annulus has not been certified. Strong admissibility is not claimed.

Implementation and receipts are under `experiments/root_st073/`:

- `lei_ren_part1_paper_interval_partial_five_bump_map.py/.json`
- `lei_ren_part1_paper_interval_repaired_reference_field.py/.json`
- `lei_ren_part1_paper_interval_repaired_reference_field_check.py/.json`

Run the corresponding Python modules from the repository root. The production module loads the completed core and prior immutable certificates once, checks dependency hashes and writes a fresh field receipt.

## Next work

1. Enclose the complete radial `[1,2]` annulus, including nonzero radial bump slopes, and certify the appropriate cone branch on each full cell. Keep point certificates separate from cell certificates.
2. Materialize the actual reshape/restoration field between `R=110` and `Rm`; an enclosed terminal defect integral is not a callable connecting field.
3. Extend the implicit inverse to C2 and higher axial derivatives so `Ur_Z` and higher smoothness are available; cover the required whole-axis domain.
4. Implement the paper's shear modulation, flatten, heat collar and exact heat exterior while preserving the same analytic pressure and functional terminal moments.
5. After leading profile matching, implement genuine time-order dependent recovery and moment repairs, then oscillatory stress cancellation and independent full Cartesian NS residual validation.

Original parameter remainders, whole-profile matching, heat exterior, finite energy of the complete field, flat remainder, temporal recursion and oscillatory correction remain incomplete.
