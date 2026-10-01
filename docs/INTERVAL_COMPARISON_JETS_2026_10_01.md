# Fresh interval-family comparison adapter

The explicit-center core for Z in [0.49,0.51] has advanced from degree 15 to 28 of 124. The maximum normalized mixed C3 analytic Taylor-tail upper bound is 230255768972826965691194877441343514900992292969874275.703552589757537765763959695687472891; its 1e-12 target remains unmet.

IntervalComparisonJets in experiments/root_st073/lei_ren_part1_paper_interval_comparison_jets.py consumes the exact checkpoint directly. It validates construction parameters, accepted pressure identity, dependency hashes and retained coefficient depths. It restores A and Uz coefficients divided by Lambda**n, the exact amplitude, physical axis pressure and amplitude derivative data. It reuses the established comparison algebra without loading the old Z=0.3 core or changing that adapter's center guard.

The evaluator returns the saved family of local axial Taylor jets, not a scalar slice or a Taylor series about the interval midpoint. Scalar Z queries are rejected. Normal use requires degree 124 and a passing analytic tail target; diagnostic_partial=True explicitly enables incomplete-core diagnostics. Neither mode encloses ODE discretization error or propagates radial tail errors through the comparison equations. Positive amplitude, coordinate denominator, inlet phi and radial slope denominators are checked before division.

The reproducible check computes the actual finite-core inlet at y=0, the switched comparison endpoint y=0.01 and the A/B driver jets using four RK4 steps at degree 28. It confirms rejection of incomplete-core normal loading and scalar-slice queries. Its receipt records the exact source checkpoint hash and finite interval outputs. Four steps serve as an interface diagnostic, not an accuracy certificate.

## Remaining tasks

1. Resume the production core to 124 and meet the analytic tail target; keep saved source identities fixed.
2. Independently assess family widths and generate narrower axial partitions if ratios remain too wide.
3. Compare independently generated scalar-center core/driver results against family enclosures; add pressure-row consistency evidence.
4. Propagate analytic tail and controlled ODE error through divisions and the switch. Current finite arithmetic intervals alone do not provide this guarantee.
5. Connect the new driver to exit continuation and transition equations with analytic axial jets, then calculate functional defects and five-moment repairs.
6. Extend controlled coverage to the full axial domain before promoting this to functional terminal closure.

No completed whole-axis transition, exact heat exterior, final stress cone, temporal n-dependent recursion or full Cartesian residual is claimed.