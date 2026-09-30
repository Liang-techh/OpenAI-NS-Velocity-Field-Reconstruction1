# Joint pressure/width prescribed exit — 2026-09-30

## Construction advance

The actual Section 9.25–9.26 base exit ODE now carries both the analytic preheat pressure-tail parameter and the physical collar width as separate MP atoms. Integration uses the dimensionless coordinate s=y/h_b. The auxiliary Section 9.23 comparison uses the same core, analytic Z/mixed core jets, pressure and five-moment ODEs in this coordinate.

This removes the earlier single-atom loss of radius/velocity/moment increments when h_b=exp(-100-1e154). R=Ra*exp(s*W), W=h_b*width_parameter, is evaluated in the bivariate ring, not by adding a tiny MP change to a nominal radius. The cutoff multiplier also retains the actual source choice epsilon=h_b; this parameter relation is explicit in metadata.

`PressureWidthExitBridge` recovers angular and axial fields from the prescribed driver equations and carries all five cumulative moments in the same state. It does not separately fit velocities or set terminal moments to target values. P0 comes from the original analytic pressure jet, avoiding subtraction of an accumulated integral from total pressure.

## Evidence

- Sparse rectangular Taylor algebra retains pressure, width and mixed atoms, including exp(-1e28) and exp(-1e154). Independent nested Taylor checks for inverse/exp/log/sqrt pass near 1e-112 scaled error.
- Resolved family, pressure order 9 and width order 2: fields, pressure, radial slopes and five moments at the boundary/end match the original scalar ExitBridge with maximum scaled difference 4.11478349528e-18. The width is 1e-6 here so the truncated-width comparison is numerically resolvable.
- Actual source: delta=1e-200, logC=5e151, Lambda=1e36, core radial degree18, coherent waiting and complete preheat datum. The source collar width is exp(-100-1e154), distinct from schedule delta and angular flattening width.
- The actual endpoint retains nonzero pressure-tail and axial-response atoms, nonzero width-driven log(F/Fa), and nonzero first-width increments for all five cumulative moments. The five increments being nonzero is not terminal five-moment closure.

Reproduce `python experiments/root_st073/lei_ren_part1_paper_pressure_width_bridge_check.py`. Main receipt uses pressure order9/width order2. A separate order3 pressure pilot is retained for overlap replay; all common serialized atoms agree exactly. This is not a remainder enclosure and the pilot is not substituted for the main receipt.

## Remaining requirements

This is a finite bivariate base exit ODE. Pressure/width truncation, finite radial/Z jets and ODE errors are not enclosed. The collar constants remain numerical inputs, not certified source contraction constants. In particular:

1. Carry Z tangents of the normalized prescribed driver and the actual state, using the same atom data. The auxiliary comparison has analytic core/moment Z jets, but the prescribed bridge currently does not expose full moment Z tangents; it is not yet a divergence-certified joined velocity field.
2. Propagate these fields/tangents through switching, continuation, long reshape and inner moment repair.
3. Regenerate angular and axial corrections/heat inputs from the updated common source and close all five terminal identities as functions of Z.
4. Establish pressure compatibility, finite radial energy and controlled exact heat exterior before declaring a matched background.
5. Establish admissible cone margins and flat remainder, then genuine n-dependent coefficient recursion and oscillatory stress cancellation.

The Taylor parameters here organize pressure and collar width; they are **not** temporal coefficient-recursion orders. The full residual 1e-3 gate remains after realizable oscillatory correction.
