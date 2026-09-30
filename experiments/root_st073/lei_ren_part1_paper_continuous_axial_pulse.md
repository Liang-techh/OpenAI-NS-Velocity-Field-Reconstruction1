# Shared continuous axial pulse integrals

The MP provider uses the same smooth startup/cutoff pulse for point values and analytic derivatives, weighted full-row integrals, and analytic plateau partial integrals and partial integrals in the resolved saddle window. Endpoint support values are exact by definition. Startup and off-window cutoff queries raise an error rather than replacing a nonzero contribution by zero.

For the actual shared candidate, 160/200-digit full-row log refinement is about 1e-133. Old float-centered inputs differ by about 6.68e-15 in logarithm. Positive omitted-piece bounds have relative logarithm approximately -1108.70. These bound omitted pieces, not numerical quadrature error.

The saddle-window primitive derivative has an independent two-step centered-difference diagnostic at mu=1e-12 recorded in JSON. This moderate-mu diagnostic does not establish a uniform bound over the candidate parameters. Full-row evaluation still uses the actual candidate mu.

Run `python experiments/root_st073/lei_ren_part1_paper_continuous_axial_pulse.py`.

Next implement startup/off-window cutoff partial integrals and shared incoming primitives/Z jets. The continuous_pulse_energy and continuous_axial_solve modules now provide the energy atom and a numerical continuous correction solve; then solve/install all components with the same continuous definitions. No global mean cancellation or finite-energy certificate is asserted.

## Cutoff partial coverage

Cutoff partial queries now use exponential endpoint quadrature on the first half, convex-phase endpoint quadrature before the saddle, truncated saddle quadrature near it, and the full central-window atom with bounded positive remainder after it. The omitted-piece bound follows gp<=11 or convexity of ku+1/u^2; it does not enclose quadrature error. The same pointwise pulse supplies the primitive derivative.

At mu=1e-12 the cutoff endpoint centered-difference relative errors are 3.19e-12 and 7.98e-13 for steps 1e-17 and 5e-18. Actual candidate receipts cover xi=10.25,10.75 and two near-endpoint positions; these are samples, not a uniform certificate. Startup weighted primitives (0<xi<=.02) still raise explicitly. Near/after the saddle a constant nominal central-window value does not mean the true primitive is constant; the retained remainder and nonzero derivative remain explicit.
