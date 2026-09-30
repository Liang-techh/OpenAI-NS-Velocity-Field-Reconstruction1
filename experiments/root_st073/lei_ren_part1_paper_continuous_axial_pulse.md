# Shared continuous axial pulse integrals

The MP provider uses the same smooth startup/cutoff pulse for point values and analytic derivatives, weighted full-row integrals, and analytic startup, plateau, and resolved-saddle partial integrals. Endpoint support values are exact by definition. A finite-band branch reports its omitted contribution or signed integration-by-parts correction explicitly; a source pulse that falls below active MP precision raises rather than being replaced by zero.

For the actual shared candidate, 160/200-digit full-row log refinement is about 1e-133. Old float-centered inputs differ by about 6.68e-15 in logarithm. Positive omitted-piece bounds have relative logarithm approximately -1108.70. These bound omitted pieces, not numerical quadrature error.

The saddle-window primitive derivative has an independent two-step centered-difference diagnostic at mu=1e-12 recorded in JSON. This moderate-mu diagnostic does not establish a uniform bound over the candidate parameters. Full-row evaluation still uses the actual candidate mu.

## Startup weighted primitive

For 0 < xi <= .02, `partial_row` uses the endpoint coordinate s = k (xi - v), with k = lambda / mu, and integrates by parts once:

`H = gp(xi)/lambda - mu/lambda^2 * integral exp(-s) sigma(50 (xi - s/k)) ds`.

The returned log primitive is `-k (13 - xi) + log(H)` over the retained band. If the band ends before k xi, the truncated sigma integral is subtracted, so the missing correction has sign -1. Monotonicity gives the absolute bound `|Delta J| <= exp(-k (13 - xi)) mu sigma(50 xi) exp(-band) / lambda^2`. The receipt records its log relative size as `log_relative_omitted_absolute_bound` and records `omitted_correction_sign = -1`; this bound does not enclose quadrature error.

The focused receipt samples xi = .005, .01, .015, .02 for both rows at mu = 1e-6 and at the actual candidate mu (about 4.78089e-28). All samples use the endpoint-scaled branch and retain 10^28-scale signed logs. With band 16, the sampled relative absolute correction logarithms range from about -20.23 to -24.52 at mu = 1e-6 and from about -69.34 to -73.61 at the actual candidate. A scaled centered derivative refinement at xi = .01 uses steps 0.05/k and 0.025/k; errors fall from about 4.17e-4 to 1.04e-4 for both the moderate and actual mu cases. This is a bounded diagnostic only: quadrature is uncertified and the samples do not establish a uniform startup result. For q = 50 xi <= .5, an analytic upper bound on the flat-endpoint integral forces an explicit error when the pulse is below active MP resolution; extremely tiny positive xi therefore raise rather than returning numerical zero.

Run `python experiments/root_st073/lei_ren_part1_paper_continuous_axial_pulse.py`.

The continuous_pulse_energy and continuous_axial_solve modules provide the energy atom and a numerical continuous correction solve; then solve/install all components with the same continuous definitions. No global mean cancellation or finite-energy certificate is asserted.

## Cutoff partial coverage

Cutoff partial queries now use exponential endpoint quadrature on the first half, convex-phase endpoint quadrature before the saddle, truncated saddle quadrature near it, and the full central-window atom with bounded positive remainder after it. The omitted-piece bound follows gp<=11 or convexity of ku+1/u^2; it does not enclose quadrature error. The same pointwise pulse supplies the primitive derivative.

At mu=1e-12 the cutoff endpoint centered-difference relative errors are 3.19e-12 and 7.98e-13 for steps 1e-17 and 5e-18. Actual candidate receipts cover xi=10.25,10.75 and two near-endpoint positions; these are samples, not a uniform certificate. Near/after the saddle a constant nominal central-window value does not mean the true primitive is constant; the retained remainder and nonzero derivative remain explicit.
