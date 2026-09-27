# Momentum cost reduction subject to resolved sampled moments

The k=11 enriched `(0,1,2)` axial family was optimized with 48-point
split radial quadrature. A minimum-coefficient-norm SLSQP initializer
reached its iteration limit; its output is only an initializer. The
subsequent full sampled momentum objective, constrained by four exact
quadratic moment equations, converged. The objective is a radial-weighted
mean of scaled squared residual plus 0.1 times its square, not a reported
three-dimensional volume L2. Its analytic gradient passed a directional
finite-difference check.

On the common 48-point grid, the previous candidate's full momentum
peak was 2.17839e7. The optimized peak is 1.38659e7, about 36.35% lower.
The minimum-coefficient-norm initializer instead increased the peak to
2.82285e7; smaller coefficients alone were not a useful momentum proxy.

| Independent radial order | Max sampled moment | Sampled momentum peak |
| --- | ---: | ---: |
| 96 | 5.99519e-6 | 1.39407e7 |
| 128 | 3.00166e-5 | 1.39285e7 |

The moment values are small at both independent orders but not monotone
in order, so these runs do not certify continuum error bounds. A held-out
radial grid at eta=-0.3,-0.1,0,0.1,0.3 gives peak 1.95827e7, versus
2.64073e6 for the baseline on the same points (7.42x). The 96-point
physical stress-cone replay still passes only 3/9 nodes. The candidate
is rejected, and neither full momentum target nor scale recursion is met.

Next cache all enriched mode responses on physical support nodes and
include the actual cone inequalities with resolved moments. Inspect
component-wise full residual costs before choosing further velocity or
pressure corrections. The OpenAI Section 8 mean/rank correction route
preserves moments while repairing remaining debt; this finite physical
fit is only a diagnostic analogue, not that coupled construction.
Keep the independent time-interval and recursive contraction requirements.

Reproduce: `python experiments/root_st073/midplane_quadratic_axial_cost.py`.
