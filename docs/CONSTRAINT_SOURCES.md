# CR001: source constraints and preregistered experiment

Configuration: `configs/constraints.json`. Status: preregistered, **not a successful candidate**.

Primary source inspected on 2026-09-16: [paper PDF](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf), printed pp. 1, 3–4 (Theorem 1.1, equation 1.1, section 2). This inspection establishes what the document states, not independent verification of its theorem.

| Item | Source statement / scope | Experiment treatment |
| --- | --- | --- |
| Incompressible viscous momentum equation | Eq. 1.1, whole space | Same differential equation; nu=0.01 is our choice |
| Compact spatial support and smooth force | Theorem 1.1 | Retain smooth compact representation; restrict force to two scalar coefficients of a fixed curl basis |
| Axisymmetric swirling inner flow | Section 2, leading core only | Chosen for the entire candidate family, an additional design restriction |
| Radial and axial contraction | Section 2, asymptotic description | Test prescribed scales on a finite window; not an asymptotic claim |
| Inward core and axial outflow | Section 2, central core | Sign probes and scaled-profile diagnostics; no claim every point flows inward |
| Zero initial data and singular endpoint | Theorem 1.1 | Not imposed on this finite-window experiment; not claimed as reproduced |
| Energy | Theorem 1.1 bounded kinetic energy | Finite-window energy bounds and quadrature only; values 0.1–10 are autonomous |

The window [0.25,0.75], support cylinder r<2 and |z|<2, unit reference energy, pressure support, force basis, seeds, budgets and thresholds are **autonomous choices**, not recovered source parameters. The restriction on pressure is substantial: an admissible velocity need not have compactly supported compatible pressure. Failure is possible and must be reported, not hidden by choosing the force equal to its residual.

The force potential is defined before candidate optimization. Its curl supplies a smooth axisymmetric strain/swirl input with only two adjustable coefficients. Their bounded optimization is allowed; introducing pointwise force parameters or fitting a new force basis to a completed residual is not. The bump is smooth across its outer edge, and its polynomial arguments avoid square-root singularities on the axis. Symbolic checks and force norms remain to be implemented.

All quantities are dimensionless. Report the unnormalized residual as well as reference-scaled values (reference scale equals one); do not divide by the candidate's shrinking amplitude. Spatial L2 norms include volume weights and are reported separately at each validation time. Initial energy normalization removes the zero optimum. The three derivative spacings hold candidate and samples fixed; three quadrature orders separately hold derivatives and parameters fixed. Boundary and axis strata supplement uniform sampling. A sampled maximum is explicitly a sampled maximum.

Preregistration is complete at the problem level; candidate coefficient bounds are supplied by CR002 and must be frozen before CR005. No optimization has been run against these thresholds. If feasibility fails, preserve the failure and issue a new version with a documented reason; do not mark this experiment passed.
