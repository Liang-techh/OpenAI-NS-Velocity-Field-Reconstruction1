# Long swirl reshaping and reference continuation

Run `python experiments/root_st073/lei_ren_part1_paper_long_reshape.py`.
Source: locally pinned arXiv:2609.35406v1, Section 9.5, equations 9.30–9.33.
Uses the shared candidate from shared_candidate_1.py and its actual R=110
velocity, pressure, five moments and Z jets. This is radial construction;
it does not establish temporal scale recursion or a global solution.

For y=log(R/110), T=400A, A=1e150, prescribe

```
log u = log u1 + y/10 - sigma(y/T) B
B = log u1 + logCstar + log(1+Z²)
V = v1
zeta = (1-sigma) zeta1 - sigma 2Z/(1+Z²)
a = .8 + 2 sigma'(y/T) B/T, b=0
```

`u` denotes angular velocity Utheta; F=u/sqrt(2R). The module transports
the actual axial primitive and mixed moment, rather than replacing them
with reference targets. ExitSwitches now exposes its separate actual
axial-square and swirl-square primitives to avoid subtracting away the
tiny swirl contribution when extracting it from z_theta.

In ell=log(R/current integration radius), define rho=u(y-ell)/u(y).
Shaping increments have endpoint-normalized integrands

```
Mtheta / (sqrt(2) R^(3/2) u): exp(-1.5 ell) rho
Mp / (u²/2): rho²
swirl-square primitive / (R u²/2): exp(-ell) rho²
```

These displayed ratios concern increments from R=110. All inherited
moments are added explicitly. Integrals use piecewise MP Gauss quadrature
over the final logarithmic units. If |B|<=2A and sup sigma'<=8,
d_y logu>=.06; using the weaker .05 bound gives omitted-tail envelopes
exp(-1.55 L)/1.55, 10 exp(-.1 L), and exp(-1.1 L)/1.1.
The local B/BZ guards do not establish a uniform mixed core norm.

The recorded Z tail bounds bound the omitted actual Z derivative divided
by the endpoint scale; they do not differentiate that scale. They multiply
the value bounds by max(|zeta1|,|2Z/(1+Z²)|), twice for quadratic integrals.
For normalized-integral derivatives one must also include the derivative
of the endpoint normalization. Finite quadrature error is separately
estimated by 16/32 order refinement, not certified by the tail estimate.

Past Rsh=110 exp(T), the exact reference angular power law continues to
Rz=exp(-8) Rref while V=v1. Moment increments are analytic with powers
1.6 (theta), .2 (pressure), 1.2 (swirl-square), and 1 (axial).
Pressure and Ur use the same accumulated moments and their Z derivatives.
Receipts for this region label inherited quadrature/tail quantities as
Rsh quantities; they are not bounds normalized at Rz.

At Z=.3, five probes (R=110, one logarithmic unit, shaping midpoint/end,
and Rz) pass the sampled relaxed cone. The stricter admissible predicate
remains false. R=110 field matching is retained and the terminal angular
reference log-value difference is about 5.59e-263. The largest recorded
16/32 normalized-integral difference is 6.57e-25. Endpoint normalized
shaping tail bounds are 6.46e-1241, 1e-79, and 9.10e-881; the pressure
Z-tail bound is 1.62e-43 under the stated local assumptions.

Next implement V restoration on [Rz,e Rz], integrating actual moments and
Z jets against the exact reference swirl. Then continue to Rh=exp(-5)Rref,
record five actual repair defects, and implement Section 10 repair before
matching the exterior. Full complex and mixed source constants, pressure
tail derivative bounds, finite energy, stress admissibility and temporal
scale recursion remain open.
