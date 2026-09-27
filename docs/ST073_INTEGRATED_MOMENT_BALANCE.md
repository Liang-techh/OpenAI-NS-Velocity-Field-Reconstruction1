# Stable integrated moment balance and local refinement

For an axisymmetric divergence-free field, write v=u_theta, w=u_z and
integrate over a fixed radius R while taking z and remaining-time tau
derivatives. With regular-axis boundary terms zero, the physical
momentum moments satisfy

```text
I_theta = -d_tau integral(r^2 v)
          + d_z integral(r^2 w v) + R^2 u_r(R) v(R)
          - nu [R^2 d_r v(R) - R v(R) + d_zz integral(r^2 v)]
I_z     = -d_tau integral(r w)
          + d_z integral(r w^2) + R u_r(R) w(R)
          + d_z integral(r p)
          - nu [R d_r w(R) + d_zz integral(r w)]
M_theta = -I_theta / R^2; M_z = -I_z / R.
```

The implementation keeps R fixed under differentiation, splits radial
quadrature at the moving bump supports, and uses an axial-scale step
for derivatives of integrated quantities. It avoids summing Cartesian
second-radial-derivative cancellation across the domain. Static and
time-dependent polynomial velocity/pressure fields pass analytic
identity checks, including physical-time versus remaining-time sign.
Regularity at the axis and axisymmetry remain assumptions of this helper;
it must not be used unchanged for nonaxisymmetric corrections.

For the previous k=13 local trajectory, independent radial order and
axial-step checks give moment maxima 3.23e-5 to 3.26e-5. At k=17 they
give 1.35e-3 to 1.42e-3, confirming a small actual remaining defect in
this numerical formulation rather than a reliably closed moment.

A two-direction outer-poloidal value correction followed by the
rank-two outer-swirl slope solve repairs the k=17 defect. Value changes
are -1.45351e-6 and -2.97845e-7; slope change norm is 6.63405e-7.

| Radial order | Axial step factor | Max refined sampled moment |
| --- | ---: | ---: |
| 48 | 0.002 | 3.89673e-5 |
| 96 (fit rule) | 0.002 | 1.50783e-9 |
| 96 | 0.001 | 7.96684e-5 |

The altered radial rule and halved axial step both remain below 1e-3,
while exposing the difference between the fitted 1e-9 value and an
independent numerical estimate. Direct stress replay still passes 9/9
inner cone nodes. No full momentum improvement or interval certificate
is inferred from these integrated-moment checks.

Next integrate the coefficient differential-algebraic continuation using
this moment formulation, retaining the axial algebraic constraints and
angular transport. Repeated independent local solves do not make a
single path; verify the resulting callable field between integration
nodes, preserve inner support separation, and evaluate full residuals.
The full momentum and volume-L2 target, finite-energy closure and scale
recursion remain unresolved. All candidates remain accepted=false.

Reproduce:
```text
python experiments/root_st073/midplane_integrated_moment_balance.py
python experiments/root_st073/midplane_integrated_moment_refine.py
```
