# Implicit time integration does not rescue the transferred three-template ansatz

`midplane_wave_implicit_template.py` restricts the evolving exact-curl
potential to the two `k=11 → k=19` transferred harmonic slopes and
their transferred mean slope. At each implicit-midpoint state it
recomputes the **full nonlinear frozen-state residual** and solves for
the three pulse-time derivatives together with compact harmonic/mean
pressure coefficients. A three-variable root solve closes each
midpoint equation. The direct full momentum check uses nine held-out
radial/axial nodes and eight angles shifted relative to the fit angles.
The physical time interval is pulse fractions `[-0.1,+0.1]`, split in
two.

| Spatial training nodes | Implicit equation defect, intervals 1/2 | Corrected / frozen held-out midpoint max, intervals 1/2 |
| ---: | :--- | :--- |
| 4 | `3.50e-12`, `1.60e-6` | `52.1`, `1070.6` |
| 9 | `2.98e-11`, `8.16e-8` | `12.8`, `315.6` |

With nine training nodes, the relative fit residual at the **training**
points is already `0.191` and `0.232`; the three endpoint template
coordinates move from zero to approximately `(-1.21,-2.43,4.99)` and
then `(25.7,99.9,-647)`. Thus a solved implicit time equation does not
imply a controlled physical momentum field. More spatial nodes reduce
overfitting, but the selected low-rank space remains inadequate. The
four-node second implicit defect is also too large to call a precise
root, although its direct momentum failure is unambiguous.

This rejects **this transferred three-template ansatz**, not the
paper's [Proposition 7.2 pulse inverse](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf).
The next construction must allow the transverse amplitude to vary
across the supported spatial patch and evolve along the moving pulse
path, with pressure constrained by the normal identity. Then the
Section 8 mean and Section 9 stress/moment operations need to be
included before any scale recursion claim. Further implicit solver
tuning inside these three fixed spatial templates is not warranted by
the present held-out results.
