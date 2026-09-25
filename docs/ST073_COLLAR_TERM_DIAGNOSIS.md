# ST073 axial-collar momentum terms

The accepted-state gate remains false. This diagnostic uses the same robust candidate and the same eight moving collar points per time as `ST073_FORCING_REGULARITY.md`. At each registered time it decomposes the full physical momentum residual at the worst of those eight points into backward-time derivative, convection, pressure gradient, and viscosity. The selected point may differ between times; these data are neither global maxima nor an asymptotic limit.

| `tau` | Radial residual | Radial `-nu Delta u` | Axial part of radial `-nu Delta u` | Azimuthal residual |
| ---: | ---: | ---: | ---: | ---: |
| 0.024 | 4.8477e4 | 4.6045e4 | 4.5251e4 | 1.0732e4 |
| 0.016 | 1.1436e5 | 1.1060e5 | 1.0910e5 | 1.9661e4 |
| 0.012 | 2.1524e5 | 2.1013e5 | 2.0778e5 | 3.0210e4 |
| 0.010 | 3.2414e5 | 3.1799e5 | 3.1484e5 | 3.9663e4 |
| 0.0084 | 4.8210e5 | 4.7484e5 | 4.7069e5 | 5.1461e4 |

At `tau=0.0084`, axial second derivatives account for about 99.1% of the radial viscous contribution. The radial pressure-gradient contribution is only `+2.797e3`; the backward-time and convective contributions are `+2.351e3` and `+2.113e3`. Existing pressure therefore does not cancel this steep axial-collar shear. The azimuthal residual cannot be canceled by an axisymmetric scalar pressure at this point.

An ablation at these **same corrected-field-selected points** shows that the local poloidal correction of strength `0.1` is not the source of the large residual: at `tau=0.0084`, baseline norm `4.84351e5` changes to `4.84844e5`, and at `tau=0.024`, `5.00418e4` changes to `5.00627e4`. This is not a comparison of independently maximized fields.

The next numerical construction should first target the leading axial diffusion in the cutoff collar, including the azimuthal equation, with a time-dependent transition profile. A pressure-only or small static poloidal fit cannot address both components. Derive the leading scaled collar equation and its boundary data from the inner and outer fields, then test a new solenoidal transition against the complete operator over multiple times. The paper's pulse/mean-correction mechanism remains relevant for the later stress and smooth-force conditions; this diagnostic alone does not supply it.

Reproduce with `python experiments/root_st073/collar_residual_terms.py` and `python experiments/root_st073/collar_mode_ablation.py`. Reports are `experiments/root_st073/compact_potential/collar_residual_terms.json` and `collar_mode_ablation.json` in the same directory. The kinematics helper returns separate terms only when `return_terms=True`; existing callers retain the original three-value return.
