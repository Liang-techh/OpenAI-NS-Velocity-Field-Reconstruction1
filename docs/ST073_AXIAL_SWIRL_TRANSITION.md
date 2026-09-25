# Multi-time axial-collar swirl trial

The previous term decomposition found an azimuthal momentum residual of about `5.15e4` at `tau=0.0084`, alongside a still larger meridional diffusion error. An axisymmetric pure-swirl increment is solenoidal and can target the azimuthal equation without changing meridional velocity directly. This trial adds a compact profile

\[
\delta u_\theta=\sum_{i,j} a_{ij}(\tau_0/\tau)^{1/2}
  y(1-y^2)^5 P_i(2y^2-1)
  1024s^5(1-s)^5P_j(2s-1),
\]

where `y=r/R(tau)`, `s=(|z|-z_flat)/(z_support-z_flat)`, `tau_0=0.0084`, and the profile is zero outside `0<y<1`, `0<s<1`. The Cartesian increment is `(-delta_u_theta*y_cart/r, delta_u_theta*x/r, 0)`, regular on the axis. The factors give fourth-order matching at the support boundaries, and compact spatial support on each registered time slice. This is an autonomous diagnostic basis, not the paper's oscillatory pulse construction.

The azimuthal equation is linear in this increment when the meridional velocity is fixed. The six coefficients are fitted across five registered times, then the **complete** momentum residual is evaluated with centrifugal feedback. Sparse fitting at two axial fractions (`s=.35,.65`) gave an apparent ~200-fold reduction at its own points but increased the angular residual about 16-fold at unseen points. This version is rejected as overfit. A one-mode constant axial basis barely changed the fitted residual; a two-mode basis improved training but worsened holdout by about 1.5-fold. The baseline angular residual changes sign between `s=.35` and `.65`, making both axial shape and spatial coverage necessary.

The dense variant fits at `s=.2,.35,.5,.65,.8` and `r/R=.25,.4,.55`. At disjoint points `s=.275,.425,.575,.725`, `r/R=.3,.5`, its angular maximum changes as follows:

| `tau` | Baseline | Dense swirl correction | Full residual max, baseline → corrected |
| ---: | ---: | ---: | ---: |
| 0.0084 | 5.2854e4 | 3.8720e4 | 5.9218e5 → 5.9139e5 |
| 0.012 | 3.0741e4 | 2.2436e4 | 3.1380e5 → 3.1327e5 |
| 0.024 | 1.0773e4 | 7.8172e3 | 9.6164e4 → 9.5934e4 |

This is a repeatable local improvement of about 27% in the angular residual, but the full momentum remains dominated by axial diffusion in the meridional equation. No global max or physical-volume L2 certificate has been run for the new candidate, and `accepted` remains false. It does not establish a smooth force at the critical time.

Reproduce the dense fit with `python experiments/root_st073/axial_swirl_multitime_screen.py --dense`, then its independent points with `python experiments/root_st073/axial_swirl_holdout.py --dense`. `load_dense_candidate()` in the first script reconstructs the candidate from `compact_potential/axial_swirl_dense.json`. The next construction step is a dynamic solenoidal meridional transition that cancels the leading axial-diffusion term while preserving the inner/outer flux and stress data. The [OpenAI paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf) uses radial moment repair, a heat exterior, and later oscillatory stress/mean corrections; this local swirl fit substitutes for none of those stages.
