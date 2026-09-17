# Kokuno oscillatory complete-curl provenance

## Source snapshot

This lane uses **KokunoYumeto, _Forced Navier–Stokes blowup: reconstruction and validation reader (corrected 208-page edition)_**, Zenodo DOI `10.5281/zenodo.22678406`, published 2026-09-09, version `2026.09.09-consolidated`. The public workbench snapshot used to re-check the displayed formulas is `KokunoYumeto/yang-mills-interacting-workbench@143f6773feb424ad9ed3a8d116653200f20346b7`, `navier-stokes/navier_stokes_workbench.tex`.

The corrected reader is an independent reconstruction/validation reader, not the original discovery transcript. This repository therefore treats it as a mathematical source to be checked, not as evidence that any field here is the hidden OpenAI field.

## Source structures used in K2-OSC-001 / K2-OSC-002

The reader's oscillatory-field section writes the normalized right-handed cylindrical curl explicitly and, for a transverse supported amplitude `t_m`, introduces

```text
C_m = i (n_Phi x t_m) / (k m |n_Phi|^2),
A_m = C_m exp(i k m Phi).
```

Differentiating the exponential gives the transverse leading amplitude `t_m`; differentiating `C_m` supplies the retained complete-curl remainder. The source writes that remainder as

```text
r_m = (-D_z (C_m)_theta,
       D_z (C_m)_r - D_r (C_m)_z,
       (D_r + R^(-1)) (C_m)_theta).
```

Thus the complete harmonic amplitude is `t_m + r_m`, not `t_m` alone. The reader explicitly warns that the full amplitude is not purely transverse: its longitudinal remainder participates in the exact divergence cancellation. Real waves are assembled from conjugate harmonic pairs.

The corrected edition also records the source phase/covector as

```text
Phi = p theta + p_z Z/epsilon + x_0 R - v H_Phi,
H_Phi = p F + p_z G,
n_Phi = (x_0 - v (H_Phi)_R,
         p/R,
         p_z - epsilon v (H_Phi)_Z).
```

Those formulas require the source cylindrical/auxiliary frame, base fields `F,G`, pulse variable, label data and scaling hierarchy. They are recorded here as provenance but are **not** silently substituted into the current Cartesian surrogate. Doing that before the required base/frame inputs exist would mix coordinate conventions.

For localization, identity (A14) expands

```text
curl(chi(a q) A)
  = chi(a q) curl(A) + a chi'(a q) grad(q) x A.
```

This is the key contract retained here: localization is applied to a vector potential, and the cutoff-gradient curl term is kept rather than multiplying a divergence-free velocity by a cutoff afterward.

## Autonomous reconstruction choices in the executable kernel

`openai_ns_reconstruction.kokuno_complete_curl.KokunoCompleteCurlCorrection` intentionally uses a simpler Cartesian local model:

```text
psi = n . (x - center) + omega t + phase,
t = normalized projection of polarization onto n-perp,
q = (n x t) / |n|^2,
A = - amplitude B(x) q sin(psi),
u_osc = curl(A)
      = amplitude [B t cos(psi) - (grad B x q) sin(psi)].
```

`B` is a separable compact C4 box window using `(1-s^2)^5` inside `|s|<1` and zero outside. The Cartesian coordinates make this first kernel axis/singularity safe. Amplitude, temporal frequency, phase, wave-vector components and half-widths are bounded in code.

These are **autonomous numerical/modeling choices**, not formulas copied from the source. In particular, this lane does not yet claim to reproduce the source's annular support, torus variables, exact `n_Phi`, harmonic scaling `k m`, physical `Q` scaling, or full phase/frame hierarchy.

## Analytic derivative increment K2-OSC-002

K2-OSC-002 changes no candidate parameters and adds no new oscillatory degree of freedom. It differentiates the same surrogate analytically and exposes:

- `time_derivative(...) = partial_t u_osc`,
- `spatial_jacobian(...) = [partial_j u_i]`,
- `vorticity(...) = curl u_osc`,
- `laplacian(...) = Delta u_osc`,
- `divergence(...) = trace(spatial_jacobian)`.

The C4 box is differentiated through third order so the componentwise Laplacian of the complete curl includes the derivative of the cutoff remainder. These derivatives are implementation identities for the autonomous surrogate; they are **not additional Kokuno source formulas** and do not promote the field to paper-exact status. Their purpose is to let later full-momentum residual code consume `u_t`, `grad u` and `Delta u` without finite-differencing the oscillatory correction itself.

## Relationship to existing repository lanes

The pre-existing `actual_signed_common_curl.py` and `curl_realization_algebra.py` implement coefficient/sign/scale algebra used by older symbolic witnesses. K2-OSC-001/002 do not replace or reinterpret them. Kokuno Agent 1 owns leading/base-profile and native-coordinate reconstruction, Agent 3 owns mean/radial defects and corrections, Agent 4 owns independent validation, and Agent 5 owns integration/packaging receipts.

## Validation boundary

K2-OSC-001 independently finite-differences the vector potential and verifies convergence to the analytic curl at three spatial resolutions; it also checks divergence, exact compact support, batch stability and composition with the public `velocity(x,y,z,t)->[...,3]` interface.

K2-OSC-002 independently finite-differences the public velocity on the same three-level ladder to check the analytic time derivative, spatial Jacobian and Laplacian, derives vorticity from the analytic Jacobian and checks the analytic divergence trace directly. Exact support is also required for all derivative outputs.

These increments perform **no amplitude optimization, no pressure/forcing fit, no mean/radial correction, no canonical candidate promotion, and no held-out Navier–Stokes residual claim**. A later residual screen must freeze a bounded oscillatory parameter set before evaluating held-out PDE metrics. Any nonlinear mean defect belongs to the typed handoff to Kokuno Agent 3 rather than being silently absorbed here.
