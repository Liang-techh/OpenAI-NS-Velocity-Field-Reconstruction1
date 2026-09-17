# Kokuno oscillatory complete-curl provenance

## Source snapshot

This lane uses **KokunoYumeto, _Forced Navier–Stokes blowup: reconstruction and validation reader (corrected 208-page edition)_**, Zenodo DOI `10.5281/zenodo.22678406`, published 2026-09-09, version `2026.09.09-consolidated`.

The corrected reader is an independent reconstruction/validation reader, not the original discovery transcript. This repository therefore treats it as a mathematical source to be checked, not as evidence that any field here is the hidden OpenAI field.

## Source structures used in K2-OSC-001

The reader's oscillatory-field section writes the normalized cylindrical curl explicitly and, for a transverse supported amplitude `t_m`, introduces

```text
C_m = i (n_Phi x t_m) / (k m |n_Phi|^2),
A_m = C_m exp(i k m Phi).
```

Differentiating the exponential gives the transverse leading amplitude `t_m`; differentiating `C_m` supplies the retained complete-curl remainder. The reader explicitly warns that the full amplitude is not purely transverse: its longitudinal remainder participates in the exact divergence cancellation.

The corrected edition also records the physical phase as containing the displayed `p_z Z/epsilon`, `v H_Phi`, `p theta`, and `x_0 R` pieces and uses cylindrical frame factors. Those source-specific phase/frame terms are **not implemented in this first executable kernel**.

For localization, identity (A14) expands

```text
curl(chi(a q) A)
  = chi(a q) curl(A) + a chi'(a q) grad(q) x A.
```

This is the key contract retained here: localization is applied to a vector potential, and the cutoff-gradient curl term is kept rather than multiplying a divergence-free velocity by a cutoff afterward.

## Autonomous reconstruction choices in this increment

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

These are **autonomous numerical/modeling choices**, not formulas copied from the source. In particular, this increment does not claim to reproduce the source's annular support, torus variables, exact `n_Phi`, harmonic scaling `k m`, physical `Q` scaling, or phase/frame hierarchy.

## Relationship to existing repository curl code

The pre-existing `actual_signed_common_curl.py` and `curl_realization_algebra.py` implement coefficient/sign/scale algebra used by older symbolic witnesses. K2-OSC-001 does not replace or reinterpret them. Its only new responsibility is an executable spatial vector potential and its complete curl that can be added to the current callable 3D velocity delivery.

## Validation boundary for K2-OSC-001

The tests independently finite-difference the vector potential and verify convergence to the analytic curl at three spatial resolutions. They independently finite-difference the analytic velocity and verify divergence converges to zero, test exact compact support and batch stability, and compose the correction with the current support-connected Eq. (4.5) `velocity(x,y,z,t)->[...,3]` interface.

This increment performs **no amplitude optimization, no pressure/forcing fit, no canonical candidate promotion, and no held-out Navier–Stokes residual claim**. The first residual screen must freeze a bounded oscillatory parameter set before evaluating held-out PDE metrics. Any nonlinear mean defect belongs to the typed handoff to Kokuno Agent 3 rather than being silently absorbed here.
