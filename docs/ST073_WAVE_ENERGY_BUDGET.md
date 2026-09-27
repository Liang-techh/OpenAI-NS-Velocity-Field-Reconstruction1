# Supported-wave production and dissipation

## Why this calculation changes the construction decision

The implicit spatial Fourier trajectory completed, but its oscillatory
energy proxy fell to 1.44% over the saved interval. This establishes decay
of that trial, not an impossibility result for pulses or scale recursion.
The [reference paper, Section 7, equations (7.2), (7.5), and (7.12)](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
selects a pulse with amplification followed by viscous decay; small tails
are needed for temporal localization. Therefore decay alone is not the
failure criterion. The missing evidence is a growth phase and a coupled
stress correction that reduces full momentum across scales.

## Fixed-support energy screen

For a compactly supported divergence-free perturbation w of a background
U, the homogeneous linearized equation gives, with E = integral |w|^2 / 2,

```text
dE/dt = - integral w . sym(grad U) w - nu integral |grad w|^2
        + (1/2) integral (div U) |w|^2.
```

The last term vanishes for a divergence-free background. Advection by U
and pressure contribute only the integration-by-parts terms displayed
above. This identity concerns the homogeneous linearized problem; it
does not include the background momentum source, external forcing, or
energy exchanged with an independently evolving mean correction.

For the retained exact-curl Fourier basis V, assemble velocity mass M,
gradient Gram matrix K, and strain form S using cylindrical volume
weights. A weighted velocity SVD removes the numerical gauge nullspace
and produces T with T* M T = I. Eigenvalues of

```text
H = T* (-nu K - S) T
```

bound instantaneous amplitude growth in that retained space; the
corresponding logarithmic energy growth rates are twice these values.
This Hermitian energy calculation is different from the eigenvalues of
the nonnormal evolution operator. Negative values rule out instantaneous
homogeneous energy growth in the tested subspace and frozen background,
not growth in every possible spatial field or at every future time.

The geometry screen keeps the center, physical carrier wave numbers and
degree fixed while widening the radial and axial support. It must not
be interpreted as an accepted wave: wider supports can leave the stress
cone region or damage the inner/outer matching. Positive growth, if found,
would only justify a subsequent coupled-field construction.

## Reproduction and decision

Run `python experiments/root_st073/fourier_patch_energy_budget.py` from
the repository root. Its JSON records the actual sampled geometry,
quadrature and spectra. The 16-by-16 results, maximized separately over
retained modes 1 through 8, are:

| Support halfwidth multiplier | Largest strain-only rate | Least negative diffusion rate | Largest combined rate |
| --- | ---: | ---: | ---: |
| 1 | 8.1561e3 | -1.6766e7 | -1.6762e7 |
| 2 | 9.1051e3 | -4.1726e6 | -4.1692e6 |
| 4 | 1.1380e4 | -1.0434e6 | -1.0419e6 |
| 8 | 1.4192e4 | -2.6523e5 | -2.6326e5 |

These are physical inverse-time amplitude rates. The extrema in the
different columns need not have the same eigenvector, so their sum is
not the combined rate. Every tested direction decays. Background
divergence is at most 3.72e-6 on these nodes and is negligible on the
displayed rate scale; the artifact also includes its energy correction.

The completed 24-by-24 refinement at multipliers 1 and 8 preserves the
negative signs. The largest relative change of any mode's largest
combined eigenvalue is 2.27e-9 and 2.55e-5 respectively. The refined
best rate at multiplier 8 is -263261.4685. This is a quadrature check
within the fixed degree-2 subspace, not spatial-basis convergence.

At multiplier 8 even the most optimistic separately maximized strain
is less than one eighteenth of the weakest viscous loss. At multiplier
1 the gap is over a factor of 2000. This explains why stable integration
alone does not yield an amplifying stress pulse in the current subspace.
It does not prove that enlarging the approximation space, changing the
mean, or changing the support geometry cannot work.

All four patches remain off-axis and inside the registered axial/time
slab. This domain check does not establish the stress cone over their
entire supports. Initial modes 1 and 4 have axial carrier phase
halfwidths only 0.0348 and 0.1391; even at multiplier 8 these are about
0.278 and 1.113. The envelope is therefore not verified to vary slowly
relative to the carrier, as the principal-wave approximation would need.

The next construction should co-design the background shear and supported
wave space, with this energy budget as an early gate. Keep the corrected
mean moment/cone constraints and require an actual full-momentum
improvement afterward. Do not rerun the same degree-2 patches expecting
time-step changes to produce growth. Neither this screen nor the
preceding BDF trajectory satisfies the full momentum max/volume-L2 or
recursion gates.
