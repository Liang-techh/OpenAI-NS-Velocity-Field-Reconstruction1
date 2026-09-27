# Wave growth after mean evolution

`evolved_wave_energy.py/.json` compares instantaneous perturbation energy
production minus viscous dissipation on the same fixed physical Fourier
patch. It uses the actual evolved mean gradient, complete exact-curl basis,
mass whitening, and explicit support checks. The operator is -nu K - S;
its largest generalized eigenvalue lambda corresponds to energy rate
2 lambda. A separate numerical divergence term is reported.

Unsplit tensor quadrature gives the following mode-1 energy rates:

| Mean state | Energy rate 2 lambda |
|---|---:|
| Initial | +20000.0 |
| Unconstrained delta k=0.01, two steps | -446060.0 |
| Unconstrained delta k=0.01, four steps | -429090.5 |
| Constrained delta k=0.0001 | +13955.5 |

All eight tested modes decay in both unconstrained endpoints. These are
finite-basis diagnostics; the tensor grid did not split all evolved mean
cutoff boundaries. They corroborate the observed weakening and do not prove
absence of growing waves outside this basis or patch.

`evolved_wave_energy_split.py/.json` separately compares mode 1 initially and
after the constrained step, using split physical radial panels including
the meridional onset y=0.62:

| Mean state | z order 12, panel order 6 | z order 16, panel order 9 |
|---|---:|---:|
| Initial | +19026.6 | +20437.9 |
| Constrained endpoint | +13161.2 | +14500.8 |

Both time slices remain supported and both tested quadratures retain one
positive direction among 27. The relative rate changes are 6.91% and 9.24%;
this is sensitivity evidence, not quadrature convergence. Split results
export the dominant potential coefficients normalized in the physical
mass matrix. Do not compare their raw coefficient norm with coefficients
normalized using unit-sum quadrature weights without rescaling.

Thus the short constrained step preserves a tested instantaneous growth
direction, whereas the examined unconstrained trajectories do not. No wave
has been integrated in these screens. Positive energy growth also does not
guarantee the required radial momentum-flux covariance. The next construction
optimizes actual wave polarization and amplitude under a positive-growth
constraint, then must replay complete mean-plus-wave momentum.
