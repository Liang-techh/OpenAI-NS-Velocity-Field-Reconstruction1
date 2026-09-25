# Dyadic transfer of the ST073-V inner kernel

`experiments/root_st073/core_dyadic_transfer.py` evaluates the actual
ST073-V field on fixed similarity coordinates `0.001 <= X <= 1/64`,
`|eta| <= .3`, at `tau=.5*2**(-k)` for `k=0..6`. It normalizes radial
velocity by `q**(1/2)` and angular/axial velocity by `q**A`, with
`q=tau/(1-eta**2)`, before comparing each pair of adjacent scales.

The normalized full-vector profile changes by relative RMS
`0.000153..0.000159` per dyadic step on this finite grid. The
similarity-normalized complete residual maximum decreases from about
`2.14e-9` at `k=0` to `1.57e-9` at `k=6`; sampled physical residuals
rise as the core concentrates but remain below `2.1e-7` at `k=6`.
Sampled azimuthal speed grows by about `1.419` per halving of `tau`.

The exact coordinate exponents of this candidate are `h=.005`,
`A=.505`, and `D=.495`. They give radial length ratio `0.7071`,
axial length ratio `0.7096`, and axial-to-radial aspect growth
`2**h=1.00347` per dyadic step, only `1.021` over the six checked
steps. Thus the direction of relative axial slenderization is built
into the local scaling, but it is quantitatively weak over this window.

This establishes a reproducible *local* profile-transfer diagnostic,
not a recursive whole-space solution. It says nothing about the
transition annulus, exterior heat/pressure compatibility, nonaxisymmetric
stress, total energy, or a controlled limit as `tau -> 0`. The full
momentum max and spatial-volume L2 acceptance gates must still be
evaluated on a matched physical field.
