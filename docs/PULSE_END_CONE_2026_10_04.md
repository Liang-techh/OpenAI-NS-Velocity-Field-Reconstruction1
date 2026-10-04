# Continuous original pulse-end cone — 2026-10-04

The complete original end region s in [-4,0], Z in [-1,1] now satisfies the
paper's regional two-vector admissible stress criterion. Both beta supports,
their complement, signed angular memory, full meridional transport, actual
absolute pressure and all ten stress sectors remain. No field is replaced
by a numerical cap.

## Full shear and continuous bounds

The original equations give

`S=A*(-1,sigma)`, `A=(2+2mu)*C*B/sqrt(2R)>0`,
`sigma=D*(Bhat_s-(.5+mu)*Bhat)/(1+mu)`.

Thus `kappa-2=2mu+(2+2mu)*sigma^2`. Axial shear remains nonzero inside
the supports. The proof keeps this positive difference symbolically instead
of subtracting a rounded kappa from 2.

Normalize the actual stress by `Qbase=sqrt(R/2)*B>0`. The angular equilibrium
is exactly `[C*g+2*b0*Z^2/(r*(1+Z^2)^2)]/L`, with
`g=(mu-delta/2)/r>0`, `r=1-mu`, `b0=(1-delta)/2`,
`L=1-delta*Z^2`. It is at least g/2 on the whole original axial interval.
Each of the three angular correction sectors is bounded by g_lower/32.
Their signed sum therefore admits the conservative lower theta_min=g_lower/4.

All six axial sectors are bounded by theta_min*exp(-1000)/6, and the actual
shear ratio by theta_min/theta_max*exp(-1000). These are bounds only, obtained
from the exact original B/D/H/R log recipes and current whole-source
coefficient enclosures. Modes are canceled symbolically before enclosure;
enormous nearly equal log boxes are never subtracted.

The checker explicitly establishes

`dot_min=theta_min-axial_max*sigma_max>0`

and

`2*dot_min^2-(kappa-2)_max*(theta_max*sigma_max+axial_max)^2>0`.

These prove T dot S < 0 and the strict original directional inequality on
the entire region. The support receipt supplies normalized coefficient
difference bounds; its absolute positive factors are restored through the
accepted parent sectors and source log recipes. Whole-source history
enclosures remain bounds for the edge functions, not chosen edge values.
The selected C5 source checker is directly consumed here.

## Composition and physical units

The actual production radius is replayed from source. Pulse s=0 equals
flatten t=0; pulse s=-4 lies four log-radius units earlier. The newly admitted
regional nonzero tail is

`Rtail*exp(-wait-Ts-106-Lrel) <= R < Rtail*exp(3)`.

Gamma stress is exactly zero beyond the accepted outer endpoint. Current
pulse/flatten full moment, pressure, physical and support interfaces are
consumed; downstream accepted proofs are not rerun.

Both physical stress components and shear components acquire the common
positive factor nu*lambda^(-2-delta). For fixed nu>0 the physical strength is
`kappa=-|Sphys|^2/(nu*Fphys*Sphys_theta)`; this preserves the dimensionless
source criterion. The unchanged physical evaluator is exposed through
CertifiedPulseEndPhysical, which only adds current regional admission.

## Evidence and remaining scope

Run the bounded stage:

`python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage pulseendcone`

Focused checker PASS: 419 current input hashes, 22 exact source/cone/radius
identities, 10 positive logarithmic source margins and 10 positive algebraic
cone margins. Four independent original-evaluator fixtures retain both
axial shear signs and exercise nu=.01/.7: 36 comparisons, tolerance 1e-60,
maximum absolute error 4.9091e-91. Those fixtures check normalization and
units; the actual cone proof uses continuous whole-source enclosures.

This admission covers the original terminal end region and the previously
accepted downstream tail. Gap/main/entrance coverage and their adjacent
interfaces remain open. Upstream finite-width bridge feedback, completed
full-tensor/global admissibility, global temporal flatness, physical-volume
norms, required-domain kinetic energy, actual n-dependent coefficient
recursion, oscillatory correction and corrected Cartesian residual/dynamics
remain unfinished. A regional two-vector cone is not a completed global
tensor cone.
