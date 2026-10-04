# First micro-switch leading correction - corrected 2026-10-04

The first chart is R=100*exp(hb*s), 0<=s<=1. The original prescription
keeps angular shear a=hb*Dbar and closes only axial shear with
1-sigma(s). Therefore d_s logF=-hb^2*Dbar/2, while
d_s V=-hb^2*(1-sigma(s))*(phi_actual/barphi)*G.

The angular hb^2 coefficient is -50*D_over_R(R100). The previous
-25 coefficient incorrectly used the axial half-weight for the angular
equation and has been corrected in both producer and checker. The
physical switch providers already used the correct source equations.
An AST bridge now binds the original controls and postpower source.

The axial pulse has exact integral1/2 by reflection symmetry. Its hydro,
pressure and swirl signed jets retain their original positive scale logs
at source points0, .5 and the exact shared root. The focused checker
passes with current hashes.

The full R100-to-R110 signed source-integral enclosures and composed
width coefficients through order2 are now documented in
[SWITCH_SIGNED_INTEGRALS_2026_10_04.md](SWITCH_SIGNED_INTEGRALS_2026_10_04.md).
They preserve the actual incoming bridge field. Higher actual Ra-to-R100
orders, nonlinear point recovery, global stress/flat remainder and
genuine recursion remain open.

Reproduce this focused first-chart packet with stage firstswitchleading,
or both corrected packet and full switch companion with stage switchintegrals.
