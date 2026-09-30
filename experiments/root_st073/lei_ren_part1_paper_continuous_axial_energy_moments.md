# Actual complete axial-square and exterior quadratic moments

ContinuousIncomingProfile.quadratic_moments_jet currently supports Rv through Rtail. It combines actual inner axial energy, continuous incoming energy, and the live solve's complete pulse/end energy atoms. Pulse and end supports are disjoint: the end band has xi = 13 + mu*t_v > 11. There is no cross term. The end Gram atoms already contain exp(-26).

The Z jet differentiates coefficients and angular amplitude. The provider combines axial_square with the actual angular swirl-square primitive as z_theta = axial_square - swirl_energy/2. Axial-square value and Z derivative remain constant after axial support. Quadratic forward integrand relative differences are 7.34e-15 and 1.88e-14 in the flattening region.

Partial axial energy before Rv, heat continuation, enclosed quadrature/input errors and terminal energy matching remain open. All full-moment, finite-energy and recursive certification flags are false. Run the matching Python module to regenerate its JSON.
