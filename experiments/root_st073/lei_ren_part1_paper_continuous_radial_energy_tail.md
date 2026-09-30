# Materialized terminal radial-energy audit

After axial support, define C(Z)=[(1-delta) Z Mz+(1-Z^2) Mz_Z]/(1-delta Z^2). The physical mapping q=tau/(1-Z^2), r=sqrt(2 nu q R), z=sqrt(nu) q^((1-delta)/2) Z gives ur=-nu C/r. No evaluated terminal value is replaced by zero.

At fixed tau, dz/dZ=sqrt(nu) tau^((1-delta)/2) (1-delta Z^2) (1-Z^2)^(-(3-delta)/2). Hence kinetic radial-energy density per dZ dlogR is pi*nu^2/2 * (dz/dZ) * C^2. A nonzero C persisting on a Z interval yields logarithmic radial growth in that numerical field. This module reports only pointwise materialized coefficients and density; it does not integrate a Z interval or certify input uncertainty.

The audit compares current post-Rv and heat values of Mz, Mz_Z, C and energy density, and checks the physical ur*r=-nu*C scaling. It uses the latest regenerated incoming swirl/energy solve. Nonzero materialized residuals are an obstruction for the current numerical representation, not a proof that the exact analytic construction fails. Exact coefficient/integral identities and input enclosures remain necessary to distinguish true closure from roundoff defects.

Do not infer total finite energy from a tiny normalized terminal residual, or impose zero in this diagnostic. The report does not certify finite energy or recursion.
