# Five-bump field and cumulative-moment recovery

The correction now has a field adapter for the same Section10 coefficients. It recovers u_theta, Uz, F, all five cumulative moment changes, pressure changes, and the value of Ur from the same axial moment and first-Z derivative. P0 stays fixed. The logarithmic derivative uses x*d/dx, and Rm must remain independent of Z. The reference join updates the stress/shear inputs instead of keeping baseline derivatives.

For normalized partial changes a1..a5, the physical moment increments are

    Delta Mz = Rm*a1
    Delta Mtheta = sqrt(2)*Rm^(3/2)*Am*a3
    Delta Mtheta_z = sqrt(2)*Rm^(3/2)*Am*(a2+4Z*a3)
    Delta Mztheta = Rm*Am^2*a4 + 8Z*Rm*a1
    Delta Mp = Am^2*a5

Pressure changes by Delta Mp, without changing P0. Radial velocity changes by

    Delta Ur = [2ZR*Delta Uz -(1-delta)Z*Delta Mz
                -(1-Z^2)*Delta Mz_Z] / [(1-delta Z^2)*sqrt(2R)].

First-Z of Ur needs second-Z moment/coefficient input and is explicitly unavailable.

## Separate defect responses and all quadratic terms

The response-field adapter preserves each defect monomial when computing velocity and moments. A degree N velocity correction produces quadratic moment terms through degree 2N. Those terms remain in the result, including terms above the formal inverse solve degree. Dropping them would make the reported moments inconsistent with the finite corrected velocity.

The actual serialized Z=.3 P9/W2 check at x=1.25 completed with 55 velocity response terms, 461 moment terms and 406 nonzero terms above degree 3. Direct d3/d5 velocity responses stay nonzero. The receipt reports per-monomial baseline/first-Z logs and finite-ring atom counts; the calculation retains all declared P9/W2 atoms. This is a correction-only evaluation, not an installation on the complete actual baseline.

The independent resolved fixture at x=1.25 and 2.2 directly integrates all five corrected velocity-density differences and checks radial recovery with an independent Z stencil. Maximum moment relative error is 4.34e-16 and radial error is below 2.64e-80. After the supports, velocity corrections vanish while cumulative moment increments remain. A separate reference-join fixture confirms P0 preservation, pressure/moment additions, logarithmic derivatives and the radial identity.

## Conditional convergence bound

five_bump_majorant computes CA and a quadratic C1 operator bound from the fixed finite weights, using the paper estimate norm(Am^-2)<=40/Pstar^2. Given an externally established uniform defect bound e, it computes the contraction condition e<=1/(8 CA^2 CQ) and a positive Catalan-series tail bound without subtracting nearly equal roots. Its resolved scalar fixture passes zero/boundary/failure cases and bounds the independently observed response error.

The finite quadrature weights are unenclosed and actual uniform e has not been verified. These majorants are conditional algebraic bounds, not a certificate of actual five-moment closure, admissible stress, exact heat, finite energy, or temporal recursion. Next priority is joining the same actual reference provider and then controlling the input/convergence/remainder and cone conditions uniformly in Z.

## Common-source integration convention

The repository raw swirl primitive is integral(R F^2) dR = integral(u_theta^2/2) dR. The field adapter now uses the same factor, so Mztheta = raw_axial - raw_swirl. Joined moment/raw part maps explicitly include the five_bump_correction label; downstream consumers must not use an inherited part map that omits the correction. The resolved join fixture checks this convention and the corrected part sums.

Old materialized long-reshape/restore moment totals can round away the flat d3/d5 terms. On Rm..Rh, source and power reference velocities agree, so each original five-moment difference is constant in R. A same-source reference provider can therefore reconstruct its moments as the exact power reference plus the separately retained centered source defects, using the original P0. Separate raw energies require their own source primitives and cannot be inferred from d4 alone.
