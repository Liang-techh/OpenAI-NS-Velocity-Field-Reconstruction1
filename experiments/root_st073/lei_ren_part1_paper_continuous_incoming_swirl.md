# Shared continuous incoming swirl energy

ContinuousIncomingSwirl regenerates reference I_swirl = (1/(2 Rref)) integral_0^Rp Utheta^2 dR from the same normalized angular propagation used by actual cumulative moments. It reuses continuous schedule primitives and analytic constant-stage exponential transfers; variable stages use declared quadrature nodes. No inherited float swirl row is reused.

The incoming reference is separated from actual inner offsets. regenerate_incoming replaces I_swirl before reapplying measured inner z_theta offsets once and solving the live energy/linear equations. The swirl reference divided by Ep^2 remains Z independent, so the existing analytic reference prior-energy tangent cancellation remains applicable. Actual inner and live future-energy jets still have their declared uncertainty.

The source reference is memoized on the shared source object and uses the same installed angular correction/schedule definitions. Nominal quadrature and working precision do not certify input accuracy. Complete five moments, finite energy and scale recursion remain unclosed.

Run the module for a shared-field comparison of regenerated reference plus actual inner swirl offset with the angular cumulative primitive at Rp, together with reference Z factorization and live tangent evaluation.

Shared-field result: Rp value/Z matching about 1.17e-84, inherited input change about 7.92e-19, updated linear replay maximum about 1.09e-173. These are nominal consistency checks.
