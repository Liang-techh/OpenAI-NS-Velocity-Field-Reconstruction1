# Pressure-weighted heat primitives from the shared heat owner

ContinuousHeatMoments.pressure_increments(t,Z) returns the heat contribution to Mp=integral Utheta^2/(2R)dR. It uses the same K0, retained heat polynomial correction and Z jet as swirl energy, with normalized radial exponent -(1+delta) and physical normalization Ustar_tail^2/2. Returned reference, heat_correction and heat_correction_Z stay separate; pressure_increment is their nominal sum.

complete_pressure_heat_integral(Z) integrates from Rtail to infinity, retaining separate reference/correction/Z atoms and an analytic heat-polynomial truncation bound. The pure-power exterior term is exp(-3*(1+delta))/(1+delta). The bound is physical_scale*c3*xi0^3/[3*(4+delta)]. Quadrature and finite arithmetic remain unenclosed.

Shared-field checks at t=.6 and 4 agree with independent radial integrands to maximum relative difference about 9.50e-14. The infinite correction Z jet agrees nominally with centered differences at working precision. All three increments are zero at Rtail.

The actual inner-seeded full pressure moment is not installed by this helper. A bounded worker is building the preheat/inner pressure provider and will consume these shared heat primitives. Full five moments, pressure matching, stress, finite energy and recursion remain open.
