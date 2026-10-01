# Controlled fresh-family switched comparison and frozen extension

The fresh degree124 core for axial centers in [0.49,0.51] now feeds a directed integral enclosure of the Section9.23 switched comparison equations on y in [0.005,0.01]. The accepted fourteen-stage pressure/amplitude data remain unchanged. This advances the uncertainty chain beyond the earlier finite RK4 diagnostic.

## Construction and error scope

Write g=s alpha Phi_s/Phi and j=s alpha U_s, with s=4 exp(y). The comparison equations are phi'=g phi and U'=j. They have integrating-factor and additive integral solutions. Each of32 cells encloses g and j over the entire cell using the finite source polynomials plus their analytic radial tails. The cutoff alpha is bounded by its monotone endpoint values, with exact limiting values at the two ends. Directed cell integrals enclose the full G/J integrals; partial-cell lengths in [0,h] enclose fields throughout each cell. The resulting field ranges bound each moment integrand before integration.

The radial-order1, axial-order0..2 error budget is within the existing mixed C3 certificate. Hence field and moment jets retain axial order2, while the D/I_z drivers retain axial order1 after axial differentiation. No unsupported third-order driver certificate is claimed.

The initial Phi(4) used in A/B includes its radial tail. The terminal driver formulas are A=-D/2 and B=-sqrt(s epsilon/2) Phi(4)/phi(s) I_z, with chi=1 explicitly stated. Physical pressure remains P0+S epsilon m_p, with S=F0 squared and the same physical P0. No pressure fit, amplitude reset or old-center tensor extrapolation is introduced.

At the switched endpoint, the normalized angular field enclosure is approximately[0.04044130,0.51696617], and D is approximately[1.06413,35.82671]. Positive source denominators were checked in every integration cell. These are bounds over the entire center family; their width includes source-family dependency and integral overestimation. They are not physical NS residual values.

After alpha vanishes, comparison phi and U are frozen. Their six moments have exact polynomial radial increments. The receipt includes this frozen comparison continuation to physical R=110 (scaled radius110 Lambda), preserving the enclosed switched endpoint. It does not extrapolate the analytic inner core beyond scaled4.1. This is the auxiliary comparison field, not a solution of the physical exit continuation ODE.

The receipt reports comparison discretization enclosed by interval cell integration, and analytic core tail propagated. It does not certify original construction-parameter errors, physical exit ODE error, final stress cone, whole-axis transition or terminal five-moment closure. Existing finite RK4 receipts remain diagnostic only.

## Independent evidence

The enclosure fixtures check24 coefficients against closed-form constant-in-radius fields with nonzero axial derivatives, another24 coefficients for their exact frozen continuation,24 nonzero-slope field coefficients against a256-step RK4 diagnostic, and6 corresponding driver coefficients. The numerical reference is supporting accuracy evidence rather than the basis of the enclosure proof. A separate read-only mathematical review confirmed the integrating-factor, partial-cell, cutoff and radial-weight arguments.

The earlier deep consistency job has also completed. It captured a coherent historical degree38 state (SHA2cf064d95c53d64c15e90028c9a492aeb1e95af2449ba76c886c3b3c453cdf42), before the live production state advanced. All4251 physical pressure coefficients satisfy interval overlap with the pressure recurrence. A fresh scalar center Z=0.5 rebuild with initial axis length128 is contained by all12753 saved A/Uz/P family coefficients. These are degree38 checks, not a new independent degree124 scalar computation; the receipt states its captured degree and hash.

## Continue

1. Integrate the physical exit equations from these controlled comparison drivers with interval errors, retaining analytic axial derivatives and actual moment data. The frozen comparison alone does not accomplish this.
2. Reduce axial-family width or exploit correlated quantities when ratios remain too broad for cone or moment repair bounds. The tiny radial tail does not remove axial dependency width.
3. Provide family-aware bridge interfaces that preserve these directed intervals instead of silently projecting to a midpoint.
4. Propagate controlled exit/switch errors into functional defects and five-moment repair; then extend the axial atlas.
5. Complete heat exterior compatibility, final admissible stress and temporal n-dependent recursion before oscillatory stress correction and full Cartesian residual validation.