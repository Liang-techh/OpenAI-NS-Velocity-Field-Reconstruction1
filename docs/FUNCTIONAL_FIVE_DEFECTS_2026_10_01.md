# Controlled fixed-parameter five-row defects on the fresh axial family

The controlled actual R110 endpoint now supplies five terminal defect rows and their first axial derivatives over every center in [0.49,0.51]. This replaces reliance on the old scalar Z=0.3 projection for this part of the construction. It computes defects; it does not yet solve five-bump controls or prove terminal closure.

## Common construction data and centered transport

The row assembly follows the physical fixed-parameter evaluation of the labeled terms in centered_component_defects: switch-endpoint seeds, direct source terms, flat-shape kernel terms and first-unit-phase axial restoration. It preserves the accepted P0, its derivative and the same amplitude/pressure scales. It does not turn fixed parameters into formal pressure/width atoms or certify their parameter remainders.

The small axial displacement g=Uz110-4Z is reconstructed algebraically. Its terms are j, nonzero radial core coefficients, the analytic core tail, the controlled initial exit U increment, the post-exit continuation increment and a full-range first-switch source integral. The second switch and constant-power tail leave U fixed. This avoids subtracting two independent copies of 4Z. The enclosed g value remains approximately 1e-14 rather than the spurious width [-0.08,0.08] from direct family interval subtraction.

The first two normalized defect rows are approximately:

- d1: [2.2313324661e-15, 2.2374172161e-15].
- d2: [5.6760876119e-16, 5.7036879397e-16].
- d4: [7.6459757764e-41, 7.9158335513e-41].

The d3 and d5 bounds are negative and extremely small under the prescribed exponential scales; their exact MP intervals are retained in the receipt. These normalized rows are neither full NS residuals nor a percentage of project completion.

## Full negative-amplitude flat kernel enclosure

For I=integral_0^T exp(-k*l) expm1(m B sigma(l/T)) dl, the enclosed B is strictly negative. Therefore abs(expm1(m B sigma)) <= m abs(B) sigma and dI/dB <= m integral exp(-k*l) sigma dl. Split at a positive L<=T/2. On the first piece sigma <= exp(-(T_lower/L)^2 + 1/(1-L/T_lower)^2); on the second piece use sigma<=1 and the exponential tail integral. The resulting directed bound encloses the entire kernel and its first axial derivative through B_Z. No saddle window or omitted quadrature tail is treated as exact.

B_Z is recovered as ell + phi_Z/phi + 2Z/(1+Z^2), cancelling the common F0 before division. This avoids exponential interval growth from F_Z/F when the amplitude family is very small.

The three restoration integrals K1, K16 and K2 use 512 directed cells on their full [0,1] support. Their uncertainty is included in the defect rows.

## Evidence and next work

Independent quadrature fixtures contain 18 kernel value/axial coefficients, three restoration constants and two centered transport coefficients. The fixtures do not independently validate the complete five-row assembly. A read-only formula review checks the row assembly against its existing source algebra.

Next construct a functional five-bump inverse/repair using these family rows and controlled integral coefficients. Retain the same pressure/amplitude data, propagate these row uncertainties, and test terminal identities as functions. Higher axial orders, whole-axis coverage, source parameter errors, heat-exterior matching and genuine temporal coefficient recovery remain incomplete.
