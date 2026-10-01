# Analytic second-width collar integral coefficients — 2026-09-30

Six of the eight normalized collar state equations are accumulated integrals. Their second-width coefficients need only the first-width dynamic solutions. The two remaining coefficients also reduce to core-inlet radial derivatives because the leading switch vanishes beyond the direct-core comparison region.

Write chi0(s)=sigma(1-s) for s<=1, zero afterwards, J(s)=integral_0^s chi0, D=D0, B=sqrt(Ra/2) Iz0, and u=u0 at the same finite core inlet. Then g1(s)=-D J(s)/2 and u1(s)=-B J(s), where these coefficients are divided by h_b. Fubini and switch symmetry give

K=integral_0^2 J(s) ds = 1/2 + integral_0^1 x sigma(x) dx.

Directed Taylor integration with 128 panels and order12, using separate flat-end bounds, gives K=0.873625407102709995 within an interval of width 5.039806934e-15. Independent MP quadrature of the weighted switch integral is contained. Exact interval endpoints are in the committed receipt.

The endpoint second-width coefficients divided by h_b^2 are:

| Normalized state | Second-width coefficient / h_b^2 |
| --- | --- |
| theta | 8-D K |
| mz | 2u-B K |
| mixed | 8u-u D K-2B K |
| axial | 2u^2-2u B K |
| swirl | 4-D K |
| p | 2-D K |

These follow by differentiating the six original right-hand sides once with respect to width before the outer width factor: 2exp(2sW+g), exp(sW)u, 2exp(2sW+g)u, exp(sW)u^2, exp(2sW+2g), and exp(sW+2g). Width-tied epsilon first appears in the second-width dynamic equations, which do not enter these six formulas.

## The two dynamic coefficients

Let M=1-K, a=Ra F_R/F, and derivatives D_R and Iz_R be evaluated from the same direct core and analytic pressure at Ra. On the support of chi0, the comparison is the finite core, so its first-width driver coefficients are s Ra D_R and s Ra Iz_R. Moreover, integral chi0(s)J(s) ds=1/8, integral s chi0(s) ds=M, and integral_0^2 (1-chi0(s)) ds=3/2. Therefore:

g2/h_b^2=-(Ra D_R M + 3D/2)/2.

u2/h_b^2=-sqrt(Ra/2)[Ra Iz_R M + Iz((1/2-a)M-D/16+3/2)].

These formulas include width-tied epsilon's second-order contribution. `second_width_coefficients` exposes all eight jets. It requires directed D_R and Iz_R input; it does not approximate them from nearby scalar samples or rerun the transition ODE.

An independent resolved fixture integrates the eight original nonlinear collar equations at both signs of a small width, using the actual smooth switch and prescribed analytic core drivers. Symmetric width differences recover the first and second coefficients. Increasing RK steps from128 to256 reduces the maximum second-coefficient discrepancy by31.84, to1.93944384262e-12. This verifies the formulas against an independently integrated resolved model; it is not an actual-source calculation or a certified RK remainder bound.

`second_width_integral_coefficients` accepts directed second-Z inlet jets. Products therefore propagate actual first and second Z derivatives, independent pressure powers, nominal coefficient enclosures, and caller-supplied perturbation enclosures. It never materializes h_b or combines the pressure atoms.

This supplies all eight analytic endpoint formulas and a directed scalar integral bound. It does not by itself establish inlet/source accuracy, evaluate the actual inlet, replace the old cache, or prove the higher-width/full nonlinear remainder. Pressure truncation, global axial coverage, stress admissibility, and true temporal recursion remain open.
