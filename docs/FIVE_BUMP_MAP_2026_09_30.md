# Section10 five-bump nonlinear moment map

The next velocity correction uses the reference interval x=R/Rm in (1,2). The fixed normalized beta_1/40 bumps are centered at 5/4, 3/2 and 7/4. Axial coefficients c1,c2 use the first and third supports; angular coefficients xi1,xi2,xi3 use all three. Set f=sum xi_j gamma_j and g=c1 gamma1+c2 gamma3. The corrected velocities are u=Am*(x^(1/10)+f), V=4Z+g, with Am=exp(logPstar-.6)/(1+Z^2).

The centered moment changes from Eq.(10.8) are integrals over (1,2):

1. g
2. x^(3/5)*g + sqrt(x)*f*g
3. sqrt(x)*f
4. -x^(1/10)*f + Am^(-2)*g^2 - f^2/2
5. x^(-9/10)*f + f^2/(2*x)

Thus the equation to solve is A*h+Q_Z(h,h)=-d, with d supplied by CenteredComponentDefects. The exact paper branch is the contraction h -> -A^(-1)*(d+Q_Z(h,h)) under its uniform C1 smallness hypotheses. Compact supports preserve endpoint matching, but cumulative moment changes persist after the supports end. The pressure datum P0 remains fixed; pressure changes through the corrected cumulative pressure moment.

The fixed linear matrix has two separated axial supports and three separated angular supports. It has no input-dependent coalescing exponents. Quadratic crossproducts occur only on matching supports. A generic float least-squares repair elsewhere in the repository is not an implementation of this coupled map.

## Numerical scope and next acceptance conditions

The finite moment map alone does not solve the actual defects or certify cone admissibility. Its quadrature is unenclosed. A resolved independent replay must integrate corrected field densities, test the inverse, and differentiate those densities separately for first-Z validation.

Actual d3/d5 log magnitudes are of order -1e152; other induced coefficients can be of order 1e-41. A fixed-precision scalar contraction can therefore lose tiny row responses through cancellation even while showing a small absolute residual. Preserve defect-source responses or a formal convergent defect expansion, and report any finite truncation remainder. A finite polynomial expansion is not an exact moment-closure certificate.

Uniform C1 input smallness, inherited source bounds and second-Z control remain necessary. After a valid correction, use partial bump integrals to recover radial velocity, pressure and stresses, check the relaxed cone throughout the supports, then append the same outer construction and assess its energy/heat conditions. Temporal coefficient recursion remains subsequent work.

## Resolved finite-map evidence

The independent five_bump_map_fixture integrates the corrected normalized velocity densities with mp.quad, including all quadratic terms. At Gauss order96 and precision90, maximum relative row, first-Z and inverse replay errors are below 4.41e-16. The first-Z replay uses a separate four-point scalar stencil. These results cover the finite map and declared bump quadrature, not the actual common-source defect solve or a uniform functional certificate.
