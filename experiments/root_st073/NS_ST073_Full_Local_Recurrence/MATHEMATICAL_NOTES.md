# All-momentum radial recurrence: local and finite, not global admission

## 1. Problem and comparison identities

This experiment uses the **local unforced** operator with physical viscosity
nu=0.01, on X=r_source^2/(2q)<=1/64, |eta|<=1/2,
q-z_source^2 q^(2h)=tau, h=.005, tau in [.5/64,.5]. The physical coordinates
are x_phys=sqrt(nu)*x_source, and physical u,p are sqrt(nu)*u_source, nu*p_source.
The source operator has viscosity one. It is not the older globally compact,
initial-energy-normalized, two-parameter fixed-force problem.

ST073-F has exactly the ST068 axis traces:
B(0,z,t)=q^(-1-h), C(0,z,t)=q^(-1/2-h)*(4 eta+.02),
p(0,z,t)=q^(-1-2h)*(-1+eta^2/2).
The generated off-axis field differs because **full**, rather than only leading,
axisymmetric momentum is imposed radially. This is not a refit of the sharp ST072-N.
ST073-V is a second registered control that changes axis swirl from 1 to 4;
its comparisons use a newly generated leading-only field with that SAME axis datum.
Changing this amplitude is explicit. No global energy conservation/normalization
or source-exactness is asserted.

## 2. Derivation at source viscosity one

Let s=x^2+y^2, u=(x A-y B,y A+x B,C), and write all four quantities as series in s.
Incompressibility gives A_n=-(partial_z C_n)/(2(n+1)). In scalar form,

    RA=A_t+A^2+2s A A_s+C A_z-B^2+2p_s-(8A_s+4s A_ss+A_zz),
    RB=B_t+2AB+2s A B_s+C B_z-(8B_s+4s B_ss+B_zz),
    RC=C_t+2s A C_s+C C_z+p_z-(4C_s+4s C_ss+C_zz).

The complete Cartesian vector is (x RA-y RB,y RA+x RB,RC).
Equating coefficient n produces a triangular recurrence:

    4(n+1)(n+2) B_(n+1)
      = dt B_n + sum_(i+j=n)[2(j+1) A_i B_j+C_i dz B_j] - dzz B_n,

    4(n+1)^2 C_(n+1)
      = dt C_n + sum_(i+j=n)[2j A_i C_j+C_i dz C_j] + dz p_n - dzz C_n,

    A_(n+1) = -dz C_(n+1)/(2(n+2)),

    2(n+1) p_(n+1)
      = -dt A_n - sum_(i+j=n)[(1+2j)A_i A_j+C_i dz A_j-B_i B_j]
        + dzz A_n + 4(n+1)(n+2)A_(n+1).

In particular, pressure is **not** fixed to the leading centrifugal relation
Pi_X=F^2. It is solved along with the two velocity series. All axial diffusion,
radial inertia and time derivatives appear before the next radial coefficient
is generated. This is an autonomous local Taylor construction, not a claim to
have implemented the source's global correction and wave-matching scheme.

For numerical conditioning, each local calculation uses s=S X with fixed
S=2q_center. Each right side above is multiplied by S where a radial derivative
is inverted. Time and z derivatives hold this numerical S constant; the actual
coefficient functions in s are evaluated anew at each (z,t). This avoids omitting
coordinate drift when calculating physical derivatives.

## 3. Independent two-variable axis jets

In local scaled coordinates z=z_center+q_center^(1/2-h)*zeta and
 tau=tau_center+q_center*theta, write q=q_center*Q. Then

    Q-(eta_center+zeta)^2 Q^(2h)=1-eta_center^2+theta.

Lagrange inversion gives

    Q^beta = sum_(m>=0) c_m (eta_center+zeta)^(2m)
                            (1-eta_center^2+theta)^(beta+(2h-1)m),
    c_0=1,
    c_m=beta/m! * product_(j=1)^(m-1)(beta+2h*m-j).

This supplies time AND axial Taylor coefficients directly, rather than truncating
the global eta power series used in ST068. The finite sum uses 120 terms and is
checked against 160 terms and independent implicit-coordinate derivatives. At
radial order N=8, the required mixed jets are carried with extra derivative rows.
The radial series and Lagrange sum are finite approximations. This file does not
prove uniform infinite-series convergence, blow-up, or continuation through tau=0.

## 4. Why the leading-only field was insufficient

For ST068-style axis data, the first angular coefficient at eta=0 is

    F_X(0)/g = [h-3+2(1+h)tau^(2h)]/4

in the full recurrence, versus (h-3)/4 in the leading-only recurrence.
At the current h and tau window, the extra axial-viscosity term is not small.
The full radial-pressure equation also changes the pressure gradient. With axis
swirl 1 it points outward over part of the tested core; simply obtaining small
momentum residual is not sufficient for the desired geometry. The independent
swirl-strength control is explicitly distinguished from that fixed-axis test.

## 5. What the residual results can and cannot say

The residual is computed from finite stored/generative fields and all physical
terms; an independent fourth-order Cartesian/time operator is compared at fixed
physical points. Sampled maxima and quadrature L2 values are NOT continuum bounds.
Even if the finite local values are below .001, pde_validated remains false for
the original global project. Local boundary jets, global support, total energy,
outer pressure, non-axisymmetric waves and stress realization are unresolved.

Correcting this tiny local field also changes its interface traces. Old ST068/69
moment numbers and old leading stress-cone formulas must not be reused as if
those traces or the high-order pressure were unchanged. The prior v>2 diagnostic
is still reported to expose the missing interface, not to claim a universal
necessary condition for all local NS solutions.

## 6. Provenance and source scope

Inputs are byte-identical files from NS_ST072_Normalized_Core.zip, SHA256
bdc40184959fdb67a269622df499e1c857500dde2cc6d51052cf3714cfe3e205.
The original reference's Section 5 separately corrects axial viscosity and radial
terms left by leading equations. Appendix B.16 requires a complex-neighborhood
bound, not merely a real-axis amplitude bound. Those statements motivated this
increment but are not independently proved here. Reference:
https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf
The full radial recurrence above was derived directly from the stated local
operator, not copied from or claimed identical to that global construction.
