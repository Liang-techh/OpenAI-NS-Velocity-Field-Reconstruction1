# Selected pulse: C4 partial primitives and C3 radial velocity

The admitted high-order amplitude and end coefficients now feed every O.4
entrance/main/inactive-gap/end chart. All five normalized partial primitives
and axial/swirl coefficients have axial Taylor jets through order four.
The original analytic preheat pressure P0 remains included through order four.

Paper (3.9) now supplies the normalized radial velocity through order three,
its log-radius derivative through order three, and their physical axial
derivatives through order two. The construction keeps the same Mz primitive
for divergence recovery; radial velocity is not fitted independently.

This is a partial implementation of F38-D3c1. Full velocity C4 is deliberately
not claimed: the radial formula consumes one axial moment derivative. The
next required input is Mz through order five, followed by higher radial
derivatives and uniform interface bounds. Full post-pulse C4, stress cone,
admissible stress, flat remainder and temporal recursion remain incomplete.

## Files and reproduction

Under experiments/root_st073, with prefix lei_ren_part1_paper_:

- compliant_pulse_high_jets.py/.json;
- compliant_pulse_high_jets_check.py/.json;
- compliant_reconstruction.py --stage pulsejets.

The ordered pipeline has 51 modules. Existing C1 pulse and C4 coefficient
receipts remain unchanged. The new provider checks their transitive hashes
and common compliant .001 source/five-defect family. Every scalar, cap and
normalization is initialized from its high-order source graph in one interval
context. The old segmented pulse formulas are reused as the actual source,
not reconstructed from a midpoint or nominal ap value.

## Incoming histories and pressure

The high-order inlet retains the exact functional shapes

    u=U/(1+Z^2), Mz/R=M*Z,
    Mtheta_z/(sqrt(2)*R^1.5*Pstar)=K*Z/(1+Z^2),
    Mztheta/(R*Pstar^2)=E_Z*Z^2+E_Q/(1+Z^2)^2,
    Mtheta/(sqrt(2)*R^1.5*Pstar)=H/(1+Z^2),
    Mp/Pstar^2=Pm/(1+Z^2)^2.

The Z-independent H and Pm are isolated from the existing source packet at
Z=0; U/M/K/E_Z/E_Q come from the admitted actual axial high-jet source.
Analytic P0 uses the original fourteen-stage normalized pressure datum API
through order four, including its enclosed variable-flattening derivatives.
The local pressure remains P0+Mp. No tail adjustment or source replacement
is introduced. Its log-radius derivative follows the original radial balance.

All radial kernel factors, support primitives and formal tiny scales remain
independent of Z. Higher jets therefore propagate through the same actual
partial pulse equations. Forward energy is used on the ordinary main chart;
the selected backward energy identity avoids large cancelling boxes near its
end. The inactive gap preserves moment histories even though axial velocity
is zero. At Rv, empty future supports enforce exact terminal linear and
radial zeros while the energy is half the positive complete future swirl
energy, including every axial derivative of that identity.

## Radial and mixed derivatives

Write y=log(R), b=Uz/Utheta, m=Mz/(R*Utheta), and define

    A = [2Z*b-(1-delta)Z*m
         -(1-Z^2)*(m_Z-2Z*m/(1+Z^2))]/(1-delta*Z^2).

Then Ur=sqrt(R/2)*Utheta*A. The normalized moment ODE gives

    m_y=b-(.5-mu)*m,
    b_y=Uz_y/Utheta+(.5+mu)*b,
    Ur_y/(sqrt(R/2)*Utheta)=A[b_y,m_y]-mu*A[b,m].

The notation A[b,m] means the same linear paper recovery operator. The
physical axial derivatives include the swirl derivative factor:

    Ur_Z/(sqrt(R/2)*Utheta)=A_Z+(Utheta_Z/Utheta)*A,
    Ur_yZ/(sqrt(R/2)*Utheta)=(A_y-mu*A)_Z
                             +(Utheta_Z/Utheta)*(A_y-mu*A).

Utheta_Z/Utheta=-2Z/(1+Z^2). Omitting this product-rule term would
differentiate only the normalized radial field, rather than the velocity.
Ordinary Taylor differentiation multiplies coefficient n+1 by n+1. Thus
order-four m yields order-three Ur and order-two Ur_Z/Ur_yZ. The exposed
factored velocity jets retain the common exact positive radial exponential
separately from their inlet/axial Taylor factors.

## Verification scope

The checker binds the production radial numerator to paper (3.9), derives
the log-radius equation independently from (3.8), and checks physical axial
product rules and the exact pressure-shape coefficients. Independent direct
physical Mz differentiation checks Ur, Ur_y, Ur_Z and Ur_yZ at the axis
endpoints, zero and an interior axial center. This finite-parameter fixture
is marked separately from the actual Md=40 source admission.

Actual checks cover C4/C3/C2 array orders, old-source C0/C1 agreement, all
three pulse-chart overlaps, P0+Mp restoration, whole-Z terminal derivative
zeros and the positive future-energy target. Whole-Z end-support and
support-crossing packets are retained. These checks do not certify all
radial mixed orders, full-field C4 or a fifth-order Taylor remainder.

## Next bounded tasks

Recover actual amplitude/end/moment fifth-order axial input so radial C4
can be constructed. Derive higher log-radius and mixed derivatives of the
true gp/beta functions with support-boundary enclosures. Prove quantitative
flat interface bounds, then propagate the admitted jets through the complete
post-pulse flatten/angular/steep/waiting/collar/Gamma exterior. Only after
those derivatives are available should the whole outer stress cone be tested.
