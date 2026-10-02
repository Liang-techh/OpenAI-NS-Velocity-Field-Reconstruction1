# Fifth-order leading source and radial axial C4

The same implicit leading construction now supplies angular correction,
complete future swirl energy, actual selected ap and both end coefficients
through axial order five. These functions feed all O.4 charts and recover
physical radial velocity and its log-radius derivative through axial order
four. Physical Ur_Z and Ur_yZ are available through axial order three.

This resolves the missing fifth-moment input from the previous checkpoint.
It does not certify every radial/axial mixed derivative or full two-variable
pulse C4. Quantitative flat interface bounds, complete post-pulse C4,
physical kinetic energy, outer stress cone, admissible stress/flat remainder
and actual n-dependent coefficient recursion remain unfinished.

## Files and reproduction

Under experiments/root_st073, with prefix lei_ren_part1_paper_:

- compliant_fifth_axial_jets.py/.json and _check.py/.json;
- compliant_pulse_radial_C4.py/.json and _check.py/.json;
- compliant_reconstruction.py --stage radialjets.

The ordered pipeline has 55 modules. No existing C4 admission guard is
bypassed, and no old source/receipt is rewritten. A separate fifth-order
layer inherits the admitted C4 coefficients and recovers only coefficient
five. Every layer checks transitive source/family hashes. All physical
scales remain exact formal positive sources with directed caps; caps do
not replace those scales or change the source functions.

## Fifth coefficients of the actual implicit branches

For ordinary Taylor coefficients, the two angular functions satisfy

    p*x+y=b1(Z),
    x+q*y+N*(x^2+q*y^2)=b2(Z),
    J=[[p,1],[1+2N*x0,q*(1+2N*y0)]].

The existing C0-C4 functions supply every term of the fifth convolution:

    J*(x5,y5)=(b1_5,b2_5-N*sum_i=1..4(x_i*x_(5-i)+q*y_i*y_(5-i))).

The same nonzero Jacobian and admitted branch are retained. No extra
binomial coefficient is required with factorial-normalized Taylor data.
The physical fifth derivative sum is uniformly below mu/100 over Z in
[-1,1]. Odd fifth angular coefficients vanish at Z=0 by exact parity;
their nonzero heat-driven C0 values are retained.

The fifth RHS includes the true flatten integral and complete Gamma
collar/infinite tails. For m=5 the exact Gamma expectation gives

    |d_d^5 Dhat / 5!| <= (4/15)*(1+a)_4*(1+a)_5*S^4*exp(-5t),
    a=delta/2, d=1-Z^2, S=1/Rtail>0.

This is a finite derivative bound for the true Gamma integral, not an
infinite S-series definition. Composition uses delta_d=-2Z*h-h^2 and all
contributing powers. Infinite angular/pressure/energy tails are integrated
with their actual positive exponential denominators.

## Complete energy and selected ap

Every fifth energy contribution is retained: flatten, power, both signed
angular bumps, steep/waiting, epsilon atoms and the complete Gamma tail.
The q^2 radial prefactors have zero fifth coefficient, but the flatten and
actual bump/Gamma terms generally do not. Section 7.34 keeps the same
mu*exp(-26)/2 normalization.

Exact incoming shapes retain g=Z+Z^3 and h=Z^2+2Z^4+Z^6. Their fifth
ordinary coefficients are zero and 6Z, respectively. Every incoming row
uses the same fixed source factor for all orders. Positive incoming axial
energy and the negative swirl contribution remain separate.

For the actual selected scalar equation A2*ap^2+A1(Z)*ap+A0(Z)=0,

    D=2*A2*ap0+A1_0>0,
    ap5=-(A2*sum_i=1..4 ap_i*ap_(5-i)
          +sum_i=1..5 A1_i*ap_(5-i)+A0_5)/D.

The actual affine end functions remain
cj=exp(log_end_scale)*(uj+vj*ap), with their admitted nonzero formal scale.
There is no nominal amplitude substitution and no imposed even parity for
the axial branch. All admitted C0-C4 prefixes are preserved.

## Radial velocity order bookkeeping

The new inlet and all five normalized partial primitives have axial jets
through five. Original analytic preheat P0 and local P0+Mp pressure are
also retained through five. Kernels and radial supports are unchanged.

Paper (3.9) differentiates m=Mz/(R*Utheta) once in Z. Its recovery operator
therefore maps m5 and b4 to radial velocity axial C4. The normalized moment
ODE m_y=b-(.5-mu)*m, the actual b_y and the physical factor
sqrt(R/2)*Utheta give Ur_y through axial C4. Physical axial derivatives
include Utheta_Z/Utheta=-2Z/(1+Z^2), leaving Ur_Z and Ur_yZ through C3.

At Rv all terminal linear moments and radial/mixed derivatives are exactly
zero because future supports are empty. Energy and its first five axial
coefficients equal half the positive complete future energy. The gap
continues to retain accumulated histories even when axial velocity is zero.

## Verification and limitations

Six fifth-equation identities, six independent stable branch/root fifth
checks, three independent power/log/flatten fifth checks and four true
Gamma fifth-deficit/infinite-tail checks pass. The latter use an
independently bounded omitted infinite tail. Actual C4 prefix inclusion
is checked in 225 comparisons, with the actual angular inverse, selected
positive denominator, end signs, whole-Z fifth smallness and parity gates.

Seventy-two independent direct physical derivatives check Ur, Ur_y, Ur_Z
and Ur_yZ at Z=-1,0,.3,1. They differentiate physical theta*Mz and paper
(3.8), rather than copying a normalized-jet calculation. The actual chart
checker additionally passes 792 lower-order source comparisons, three
high-order/pressure interfaces and whole-Z terminal derivative/energy gates.
Finite-parameter fixtures are explicitly distinct from actual Md=40 source
admission; source-bound continuum enclosures supply that admission.

The reported C4 is AXIAL radial-velocity regularity. Higher y derivatives
of gp/beta and all mixed orders required by stress/C4 have not been
installed. No cell Taylor remainder or whole-field C4 acceptance is claimed.
These fifth derivatives belong to the leading profile; they are not the
n-dependent temporal recovery equations.

## Next dependency

Derive directed higher radial/log-radial gp, sigma and beta derivatives,
including endpoint-crossing boxes and exact flat support limits. Recover
the required Mz and velocity mixed derivatives from source ODEs. Establish
uniform pulse bounds and quantitative flat interface compatibility, then
transport the same data through flatten/angular/steep/waiting/collar/Gamma
regions. Only after these dependencies are available can the full outer
stress cone and stress/flat-remainder lift be admitted.
