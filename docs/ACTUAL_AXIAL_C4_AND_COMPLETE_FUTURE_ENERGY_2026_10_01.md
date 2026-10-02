# Actual axial coefficients and complete future energy through order four

The same compliant .001 leading family now supplies axial Taylor coefficients
through order four for the complete corrected future swirl energy, actual
pulse-inlet moment/energy functions, selected positive amplitude ap(Z), and
both actual affine end coefficients c1(Z), c2(Z). Angular coefficients and
true Gamma derivatives were already admitted in the previous checkpoint.

This completes F38-D3b1 through D3b3. These derivatives belong to the leading
profile. They are not the n-dependent temporal coefficient recursion.
Full pulse/outer C4, mixed velocity derivatives, physical kinetic energy,
the outer stress cone, global admissible stress and flat remainder remain
unfinished. No fifth-derivative Taylor remainder is provided.

## Reproduction and source identity

Under experiments/root_st073, with prefix lei_ren_part1_paper_:

- compliant_future_energy_high_jets.py/.json and _check.py/.json;
- compliant_axial_high_jets.py/.json and _check.py/.json;
- compliant_reconstruction.py --stage axialjets.

The ordered pipeline now has 49 modules. Previous receipts remain unchanged.
Both new layers check prerequisite hashes and the same five-defect family.
Each returned coefficient is derivative divided by factorial. A whole-Z
interval encloses derivatives at every real center, rather than serving as
a Taylor expansion about an interval or as a certified cell approximation.

## Complete future swirl energy

Every radial contribution is included: flatten, power buffer, both signed
angular bumps, steep entry/power/exit, waiting, epsilon collar and the entire
Gamma heat tail. The two quadratic swirl prefactors retain their full
(1+Z^2)^2 dependence. The final Section 7.34 weight is mu*exp(-26)/2 in
Rv*Utheta(Rv,Z)^2 units. C0/C1 are intersected with prior enclosures of the
same source. Odd coefficients at Z=0 remain zero for this even energy.

Six source/normalization identities and ten independent integrand/physical
bump derivative checks pass, together with actual C1 inclusion and parity.
The independent finite-parameter fixture is explicitly distinct from the
actual Md=40 source admission. This is remaining swirl energy; it does not
claim the full physical kinetic energy including radial and axial velocity.

## Exact incoming functional shapes

Let r(Z)=1/(1+Z^2). The reference, slope, axial turnoff and both buffers
preserve these shapes in their existing normalized units:

    u=U*r, m=M*Z, k=K*Z*r, e_raw=E_Z*Z^2+E_Q*r^2.

U, M, K and E_Z are positive Z-independent source constants; E_Q is negative.
The actual pulse-inlet normalized functions are therefore

    m1=C1*(Z+Z^3), m2=C2*(Z+Z^3),
    e=C0+C_E*(Z^2+2Z^4+Z^6),
    C1=exp(-logPstar)*M/U, C2=exp(-logPstar)*K/U^2,
    C0=E_Q/U^2, C_E=E_Z/U^2.

There is no additional inverse pressure factor in normalized incoming energy.
U, M, K and E_Q are isolated from the existing packet at Z=0. E_Z is
accumulated from its positive source, avoiding cancellation against E_Q:

    yd=exp(Md)+11, td=yd-1,
    K_B2(yd)=integral_1^yd exp(s-yd)*B(log(s)/Md)^2 ds,
    E_Z=16*Pstar^-2*(exp(-td)+K_B2(yd))*exp(-1-Tw).

The kernel already includes the eleven-unit zero-axial buffer attenuation.
Its omitted positive tail remains enclosed, not zeroed. The factor sixteen
comes from (4ZB)^2. Every ordinary coefficient of each incoming moment is
multiplied by the same inherited fixed positive row factor enclosure for
exp(-13*lambda_i/mu-common_logpref), including orders two through four.

Nine independent functional identities check separated energy propagation,
normalizations and polynomial coefficients. Actual C0/C1 packet agreement,
constant signs and source hashes are also checked.

## Differentiating the selected amplitude

The exact affine end functions remain c_j=exp(log_end_scale)*(u_j+v_j*ap).
The positive scale is kept in formal logarithmic form. The selected amplitude solves

    A2*ap^2+A1(Z)*ap+A0(Z)=0,
    A2=Kpulse+sum(nu_j*v_j^2),
    A1=2*sum(nu_j*u_j*v_j),
    A0=sum(nu_j*u_j^2)-target(Z).

Here nu_j=mu*K_j*exp(2*log_end_scale) is exactly positive; finite caps only
enclose it. A2 is Z-independent. At each higher ordinary coefficient n,

    D=2*A2*ap_0+A1_0>0,
    ap_n=-(A2*sum_i=1..n-1 ap_i*ap_(n-i)
           +sum_i=1..n A1_i*ap_(n-i)+A0_n)/D.

The already admitted actual C0 root and C1 enclosure are retained. No nominal
amplitude, midpoint selection or derivative of an interval root iteration
is substituted. Incoming affine terms can be odd in Z; no even-parity reset
is imposed on ap or the end coefficients. Fifteen independent direct
positive-root/end-function derivative checks pass through order four, along
with the actual whole-Z positive denominator and c1<0<c2 signs.

## Next dependency

Install these selected jets in every pulse chart and all five partial
primitives. Recover Ur_Z and radial/axial mixed derivatives from the same
Mz primitive and paper (3.9). Propagate through the post-pulse field, bound
all required interfaces, and only then certify whole-field C4 and the outer
stress cone. Continue to admissible stress/flat remainder and actual n=1,
n>=2 recovery equations before claiming temporal scale recursion.
