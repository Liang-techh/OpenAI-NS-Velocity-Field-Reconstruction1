# Cartesian velocity functions: paper source and executable candidate

The immediate user priority is functions of (x,y,z,t), not animation.

## Paper leading field: explicit Cartesian coordinate formula

Source: https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf
Printed page 25 (PDF page index 24), equations (4.3)-(4.5).

For t<1, tau=1-t, A=1/2+h, D=1/2-h, and 0<h<1/100, solve

    q - z^2 q^(2h) = tau, q>0
    eta = z/q^D
    X = (x^2+y^2)/(2q).

The physical root satisfies q>|z|^(1/D). On that interval the derivative
is 1-2h*z^2*q^(2h-1)>=1-2h>0, establishing uniqueness. The leading velocity is

    u(x,y,z,t) = x*v0(X,eta)/(2q) - y*q^(-1-h)*F(X,eta)
    v(x,y,z,t) = y*v0(X,eta)/(2q) + x*q^(-1-h)*F(X,eta)
    w(x,y,z,t) = q^(-1/2-h)*U(X,eta).

Here E=sqrt(2X)*F, V0=X*v0, and F=phi/C. These formulas avoid division
by radius and extend regularly to the axis if the profiles have the paper's
smoothness: at x=y=0, u=v=0 and w=q^(-A)*U(0,eta).

This identifies the coordinate formula, not numerical values of the profiles.
The full constructed field includes further corrections; the leading field
alone must not be called the complete theorem solution. Profile construction
and concrete parameter choices remain the next implementation task.

## Already executable independent candidate

The saved coupled candidate has all coefficients bundled and is callable now:

```python
from openai_ns_reconstruction.velocity_components import velocity
u, v, w = velocity(0.1, 0.0, 0.1, 0.5)
# (-0.017550887073245407, 0.07337970151893687, 0.034940231516750846)
```

Its specification and parameter provenance are in VELOCITY_API.md. It is an
independent approximation, not the identified numerical OpenAI profile.
It remains outside the configured PDE acceptance thresholds.

## Computable paper-based near-axis reference (FUN002 partial delivery)

```python
from openai_ns_reconstruction.paper_core_reference import PaperCoreReference
field = PaperCoreReference()
u, v, w = field.velocity(0.1, 0, 0.1, 0.25)
# [-0.25855442999451683, 0.0015492527377398112, 0.560607549191811]
```

This reuses coordinates.py, profiles.py and velocity.py. It implements the
reference term in B.13, not the nonlinear fixed point:

    Ustar=4*eta+j; Hstar=D*eta+(1-eta^2)*Ustar
    chi=Hstar^2/(Hstar^2+sigma^2)
    g=exp(Lambda*integral_0^eta[-L*Hstar/(Hstar^2+sigma^2)]de)/C
    f0(s)=sum_n [(-s/2)^n/(n!*(n+1)!)]
    F(X,eta)=g(eta)*f0(Lambda*X*chi(eta))
    U(X,eta)=Ustar(eta)-X*Zstar(eta)/(2*L(eta))

Zstar follows B.1 with the independently chosen pressure datum
Pi0=-pressure_scale^2/(1+eta^2)^2. Defaults h=.005, j=.02, sigma=.1,
Lambda=10, C=2, pressure_scale=1 are explicit autonomous finite choices.
They are NOT demonstrated to satisfy the theorem's existential thresholds
or matched to the constructed exterior pressure. All three velocity components
then use (4.5) and (4.7), with exact radial averages of affine U and an analytic
eta derivative. The domain is 0<=t<1 and 0<=Lambda*X<=4.1. Calls outside
that radial domain raise an error; no exterior is silently filled with zero.
The initial field is nonzero, so this is not the theorem's zero-initial-data field.

Samples at three times and both axial sides, plus parameter metadata, are in
artifacts/function_first/core_reference/samples.json. Two focused tests pass
for the B.11 series, correct q coordinate identity, axis limits, axial derivatives,
and numerical divergence. This does not establish momentum-equation acceptance.

## What the source leaves to construction

The paper defines E and U through axis construction, continuation and moment
matching; it supplies no final numerical coefficient table. V0 is determined
from U by (4.7), so it is not a third independent profile to fit. Appendix B
prescribes axis data and a fixed-point construction, with thresholds Lambda0
and C0(Lambda) rather than numerical final choices. Theorem 4.6 and Appendix A
also require joining to an exterior. Sections 5 and 7-10 add further corrections
to obtain the complete field. Next work must solve the nonlinear profile equations
and compare against this reference, then address the exterior; the present
reference must not be relabeled as those completed constructions.

Correction: the prior documentation mistakenly wrote q^(-2h) in the coordinate
equation. Equation (4.1) requires q^(2h); the inherited numerical solver already
used that correct positive exponent.
