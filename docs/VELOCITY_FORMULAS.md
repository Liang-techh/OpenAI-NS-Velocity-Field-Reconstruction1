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

## Nonlinear profile equations and baseline error

Let avgU be the radial average of U, d=1-eta^2, L=1-2h*eta^2.
Equations (4.8), (4.9), (4.13) reduce to

    W = 1-2D*eta*avgU-d*partial_eta(avgU)
    Hc = D*eta+d*U
    Sq = -W*(1+X*F_X/F)-h*(1-2*eta*U)-Hc*F_eta/F
    Sn = -W*X*U_X-A*(1-2*eta*U)*U-Hc*U_eta
         -d*Pi_eta+4*A*eta*Pi+2*eta*X*Pi_X
    -2L*(X*F_XX+2*F_X)/F = Sq
    -2L*(X*U_XX+U_X) = Sn
    Pi_X = F^2.

The new independent finite-difference evaluator paper_profile_residual.py
checks these equations, not the full NS system. The B.13 reference with the
recorded defaults gives maximum absolute angular/axial equation errors
9.08389245 / 4.40968923 on 20 interior points. Halving difference step from
1e-4 to 5e-5 changes them by less than 3e-7. Its radial pressure relation
has error below 1.4e-12. Thus the reference is not already a nonlinear solution;
the remaining defect is substantive, not a derivative-resolution artifact.
Raw points and both steps are saved in core_reference/residual.json.

## Nonlinear radial-series delivery

```python
from openai_ns_reconstruction.paper_core_reference import PaperCoreReference
from openai_ns_reconstruction.paper_core_series import PaperCoreSeries
field = PaperCoreSeries(PaperCoreReference(sigma=.3), maxdegree=12, eta_nodes=513)
u, v, w = field.velocity(.1, 0, .1, .25)
```

The new module solves the triangular X-series recurrence from (4.9)/(4.13),
retaining normalized F/g and U through radial degree N. Eta differentiation
uses Chebyshev nodes; g and its log derivative are evaluated from the reference
formula at actual requested coordinates. Pressure integrates the entire square
of the retained F polynomial, including off-grid points, so it has degree 2N+1.
Four focused tests pass for eta differentiation, first-row identities,
vectorized evaluation and pressure/axis consistency. The parent reviewed and
corrected an extra source factor and off-grid pressure integration before
running the recorded experiments.

The default sigma=.1 has poor eta resolution and unstable high radial orders.
The explicit alternative sigma=.3, degree12, 513 eta nodes gives the following
maximum leading-equation errors over 39 off-grid eta points in [-.95,.95]:

| X | Angular | Axial |
|---|---:|---:|
| .01 | 4.13e-9 | 5.67e-9 |
| .1 | 6.43e-8 | 1.10e-7 |
| .2 | 1.50e-5 | 1.92e-4 |
| .3 | .00261 | .0168 |
| .4 | .122 | .403 |

At X=.1 the same-parameter uncorrected reference has errors 5.54 and 3.81.
Two derivative steps confirm the improved order of magnitude, with roundoff
visible near 1e-7. These are development profile-equation samples, NOT full
physical NS acceptance, a uniform convergence proof, exterior matching, or a
verified numerical match to the announcement image. Higher order16 can be much
worse and even lose positive swirl in some parameter runs; failed runs remain
in sigma_study.json and roundoff_study.json. Do not simply increase degree.

Artifacts under artifacts/function_first/core_series include all degree/node
runs, independent parameter exploration, and selected/coefficients.npz,
selected/samples.json, selected/probe.json, selected/comparison.json.
Reproduce the selected spatial probe from the repository root:

```sh
python -m openai_ns_reconstruction.paper_core_experiment --sigma .3 --degrees 12 --nodes 513 --radii .001 .01 .05 .1 .2 .3 .4 --eta-count 39 --output artifacts/function_first/core_series/reproduced
```

Set PYTHONPATH=src as in the API instructions. All model parameters remain
independent choices; no theorem threshold or force restriction was relaxed.
Next work: stabilize radial continuation beyond the reliable inner region,
then construct the outer connection while preserving incompressibility.

## Extended inner range: selected independent parameters

A bounded study of sigma=.5,.7,1 and radial degrees8,12,14 improves the
usable inner profile without changing the equations. The selected new seed is

```python
field = PaperCoreSeries(PaperCoreReference(sigma=.5), maxdegree=14, eta_nodes=257)
u,v,w = field.velocity(.1,0,.1,.25)
```

At 39 eta points in [-.95,.95], two derivative steps give consistent leading
residuals. With step5e-5, X=.3 gives2.51e-6/4.05e-7, X=.4 gives
1.66e-4/2.14e-5, and X=.409 gives2.30e-4/2.91e-5 (angular/axial).
All312 sampled points across8 radial levels are finite and have positive F.
These are selected development samples, not uniform or full NS acceptance.
The spatial graph correspondence of the changed parameter is not yet verified.
The new coefficient/sample files are under
artifacts/function_first/core_series/smooth_parameters/selected.

Reproduce:

```sh
python -m openai_ns_reconstruction.paper_core_experiment --sigma .5 --degrees 14 --nodes 257 --radii .001 .01 .05 .1 .2 .3 .4 .409 --eta-count 39 --output artifacts/function_first/core_series/smooth_parameters/reproduced
```

Two other continuation approaches were implemented and rejected as the default:

- paper_core_continuation.py integrates a first-order radial system using
  DOP853 from X=.03. With sigma=.3, increasing eta resolution129->257 worsens
  outer residuals and causes negative swirl at some points.513 nodes stops with
  a step-size-underflow error. The seed-to-integrator switch can also have an
  off-grid pressure interpolation mismatch, so this exploratory class is not
  a validated piecewise field. No damping or external force was introduced.
  Reproduce with python -m openai_ns_reconstruction.paper_continuation_study;
  failures are saved incrementally in core_continuation/resolution.json.
- paper_core_pade.py evaluates rational approximants and rejects denominator
  poles on the entire radial integration path. Degrees6/8/10/12 all reject
  some tested points; smaller residuals among surviving points are not success.
  core_pade/comparison.json records both residuals and rejected points.

Next: use the selected series as an inner seed and construct a smooth exterior
through a streamfunction/potential, preserving incompressibility. Exterior
pressure, restricted-force compatibility, full NS residual and visual-source
correspondence still need separate evidence.
