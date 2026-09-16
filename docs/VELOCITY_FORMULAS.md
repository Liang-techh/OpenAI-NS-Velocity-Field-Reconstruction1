# Cartesian velocity functions: paper source and executable candidate

The immediate user priority is functions of (x,y,z,t), not animation.

## Paper leading field: explicit Cartesian coordinate formula

Source: https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf
Printed page 25 (PDF page index 24), equations (4.3)-(4.5).

For t<1, tau=1-t, A=1/2+h, D=1/2-h, and 0<h<1/100, solve

    q - z^2 q^(-2h) = tau, q>0
    eta = z/q^D
    X = (x^2+y^2)/(2q).

The positive root is unique: the left-hand side has positive derivative
1+2h*z^2*q^(-2h-1), and ranges across tau. The leading velocity is

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
