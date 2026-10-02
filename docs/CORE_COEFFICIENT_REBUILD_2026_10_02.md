# Fresh compliant core coefficients with controlled radial remainder

The coupled radial core coefficients are now recomputed directly from the
current `.001` analytic preheat pressure datum and the admitted selected
Cstar. The new implementation does not read, shift or relabel old finite
coefficient rows. It provides finite source coefficient enclosures and
controlled infinite radial remainders. It does not select exact parameter
representatives, materialize F0, or implement temporal coefficient recursion.

## Source seed and units

With `rho=Lambda*R` and `epsilon_core=1/Lambda`, the recurrence uses

```
U0=4Z+j; d=1-Z^2; L=1-delta*Z^2
H=d*U0+(1-delta)*Z/2; sigma=j/500
G_Z=L*H/(H^2+sigma^2)
ell=epsilon_core*(F0_Z/F0)=-G_Z
S=epsilon_core^2*F0^2
P0_scaled=exp(2logPstar-logLambda)*p0_normalized(Z)
```

The full ordinary axial Taylor jets are rebuilt to `radial_degree+depth+1`
coefficients. The pressure jets come from CompliantPressureDatum and retain
its fourteen source atoms. The pressure source's `.001*delta` epsilon is
distinct from `epsilon_core`; it is never passed as the core recurrence scale.
Likewise, the physical logarithmic derivative `-Lambda*G_Z` is not passed to
the scaled recurrence: that would introduce an erroneous Lambda factor.

Selected Cstar satisfies the admitted complex amplitude guard. This gives
`|S|<=exp(-6logLambda-2000)` on the same complex axial tube. Cauchy bounds
retain S and its derivatives as nonzero source enclosures. Their numerical
intervals include zero, but the underlying implicit F0 is strictly positive;
the cap is not a selected amplitude value.

The appended rows have units

```
A[n] = Phi_n in rho
Uz[n] = physical axial-profile coefficient in rho
P[n] = epsilon_core * physical-pressure coefficient in rho
```

All frozen nonlinear, pressure, angular and axial couplings in the original
radial equations are retained. Initial degree24 runs at Z=.3,.5,-.5,0 produce
5,400 finite scalar coefficient enclosures. The evaluator accepts other real
axial centers and supported derivative orders; these initial runs are not
an exhaustive whole-axis numerical grid.

## Controlled infinite remainder

For mixed orders through total5, the Phi remainder is the explicit Bessel
model's factorial tail plus the admitted nonlinear fixed-point correction
Xh tail. The Psi model is linear in rho and is fully present for N>=1.
Consequently the physical Uz and its radial average have tails bounded by
`epsilon_core * correction_Xh_norm * tail_factor`. The raw Psi tail must not
be used without that outer epsilon_core. Higher requested mixed orders use
conservative admitted full analytic norms.

The radial average uses the same axial rows divided by n+1. The original
radial recovery then gives

```
Q = [2Z*Uz-(1-delta)Z*(Mz/R)-(1-Z^2)*partial_Z(Mz/R)]/L
Ur=sqrt(R/2)*Q
F=F0*Phi; Utheta=sqrt(2R)*F
```

Pressure is integrated from that same swirl. From the coefficient norm,
`|(Phi^2)_n|<=Bphi^2*(n+1)/20^n`. Therefore the pressure remainder after
degree N has the conservative bound

```
|tail(P_scaled)| <= S_bound*Bphi^2*rho^(N+1)
                   /[20^N*(1-rho/20)]
```

The known axis pressure is retained; no new pressure tail is appended to
reduce a residual. Finite polynomials and every remainder bound are returned
separately even when a bound is below the displayed numerical precision.

## Independent checks and reproduction

A finite nonzero fixture solves the original unscaled coupled radial
equations independently. Ninety coefficients verify the scaled recurrence,
F0 factorization and epsilon pressure units. Independent coefficient-norm
sums check sixty mixed tail cases; a saturated sequence checks the pressure
convolution/integration bound. The fresh-source check also verifies 48
original first-radial axis coefficients and reproducibly re-evaluates sixteen
core value packets. Interval overlap on the first-row axis identities is
diagnostic; it is not used as a functional source proof.

```python
from lei_ren_part1_paper_compliant_core_coefficient_rebuild import (
    CompliantCoreCoefficientRebuild,
)
core=CompliantCoreCoefficientRebuild()
rows=core.rebuild(Z='.3', degree=24, depth=5)
value=core.values(rows, rho='4')
mixed=core.profile(rows, rho='2', radial_order=1, axial_order=2)
```

Focused pipeline stage, now part of 109 ordered modules:

```
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage freshcore
```

The next bottleneck is usable source values throughout the chain: implicit
F0/pressure uncertainty, signed bridge and switch integrals, and the admitted
five-bump branch. After these are resolved, compose a physical-point evaluator
and measured dynamics. Required-domain energy, admissible stress, independent
flat remainder, temporal coefficient recursion and oscillatory correction
remain unfinished. The unlocalized whole-space energy remains infinite.
