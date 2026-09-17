# CR-ROOT-ST001 / ST002 / ST003: full-momentum research

**The requested 1e-3 full NS candidate was NOT obtained. Do not promote this field.**

Issue #205; isolated branch `research/root-st001-full-momentum`, based on active
integration commit `b4e72304a87f596177786ddeebcdfc433f338c0b`.

## Unchanged physical gates

The reference is `../../configs/constraints.json`: nu=0.01, time [0.25,0.75],
physical R3, evaluation box [-2,2]^3, smooth compact velocity AND pressure in
r<2 and |z|<2, E(0.25)=1, original two-parameter compact divergence-free force,
a,c in [0,10], original 1e-3 momentum max and volume-L2 gates. The original
core-sign, 0.05 scaled-core drift, energy, support and divergence gates remain.
No residual-defined force, reduced amplitude, changed viscosity, or relaxed
threshold is used. This is a new explicitly autonomous representation, NOT
a replacement of the Eq45 canonical field or its coefficients.

## Representation and actual optimization

With s=r^2, psi=r^2 F(s,z,t), u_theta=r B(s,z,t):

```
u_r = -r F_z
u_z = 2 F + 2 s F_s
A = -F_z; C = 2 F + 2 s F_s
[u,v,w] = [x A-y B, y A+x B, C]
```

F is odd in z; B and p are even. A common C-infinity bump enforces support.
Spatial tensor Legendre columns are energy-orthonormalized; time uses
Chebyshev polynomials. Initial energy determines amplitude. Largest size is
7x7x8, with 1178 stored parameters including pressure and force parameters;
normalization leaves parameter redundancy, so this is not an independent-DOF count.

ST001 uses three capacity levels (3x3x4,4x4x5,5x5x6), each with a hard 2000
function-evaluation budget. ST002 adds weak anisotropy/torque constraints and
endpoint momentum strata, using 5x5x6 and 7x7x8, also 2000 evaluations each.
ST003 uses a bounded trust-region least-squares solver and a checked analytic
Jacobian: 100 function evaluations, 81 Jacobian evaluations. All runs stop
at their budgets, not a demonstrated local/global minimum.

Training/validation seeds are 9172621/9172622, 9172623/9172624,
and 9172625/9172626. Registrations were written locally before their fits and
uploaded afterward; this is not a claim of an external preregistration service.
Earlier holdouts become diagnostics for subsequent experiments, which use
fresh validation samples. Three early pilot runs reached 2001 evaluations
because SciPy maxfun is soft; a hard guard was added and all three production
baseline fits were rerun. Pilot records are retained only in the audit bundle.

## Final frozen ST003 results: FAIL

4096 fresh uniform Cartesian points, seed9172626; independent fourth-order
Cartesian space/time finite differences. L2=sqrt(64*mean(|R|^2)), not RMS.
The full vector is R=u_t+(u.grad)u+grad(p)-nu*Laplacian(u)-f.

| t | full momentum max | volume L2 |
|---|---:|---:|
| 0.2500 | 0.0450368088 | 0.0916426041 |
| 0.3125 | 0.0572710711 | 0.0904663301 |
| 0.4375 | 0.0354286651 | 0.0655667322 |
| 0.5625 | 0.0432841496 | 0.0665025701 |
| 0.6875 | 0.1316803781 | 0.1542619793 |
| 0.7500 | 0.0684446480 | 0.1015494740 |

Force a=0.19505910862247292, c=0.2753350680165956.
E(0.25)=1.0000000000009126; energies at the six times are in [1,1.20048].
Worst divergence max=1.81518e-6, L2=1.10002e-6; core scaled drift=0.0463000.
Boundary velocity and pressure samples are exactly zero by representation.
All reported non-momentum gates pass; both momentum gates FAIL.

Space steps .02/.01/.005 (time step fixed .0025) give worst momentum maxima
.1316796831/.1316803371/.1316803781: a clear numerical plateau. The finest
analytic-vs-FD discrepancy is at most 7.64e-6 on these samples. Time refinement
is separate at fixed space step .005; energy quadrature uses 24/48/96.
48-to-96 energy relative change is 5.63e-9. These samples are not a uniform
continuum space-time residual certificate.

ST002 had max .1375902173 and L2 .1525902202 on its own seed. ST003 has slightly
lower sampled max but slightly higher L2; different validation seeds mean
this is not a controlled same-sample improvement percentage.

Original final candidate JSON SHA256:
`11624f669120e9e226db61dcdf1bec77dde5a06568c479d31724fe14cc70e4c4`.
The frozen JSON, all production/pilot histories and full validation JSONs are
in the accompanying user-delivery research archive. This PR tracks the code,
registrations, results and complete regeneration command; the bulky frozen
run directories are not tracked here. `reproduce.py` regenerates those paths.

## New pressure-independent obstruction

For any harmonic polynomial H and compact smooth divergence-free u, compact
p, and the prescribed divergence-free force, integration by parts gives

```
 integral grad(H).R = -integral u^T Hess(H) u.
```

The time, pressure, viscosity and force terms vanish in this identity.
For H=(x^2+y^2)/4-z^2/2, phi=grad(H) and
D=integral(uz^2-(ur^2+u_theta^2)/2), it follows that

```
 ||R||_L2 >= |D| / sqrt(88*pi/3).
```

The inequality is exact under its hypotheses. Numerical D values are
quadrature estimates, not interval-certified numerical lower bounds.
ST001 at t=.25 had a lower-bound estimate .05394. Adding this necessary
condition in ST002 substantially repaired the degree-2 imbalance, but not
higher-degree pressure compatibility.

For harmonic degrees 2,4,6,8, define D_i=-integral u^T Hess(H_i)u and
G_ij=integral_support grad(H_i).grad(H_j). Then

```
 ||R||_L2 >= sqrt(D^T G^{-1} D).
```

For the final ST003 velocity at t=.25, this joint lower-bound estimate is
**0.0247955861**. Across the six times it is .0167381 to .0247956. The largest
48-to-96 quadrature change is 4.81e-10; normalized Gram condition number=5.42.
Thus fitting only pressure, or changing only a,c within the original force
family, cannot repair the current frozen velocity to 1e-3. This is not a
proof that every allowed velocity family is infeasible. Interval quadrature
would be required to certify the numerical lower-bound values rigorously.

The original force also has total axial torque -c*g(t)*integral r^2*b dx;
positive central swirl force does not imply positive net torque.

## Actual checks and replay

Local focused regression run: **7 passed in 4.28s**, warnings treated as errors.
Five symbolic structural checks pass; harmonicity for degrees2/4/6/8 is also
checked symbolically. The analytic Jacobian passes a directional FD check.
The acceptance gate exits1 and prints `FAIL: momentum_max, momentum_L2`.
No full inherited repository suite, Lean build, or continuum NS proof was run.
CI status must be checked independently; a green software test is not PDE acceptance.
The torch-dependent Jacobian regression is skipped when torch is unavailable.

From this directory:

```bash
python -m pip install -r requirements.txt
python -m pytest -q -W error test_spacetime.py
python symbolic_checks.py
python reproduce.py --out rerun
python check_acceptance.py rerun/gauss_7x7x8/validation.json
```

The last command should fail the scientific gates. Never lower thresholds to
make it green. `reproduce.py` refuses to overwrite a nonempty output directory.
To evaluate a frozen or regenerated candidate without the optimizer:

```python
import numpy as np
from spacetime import Family, force
family, params = Family.load('rerun/gauss_7x7x8/candidate.json')
points = np.array([[0.1,0.,0.1], [0.,0.,0.2]])
u, p = family.fields(params, points, 0.5)
f = force(points, 0.5, *params[-2:])
```

NumPy and SciPy suffice for evaluation/FD validation. Torch is only for fitting,
SymPy for symbolic checks. The supplied requirements record the run environment,
not a promise of bitwise reproducibility across BLAS/platforms.

## Remaining work / truth boundary

Construction must address higher harmonic moment compatibility jointly with
full momentum, energy, torque and core structure. Allowing a noncompact
Poisson pressure tail would be a NEW physical experiment, not a pass of this
compact-pressure contract. No canonical field or existing truth flag is
changed. `pde_validated`, `paper_exact`, and `blowup_proved` remain false.
Scientific target: NOT ACHIEVED. Review acceptance: pending.
