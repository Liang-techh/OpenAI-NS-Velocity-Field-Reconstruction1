# Original steep exit: full-future stress and pressure — 2026-10-04

The original steep_out transition now has a companion that recovers the same complete angular, energy and absolute-pressure moments, gives signed similarity stresses through mixed order 3, and restores pressure through mixed order 4. It joins the accepted waiting companion by exact source-functional identities. Physical stress/remainder transfer and the steep-exit directional cone remain separate unfinished steps.

## Source and normalization

Domain: original steep_out phase t in [0,1], source Z in [-1,1]. Physical Z endpoints are infinity limits, as in the accepted source map. The reference is unchanged:

```
q = log(R/Rtail) = -wait-1+t
a = delta/2, k = 1-a, p = 1+delta, bh = 1/2+a
K0 = 1-epsilon, epsilon = .001*delta
f(t) = integral_t^1 (1-sigma(v))dv = J(t)-t+1/2
K(t) = K0*exp(k*f(t))
B = Ev0*theta_base*exp(-bh*q)
A = K*X_original
E = 2*K^2*energy_original
P = -pressure_absolute/(C*exp(-p*q))
C = (Ev0^2/Pstar^2)*theta_base^2
```

The exact Ev0 source is retained; the runtime Ev2 interval encloses its squared scale. The bridge explicitly consumes the accepted actual inlet/amplitude, pressure-scale, absolute-datum and complete defining/history identities. It executes the original steep_out theta/X/energy/pressure expression ASTs and their FTC equations, and replays the production radius branch against the original waiting tail.

The same actual collar endpoint supplies the waiting left moment defects at q=-wait. No collar source is evaluated at negative q. Backward finite integrals extend that already complete sigma/phi/Gamma future:

```
Ad(t) = exp(k*(1-t))*Ad1
        - integral_t^1 exp(k*(v-t))*(K(v)-1)dv
Ed(t) = exp(-delta*(1-t))*Ed1
        + integral_t^1 exp(-delta*(v-t))*(K(v)^2-1)dv
Pd(t) = exp(-p*(1-t))*Pd1
        + 1/2*integral_t^1 exp(-p*(v-t))*(K(v)^2-1)dv
A = 1/k+Ad, E = 1/delta+Ed, P = 1/(2*p)+Pd
```

Signed integrands use expm1 to retain small defects before cancellation. Reverse positive-cell accumulation encloses f and the three integrals. Cell lengths are (1-t)/cells and v-t uses the correlated distance (1-t)*u. Exact sigma reflection gives f(0)=1/2 and 0<=f<=1/2. Endpoint t=1 has exactly zero remaining integrals.

Original forward X, energy, pressure and their derivative histories remain accessible under original_forward aliases. Equivalent full-future representations replace the corresponding companion rows, including X=A/K and energy=E/(2*K^2) derivatives. The original velocity field is unchanged.

## Stress and joins

The accepted general-K (3.16)-(3.18) recovery is applied with the same positive factors Qtheta=sqrt(R/2)*B and Qz=sqrt(R/2)*B^2. Exact unit baselines cancel before enclosure. K ordinary derivatives through 4 and common moment FTC rows supply the stress mixed3 grids; pressure uses the full-datum FTC through mixed4.

At t=1, the original sigma is flat, K=K0, and every positive K derivative through 4 vanishes. Actual production formula ASTs, executed on arbitrary terminal axial functions, establish 20 stress and 15 pressure mixed joins to waiting. This is a functional join, not sampled interval overlap.

The source shear sign follows from exact S=exp(-actual finite logRtail)>0 and K>0:

```
K_y/K = k*(sigma-1)
kappa-2 = delta+2*k*(1-sigma) >= delta > 0
K_y-(1+a)*K < 0
```

The numerical S box can include zero. Its lower endpoint does not define the source. Directional stress-cone inequalities for steep exit have not yet been proved. The previous waiting/collar cone domain is unchanged.

## Reproduction and evidence

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage steepexitstress
```

The ordered producer/checker passed with current sources and 239 checker input hashes:

- 80 finite actual signed stress coefficients and 60 pressure coefficients, including whole t/Z enclosures.
- 96 exact meridional primitive/velocity zeros, consuming actual selected terminal histories and the original zero-density FTC route.
- 24 flat endpoint K coefficient checks.
- Functional source normalization, backward FTC, exact amplitude/datum/radius links, and actual formula joins.
- One bounded independent moderate fixture: 8 original-sigma signed-integral checks, 40 unnormalized original stress mixed3 checks, 30 pressure mixed4 checks, 30 complete-future moment derivative checks and 20 original normalized X/half-energy derivative checks.

The independent fixture uses direct quadrature of the original sigma transition and a convergent analytic future beyond a constant waiting segment. It differentiates the full original moment/stress formulas using local scalar primitives. Its numeric comparison tolerance is 2e-7; its interval-containment error is not a global NS accuracy claim. Actual full Gamma history is consumed separately through accepted current-source receipts.

## Ready next work

1. Implement steep_exit_physical_C2 with an original steep_exit phase selector and the accepted general-K physical mapper. Preserve q=-wait-1+t and exact nu/lambda units. Complete the symmetric tensor diagonal Ttheta_theta=r*partial_z(Tz).
2. Retain Etheta=-nu*partial_zz(utheta), generally nonzero here. Differentiate the actual source phase/radius map; do not extend waiting's pure-radial-power zero-remainder identity.
3. Add one bounded independent native Cartesian decomposition with the complete tensor, pressure datum, variable K and positive viscosities. Admit only the regional physical identity after it passes.
4. Prove the steep-exit whole-domain two-vector cone with the actual variable kappa-2 and B source logs, then compose it with the accepted waiting/collar tail.
5. Continue the same full moments/stresses back through steep power and entry, and then through angular/power/flatten, with exact endpoint functions and all source joins.
6. Independently close actual upstream finite-width bridge feedback, required-domain energy and global flat/volume bounds before common leading-input/global completion.

Global stress, global flat remainder, physical energy, genuine n=1/n>=2 recursion, oscillatory correction and full corrected NS validation remain unfinished.
