# Original steep-entry full-future stress and pressure — 2026-10-04

Update: the regional physical decomposition and native physical connections are now accepted in docs/STEEP_ENTRY_PHYSICAL_2026_10_04.md. The similarity receipt below retains its original scoped flags; the entry cone and upstream angular stress remain open.

The original sigmoid entry now has a callable similarity-stress companion on
the entire original t in [0,1], Z in [-1,1]. It transports the SAME complete
power/exit/waiting/collar/Gamma future, preserving the original velocities,
repair coefficients, source amplitudes and analytic absolute-pressure datum.
This closes the entry similarity stress/pressure increment. Its regional
physical decomposition and cone are still open.

## Implementation and reproduction

- Producer: experiments/root_st073/lei_ren_part1_paper_compliant_steep_entry_stress_C3.py
- Checker: experiments/root_st073/lei_ren_part1_paper_compliant_steep_entry_stress_C3_check.py
- Class: CompliantSteepEntryStressC3; method steep_in(Z,t).
- Ordered stage: python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage steepentrystress
- The corresponding .json files retain source hashes and exact scope flags.

## Original chart and common future

Let a=delta/2, r=1-mu, k=1-a and p=1+delta. The actual physical-assembly
radius branch gives q=log(R/Rtail)=-wait-Ts-2+t. Ordinary q derivatives
therefore equal t derivatives, with no Ts phase multiplier. The right endpoint
is original steep_power phase zero, not the later steep_exit inlet.

The original sigmoid primitive is J(t)=integral_0^t sigma(v)dv, J(1)=1/2.
Entry uses f(t)=integral_t^1 sigma(v)dv=1/2-J(t). The source bridge binds
both actual accumulators, their common sigmoid callable and their actual
half-integral endpoint guards. The previously accepted exit primitive
integral_t^1(1-sigma)=J-t+1/2 is a different function.

At the original power inlet,

    KS=(1-epsilon)*exp(k*(Ts+1/2)),
    K(t)=KS*exp((mu-a)*(1-t)+r*f(t)),
    K_q/K=a-mu-r*sigma(t).

Original Ts=4*log(2/delta), epsilon=.001*delta, thetaR and thetaS are
explicitly source-bound. Exact source shear has
kappa-2=2*mu+2*(1-mu)*sigma in [2*mu,2], so it is positive.
The original strictly negative angular shear follows with K>0 and exact S>0.
This sign alone is not an admissible cone certificate.

The same full moments obey

    A_q=K-k*A, E_q=delta*E-K^2, P_q=p*P-K^2/2.

For ell=1-t and defects Ad=A-1/k, Ed=E-1/delta, Pd=P-1/(2*p),

    Ad(t)=exp(+k*ell)*AdS-integral_t^1 exp(k*(v-t))*(K(v)-1)dv,
    Ed(t)=exp(-delta*ell)*EdS+integral_t^1 exp(-delta*(v-t))*(K(v)^2-1)dv,
    Pd(t)=exp(-p*ell)*PdS+1/2*integral_t^1 exp(-p*(v-t))*(K(v)^2-1)dv.

Angular evaluation uses the equivalent native correlation A=K*X, retaining
the original positive forward X history instead of subtracting large
independent backward histories. Native X_q=1-r*(1-sigma)*X, the actual
XS endpoint and the admitted power angular history prove it is the same
complete future. Likewise native half-energy e has E=2*K^2*e; its source
ODE and common right endpoint identify the same full energy.

Quadratic defects use directed signed remaining integrals, expm1 and
correlated positive lengths. Stress derivatives retain the shared K in
K*(k*X-b*Z*X_Z-1) before enclosure. The meridional primitive/velocity
routes retain their admitted exact zeros.

## Absolute pressure and connections

Absolute pressure is -pressure_scale*exp(-p*q)*P. The exact Ev0 amplitude,
Rtail units and analytic datum are consumed from the accepted history bridge;
runtime caps do not define a replacement pressure or source field.

The bridge replays the actual native PS/PQ assignments and power pressure
rows. Backward recovery of the admitted full PQ cancels the same original
power pressure integral and gives native PS at the entry right endpoint.
Both packet pressure routes (mixed bounds and Taylor rows) are updated to
this equivalent complete-future representation, with the forward history
retained under explicit original_forward aliases.

Actual source-row AST replay proves 55 entry/power identities on arbitrary
terminal axial functions: 5 K derivatives, 15 full-moment derivatives,
20 stress mixed3 rows and 15 absolute-pressure mixed4 rows. The original
native angular/entry left velocity, angular, energy and pressure joins are
also consumed from the current C4 source receipt. An upstream angular stress
companion and its physical/cone join are not constructed by this increment.

## Accepted checks and limits

The focused producer/checker passes with 247 current checker input hashes:
80 finite signed stress rows, 60 pressure rows and 96 exact meridional zero
rows, including the entire original t/Z box. The source bridge consumes 187
current native source/field/history identities.

One independent moderate-parameter fixture uses the original sigmoid and
the complete convergent analytic future
K_future(v)=1+(KS-1+g(Z)*v^2)*exp(-rho*v). It checks 8 primitive/signed
integrals, 30 full-moment derivatives, 20 normalized quotient derivatives,
40 original unnormalized stress mixed3 derivatives and 30 pressure mixed4
derivatives. Numeric enclosure comparison tolerance is 2e-7. The recorded
maximum positive enclosure miss is about 1.65e-83; this is neither the
numerical solution error nor a full NS residual accuracy claim. The actual
Gamma history is admitted separately through source identities.

Entry physical stress/remainder, entry cone, upstream angular stress,
global admissible stress, physical energy and true temporal coefficient
recursion flags remain false. No existing cone domain is extended here.

## Next dependency

Implement an original-entry adapter to the unchanged general-K physical
map, with full tensor diagonal, exact nu/lambda units and retained
Etheta=-nu*partial_zz(utheta). Bind native packet and physical radius/field
connections, derive the variable-K divergence cancellations, and cover only
new operators with one independent native Cartesian fixture. Then prove
the entire original entry cone using its actual variable shear and source
B logs before composing it with the accepted power tail. Continue through
the upstream angular repair, power and flatten regions afterward.
