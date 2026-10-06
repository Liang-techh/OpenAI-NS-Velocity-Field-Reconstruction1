# Current whole O3 power cone

Implementation and scoped current O3 power cone receipt: commit [3e9a4064](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/3e9a406471b694431f3f8d0fabd62f07691e6ce1).

## Constructed result

**F57C-cone1b-O3-power is complete.** The actual O3 power satisfies the original strict two-vector cone continuously on phase [0,1], all Z [-1,1], including both ends and Z=0. The directional expression has a conservative upper bound below **5.076e-435**; its dimensionless directional bracket exceeds **1.999999**. Fifteen current nonzero regions now have whole-domain cone proofs. The exact zero heat edge and unbounded exterior are a separate case. Seventeen registry regions remain without current cone admission.

The original ten O3 theta/axial stress sectors, full incoming M/H/K histories, full energy, same analytic absolute pressure, radial velocity, ordinary log-radius derivatives, physical tensor/divergence/diagonal and remainder records are retained. The axial velocity is exactly zero in this chart. Five axial stress sectors vanish by the original operator and canonical source identities; the full energy/pressure sector is retained. The variable O3 slope-to-mu transition is **still open**.

The source inventory is unchanged: 33 regions, 32 adjacent and 14 internal completed tensor traces; primitive atlas 14 / 8. This milestone does not implement higher-order coefficient recursion, actual oscillatory waves, a completed global admissible tensor/lift, temporal flatness, corrected Cartesian NS/energy, or resolved physical point u/v/w.

## Exact actual history correlation

Write t=Tw*phase, r=1-mu, C=1/(1+Z^2), L=1-delta*Z^2, k=1-delta/2. Raw canonical coefficients are U0, M0, H0, K0 at actual power phase0. The common theta stress factor is nu*lambda^(-2-delta)*sqrt(R/2)*B, with B=Pstar*U(t). The energy/pressure axial sector has B^2. Raw M/K normalization is m=raw_m/(Pstar*u) and n=raw_k/(Pstar*u^2).

The current actual production proof independently composes the original slope, axial, slope_mu and power programs. Each composed stage is verified to have canonical whole-Z shapes: u=U*C, m=M*Z, h=H*C, k=K*Z*C, e=EZ*Z^2+EQ*C^2, p=Pin*C^2. Scalar endpoint enclosures bound these identified coefficient functions; Z=0 is not used to fit a profile.

At the O2 slope endpoint let Xs=Hs/Us and tau=logPstar-1=exp(Md)+10. The same O2 turnoff mass cancels exactly in K/U-M=4*(X-1). The original transition kernel has

`KT(1)=exp(1-mu/2)-1+mu*Q`,

`Q=integral_0^1 sigma(v)*exp(v-mu*J(v))dv >= 0`, `J(1)=1/2`.

The defining sigma is the same positive-exponential flat cutoff, with reflection sigma(v)+sigma(1-v)=1. The actual pre-transition callable uses the same transition kernel implementation. Consequently

`D=(1-Xs)*exp(-tau-1+mu/2)>0`, `V=Q*exp(-1+mu/2)>=0`,

`X0=1-D+mu*V`, `K0/U0=M0*exp(mu/2)-4D`.

Both relations are replayed on the actual composed production functions as well as the generic original recurrence. Independent K/H enclosure boxes are retained but are not subtracted to define the tiny deficit.

## Continuous theta bound

Define A=k*C+(1-delta)*Z^2*C^2 and

`N=((2*delta*Z^2-1)*C+2*(1-Z^2)*Z^2*C^2)/L`.

The complete theta identity includes equilibrium, signed angular memory and incoming meridional transport:

`Theta=(A*X(t)-C)/L+C*M0*exp(-t)+N*(K0/U0)*exp(-r*t)-2*(1+mu)*C/R`.

After inserting the actual correlated histories, the retained nonnegative terms include mu*A*V*exp(-r*t)/L and (C+N)*M0*exp(-t). The sole signed drift is

`N*M0*exp(-r*t)*(exp(mu/2)-exp(-mu*t))`.

For x=Z^2 in [0,1], C+N is nonnegative, N is negative, and abs(N)<=2/(1-delta). The positive deficit coefficient has numerator

`8*(x-3/8)^2+15/8+delta*(1/2-13*x/2-8*x^2)`.

Thus Dshape>=Dmin=(15/8-14*delta)/4. Let cD=Dmin/2, cX=(1-delta)/4 and Eqmin=(mu-delta/2)/(2*r). The equilibrium is at least cX*Z^2+Eqmin; the finite-time non-equilibrium shape is at least cX*Z^2-delta/(2*(1-delta)). The actual deficit reserve satisfies

`cD*D > delta/(2*(1-delta))+Eqmin`.

The exact integral for exp(mu/2)-exp(-mu*t) is bounded by mu*(t+1/2)*exp(mu/2). Using t*exp(-r*t)<=1/r gives the uniform drift bound

`drift <= 2*M0*mu*exp(mu/2)*(1/r+1/2)/(1-delta) < Eqmin/4`.

The inverse-R shear is below Eqmin*exp(-1000), by its exact radius log source. Its defining field is not replaced by that cap. Therefore, throughout the whole domain,

`Theta >= cX*Z^2+cD*D*exp(-r*t)+theta_floor > 0`,

`theta_floor=Eqmin-drift-Eqmin*exp(-1000) > Eqmin/2`.

## Full energy and same absolute pressure

The original normalized power source obeys E_t=2*mu*E-1/2 and P_t=(1+2*mu)*P+C^2/2. Both recurrences, the raw pressure normalization, and the actual ordinary pressure rows are consumed. Put p=1+2*mu. The same current signed Rv datum, transported backwards from the actual phase1 inlet, gives

`P=-C^2/(2*p)+(Pv+C^2/(2*p))*exp(-p*(13/mu+Tw-t))`.

The actual main chart baseline, Ptilde/Pmemory, Q exponent and right-hand P0 publication are source-bound. The phase1 formula is proved to be the actual pulse inlet; no added pressure tail is introduced.

The full energy coefficient is CE=Z^2*AZ+C^2*AQ, with

`AZ=(EZ0/U0^2)*exp(2*mu*t)`,

`AQ=(EQ0/U0^2)*exp(2*mu*t)-(exp(2*mu*t)-1)/(4*mu)`.

Replaying the actual pulse_coefficients with full energy and absolute pressure proves that energy plus baseline pressure is Z*F(Z,t). The bound Lz for abs(F) retains both AZ and AQ and the full baseline pressure. The other five axial sectors are proved zero by the original operator, including the linear moment identity Cm-Z*(Cm)_Z=0.

AM-GM gives Theta>=2*abs(Z)*sqrt(cX*cD*D*exp(-r*t)). Since B(t)^2/exp(-r*t)=Pstar^2*U0^2*exp(-3*mu*t),

`2*mu*w_energy^2 <= mu*Pstar^2*U0^2*Lz^2/(2*cD*cX*D)`.

The logarithm of this actual source bound is below -4.7077e17. The full remaining pressure-memory ratio has a source log bound below -5.141e408906090034569137. These are relative-error bounds, not physical velocity values. Combining both terms gives the directional expression below 5.076e-435. Positive viscosity and lambda factors cancel from the original cone ratio.

## API, artifacts and focused verification

Prefix: `experiments/root_st073/lei_ren_part1_paper_compliant_`.

- `current_O3_power_cone_operator.py`: exact current history/operator identities and directed continuous bounds.
- `current_O3_power_cone.py`: checked `CurrentO3PowerCone(entrancecone=checked_current_entrance_cone)`; `native('O3_power',Z,phase,log_tau,theta,viscosity)` returns the original complete view and regional proof.
- `current_O3_power_cone_check.py`: focused proof/source/scope check.
- `current_O3_power_cone.json`, `_check.json`: producer and accepted scoped receipt.
- `current_O3_power_cone_views.json.gz`: complete whole source view, actual slope endpoint and phase0 source; deterministic gzip.
- Controller stage: `currento3powercone`; it also propagates inherited regional construction gates.

Verification: 50 new original history/axial identities, 12 full-theta identities, 24 canonical whole-Z identities, 3 actual composed correlation/pressure endpoint identities, 112 consumed production identities, and 25 strictly positive directed bounds. Foreign source, partial domain, missing memory/history/energy-pressure, altered mode and changed physical row are rejected. Checked endpoints/interior API, focused controller, inherited entrance and zero exterior scope, compilation and publication hashes are verified before GitHub publication. Working/index source audit: 813 dependency files matched. Read-only mathematical reviewer: GPT-5.6 Luna / max.

## Detailed next-agent tasks

- [x] **F57C-cone1b-O3-power:** original complete theta correlation, full energy and same absolute pressure, continuous cone, both source joins, checked scoped API and controller.
- [ ] **F57C-cone1b-O3-transition-1:** start at `registry.owners['o3']`. Preserve variable logU derivatives from actual slope_mu, not the fixed power slope. Consume the checked O2/transition and transition/power source jets before bounding tensors.
- [ ] **O3-transition-2:** derive actual offset-dependent O2 deficit and nonnegative transition kernel identities for every offset [0,1]. Preserve the common B_mass in M/K; do not box these histories independently before cancellation. Use the same positive cutoff and exact reflection kernel.
- [ ] **O3-transition-3:** reconstruct the original full theta/axial stress on variable amplitude, with all five histories, analytic P0, pressure and radial velocity. Prove uniform sign/ratio bounds with complete energy; endpoint power admission cannot be extended across transition by continuity alone.
- [ ] **O3-transition-4:** add the smallest independent operator/producer/checker/API stage. Preserve OPEN gates. Publish the complete source views and mark only the actual transition row complete; if admitted, count becomes 16 nonzero plus the separate zero exterior, 16 regions remain.
- [ ] **F57C-cone1b-O2-1:** establish Rh_reference and O2_slope original shear-loop bounds. Keep full axial velocity, derivative shear, M/K/E/P and inverse-R terms. Identify which signed cancellations require finite uniform N or refined source intervals.
- [ ] **O2-2:** handle O2_axial and O2_buffer as distinct charts, with the actual turnoff profile, original ordinary-y derivatives, positive omitted kernel tails and all histories after Uz=0. Use common kernel correlations rather than sample fits.
- [ ] **F57C-cone1b-inner-1:** work inward through patch, restoration, reshape, switch, bridge and core. Use each registry owner and actual full stress operator. Keep independent moment repair and all interface functions. Do not promote old-family receipts.
- [ ] **F57C-cone1b-global:** after every required local chart and axis is admitted, construct the smooth completed admissible tensor/lift with original overlap/flat edge weights, finite uniform bounds and full interface remainder control. A registry coverage inventory alone does not complete this task.
- [ ] **F57C-cone1c-pulses:** construct two actual homogeneous oscillatory velocity families and the mean correction. Verify positive amplitude, original frequency/phase/flat edge weights, divergence structure, signed lift, covariance and finite errors, then averaged quadratic stress cancellation. Covariance algebra alone is insufficient.
- [ ] **Recursion n=1:** implement the distinct first coefficient recovery equations and independent five-moment repair on the common inner interval. Cancel the known nonflat leading origin remainder. Report source functions, finite-order remainder and compatibility.
- [ ] **Recursion n>=2:** implement the genuine n-dependent operators, moment conditions, shared support and smooth summation estimates. Scaling an existing profile is not this recursion.
- [ ] **Physical completion:** provide resolved multi-time u/v/w with the actual corrections, then independent Cartesian NS residual and finite energy. Measure radial/axial scale exponents, aspect ratio, vorticity and material winding across scale levels. Keep the full long-term goal active until these are achieved.

Reuse the warm checked source graph when available. Do not cold-start the inherited producer chain for a routine focused check. Preserve unrelated dirty files. Commit and push owned artifacts, update the four current handoff files, and change a task checkbox only when its defining code and evidence are complete.
