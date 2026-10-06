# Current whole pulse entrance cone

Current successor: [CURRENT_O3_POWER_CONE_2026_10_06.md](CURRENT_O3_POWER_CONE_2026_10_06.md), commit [3e9a4064](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/3e9a406471b694431f3f8d0fabd62f07691e6ce1). O3 power is now complete; variable O3 transition, O2/inner cones and actual recursion/waves remain open. Use the latest handoff for current counts and tasks.

Implementation and scoped current entrance cone receipt: commit [5cf6b3fd](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/5cf6b3fdc6679c716003faf9c7eb7165e0203aea).

## Constructed result

**F57C-cone1b-entrance-1 and entrance-2 are complete.** The actual pulse entrance satisfies the original strict two-vector cone throughout xi[0,.02], all Z[-1,1], including xi=0. The conservative dimensionless directional bracket is **greater than 1.9586218**. Fourteen current nonzero regions now have continuous source-function cone bounds; the exact zero heat edge/exterior remains separate. Eighteen upstream registry regions have no current whole-domain cone admission.

The source still contains all 15 signed theta/axial sectors, current selected C5 amplitude and controls, incoming moments and energy, current absolute pressure, full radial velocity, local/incoming cross and square terms, completed diagonal/divergence/Cartesian tensor and remainder records. The local axial velocity vanishes at xi=0, but incoming radial velocity, moments, energy and pressure are not reset. Forward energy is the same function as the backward energy/end-loss representation, not an extra contribution.

This advances the admissible-stress construction. It does not implement n-dependent coefficient recovery, actual oscillatory velocity corrections, global completed-tensor admissibility, temporal flatness, corrected NS/energy or resolved point u/v/w. The certified nonflat leading origin remainder is still present.

## Source-functional reduction

Use r=1-mu, lambda1=.5-mu, qmin=mu-delta/2, k=(1-delta)/2 and the current selected amplitude ap(Z). The original primitive is gp(xi)=integral_0^xi sigma(50a)da. Positivity, reflection and the exact flat endpoint imply

`0 <= gp <= gp(.02)=.01`, `0 <= gp_xi <= 1`, `abs(gp_xixi) <= 50*sup abs(sigma')`.

The same current main kernel, stress, logs, physical operators and velocity exporter are used on entrance, with only the reviewed coordinate guard extended. Current incoming/entrance and entrance/main tensor function joins are consumed before enclosure.

The pre-power angular history must be factored before subtracting nearly equal final enclosures:

`Xp-1/r = (h_pre/u_pre-1/r)*exp(-r*Tw)`.

Its full entrance memory is proportional to

`(h_pre/u_pre-1/r)*exp(-r*(Tw+xi/mu))`.

The source proof independently composes both actual pre and native O2/O3 recurrences, so h_pre/u_pre is an exact scalar function over all Z. Z=0 evaluates that already identified scalar function; it is not a fitted representative point. The bound uses the actual `physical.pre.power(0,0)` packet and the conservative coefficient `2*(abs(h_pre/u_pre)+1/r)/(1-delta)`. The native stress coefficient remains signed and unchanged.

Current pre and native inlet parameters are separate objects. A typed bridge binds both to `CompliantPressureDatum -> CompliantOuterParameters(Md='40', precision=160)`, their common defining pressure records, original parameter paths, and `Tw=-60*log_mu`. No cross-instance pointer equality is required. The selected future mu enclosure pointer is not used as the incoming parameter definition.

The source U used in the log amplitudes is exactly `flatten.U is pulse.high.constants['U']`. The accepted current production-function theorem binds that native constant to the actual pre-power phase1 function and its whole-Z canonical shape. The new construction also consumes the actual pre-power phase1 U packet. Neither interval overlap nor a midpoint establishes this identity.

All 12 signed error sectors have entrance-specific grouped log bounds. Angular memory has the additional -r*Tw factor. The absolute-pressure memory reaches its maximum at xi=.02; the other sectors use xi=0. Main xi=.02 estimates are not extrapolated across the entrance left edge. Individual absolute caps only bound errors relative to the signed baseline.

Two integrations by parts in the original full linear kernel retain its exact remainder. The current bounds are

`Bmax=.01*apmax`, `Dmax=apmax`, `Kmax=r*mu*k/(lambda1^2*qmin)`, `Wmax=2*Bmax+Kmax*Dmax`.

The w error keeps the radial-rate term, the two-IBP remainder and the actual ap_Z term. The b error keeps ordinary log-radius differentiation with d_y=mu*d_xi. Both errors are below 1e-6. The complete source shear has a=2+2mu and the full actual b source, including the incoming velocity rows. With the error corrections retained, the bounds give

- theta normalized lower >0;
- a-bw lower >1.9793365;
- 2 minus the full directional expression >1.9586218;
- normalized directional margin >7.8344874.

The positive common physical lambda/constant-viscosity factor cancels algebraically. Direct signed interval diagnostics remain available even if they lose source correlation; they are not the continuous-domain proof.

## Entry points and evidence

Prefix: `experiments/root_st073/lei_ren_part1_paper_compliant_`.

- `current_pulse_entrance_cone_operator.py`: original memory/endpoint-log identities, current typed parameter/U/source binding, full signed view validation and whole entrance bound.
- `current_pulse_entrance_cone.py`: `CurrentPulseEntranceCone(tailcone=checked_current_angular_tail_cone)`. It obtains the whole current entrance and both actual pre-power source packets from the same checked graph.
- `current_pulse_entrance_cone.json` / `_check.json`: 14 new memory/endpoint identities, the consumed 112 current incoming production-function identities, 43 strictly positive directed whole-domain bounds and scoped gates.
- `current_pulse_entrance_cone_views.json.gz`: deterministic gzip retaining the complete unpruned whole entrance and actual pre-Tw source packet. Decode with `json.loads(gzip.decompress(path.read_bytes()))`.
- Focused controller stage: `currententrancecone`; it is also the last producer/checker pair in `all`. Earlier checked scope is inherited separately.

Producer/checker, focused controller, fresh checked entrance API at xi=0/interior/.02, inherited main scope and exact-zero exterior passed. The focused checker rejects foreign source, a partial domain, omitted memory/pressure/forward energy, incomplete incoming histories, changed pressure mode, mismatched physical signed rows and oversized amplitude. All 806 dependency hashes matched the working tree and Git index. Read-only review used **GPT-5.6 Luna / max**; no new Astra child or cold inherited producer was started.

Inventory remains 33 tensor regions / 32 adjacent / 14 internal tensor traces, primitive atlas 14 / 8. The registry alias for entrance/incoming is **`registry.owners['incoming']`**; there is no `owners['entrance']` alias.

## Next executable tasks

### New O3 construction available before the next cone bound

`current_O3_theta_correlation.py` and its `.json` now provide an exact source-functional reduction of the COMPLETE order-zero O3 power theta stress. The actual `pulse_coefficients` AST and current `actual_power` inputs/normalizations are replayed; 12 symbolic identities passed. This adds a construction, without admitting the O3 cone.

Let U0, M0, K0, Xpre be the actual scalar/raw canonical coefficients at O3 power phase0, t=Tw*phase, C=1/(1+Z^2), L=1-delta*Z^2, r=1-mu, b=(1-delta)/2 and k=1-delta/2. Define

`A=k*C+2*b*Z^2*C^2`,

`N=((2*delta*Z^2-1)*C+2*(1-Z^2)*Z^2*C^2)/L`.

After dividing by the common positive `nu*lambda^(-2-delta)*sqrt(R/2)*B`, the full inertial theta numerator is

`Theta=(A*X(t)-C)/L + C*M0*exp(-t) + N*(K0/U0)*exp(-r*t)`.

The different radial-shear mode remains `-2*(1+mu)*C/R`. Incoming meridional transport has B power 2, so cancellation of the common B leaves ONE B factor. The actual normalized moment inputs are `raw_m/(Pstar*u)` and `raw_k/(Pstar*u^2)`; the latter is not `raw_k/(Pstar*u)^2`. Raw M/K constants and normalized moment hats must not be interchanged.

With `Eq=(A/r-C)/L`, `Theta0=(A*Xpre-C)/L+C*M0+N*K0/U0`, and `D=A*(Xpre-1/r)/L+N*K0/U0`, the same exact source function is

`Theta(t)=Eq+(Theta0-Eq)*exp(-t)+D*exp(-t)*expm1(mu*t)`.

This keeps the M/K-memory correlation and the small expm1 drift before interval widening. At Z=0 it reduces to `k*X(t)-1+M0*exp(-t)-(K0/U0)*exp(-r*t)`. The zero-mu/delta O2 axial identity `K/U=M+4*(X-1)` gives `3*(1-X)`, illustrating why the memory sign alone is insufficient; that O2 relation is not asserted for the later O3 phase without its transition correction. Xpre is the initial source, not the final pulse-inlet Xp.

Next work is a uniform positive bound for Theta0 and the correlated drift, retaining the separate inverse-R shear and full axial energy/pressure. The new record explicitly keeps both O3 cone gates false. It proves order-zero theta algebra only; its unused zero axial slots are not current pressure, energy or remainder values.

Use branch `codex/st073-transition-next`. Read this latest checkpoint first and retain unrelated dirty experimental files. Reuse the warm checked source graph when available. Check only affected dependencies after a change. A regional cone result must not promote global or temporal completion gates.

- [x] **F57C-cone1b-entrance-1:** recover actual current incoming/local source, full shear, energy and pressure, preserving both functional tensor joins.
- [x] **F57C-cone1b-entrance-2:** whole xi[0,.02]/Z[-1,1] strict cone with native pre-Tw memory factorization and entrance-specific correlated log bounds.
- [ ] **F57C-cone1b-O3-power-1:** start from `registry.owners['incoming'].chart('O3_power',...)` and the actual `physical.pre.power`. Derive the theta inertial numerator as an exact source function. Factor Xpre-1 and Xpre-1/(1-mu) BEFORE enclosing; a rounded Xpre~1 box cannot resolve a mu-order gap. Keep the current full five histories and absolute datum. Output the actual source identity, its complete domain and a uniform positive margin or an analytic obstruction.
- [ ] **F57C-cone1b-O3-power-2:** use the native variable radius and positive amplitude logs to bound full energy/pressure and incoming axial stress relative to theta. Bind actual full shear and ordinary log-radius derivatives. Close phase[0,1], all Z[-1,1], including the O3-power/entrance source join; do not use a finite list of phase points as the proof.
- [ ] **F57C-cone1b-O3-transition:** `registry.owners['o3']` provides the actual slope-to-mu tensor. Preserve the original transition integral and its complement, full m/h/k/e/p, and the O2/O3 plus transition/power joins. Derive the variable-slope stress/shear correlation before applying any broad caps. Publish a separate full-domain result.
- [ ] **F57C-cone1b-O2-reference:** `registry.owners['o2']` and the current pre dispatcher provide Rh reference and O2 slope/axial/buffer. Preserve nonzero axial velocity, current pressure, incoming energy, actual original kernels and all source-coordinate joins. Handle each domain separately; diagnose an unresolved sign by exact source numerator and required relation, not an interval crossing zero alone.
- [ ] **F57C-cone1b-patch/restore:** recover the original shear loop and current signed tensor correlation for actual patch and restore. Keep K_Z, complete mixed histories and pressure. Reserve independent moment-repair intervals. Any defining-source modification requires updating downstream closure, pressure, tensors and joins affected by it.
- [ ] **F57C-cone1b-reshape/switch:** treat reshape, switch-power and both microswitches with their actual source slopes and periodic shear construction. Establish a finite uniform N and independent five-moment repair where the original method requires it. Do not substitute an unrelated fitted velocity to force a sign.
- [ ] **F57C-cone1b-bridge/core:** cover first/second/macro bridges and the positive-radius core with the same current family and exact source coordinate maps. Keep the analytic core-axis extension distinct from positive-radius cone formulas. Preserve all current tensor interfaces and completed diagonal/divergence records.
- [ ] **F57C-cone1b-global:** after every remaining current region and relevant zero-edge limit is handled, assemble a smooth globally admissible completed stress/lift and derivative bounds across joins. Two-vector meridional cone bounds alone do not certify the completed diagonal tensor.
- [ ] **F57C-cone1c-pulses:** construct both actual homogeneous oscillatory pulse velocity/potential families, supports and divergence identities. Compute their actual covariance integrals, frequency parameters and finite errors. A positive reference matrix alone is not a pulse family.
- [ ] **F57C-cone1c-amplitudes:** bind squared amplitudes to current background stress, prove uniform positivity, and construct flat edge weights and smooth square roots with derivative bounds. Keep the exact heat zero edge distinct from strictly positive regions.
- [ ] **F57C-cone1c-cancellation:** implement the signed linear lift and mean correction, then actual averaged quadratic momentum-flux cancellation. Record all finite-frequency and derivative errors against the uncorrected current stress in common physical coordinates.
- [ ] **F57C-recursion-n1:** implement the actual first coefficient equations on common inner domains, with independent moment repair. Cancel or absorb the certified nonflat leading origin Ez before declaring a flat remainder.
- [ ] **F57C-recursion-nge2:** implement the distinct n-dependent recovery equations, finite-order estimates and smooth sum. Truncate streamfunction/potential before curl. Scaling n=0 in coordinates does not implement higher coefficients.
- [ ] **F57C-final-fields:** resolve u/v/w/p at multiple times with tolerances; independently evaluate corrected Cartesian divergence/NS, prescribed-domain energy/radial tail, scale recursion and true accumulated material winding. Distinguish directed source enclosures from resolved field values.

Mark tasks complete only when their source construction and focused affected checks are published. The long-term goal remains active.
