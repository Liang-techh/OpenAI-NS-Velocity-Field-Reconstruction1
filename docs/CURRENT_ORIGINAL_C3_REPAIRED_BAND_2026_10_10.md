# Actual C3 repaired-band velocity and histories (2026-10-10)

Full reconstruction **ACTIVE / INCOMPLETE**. Source [89908412](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/89908412dfde3855892b82fb660121cda327142f) installs genuine third axial derivatives of the original reserved-band velocity profiles, signed correction densities, partial and complete five histories, independent P0, first-y histories, and relative terminal identities. Directed whole-band bounds cover all four original axial cells. Every accepted C2 function handle and the same original repair integer are retained. The next task is actual radial/inertial recovery through Z2 from these repaired inputs. A complete Cartesian field, admissible stress and temporal coefficient recursion remain open.

## Installed functions and source contract

The band uses x in [1,2], y=log(x), original power coordinate t=2+log(x), and the same original radius offset Rc_offset+log(x). The strict native O3_power t<=2 guard is unchanged; the existing separate reserved-band source supplies the extension. The actual source family, live owner, context, scale basis, ledger, original weights and Banach limit are shared with accepted C3 controls and the lower C2 band.

Normalized E/V have radial y orders 0,1,2 and ordinary axial Z orders 0,1,2,3. Bumps are Z-independent. The original beta, beta', beta'' formulas supply the y rows; third Z differentiates the actual amplitude and actual controls, retaining cubic product factors 1,3,3,1. The source endpoint and whole-band rows are projected as [coefficient0, coefficient1, 2*coefficient2, 6*coefficient3], with each factorial applied once. Endpoint common histories are already normalized; whole-band raw V/m/k receive the original S^-1 factor once. E/h/e/p and the independent P0 keep their original units.

With E0 the leading power and dE/dV the repaired increments, the signed density differences are

```
m=dV; h=dE; k=(E0+dE)*dV
e=dV^2-E0*dE-dE^2/2
p=E0*dE+dE^2/2
```

These are actual signed function expressions. Magnitude caps do not replace their signs. Each partial history retains the nonzero incoming C3 memory and its own rate lambda in (1,3/2,3/2,1,0):

```
D_j(x)=x^-lambda_j * (incoming_j + integral_1^x t^(lambda_j-1)*density_j(t) dt)
```

The physical Jacobian appears once and the x endpoint is independent of Z. The p rate is zero, with kernel 1/t and no artificial memory decay. Leading histories follow the unchanged exact power semigroup from the actual original Rc endpoint. All amplitude-squared third cross terms are retained. Complete histories are leading plus partial; the first-y rows follow the original full signed FTC equations. P0 stays separate from the p history.

All five relative endpoint identities now hold through ordinary third Z derivatives on the same repaired limit:

```
m(2)=Am*R_M/(2N)
h(2)=2^-3/2*Am*R_I/N
k(2)=2^-3/2*Am^2*(R_M+mu*R_D)/N
e(2)=Am^2*R_S/(2N)
p(2)=Am^2*R_Cp/N
```

The new certificates bind the already accepted C2 relative zeros and actual third control zeros. They close the correction histories relative to the leading background at x=2. They do not assert absolute exterior five-moment or heat/Rh closure.

## Directed bounds and acceptance evidence

Each magnitude jet has four rows. Profile bounds use the genuine H3 control/iterate cap, H2, the original C1 ball radius and the same selected 1/N. The original compact bump derivative caps are e^-1, 8e^-2, 2048e^-4. Product bounds retain every cubic binomial term. The actual positive terminal amplitude source gives the power profile caps. Incoming histories use the accepted actual N-scaled C3 transport bounds divided by the same N.

Own-rate whole-band masses are (1-2^-lambda)/lambda, with log(2) for p. Independent full reserved-band leading history rows are added to the partial bounds; complete first-y bounds use the signed FTC with absolute values. Function handles, directed source witnesses and logarithmic caps remain separate.

The dedicated checker passed:

- 168 mixed profile rows from both band-x and integration-t functions (2 locations * 3 y orders * 7 quantities * 4 Z orders).
- 20 signed density rows, 20 original power semigroup rows and 20 relative endpoint RHS rows.
- 5 third inlet/FTC identities and 5 complete first-y third identities.
- 294 unchanged accepted C2 function handle aliases, including independent P0.
- 4 actual axial cells, 80 actual endpoint history rows and 16 independent pressure rows.
- Exact lower whole-band source prefixes, quiet leading V=0, same-source ledger/basis, third factorial and unit guards, and replay of the four directed range records.

The default accepted constructor succeeded using the reused live owner. The exact Git-index audit passed all 1198 receipt dependencies. The sole reused static reviewer was GPT-5.6 Luna / max; no new agent was spawned. Producer gates remain false; the separate checked receipt admits only the C3 repaired-band functions, relative identities and magnitude range scopes.

```python
from lei_ren_part1_paper_compliant_current_original_C3_repair_band import CurrentC3RepairBand

band = CurrentC3RepairBand()  # reuse a live accepted control/owner when available
profiles = band.velocity_functions()  # y0/y1/y2, each value/Z/ZZ/ZZZ
histories = band.history_functions()   # complete m/h/k/e/p, value/Z/ZZ/ZZZ
P0 = band.functions['original_leading_power']['independent_P0']
```

Installed quartet: `lei_ren_part1_paper_compliant_current_original_C3_repair_band.py`, `.json.gz`, `_check.py`, `_check.json`. The symbolic/integral provider is accepted; a numerical point-field oracle is not installed by this quartet.

## Completed bounded tasks

- [x] **CURRENT-C3-SOURCE-AND-PHASE** Same-source C3 roots and same phase inverse on all original charts.
- [x] **CURRENT-C3-TARGET-FUNCTIONS** Signed third densities, own-rate transported histories and targets.
- [x] **CURRENT-C3-TARGET-RANGES** Directed source/primitive/density/target bounds over all original cells.
- [x] **CURRENT-C3-IMPLICIT-CONTROLS** Genuine third controls and quantitative third Picard tail of the same limit.
- [x] **CURRENT-C3-REPAIR-BAND-PROFILES** Actual C3 normalized E/V with radial y0/y1/y2 on the original reserved band; every lower C2 handle retained.
- [x] **CURRENT-C3-REPAIRED-HISTORIES** Signed C3 densities, incoming memories, partial/leading/complete histories, first-y rows, independent P0, relative five terminal identities and directed whole-band bounds.

## Next tasks, in dependency order

- [ ] **CURRENT-PATCHED-RECOVERY-SOURCE** Own a new companion module/checker for actual repaired radial and inertial functions through ordinary Z2. Consume this accepted band's complete E/V, E_y/V_y, m/h/k/e/p and independent P0. Use the unchanged equations in `current_original_Rm_generic_inputs.py::recover_inputs`; dz consumes one extra axial row, so C3 histories supply Q_ZZ and inertial p1/p2_ZZ. Do not substitute leading-only high rows. Bind physical R to the actual original Rc_offset+log(x), same Pstar/S/physical scales and source family. Recover genuine signed expressions, not derivatives of magnitude bounds.
- [ ] **CURRENT-PATCHED-RECOVERY-POSITIVE-DENOMINATORS** Before admitting quotients, prove L=1-delta*Z^2>=1-delta>0 on the same actual source and Z[-1,1]. Prove whole-band E>0 and C=E-2E_y>0 from the positive leading amplitude and actual selected-N perturbation bounds. C0=(2+2mu)*E0, with corrections dE and dE-2*dE_y. Derive strict quantitative margins; leading source positivity or interval overlap is insufficient. Keep generic stress cone as a separate unpassed gate.
- [ ] **CURRENT-PATCHED-RECOVERY-C2-RANGES** Enclose Q, inertial linear/quadratic functions, and quotients a/b/t0/p1/p2 through Z2 using same-source denominator lower bounds and all actual repaired third history rows. Include complete absolute pressure P0+p, physical radius factors, product/quotient cross terms, and source-bound acceptance evidence. Do not promote this band-local result to global stress admission.
- [ ] **CURRENT-PATCHED-RADIAL-Y4-CONTRACT** Extend the corrected velocity to y0..4 and necessary mixed Z rows before using `current_generic_shear_moment_recovery.py::field_rows`, which requires five radial rows. Bind all higher beta derivatives to the original checked beta source and compact endpoint flatness. The current y0..2 block is enough only for consumers explicitly limited to those rows.
- [ ] **CURRENT-PATCHED-RECOVERY-Z3-PREREQUISITE** If a downstream result needs Q/inertial Z3, first build actual fourth-Z histories/P0, source/phase/target bounds and same-limit control rows. A third source input provides second recovered Z derivatives; it does not provide third ones. Do not raise the admitted order by appending leading-only rows.
- [ ] **CURRENT-C2-RH-CORRECTION-JOIN** Continue actual repaired second histories through Rh with the same phase/radius and nonzero memory. Prove leading and correction continuation independently and preserve original interface derivatives and units.
- [ ] **CURRENT-C2-ABSOLUTE-FUTURE-INTEGRALS** Recover source-owned absolute exterior five identities through Z2 with real heat amplitude, independent P0, future-integral FTC/Gamma tails and all original cancellations. Relative repair zeros are a prerequisite, not an absolute heat closure receipt.
- [ ] **CURRENT-C2-PRESSURE-HEAT-ASSEMBLY** Assemble actual signed pressure restoration and heat exterior, source-derived axis/interface regularity, full finite-energy radial tails and exact/controlled heat evolution.
- [ ] **CURRENT-GLOBAL-MIXED-DERIVATIVES** Supply all chart/interface mixed derivatives needed by the physical velocity and stress beyond this band-local y0..2/Z0..3 block. Record every exact derivative order and source identity.
- [ ] **CURRENT-PHYSICAL-DISPATCHER-AND-ORACLE** Assemble one core-to-heat similarity/physical dispatcher returning Cartesian u/v/w and source-derived axis limits. Combine phase/Picard/integral numerical errors and select a finite common frequency meeting physical/matching/stress/oracle bounds. The present N is repair-only.
- [ ] **CURRENT-STRESS-CONE** Recover admissible divergence-form stress and flat remainder separately; prove signed margins at exit/matching/pulse/end/flatten/collar with actual physical scales and complete absolute histories. Quantify magnitude and scale dependence.
- [ ] **CURRENT-TEMPORAL-RECURSION** Implement the original n-dependent coefficient recovery, independent moment repair on the common core domain, curl-preserving truncation, finite-order remainder and final smooth sum. Spatial C3 and Picard iteration do not constitute temporal scale recursion.
- [ ] **CURRENT-OSCILLATORY-CANCELLATION** Install both original pulse families only after admitted background stress/recursion; prove averaged quadratic momentum-flux cancellation, common frequency hierarchy and flat forcing remainder.
- [ ] **CURRENT-DYNAMICS-AND-CARTESIAN-RESIDUAL** On the assembled time-dependent field, measure radial contraction, axial aspect ratio, swirl/vorticity amplification, true material winding, scale recursion and finite energy; independently evaluate Cartesian divergence and full forced NS residual.

Read this handoff first. Continue actual repaired radial/inertial recovery through Z2 next. Completed C3 source/range/control/band tasks remain complete. Full reconstruction stays ACTIVE / INCOMPLETE.
