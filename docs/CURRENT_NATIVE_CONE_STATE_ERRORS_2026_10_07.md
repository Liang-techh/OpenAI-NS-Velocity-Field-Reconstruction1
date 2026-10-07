# Original continuous modulation state errors and active-loop cone budget

Checked source: [1eb86191](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/1eb861918c7910b629e7113c274838c81999cd25). Predecessor: [original all-N coefficient functions](CURRENT_NATIVE_ALL_N_FUNCTION_CONTROLS_2026_10_07.md). Producer 7.297s, focused checker 0.625s. Working/index audit: 1,139 matching dependency hashes. Read-only mathematical/source reviewer metadata: **gpt-5.6-luna / max**. The existing original all-N route was reused; no source ancestor constructor or route rebuild was needed.

The reconstruction now has quantitative errors for the actual modulation's normalized cone state `(a,b,p1,p2)` throughout all 24 original continuous radial cells and full Z[-1,1], with 35 conditional source branches. These errors retain the full original inertial operator, nonzero axial background, pressure memory and meridional cross terms. They supply a sufficient frequency for stability of the frozen loop **where q is nonzero**. The new active-state threshold is strictly larger than the preceding combined source/repair threshold: the missing stress comparison was a real additional sufficient requirement.

The q-flat input margin, post-modulation quiet-region margin, corrected repair-band error/margin and global joins remain open. No actual compatible whole N, solved controls, instantiated tail, terminal closure, global cone or coefficient recursion is admitted.

## Relevant paper estimate and derivative scope

The supplied Lei–Ren v2 Section 11.6, equations (11.29)–(11.32), compares the C0 state `(a,b,p1,p2)`, where p is the full inertial stress divided by the positive swirl factor F. The required inertial comparison uses profile and cumulative moment values and their first Z derivatives. The shear comparison uses the original phase-held A_y/B_y. OpenAI Appendix C, equations (C.13)–(C.16), states the same distinction.

Higher total-y derivatives may grow with N. They need finite, source-faithful bounds for later smooth reconstruction, coefficient recursion and waves; they do not all need to become small before this C0 cone comparison. The existing fast-y N^0 terms and higher growing powers remain valid sidecars. This stage does not admit new mixed4 velocity/mixed3 stress rows or replace their later tasks with a smallness claim.

## Continuous history and complete inertial differences

For a cell of true log-radius width w and each original rate r=0,1,3/2, the coefficient history at any interior point satisfies

```text
D_p(s)=exp(-r*s)*D_p(in)+integral_0^s exp(-r*(s-t))*density_p(t)dt, 0<=s<=w
|D_p(s)|<=|D_p(in)|+density_p_absolute_cap*mass_r(w).
```

The existing all-N whole-cell contribution is the same uniform source cover times its own true mass. Its magnitude therefore bounds every partial source integral. This bound includes the incoming memory, uses each distinct rate, and applies to C0 and ordinary Z coefficients because widths/endpoints are Z-independent. It is an absolute bound; no joint signed cancellation is claimed. Pressure rate zero retains memory. It covers the entire cell, rather than using a right-endpoint enclosure as an interior value.

Let D=(Dm,Dh,Dk,De,Dp) be the actual five history differences and DZ their ordinary first-Z rows; write dE,dV for the profile differences. In common S=Pstar units put

```text
L=1-delta*Z^2, d=1-Z^2
T=(1-delta)*Z*m+d*m_Z
DT=(1-delta)*Z*Dm+d*Dm_Z

delta_pressure=Dp
delta_radial=(2*Z*dV-DT)/L
delta_Itheta_linear=(-dE+(1-delta/2)*Dh-(1-delta)*Z*Dh_Z/2)/L
delta_Itheta_quadratic=((2delta-1)*Z*Dk-d*Dk_Z+E*DT+dE*T+dE*DT)/L
delta_Iz_linear=(-dV+(1-delta)*(Dm-Z*Dm_Z)/2)/L
delta_Iz_quadratic=(V*DT+dV*T+dV*DT+2delta*Z*De-d*De_Z
                   +2(1+delta)*Z*Dp-d*Dp_Z)/L.
```

`stress_increment()` implements these exact differences. The checker derives the original recovery expressions directly from `GenericMomentRecovery.field`'s AST and subtracts changed minus original. Original P0 cancels algebraically; it is never renormalized, reset or merged into a new pressure datum. Its inherited source preservation is separate from this newly checked difference identity.

The bound path uses L^-1<=2 for actual 0<=delta<=1/2 and full |Z|<=1. It carries the N^-1/N^-2 history coefficient bounds, original T and E/V caps into executable negative-power log polynomials. Products generate further negative powers. Source caps are used only for conservative bounds, never as defining functions or point coefficients.

Normalization retains the actual original radius and raw Pstar factors once. With x=mathcal_A/N and E_N=E*exp(x),

```text
delta_p=exp(-x)*(R/E*(delta_I_linear+Pstar*delta_I_quadratic)
                 -p*(exp(x)-1))
delta_a=-2*A_y_slow/N
delta_b=(exp(-x)-1)*bL+exp(-x)*2*B_y_slow/(N*E).
```

The actual positive E theorem supplies the denominator bound. The original |mathcal_A|<=5/4 and N>=160 give |x|<=1/128. Thus both exponential factors and the relative swirl difference have directed uniform covers. The errors have only negative N powers; radial fast derivatives do not enter this comparison.

## Frozen active polynomial margin and explicit tolerance

For X=(a,b,p1,p2), use the weighted polynomials

```text
G1=a
G2=a^2+b^2-2a
G3=p1*a-p2*b-a^2-b^2
G4=2a*G3^2-G2*(p1*b+p2*a)^2.
```

The actual original loop scales give, on q!=0,

```text
aL>=a_floor=2/(1+Bstar^2), v-2>=3eta/8
D=H(t)-v>=m/2, Q=2D^2-(v-2)J^2>=m^2/4.
G1>=a_floor, G2>=a_floor*3eta/8,
G3>=a_floor*m/2, G4>=a_floor^3*m^2/4.
```

G4 is a weighted lower bound `aL^3*Q`; it is not labelled a raw Q bound. All positive scales and margins remain in finite directed logs, including the inverse aL/Bstar losses. For the active state, |aL|,|bL|<=3. Set H=2+max(3,|p1|,|p2|), g equal to a positive lower bound for all four G values, and rho=min(1,g/(512H^5)). H bounds the state and its radius-one perturbed segment.

The independent gradient expansion has total absolute coefficient sums `(1,6,8,208)`, derivative degree at most five, and H>=1. Hence the gradient L1 norm is bounded by **208H^5**. Since 208/512<1/2, state error at most rho preserves G_j>=g/2. This conservative proof verifies the stated tolerance without assuming bounded t0 or inverse-a derivatives; those rational quantities are absent from the polynomial gradient and their losses are already in g.

For every error polynomial with powers p<=-1 and N>=1, `sum C_p*N^p <= (sum C_p)/N`. The recorded log-N threshold makes all four state errors at most rho, and its maximum is combined with the already accepted 17 source plus repair/positivity requirements. This is a finite active-branch sufficient lower threshold, not a global construction N. q=0, the quiet continuation and the correction band need their own strict input margin and error proof.

## Evidence and APIs

- Seven independent exact full original recovery difference identities and four weighted frozen cone identities.
- Two exact normalized stress/shear quotient identities and independent polynomial gradient coefficient bounds.
- 36 independently quadrature-integrated signed continuous partial Duhamel references, with both incoming signs and rates 0,1,3/2.
- 324 full original recovery error enclosure comparisons and 216 normalized stress/shear comparisons at N=160,257,2048, including nonzero V/P0/transport/history-Z rows and zero/nonzero primitive values. These are manufactured operator references; they are not actual original-field point evaluations.
- 240 genuine continuous-prefix C0/Z coefficient polynomials on 24 original cells and 140 negative-power normalized state error polynomials on 35 original source branches. Original true widths, positive E proofs, pressure memory and common source/datum family are checked.

Files use prefix `experiments/root_st073/lei_ren_part1_paper_compliant_` and stem `current_native_generic_cone_state_errors`: `.py`/`.json` and `_check.py`/`_check.json`.

API: `NativeGenericConeStateErrors(existing_allN_owner).compute(existing_allN_live)`, or `run(owner,live,return_live=True)`. It reuses the accepted all-N source queries and histories. No original graph ancestor, source route or actual point field is rebuilt. The result supplies complete modulation state errors and a q-active tolerance; it leaves the missing global gates false.

## Detailed continuation

- [x] **CONE1-prefix:** original continuous-cell C0/Z history coefficient bounds with true distinct masses and incoming/quiet pressure memory.
- [x] **CONE1-full-error:** exact full inertial/radial/pressure differences, normalized p quotient, phase-held a/b errors and original E/R/Pstar provenance.
- [x] **CONE1-active-margin:** weighted frozen q-active G margins, explicit polynomial Lipschitz/tolerance proof and new source+repair+active-state lower threshold.
- [ ] **CONE2-q-flat:** extract quantitative original signed D/Q (or weighted G3/G4) margins on the actual restricted set Delta>=eta. Use original correlated source proofs; a global relaxed H0-2 margin alone does not give this D/Q lower. Keep branch crossings and complete Z coverage. Apply the same actual state-error envelopes with this restricted margin.
- [ ] **CONE2-quiet:** obtain uniform original G margins from r_plus through Rc using the actual whole power cone and complete incoming M/K/energy/absolute-pressure correlations. Keep modulation defects even where the local source is zero, and compare its p errors against that margin.
- [ ] **CONE3-repair-profiles:** propagate the same-source C1 control ball through the actual fixed g_i and g_i' to bound dE,dV and shear errors on Rc..2Rc. Use exact normalized bump shapes, A/A_Z and positive mu; do not reuse a scoped old O3 repair frequency.
- [ ] **CONE3-partial-moments:** construct corrected partial five moments on the repair band from the actual functions plus incoming defects. Retain the separate pressure datum, all quadratic terms, radial recovery and first Z derivatives. An endpoint closure does not bound the band interior.
- [ ] **CONE3-band-state:** combine repair profile/partial-moment bounds with the actual full inertial operator, original Rc..2Rc G margins and positivity. Produce a same-family finite lower condition including derivative/support losses.
- [ ] **CONTROL1b-common-N:** combine source, repair, active, q-flat, quiet and repair-band conditions only after the actual source/function joins and unchanged outer admission are connected. Record one typed common integer definition, directed inequalities and domain; original fixed higher derivatives may grow with N.
- [ ] **CONTROL1b-oracle/tail:** implement directed evaluation/averaging of the actual phase-bound functions, instantiate useful C1 controls and a certified tail, and retain exact function roots. Do not use caps as coefficient values or phase samples as integral proof.
- [ ] **CONTROL2-field/closure:** install convergent band functions, recover radial velocity/pressure, independently integrate all five corrected moments/Z rows with the tail and establish terminal support/pressure/heat function joins.
- [ ] **HIGH/OUTER:** obtain the required finite higher mixed rows and changed source seams; complete the global cone, analytic heat/pressure exterior and physical finite-energy tail.
- [ ] **REC/WAVE/PHYS:** genuine n-dependent coefficients with independent per-order repairs and flat summation, the two oscillatory stress pulse families and cancellation, then corrected Cartesian uvw/NS residual and measured core contraction/elongation/material winding. Animation remains secondary.

Gate: `current_original_whole_modulation_C0_state_errors_and_active_loop_frequency_bound_certified`. The full long-term goal remains active and incomplete.
