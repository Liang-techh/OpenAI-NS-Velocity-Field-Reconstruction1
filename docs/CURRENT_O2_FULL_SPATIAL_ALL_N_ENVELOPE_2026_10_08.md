# Actual original O2 full spatial all-N envelope

Checked source [8d7ca2a0](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/8d7ca2a03ba7546b72871131d3a80f041160b787). Predecessor: [CURRENT_O2_COLLECTED_SLOW_JETS_2026_10_08.md](CURRENT_O2_COLLECTED_SLOW_JETS_2026_10_08.md), source42738b08/docs3ce577a1. The preceding goal turn was progress: actual slow source jets and five first-order native norms were accepted and pushed. The full paper-faithful reconstruction goal remains active.

## Concrete advancement

The actual original O2 five **complete changed spatial contributions** now have a conservative **all integer N>=160 N^-2 envelope**, over original y[0,1], Z[-1,1], all64 source cells and192 regular/positive/negative conditional predicates. Both C0 and ordinary-Z contributions are included. This extends beyond the accepted phase-mean component: the first-order oscillatory spatial contribution, full nonlinear bias/remainder, slow transport and every local/global endpoint are retained.

Four files under experiments/root_st073 use stem **lei_ren_part1_paper_compliant_current_original_O2_all_N_spatial_envelope**: producer, deterministic compressed report, focused checker and receipt. APIs: **OriginalO2AllNSpatialEnvelope.integrate / at_candidate / primitive_bounds / nonlinear_remainder**. The existing GPT-5.6 Luna/max worker reviewed read-only; root implemented, computed and accepted. No new child.

The new envelope is **not sharp**: original large factors stay formal in the output. It does not imply small physical gradients, a low momentum residual, solved terminal identities or selected global N. The tighter direct finite-N contribution bounds and mean-component bounds remain accepted. No improvement ratio is claimed between their different units and this new all-N spatial estimate.

## Same-source periodic primitive

The actual signed first-order functions are f1=(B,E*A,E*(V*A+B),2*V*B-E^2*A,E^2*A). Original E,V depend only on slow variables. Actual A(1-phi)=-A(phi), B(1-phi)=-B(phi) therefore imply zero true-phase mean for f1 and all fixed-phase ordinary slow C0/y/Z/yZ derivatives.

For the same source function, G_j(y,Z,phi)=integral_0^phi f1_j(y,Z,s) ds, G_j(0)=G_j(1)=0. Zero mean gives

~~~
|G_alpha(phi)| <= min(1/2, upper(phi), 1-lower(phi))
                     *sup_s |D_slow^alpha f1(s)|,
alpha=0,y,Z,yZ.
~~~

The API returns source-backed absolute upper envelopes and exact zero endpoint traces. It does **not** select or evaluate a signed G point value. Original signed density rows and accepted native source/derivative inputs remain attached; bounds are never differentiated as functions. Same-source physical |A|<=1/2 also intersects first-order C0 magnitudes without changing the signed jets or derivative rows.

## Exact full nonlinear remainder and its Z derivative

For x=A/N, define RE=N^2*E*(exp(x)-1-x), FN=N*E*(exp(x)-1), rp=N^2*E^2*(exp(2*x)-1-2*x)/2. The literal full changed densities satisfy delta=f1/N+r/N^2 with

~~~
r_m=0; r_h=RE; r_k=V*RE+FN*B; r_e=B^2-rp; r_p=rp.
RE_Z=E_Z*A^2*R2(x)+E*A*A_Z*exprel(x)
FN_Z=E_Z*A*exprel(x)+E*A_Z*exp(x)
rp_Z=4*E*E_Z*A^2*R2(2*x)+2*E^2*A*A_Z*exprel(2*x)
R2(x)=integral_0^1(1-s)*exp(s*x)ds.
~~~

The original half-phase monotone-inverse |A|<=1/2 bound extends to **all phase** using actual reflection, with slow parameters fixed. Thus rho=exp(1/320), rho2=exp(1/160) are valid for every N>=160 on the unpaired full-period remainder. Directed native and physical amplitude envelope families retain all E_Z,V_Z,A_Z,B_Z terms and intersect only bounds for the same actual functions. No source q floor or flat-limit replacement occurs.

The remainder is not zero mean. Its nonzero nonlinear mean bias is already included; the accepted separate mean contribution must not be added a second time. The positive envelope formulas coincide with the earlier H/K/P table, but have a separate exact **unpaired remainder** proof and480 independent exponential C0/Z comparisons.

## Original radius phase and all endpoints

Original logR=logRref+y, logRref=log110+10(logCstar+logPstar), r_minus=4*epsilon_core*exp(hb*s_c/2). Therefore the same literal common phase is

~~~
frac(N*(log(110/4)+14*logPstar+10*logCstar+1000+y-hb*s_c/2)).
~~~

Original phi0_Z=0. The strictly positive microscopic hb*s_c/2 origin is retained as a directed error budget, not zeroed or replaced by an approximate phase point. Each candidate query obtains all65 original radius/phase endpoints afresh for that N; no diagnostic N7 phase archive is substituted. Candidate endpoint queries inherit the accepted phase oracle's4096-bit and microscopic-origin budget checks; the analytic envelope itself is an all-N theorem. No global N admission is claimed.

For each original cell[a,b] and final-endpoint own-rate kernel K(y)=exp(-rate*(1-y)):

~~~
integral_a^b K*f1(y,Z,frac(N*y+phi0))/N dy
 = ([K*G]_a^b-integral_a^b K*(G_y_slow+rate*G)dy)/N^2.
~~~

Ordinary Z uses G_Z,G_yZ because phase0_Z=0 and the fixed y endpoints/kernel are Z independent. Both endpoint magnitudes K(a)*sup|G| and K(b)*sup|G| remain. No local/global/source-join cancellation is assumed. Each exact positive Duhamel mass multiplies transport and remainder once. Own rates are m/e:1, h/k:3/2, p:0. Branch alternatives are hulled after local normalization, not added over overlaps; cell contributions are then summed. No extra R Jacobian, angle Jacobian, period multiplier or reset incoming pressure history appears.

## Uniform full-spatial coefficient bounds

Let Lambda0=Pstar^11*Cstar^10. For all N>=160,

~~~
|full changed C0 contribution_j| <= (Lambda0/L(Z))*column2/N^2
|ordinary-Z full changed contribution_j| <= (Lambda0^2/L(Z)^3)*column3/N^2.
~~~

These are dominating formal output units, not raw derivative identities. Native G/remainder/endpoint products are formed in the source q/u atlas before division; the exact Z-only unit is restored in the final native enclosure. An ordinary-Z enclosure is not the derivative of a C0 cap.

| Moment | C0 coefficient divided by Lambda0/L | Z coefficient divided by Lambda0^2/L^3 |
|---|---:|---:|
| m | 303630.5197 | 8.222099752e+13 |
| h | 58636.09544 | 1.387936067e+13 |
| k | 239998.1871 | 6.309979843e+13 |
| e | 79426.36916 | 1.94038493e+13 |
| p | 136017.5769 | 3.544548737e+13 |

The large units matter: this table does not assert a small N160 physical contribution. It establishes useful frequency dependence with the actual source and complete spatial structure. Raising N by a factor10 reduces these envelopes by100 with source data fixed. Sharpening or admitting the required common N across all routes remains future work.

## Focused acceptance

Five exact full-density first-order/remainder identities, two exact ordinary-Z remainder cancellations and two C0/Z IBP chain identities.24 compatible full unpaired exponential fixtures at N160,257,10^9, both V signs and amplitude boundaries/zero/tiny values:480 independent C0/Z comparisons. Independent finite analytic fields with genuine rapid phase and nonzero endpoints give20 direct C0/Z IBP identities and20 complete nonlinear spatial contribution enclosure checks. Finite fixtures do not replace original source functions.

One actual source reconstruction checks **3840 periodic primitive native cap rows**, **3840 full nonlinear native remainder rows**, **7680 native endpoint/transport/remainder/total coefficient components**, **320 independent positive own-rate masses**, **640 independent kernel endpoint values**, **650 branch-union/whole-window coefficient checks**. Three candidates verify195 actual original radius/phase endpoints and native C0/Z units; exact G0/G1 traces and eight owner/N/phase guards pass. **1288 staged exact dependency hashes PASS**. Accepted ancestor suites are reused, not rerun.

Producer 48.765s; checker 72.063s. The read-only worker confirmed remainder terms, full-phase amplitude extension, source binding, masses/endpoints/units and limited scope.

## Ordered next tasks

- [x] **O2-SAME-SOURCE-PERIODIC-PRIMITIVE-BOUNDS:** actual signed f1 function definition, zero true mean, all four slow primitive bounds and exact endpoint traces; source 8d7ca2a0. Signed G point evaluation remains open and is not required to use this enclosure.
- [x] **O2-FULL-SPATIAL-ALL-N-NMINUS2-ENVELOPE:** original source/candidate phase, full nonlinear remainder and ordinary-Z terms, all own-rate/local/global endpoints, native dominating units; source 8d7ca2a0. Conservative complete O2 contribution, not sharp or all-route admission.
- [ ] **NEXT ACTUAL SOURCE-OWNED INCOMING/CONTINUOUS DRIVER:** reuse the existing genuine candidate-N1024 upstream correction bridge in current_original_O2_source_incoming_common_N.py and current_native_middle_O2_inlet_C1_histories.py, then bind its compatible five C0/Z corrections to this newly completed actual O2 full-predicate source/envelope owner. The old bridge does not by itself prove this new full driver or all-route admission. Preserve original analytic P0 and rate0 pressure memory; translate native bases only through checked same-source factor identities. Acceptance: owner/family/partition/unit-consistent actual upstream functions and covers consumed by the new real O2 driver, exact source hashes and propagated source errors; no arbitrary independent rows, zero/affine fixture substitution or reset incoming. Do not call original_full_inlet the correction inlet.
- [ ] **SHARP O2 FREQUENCY/PHASE BUDGET:** compare the conservative derivative/IBP units with original admissible N inequalities; retain paired flat endpoint suppression and necessary source correlations. If required N exceeds the accepted4096-bit endpoint oracle, extend directed modular/origin arithmetic without losing positive hb*s_c/2, or provide an exact symbolic phase contract with adequate error bounds. Acceptance: useful complete spatial C0/Z norms and no silent oracle restriction/global-N claim.
- [ ] **ACTUAL AXIAL BUFFER / ALL17-24 ROUTES:** nonconstant original a,b,q,t0 mixed roots and slow offsets, buffer/route join derivatives and continuous inverse/primitive/density/source-error coverage. Rebind the actual driver; generic cap-only graphs and constant-q references are insufficient.
- [ ] **FIVE FUNCTIONAL TERMINAL CONTROLS / ONE GLOBAL N:** join actual route histories with all signed mean/remainder/endpoint/oracle/error terms, solve five whole-Z terminal identities and satisfy original frequency inequalities using one common N. Existing symbolic transport/controls or candidate queries do not prove this.
- [ ] **JOIN/PREHEAT/HEAT/ENERGY / STRESS-FLAT / n-RECURSION / PULSES / CORRECTED UVW:** unchanged full objective: same-data matching and analytic pressure, exact heat/finite physical energy, admissible stress/cone and flat remainder, actual n-dependent recovery/moment repairs/smooth sum, both pulse families and quadratic stress cancellation, Cartesian velocity/forced NS and contraction/elongation/winding diagnostics. Actual scale recursion remains open.

Use this successor as active queue. Mark scoped tasks DONE only with evidence and a commit. Full reconstruction remains active.
