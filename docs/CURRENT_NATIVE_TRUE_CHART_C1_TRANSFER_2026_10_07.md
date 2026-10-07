# True original chart lengths and initial C1 collar history transport

Checked implementation: [45bcf343](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/45bcf343eb03513d137ab150db8170034580dc2d). All **17 original chart lengths** now have positive directed covers with their true Jacobians and source parameters. The same C1 integration backend executes on **five declared cells in different charts**. Actual inlet-to-collar C1 histories now extend from native phase s_c/2 to3s_c/4 on whole Z=[-1,1] and Z=[.49,.51]. The active bridge and the full inlet-to-O2/Rc route are still missing.

## Original log-radius lengths

Write P=log Pstar,C=log Cstar,B=h_bridge,S=h_switch,T=400*Abar,W=Tw,M=Md. These are the unchanged source constants; analytic interval constants stay intervals and selected binary source constants retain their original assignments. Widths are formed from the exact original radius identities before huge endpoint subtraction. The inlet is interior to bridge_first, not phase0.

| Original chart | Full native coordinate domain | True full log-radius length |
|---|---|---|
| bridge_first | [0,1] | B |
| bridge_second | [1,2] | B |
| bridge_macro | [0,1] | 4P+log(100/4)+1000-2B |
| switch_first | [0,1] | S |
| switch_second | [1,2] | S |
| switch_power | [0,1] | log(110/100)-2S |
| reshape | [0,1] | T |
| inner_reference | [0,1] | 10(C+P)-T-8 |
| axial_restore | [0,1] | 1 |
| restore_buffer | [-7,-6] | 1 |
| actual_patch | [1,e] | 1; subcell width log(right/left) |
| Rh_reference | [-5,0] | 5 |
| O2_slope | [0,1] | 1 |
| O2_axial | [0,1] | exp(M)-1; subcell exp(M*left)*expm1(M*(right-left)) |
| O2_buffer | [0,11] | 11 |
| O3_slope_mu | [0,1] | 1 |
| O3_power | [0,1] | W; original-power-offset cell width right_offset-left_offset |

The existing checked affine-radius theorem also supplies all16 exact neighboring radius seams. It proves geometry/phase continuity only, not all velocity or stress derivatives. The O2-buffer/O3-slope seam uses the original relation P=exp(M)+11.

Microscopic B/S and B*s_c/4 stay nonzero formal factors. Mixed macro/power widths keep separate regular and microscopic terms; strict positivity is checked on the full directed width. No rounded huge endpoint difference, radius cap or midpoint replaces a length. All widths/endpoints and the original spatial phase are Z independent.

## True-width C1 integral adapter

For the original five rates lambda=1,3/2,3/2,1,0 the adapter encloses

```text
D_lambda = exp(-lambda*w)
M_lambda = integral_0^w exp(-lambda*t)dt
I_j(Z)   = integral_left^right exp(-lambda*(right-y))*f_j(y,Z)dy
I_j_Z(Z) = integral_left^right exp(-lambda*(right-y))*f_j_Z(y,Z)dy
```

Entire native radial/Z/actual-phase signed density C0/Z covers are multiplied by the positive mass. Density rows already use ordinary y=log R, so the true-width calculation accounts for the native Jacobian once. There is no extra width/Jacobian/Pstar/N normalization.

Tiny-width masses retain the formal positive width times a directed near-one exponential-average cover. Positive-rate tiny attenuation has a directed cover and is never declared exactly1. Moderate/long masses use (1-exp(-lambda*w))/lambda with a directed nonzero exponential upper tail when necessary. Long positive decays remain formal. Rate0 has exact coefficient1 and formal massw, preserving pressure C0/Z memory.

Executed whole-Z-interval cells at candidate N=1024:

| Chart | Actual native endpoints | Scope |
|---|---|---|
| inner_reference | [.12,.15] | True very long source width, all five signed C0/Z contributions |
| O2_slope | [.12,.15] | True width.03, all five signed C0/Z contributions |
| O2_buffer | [5.33,5.34] | True width.01, all five signed C0/Z contributions |
| O3_slope_mu | [.53,.54] | True width.01, all five signed C0/Z contributions |
| O3_power | original power offsets [.53,.54] | True width.01; native phases offset/Tw, all five signed C0/Z contributions |

These are local function integrals. They are not adjacent to each other and do not supply the missing incoming corrections at their left endpoints. Signed phase range hulls remain conservative; tight oscillatory cancellation is still needed for quantitative terminal targets.

## Actual initial histories, beyond a single inlet point

The accepted two-sided left collar is native phase[s_c/4,3s_c/4]. On the whole collar and whole Z domain, the original kappa excess bound is at least2*eta. The existing original flat-branch graph then proves q and its derivatives, A/B and their derivatives, and all changed signed density C0/Z functions vanish exactly.

The true inlet-to3s_c/4 width is **B*s_c/4**, kept formal. The affine operator carries the checked correction initial condition zero to that endpoint. Original background histories, their ordinary Z rows and separate P0/P0_Z are queried at the true right endpoint and retained. Thus this segment now has actual own C1 history covers, not just an operator with an unknown incoming argument. This does not reset the original background memory, and no later bridge interval is inferred quiet.

Two native Z queries execute: whole[-1,1] and[.49,.51]. Their40 checked rows comprise20 initial correction C0/Z rows and20 actual own-history C0/Z rows. The active portion after3s_c/4 must still be integrated from these known initial rows.

Evidence:17 positive native lengths and16 reused exact radius seams;40 actual initial history C0/Z rows;50 signed true-chart integral C0/Z rows;17 independent original absolute-radius difference references;102 independent mass/attenuation comparisons;30 signed nonzero-incoming C0/Z comparisons across nonlinear patch/axial and quiet coordinates; tiny/huge retention and forbidden collar-extension guards. Read-only review: **GPT-5.6 Luna / max**. 1052 working/index hashes pass. Production 25.765s, focused checker 22.218s on the same live native source.

## Agent tasks, in dependency order

- [x] **LEFT4c3-original-inlet/P0_Z/common-basis/local-serial:** original inlet C1 memory and directed rebasing; adjacent O2[.12,.15] affine transfer. See predecessor.
- [x] **LEFT4c3-true-chart-length-backend:** all17 directed positive physical lengths; mixed microscopic terms, nonlinear patch/axial and original O3 offset geometry; original16 radius seams.
- [x] **LEFT4c3-general-true-width-C1-adapter on five declared cells:** finite/tiny/huge mass and attenuation, same signed density C0/Z functions, one Jacobian and fixed-Z differentiation.
- [x] **LEFT4c3-actual-initial-flat-collar-transport:** true inlet sc/2 ->3sc/4, wholeZ and intervalZ, original background/P0/P0_Z retained.
- [ ] **LEFT4c3-active-first-bridge:** start at3sc/4 with the now-known correction C0/Z rows. Resolve native source/cutoff/inverse covers on the remaining actual bridge interval through phase1. Subdivide by original cutoff/source branches as needed; do not initialize at an arbitrary easier phase, taper the source, or extrapolate the collar.
- [ ] **LEFT4c3-first-bridge-whole-Z:** current bridge root/cutoff boxes can straddle branches. Preserve the original positive eta/Delta and all derivative factors; return covers for the full requested Z functions, with actual partition boundaries and errors.
- [ ] **LEFT4c3-second-bridge/macro/microswitch:** use exact neighboring radius seams and true B/S/macro lengths. Carry the same five correction C0/Z functions across all active and quiet pieces; no reset at chart switches. Every missing interval must remain explicit.
- [ ] **LEFT4c3-reshape/restoration/patch-route:** integrate every intervening chart, not just the declared local cells. Keep original analytic background histories and P0/P0_Z separate. Original patch log-coordinate width and ordinary-y density normalization must be used once.
- [ ] **LEFT4c3-full-O2-coverage-and-axial-crossings:** complete Rh_reference/O2 slope/axial/buffer with their original true lengths. Resolve O2_axial whole-Z cutoff crossings and preserve buffer[0,9]/[9,11] admission distinctions.
- [ ] **LEFT4c3-inlet-to-local-incoming:** provide actual correction C0/Z functions at O2 coordinate.12, then apply the existing adjacent operator. The present five local cells and initial collar do not fill this route.
- [ ] **LEFT4c3-signed-oscillatory-cancellation:** implement phase primitive/zero-mean integration-by-parts bounds, separating quadratic averages and first-Z slow variation. Avoid allocating huge numbers of cycles or inferring signs from broad phase hulls.
- [ ] **LEFT4c3-rebase-error-control:** quantify lost radius/amplitude correlations on long cells; refine original correlated factors or certified partitions when covers are too wide, without selected field caps.
- [ ] **LEFT4c3-right-flat-collar:** prove exact supported quiet segments across the O3-slope/power seam and carry inherited pressure C0/Z memory. The current accepted right collar is not power-only.
- [ ] **LEFT4c3-global-Rc-C1-functions:** complete the actual inlet-to-Rc route including all original charts, signed integrals and initial data, then combine with endpoint background/P0_Z.
- [ ] **LEFT4d-terminal-control-functions:** derive A/A_Z and divided(J-M)/mu from the completed histories; use the reserved independent inverse with formal inverse-mu. Verify all five terminal identities as Z functions.
- [ ] **HIGH/CONT/OUTER/ENERGY/LEFT4e:** higher phase/primitive recurrence, same-function C4 seams, pressure/heat/finite energy, inherited constants and one finite N/global stress cone.
- [ ] **REC/WAVE/PHYS:** actual n-dependent recovery and independent repairs/smooth sum; quadratic stress-canceling oscillations; corrected Cartesian NS and multi-time dynamics.

Next production starts at the known3sc/4 endpoint and addresses the active bridge. Do not rerun unchanged inlet/local C1 stages. Record code/result/check commits, exact domains and remaining gaps before marking any task complete. Scoped gate: `current_original_true_chart_lengths_factored_C1_integrals_and_initial_collar_transport_executed`. Global terminal five moments, common finite N/cone, coefficient recursion and full NS remain incomplete.
