# Actual inlet pressure coefficients and n=1 source-vector identity

Tasks 14/15 remain **IN PROGRESS**. The long-term goal remains **ACTIVE / INCOMPLETE**. This stage recovers the known signed O3 background coefficient, implements a scale-correct inlet functional on the actual unknown history, evaluates its coefficients on whole axial cells, and proves the forward/canonical n=1 known-source vector identity. No finite-N error bound is selected as an inlet value, and no complete n=1 solution or terminal closure is claimed.

## Actual O3 inlet functional

The actual normalized background mean is m=4Z*mstar*exp(-12)/Pstar at O3 native0. The rational pressure recovery gives Q0=-4*mstar*exp(-12) for every Z on the original domain, with exact ordinary Z1/Z2 equal to zero. This follows algebraically from the actual background source and rate-one history; matching three sampled enclosures is only a consistency check. The saved true signed mean enclosure and its original positive-memory receipt are reused without rerunning that integral. The original Q0 source remains Pstar times background.generic.Q at native0; its numerical enclosure is separate from the function node. Q0 is strictly negative.

Write H=N*delta-m_in for the actual signed normalized incoming correction. The original incoming mean already includes its one N division. Multiplication by N defines H; the recombination below restores Pstar/N once. With D=1-Z^2 and L=1-delta*Z^2:

```
qhat = -[Z*(1-delta)*H + D*H_Z]/L,
Delta m15_O3(x) = R0^2 * [ (Pstar/N)*A + (Pstar/N)^2*B ],
A = Z*(1-delta)*qhat_Z*expm1(x)/(4L) - Q0*qhat*x/2,
B = -qhat^2*x/4,                         0<=x<=2.
```

The background cross term is linear in Pstar/N after using the physical Q0; it is not incorrectly left in a Pstar-squared sector. This is the finite-N correction increment relative to the background, not full O3 m15. Incoming pressure, both nonlinear factors, exact N and ordinary Z1 are retained. R*Delta Q_y is constant along the quiet chart and its endpoint difference is zero.

`original_O3_inlet_functional.py` requires four explicit ordinary H rows and an explicit constant-background Q0 jet. Missing H_Z3 and Z-dependent Q0 are rejected. It constructs 133 local expression nodes from the actual incoming signed source ports. H is neither zeroed nor assigned a symmetric error-center value. Pstar, N and physical R0 are not materialized.

The linear form on H0..H3 and the symmetric quadratic form are also implemented. Their ordinary-C1 coefficients are interval-evaluated at x=2 on the two full cells [-3/10,0] and [0,3/10], at160 decimal digits. This yields16 linear and64 quadratic directed coefficient enclosures. The bounded original delta parameter graph (10 scalar nodes) is evaluated; archived owners and history functions are not. These are coefficient ranges of an operator acting on unknown H, not numeric H values or completed pressure increments. This cell coverage is not global five-moment closure.

Independent symbolic controls compare the derivative of the new primitive against the original reduced-S difference for arbitrary differentiable H, check ordinary C1, zero inlet increment, viscous endpoint cancellation, universal background recovery, missing-row rejection and four coefficient reconstruction identities.

## n=1 actual source-vector identity

The accepted common fixed-point receipt binds the actual leading input functions and all finite rows/tails to the unique analytic n=0 solution. The new adapter verifies its 557 inherited input hashes and descendant source-family/hash consistency. Source ASTs bind the actual forward seed, its gradient/ell/S inputs, canonical F0 and its derivative, F0*Phi, U0+epsilon*Psi and physical Q recovery. The accepted scaled-S receipt supplies inclusion of the exact epsilon^2*F0^2 jets. The numerical zero-S proxy does not become that function.

Let Z_a(v)=[a*Z*v+(1-Z^2)*v_Z-2Z*rho*v_rho]/L, x=sqrt(rho), and epsilon=exp(-logLambda). The two sources have

```
N1theta = -Z_(-3) Z_(-2-delta)(F0*Phi),
N1z     = -Z_(-2) Z_(-1-delta)(U0+epsilon*Psi),
N1p     = -S/2 + (2Q_rho+rho*Q_rhorho)/epsilon.
```

The canonical expanded axial_second formula is independently proved to equal the weighted composition for arbitrary weight and arbitrary source function. The forward pressure S/head and canonical S/head are independently assembled from the pinned formulas. At every radial order n, the head coefficient is (n+1)(n+2)Q_(n+1); no finite-prefix computation is replayed.

With source order (pressure, angular, axial),

```
b_forward = (N1p/Pstar^2,
             2*epsilon^2*N1theta/F0,
             2*N1z - 4*epsilon*Z*rho*N1p/L),
Dscale = diag((2*epsilon*Pstar^2*x)^(-1), epsilon/F0, epsilon^(-1)),
Dscale^(-1)*b_forward = (g4,g5,g6).
```

P2 is explicitly substituted as Pstar^2 in the symbolic algebra. Source guards tie the forward P2=exp(2*logP) to the accepted common logP and tie both producers to logLambda=4*logP+1000 and epsilon=exp(-logLambda); the exponential square and reciprocal identities are checked.

In particular the raw axial pressure-feedback term is -4*epsilon^2*Z*rho*N1p/L. The extra epsilon comes from the axial entry of Dscale inverse. The proof holds for rho>0, L!=0; the regular pressure source g4 extends as zero at rho=0. Dscale itself is not evaluated at the axis.

The true P0 definition/scaling is inherited from the accepted report, including all14 sources and the flatten function. Its flags are premises, not tautological substitution tests. The exact relation F2=F0^2 identifies the true homogeneous input; forward S/(epsilon^2) boxes remain enclosures, not selected values. All 14 new source/algebra checks pass. This resolves the known source-vector seam recorded in the preceding milestone; it does not certify the complete solved n=1 numerical field.

## Actual run and preserved state

New production exits0 in 4.638 s. Report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_inlet_source_vector.json` is 4,452,506 bytes, SHA256 `2ac1269c3da29466067d1dfd50a92be321e85291b2634ffa6d2eafd679e9941e`. It pins 1973 input files and 284 implementation files. Saved-only inspection validates all2257 pins and packet invariants; compilation passes. No accepted owner, producer, constructor, fixture, core recurrence, prepare step or quadrature is rerun. New bounded symbolic algebra and scalar/coefficient interval arithmetic are the executed work.

All three Rc178 discrepancy views and their error terms remain exactly unchanged. Original raw213/reduced183 and quiet179 archives, six signed core/lower contributions, four signed leading leaves19/20/41/42, remaining O3 lower endpoint enclosure, all volume caps, external viscous packet and pressure throughRv remain. The target and quiet pair were already consumed; do not subtract them again. P1/m13 remains distinct from m15.

Three fresh fixed `luna_max_worker` children with `fork_turns:none` handled the inlet-data scan, n1 source proof and independent review. Requested routing was GPT-5.6 Luna/max; tool responses exposed task/nickname rather than independent runtime-model metadata. No Astra child was requested. Root integrated source binding, scale separation, coefficient tensors, the axial-composition proof and publication.

## Remaining work

The saved finite-N O3 inlet correction rows are unsigned cap/error records, not genuine signed numeric history values. The n1 O3 pressure-output intervals are not inlet m/q data. The new functional therefore leaves H unresolved. A new bounded evaluator must consume the true signed O2-buffer/O3-transition history integrals, including physical measure, ordinary Z0..3, original phase and selected N. It must not recycle cap centers or reset incoming history. True outer source ports are recorded in the new report; the older C3 outgoing leaves [4820,4821,12010,22503] locate the same predecessor stage but are not promoted here into a numeric oracle.

For the core response, next prove the general radial coefficient conjugacy of the actual forward _advance recurrence and the canonical regular raw system, then attach a genuine radial Taylor remainder or another direct canonical function enclosure. The existing homogeneous matrix comparisons and new known-source identity now supply its inputs. A Volterra-word remainder is not automatically a radial Taylor remainder. Complete core response equality and the q-dependent collar bridge are still required before the six signed core/lower terms can be identified or removed.

Global functional all-five/common-N terminal closure, exact heat and energy, admissible stress/flat remainder, true temporal n>=2 recursion and smooth summation, realizable mean/oscillatory correction and independent corrected Cartesian NS remain incomplete. The final corrected-NS1e-3 threshold does not apply to this source stage.
