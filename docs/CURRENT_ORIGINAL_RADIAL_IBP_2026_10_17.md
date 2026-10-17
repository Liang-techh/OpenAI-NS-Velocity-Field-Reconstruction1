# All-order n1 radial coefficients and actual C5 signed-history IBP

Tasks 14/15 remain **IN PROGRESS**. The long-term goal remains **ACTIVE / INCOMPLETE**. This stage proves the general radial coefficient rules of the actual temporal n=1 core and connects the signed mean integration-by-parts scheme to actual C5 source/slow-derivative ports. It does not compute signed inlet values, certify a radial remainder, remove any core/lower term or close the terminal moments.

## General radial coefficient rules

The actual forward state is normalized by f=epsilon*F1/F0, u=Uz1/epsilon, k=K1/epsilon and p=P1/(epsilon*Pstar^2). The canonical raw state is W=(F1,Uz1,K1,P1,F1x,Uz1x), x=sqrt(rho), D=(0,0,2,0,3,1). The previously accepted 30 homogeneous matrix comparisons, 36 mixed prefix comparisons and known-source vector identity are read from pinned reports. They are not recomputed in the published production. Source ASTs bind the actual _advance/_rhs, canonical raw equations and required leading inputs.

For every n>=0, with m=n+1, the new general rules are

```
f_m = [rho^(m-1)]R_f / (4*m*(m+1)),
u_m = [rho^(m-1)]R_u / (4*m^2),
p_m = [rho^(m-1)]R_p / m,
k_m = -m*u_m/(m+1).
```

The canonical pressure row has factor2m; its source-to-forward normalization supplies the factor2. Angular and axial rows have factors4m(m+1) and4m^2. Zero-axis data fix the constants. The pinned RHS is triangular in the unknown radial degree, so these nonzero factors determine the formal series successively. Ordinary Z differentiation commutes termwise with radial coefficient extraction. The actual source budget for coefficient rho^n at ordinary Z order k requires leading radial rows through n+1 and axial rows through k+2; no missing rows are padded with zero. The pressure head coefficient is (n+1)(n+2)Q_(n+1).

This is an all-radial-order formal identity at temporal n=1. It is not the requested temporal n>=2 recursion. The accepted directed prefix through radial degree16 remains an enclosure of coefficients, without a selected midpoint or an inferred function value. Full n1 function equality and its radial remainder remain false gates.

An exploratory worker draft repeated the old 30 matrix comparisons. Root rejected that path and replaced it with receipt consumption before production. The published runner executes only the new universal coefficient/commutation rules and source binding; no accepted numerical recurrence or quadrature runs.

## Actual C5 mean-history ports

The new source-only adapter reads the saved true outer C6/C5 graph, the accepted original zero-phase-mean receipt and pinned ASTs. For O2_buffer and O3_transition it binds the actual four ordinary Z0..Z3 rows of B, fixed-phase slow_B_y, incoming/outgoing mean, native coordinate, log-radius measure, original global phase and selected common N. The transition incoming rows are literally the buffer outgoing rows. The predecessor history remains nonzero and is retained.

Use the accepted anchored phase primitive C(y,Z,phi)=integral_0^phi B(y,Z,theta)dtheta. Original B parity gives C(0)=C(1)=0, and C_phi=B. The slow derivative is the phase primitive of actual slow_B_y, taken at fixed fast phase; it excludes the N*B fast term. C5 regularity permits ordinary Z0..Z3 and the slow derivative to commute with this integral. The derived issued phase node is retained as a source handle rather than misidentified as a bound variable.

For H=N*delta-m, rate one and log-radius length L, the exact update is

```
H_out = exp(-L)*H_in
      + (C_out-exp(-L)*C_in
         - integral exp(-(y_query-y))*(C+C_y_slow) dy) / N.
```

Both endpoint terms and the integral remain signed functions. The outer integration variable is distinct from the moving query/native coordinate and the phase variable. Physical dlogR has its Jacobian exactly once (the two actual chart Jacobians equal1). Scaling the original mean by N leaves one N division in this update; no extra N is applied. The buffer endpoint is substituted at actual native11, followed by the transition endpoint at native1.

The adapter creates 197 local expression nodes without rewriting the archived graph. These ports make the actual mixed slow derivative and exact signed transport available as functions. Three new symbolic controls verify H_y+H=B, zero-length recovery of an arbitrary nonzero inlet, and exact incoming response exp(-L). The explicit slow-integral operation is distinct from the older operation that differentiated an existing primitive. A signed numerical evaluator for the phase primitive and remaining slow integral is still required. No saved cap, interval center or finite-N error exponent is used as H. No internal phase endpoint is deleted without a same-source identity.

A bounded source audit locates the next useful seam proof. At O2 native11 and O3 native0 the actual geometry offsets are respectively6+exp(40)+11 and17+exp(40). O3 E starts from the saved buffer endpoint times exp(-x/2-mu*J(x)), with J(0)=0; sigma(0)=0 gives the same E_y=-E/2 and V=V_y=0. Rate-one Duhamel starts from the buffer endpoint with a zero-width integral. Deterministic generic recovery and roots therefore suggest the same tuple (E,V,a,b,p2,t0,Delta) as a function of ordinary Z. The prescribed loop then uses one universal function of independent issued phase, with only the bound angle renamed. Common R/Ra, N and phase_origin give the same unwrapped phase. This argument has not yet been implemented as a dedicated executable seam certificate, so the local graph retains both C endpoint nodes. The certificate must explicitly bind the lower O3 endpoint0, not its upper endpoint1, and prove equality as a function of arbitrary phase before contracting them.

## Radial remainder domain audit

The existing n1 Volterra/word bound is uniform on real rho<=4.1 and a complex Z tube. Its real-x simplex factor and Z Cauchy loss do not prove a bound on a complex rho disk. The source/n0 metadata |rho|<5 and Cauchy radius1/2 cannot upgrade that n1 bound.

The regular parity is W1..W4=O(x^2), W5,W6=O(x). A rho-holomorphic vector must use the first four fields together with F1x/x and Uz1x/x, or all six axis-weighted fields. Raw W5/W6 do not have an integer-power rho expansion at the axis. The next core task is a complex-radial regularized source/operator/Green and solution majorant on some R>4, for example R=9/2, in the existing complex Z tube. Only after a proved uniform solution norm M(R) may a Cauchy tail factor (4/R)^(N+1)/(1-4/R) be attached. The current word tail is not reused as that tail.

## Actual run and preserved state

Production exits0 in 2.629 s. Report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_radial_IBP.json` is 4,671,502 bytes, SHA256 `67945b4e9a8f2b1cdaa9a255b42a3bd9e5377f514cb12409ba41e9f17f957af0`. It pins 1974 input files and 287 implementation files. Saved-only inspection checks these pins and the packet invariants; focused compilation passes. Accepted owners, constructors, numeric core recurrence, old fixtures and quadratures are not executed by the published runner.

Rc178, quiet179 and inherited raw213/reduced183 packets are unchanged, including all active errors, six signed inner/lower terms, four signed leading leaves19/20/41/42, remaining O3 lower endpoint, volume caps, external viscous contraction and pressure throughRv. Target and quiet were already consumed once. P1/m13 remains distinct from m15.

The bounded workers used fixed luna_max_worker routing, fresh initial context (fork_turns:none), requested GPT-5.6 Luna/max. Tool responses exposed task/nickname rather than independent runtime model metadata. Root reviewed the proof scope, rejected matrix replay, integrated and published. No Astra child was requested.

NEXT implement the same-source O2/O3 primitive seam certificate, then a directed signed evaluator of the remaining actual C5 phase primitive/slow integral while retaining the older incoming history; independently establish the complex-radial regularized n1 bound. The q-dependent collar response bridge and full six-source equality remain open. Global functional all-five/common-N closure, exact heat/energy, admissible stress/flat remainder, genuine temporal n>=2 recursion/summation, realizable mean/oscillatory correction and corrected Cartesian NS remain incomplete. The final corrected-NS1e-3 gate is not applied to this stage.
