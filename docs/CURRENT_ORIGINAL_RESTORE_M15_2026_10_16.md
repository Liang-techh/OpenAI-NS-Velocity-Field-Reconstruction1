# Original changing-source cumulative pressure through Rm - 2026-10-16

`OriginalRestoreM15` continues the complete124-term reference moment through axial restoration and postrestore toRm. Each chart appends16 signed/error terms:140 after restoration,156 after postrestore. Incoming terms are retained exactly once. These are numerical signed proxies with complete directed absolute errors; selected-N signed correction values and useful total relative accuracy remain unresolved.

The full goal is **ACTIVE / INCOMPLETE**; tasks14/15 remain IN PROGRESS. This is pressure-history integration throughRm, not terminal all-five closure, global heat matching, an admissible-stress result, genuine temporal n>=2 recursion, or corrected NS validation.

## Actual changing-source mean and scalar kernels

Use the same original source, exact selected N and genuine five-jets as the accepted reference pressure adapter. With a=1-sigma(t), g=exp(-t)+Kmean(t), w=exp(-t),

```
U = 4Z+(j+epsilon*Psi_Rz)*a
M = 4Z+(j+epsilon*Psi_Rz)*g+epsilon*D_Rz*w
D_Rz = (Mpsi_Rsh-Psi_Rsh)*exp(-Lref).
```

The reference difference D_Rz is obtained by the original coefficient-before-decay routine with its complete proof/bound. It is not reconstructed by independently subtracting the final correlated Psi/Mpsi interval boxes. The changing j baseline remains in epsilon degree0; epsilon degrees1/2 are never merged into it.

Pressure recovery is linear in a,g,w. The F/G half-S1 density coefficients are therefore polynomials in these scalar variables of total degree at most2. Keys use `(epsilon_degree, scalar_code)`, with code=i+3j+9k for a^i*g^j*w^k. Since the total degree is at most2, adding codes in products has no base3 carry.

`OriginalRestoreKernels` computes all ten scalar kernels

```
I_ijk(t)=integral_0^t exp(2s)*a(s)^i*g(s)^j*exp(-k*s) ds,
i+j+k<=2.
```

The original a is decreasing. From g'=a-g and g(0)=1, one has0<=a<=g<=1 and g is decreasing. For a cell[A,B], h=B-A,

```
g(B) in exp(-h)*Gamma_A+(1-exp(-h))*[a(B),a(A)].
g([A,B]) in [lower_Gamma_B,upper_Gamma_A].
```

Every scalar monomial uses these whole-cell boxes and the exact positive exp(2s) weight `exp(2A)*expm1(2h)/2`. The uncertain decimal endpoint segment is included separately with a positive upper bound. The directed enclosure defines truth; separate RK4/Simpson point values are numerical proxies with their full errors retained. Refinement acts only on these scalar cells and never changes an original source parameter.

Pure exponential modes have exact primitives: I000=expm1(2t)/2, I001=expm1(t), I002=t. These retain their exact exponential correlation and avoid unnecessary cell accumulation. The original run uses85-digit point arithmetic and115-digit directed arithmetic; independent scalar reference quadrature and serialized-row audits use190 digits.

## Exact source identity and both pressure endpoints

The driver executes the original `pressure_functions` reduced-S1 source. Its F/G coefficients already include the one-half factor in S1/2:

```
F = Q/(2A)+Q^2/4-Z^2*U^2/A^2
G = -Q/(2A)+Z*(1-delta)*Q_Z/(4A)-Q^2/4
    +(1-Z^2)*U*Q_Z/(2A)+Z^2*(1-delta)*U^2/A^2
    +Z*(1-Z^2)*U*U_Z/A^2
A=1-delta*Z^2.
Delta m15 = [R^2 F] + integral R^2 G dy - [R Q_y].
```

No additional1/2 multiplies G. The complete varying-U identity is used rather than extending the constant-U frozen coefficients. Ordinary C1 coefficients retain U/Q through Z2 and mean through Z3. All four signed baseline/epsilon RQy endpoint terms remain. Finite-N endpoint, volume and viscous budgets keep both actual endpoints; the full-prefix uniform background and min(t,1/2) positive mass control volume errors. Finite-N pressure errors receive no extra epsilon or hb factor.

## Physical scale and postrestore

Restoration starts at logR=log110+10logC+10logP-8. Postrestore starts with offset-7 and stops at offset-6, Rm. Each local integral retains the base log radius, epsilon degree, fixed unit offset and local radius growth as separate atoms. No huge physical radius is materialized.

For postrestore p0..1, U0=4Z, Psi=Uy=0, while M0 and Mpsi decay from the actual restoration endpoint. The G coefficients use modes exp(-kp), k0,1,2, with exact finite kernels `expm1((2-k)*p)/(2-k)` or the resonant p. Both signed RQy endpoints remain even where their exact background increment cancels. Postrestore starts from all140 restoration terms; it does not reuse the original124-term inlet or reset the mean.

## Executed evidence

The focused run exited0 in272.883 seconds; saved-report inspection exited0. It evaluates26 original restoration/postrestore packets, including an uncached Z=.17 moment inlet. The exact source and each incoming124/140-term packet are preserved.

- 26 exact inlet replays,832 serialized new coefficient/error conditions and192 exact empty rows.
- 80 retained lower finite-N endpoint conditions and80 actual signed RQy rows;338 separate unit-offset conditions.
- 130 scalar-kernel conditions,78 retained correlated reference-decay proof rows, three independent original-flat-source gamma inclusions.
- All ten scalar widths contract under128-to256-cell refinement or retain the identical exact exponential primitive. The uncertain decimal endpoint segment remains included; five invalid queries are rejected.
- The driver executes the original reduced-S1 pressure identity. An independent changing-source polynomial-alpha fixture checks18 physical value/Z1 integrals,30 scalar monomial inclusions, three gamma inclusions,180 polynomial coefficients and the exact varying-U mean ODE.

The39,205,706-byte report SHA256 is `8d432e3fa0c5dbba2d7bcbcb95db46fa8b206dd69460b3853bc0f3f7e7d1cab7`. It pins1768 input sources and73 implementations. Inspect with `python scripts/run_original_restore_m15.py --inspect-report`; do not rerun accepted ancestors. The results establish the stated cumulative adapters/enclosures throughRm; they do not establish global five-moment closure or useful relative cancellation of the full physical sum.

Development attempts retained in the external checkpoint: a100-digit independent target was insufficient to resolve a160-digit fixture interval, so the independent target was raised to190 digits; a later refinement assertion exposed avoidable pure-exponential cell roundoff at the actual115-digit original interval precision. Pure exponential modes now use exact primitives, and a focused ten-kernel check passed before the final original-source run. No physical equation, original parameter or source was weakened to pass these checks.

Existing read-only `/root/n1_atlas_map`: **GPT-5.6 Luna / max**, reused without a new agent. It confirmed the monotone gamma recurrence, polynomial encoding, coefficient/endpoint orientation and source nesting, then identified the genuine next Rm derivative gap described below.

## Next actual source connection toward Rp

1. Start the active Rm patch from all156 cumulative pressure terms and `post.rows.inputs(Z,iv.mpf(1))` genuine source jets. Keep the signed moment state distinct from older truncated public velocity/six-history APIs.
2. The existing numeric Rm defects/implicit controls are only Taylor Z2 and recovered Q is only Z1. This is a genuine derivative gap for pressure C1. Read `original_Rm_defect_numeric.py`, `original_Rm_patch_numeric.py`, `original_Rm_finite_N_numeric.py` and `original_Rm_patch_velocity.py`; do not feed their truncated public rows into QZZ or zero-pad a missing meanZZZ.
3. The pinned `current_original_true_Rm_patch_C6_functions.py` and its json.gz archive define the full nonlinear same-Jacobian control recurrence through Z6, actual fields/radial Ey,Vy and generic C5 recovery. Use those source identities to build genuine Z3 defect/control numerical enclosures. The archive explicitly leaves control C6 numerical containment/mixed-norm admission unproved; loading its graph is not numerical acceptance. Extend the existing implicit interval inverse to the required third row with actual source derivatives and a quantitative containment receipt.
4. Extend finite-N history transport from Z0..2 to its genuine C3 source rows, using caps only as errors. Preserve actual By/ByZ phase/slow bounds for Uy/UyZ. The identity Q+Qy=((1+delta)ZU+2ZUy-(1-Z^2)Uz)/A can help recovery; it does not remove QZZ from the G_Z integral.
5. Build actual patch U/Q/Uy/Qy pressure rows with the correct generic derivative owner, physical x=R/Rm coordinate, original dx/x measure and complete incoming prefix. Preserve all controls, their source errors and finite-N loop histories.
6. Integrate only the new patch interval into the156-term state with the exact F/G/RQy identity, actual positive kernels and both signed/finite-N endpoints. Expose source-correlated value/Z1 data for later independent Rp comparison.
7. Continue outer reference/slope/axial/buffer/O3/compact regions from their actual source definitions and native measures. Existing numerical velocity adapters do not prove pressure moment closure; retain cumulative inlet and mean jets at every seam.
8. Independently compare the complete signed pressure history against the reserved Rp pulse source; then complete all-five histories/global mixed joins/common N, original flatten/preheat/exact heat and finite-energy tail.
9. Establish admissible stress/flat remainder, then implement true n-dependent temporal coefficient recursion and smooth summation, realizable mean/oscillatory stress cancellation and independent corrected Cartesian NS validation.
