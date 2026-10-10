# Exact original leading source kernels, C3 ODEs and seams (2026-10-10)

Full reconstruction **ACTIVE / INCOMPLETE**. Source [ad0c9c4a](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/ad0c9c4a7ad11e1e38fdf65dbb25315d96051bef) installs explicit mathematical leading functions for Rh_reference, O2_slope, O2_axial, O2_buffer, O3_transition and O3_power. These functions use the original defining integrals, actual parameters and complete predecessor memory. The independent checker establishes their ordinary Z0..3 derivatives, all five own-rate ODEs, all source state/radius seams and the true x=2 Rc endpoint. It does **not** yet identify the old opaque outer Duhamel/provider leaves with these explicit functions or close the repair-band leading inlet identity.

## Concrete new functions

The exact common units are

```
A(Z) = 1/(1+Z^2), B(Z) = 4Z/Pstar
m_y = V-m; h_y = E-3h/2; k_y = E*V-3k/2
e_y = V^2-E^2/2-e; p_y = E^2/2
```

Raw V/m/k are divided by Pstar once. Common histories do not gain an Rm multiplier. The original independent P0 remains in the accepted predecessor chain; it is not reset or fitted by this module.

| Segment | Actual native domain | Exact kernel/source construction |
| --- | --- | --- |
| Rh_reference | y in [-5,0] | E=A exp(y/10), V=B; exact original reference five histories |
| O2_slope | x in [0,1] | J(x)=integral_0^x sigma; theta/pressure/energy weights (8/5,3/5), (1/5,6/5), (6/5,6/5) |
| O2_axial | phase in [0,1] | y=exp(40 phase); K_j(y)=integral_1^y exp(u-y) beta(u)^j du, j=1,2 |
| O2_buffer | tau in [0,11] | y=exp(40)+tau; retained K_j(y)=exp(-tau)K_j(exp(40)), V=0 |
| O3_transition | x in [0,1] | weights exp(u-mu J), exp(-2mu J), exp(-u-2mu J) |
| O3_power | x in [0,2] | accepted exact power semigroup from the actual transition endpoint; x=2 is Rc |

The source sigma is zero at x<=0, one at x>=1, and

```
sigma(x)=exp(-1/x^2)/(exp(-1/x^2)+exp(-1/(1-x)^2))
beta(u)=sigma(1-log(u)/40)
```

inside its transition. Squared inverse exponents are essential. J and all weighted kernels are mathematical definite-integral graph nodes. Enclosure algorithms `transition_integrals`, `turnoff_kernels` and `transition_kernels` provide bounds; their outputs are not selected as exact function values. The axial measure is true dy, with original native Jacobian 40 exp(40 phase) applied once. The buffer retains the full nonzero axial kernel memory.

## Independent evidence and limits

- 34 defining source AST predicates bind scalar sigma, the imported alpha/sigma chain, primitive updates, integrands, exact cell weights and kernel order. The scalar identities between the direct cutoff, logistic odds and 1-alpha definitions are checked algebraically.
- 10 defining kernel/primitive identities check integrand, limits, dummy variable and nested J. All 16 true integral nodes are admitted separately from scalar covers.
- 288 profile/history/own-rate rows establish ordinary Z0..3 derivatives, including the exact A and B source definitions.
- 120 independent five-history ODE rows use the fundamental theorem of calculus. Kernel functions are abstracted only after their defining integrals are checked, with exact derivative laws retained.
- 240 state and own-y C3 rows establish the five Rh/O2/O3 source seams. Twenty terminal rows identify the explicit O3 x=2 endpoint.
- Six actual native domains, original offsets and Jacobians agree with the accepted current geometry. All five log-radius seams match.

The reused read-only reviewer is **GPT-5.6 Luna / max**. Root implements, computes and publishes. Producer flags remain false; the separate checked receipt enables only

```
current_original_outer_leading_C3_source_ODE_and_seam_functions_installed
```

The explicit source functions are now available through

```python
from lei_ren_part1_paper_compliant_current_original_outer_leading_C3_identity import CurrentOuterLeadingC3Identity

leading = CurrentOuterLeadingC3Identity()
leading.functions['charts']['O2_axial']
leading.functions['actual_source_Rc_leading_C3']
```

These are exact formal function handles, not a numeric point field. Installed quartet: `lei_ren_part1_paper_compliant_current_original_outer_leading_C3_identity.py`, `.json.gz`, `_check.py`, `_check.json`.

## Completed bounded tasks

- [x] **LEADING-SOURCE-ODE-MAP** Recover explicit original Rh/O2/O3 leading profiles and five histories from their defining kernels, with correct common units and ordinary Z0..3 derivatives.
- [x] **LEADING-SEAM-DATA** Establish the explicit source leading state/own-y C3 seams and original radial coordinates, including the actual axial Jacobian and retained buffer memory.
- [x] **LEADING-EXACT-KERNEL-SOURCE-BINDING** Bind the defining scalar sigma/alpha/beta chain and all slope/turnoff/O3 integral definitions independently; keep caps distinct from functions.

## Immediate next task: close the old leading interface

- [ ] **LEADING-PROVIDER-RECIPE-FUNCTION-BINDING** For each accepted current leading leaf, bind its quantity path and original pure frontend projection to this module's explicit function. Inspect actual Rh/slope/axial-buffer/O3 background assignments and normalized recovery. A matching AST hash or shared source family alone does not prove function identity.
- [ ] **LEADING-RH-INLET-TWENTY-ROWS** Identify the accepted original Rh leading packet at y=-5 with the explicit source inlet for all m/h/k/e/p Z0..3 rows. Preserve the independent P0 guard.
- [ ] **LEADING-DUHAMEL-UNIQUENESS** Bind actual leading densities E/V and each accepted Duhamel incoming history to the explicit source. Use the checked common own-rate ODEs and inlet data for uniqueness, segment by segment, with the true native-to-physical Jacobian. Do not zero predecessor memory.
- [ ] **LEADING-RC-TWENTY-ROWS** Identify the accepted outer O3 leading right endpoint at x=2 with this explicit source endpoint, then with the genuine compact band's leading x=1 seeds for all 20 ordinary rows. Bind the canonical source endpoint and actual amplitude definition, not an archived magnitude cover.
- [ ] **CURRENT-OUTER-LEADING-TO-BAND-SEED** Enable the leading gate only when all actual leaf/Duhamel/seed function identities above are independently proved.
- [ ] **COMPLETE-RC-INTERFACE** Combine that leading identity with the already accepted correction Q/N identity. Prove complete histories and independent pressure/P0 in the genuine Rc inlet units before enabling the complete interface.

## Subsequent tasks

- [ ] **CURRENT-RP-NATIVE-FRAME-IDENTITY** Identify the genuine compact-repaired C3 quiet Rp frame with the selected native pulse constructor, including its common Rw parent, all amplitude quotient derivatives, Mp and original P0. Prove actual consumption of this frame.
- [ ] **CURRENT-C2-EXTERIOR-SEGMENT-MEMORY** Continue selected pulse/flatten/postpulse/collar histories through Z2 with exact lower aliases, true Jacobians and complete incoming memory.
- [ ] **CURRENT-C2-ABSOLUTE-FUTURE-INTEGRALS** Close full absolute Gamma/heat terminal identities through Z2, using the final heat amplitude and correct reference power. Relative correction zero does not prove absolute matching.
- [ ] **CURRENT-SELECTED-EXTERIOR-MIXED-CONTRACT** Recover selected pulse/end/flatten/collar high-y mixed derivatives and true seam identities.
- [ ] **CURRENT-GLOBAL-MIXED-DERIVATIVES** Join all current core, transition, Rh, outer, compact repair, quiet and selected exterior functions under one compatible contract, including corrected high-y Rh data.
- [ ] **CURRENT-FREQUENCY-AND-NUMERIC-ORACLE** Build usable original-frequency, phase inverse, Picard and integral evaluators with quantified errors. Exact formal N is not a numeric point oracle.
- [ ] **CURRENT-PHYSICAL-TIME-AND-CARTESIAN-FIELD** Assemble actual time/anisotropic maps into numeric u/v/w with axis limits and error bounds.
- [ ] **CURRENT-STRESS-CONE-AND-REMAINDER** Recover signed divergence stress, flat remainder and regionwise cone margins on the full matched background.
- [ ] **CURRENT-TEMPORAL-RECURSION** Implement original n=1 and n>=2 equations, independent repairs, common inner domain, curl-preserving truncation and smooth-sum remainder. Coordinate scaling alone does not complete this task.
- [ ] **CURRENT-OSCILLATORY-CANCELLATION** Construct both original oscillatory families and averaged quadratic cancellation after stress and recursion are accepted.
- [ ] **CURRENT-DYNAMICS-AND-FULL-RESIDUAL** Measure actual contraction, aspect ratio, amplification, material winding, interscale recurrence and finite energy; independently validate final forced Cartesian NS residual.

Read this handoff first, then [accepted correction inlet bridge](CURRENT_ORIGINAL_OUTER_TO_REPAIR_C3_BRIDGE_2026_10_10.md) and [actual compact exit / quiet Rp frame](CURRENT_ORIGINAL_C3_POWER_TO_RP_2026_10_10.md). Claim one bounded task, publish source and independent evidence, then mark that task complete. Full reconstruction remains ACTIVE / INCOMPLETE.
