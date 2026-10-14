# True core inlet axial6 and generic fifth derivatives, 2026-10-13

Continuation: [true prescribed micro/macro bridge C6](CURRENT_ORIGINAL_PRESCRIBED_BRIDGE_C6_FUNCTIONS_2026_10_14.md) now propagates this inlet's background functions and six histories through R100. Quantitative C6 containment, R100-to-Rm and modified finite-N/all17 sources remain open.

Status: the same admitted nonlinear core now supplies genuine C6 source functions and directed seven-row numerical enclosures at its rho=4 inlet. Full signed Q, angular/axial inertial sectors and p2 have six genuine Taylor entries through order5. This is the core inlet; the seventeen downstream active charts still require their own higher-order construction. The full goal remains **ACTIVE / INCOMPLETE**.

Use `python scripts/run_original_core_inlet_axial6_generic_source.py` to inspect the sealed source archive. It does not rebuild the core or pressure ancestors. On the original host, reuse retained Python session55770, cache26, `public_core_lower_C3`, `public_pulse_C3_pressure`, `hier`, and the new `inlet6_owner`. Checked archive restoration is available as `CurrentOriginalCoreInletAxial6GenericSource(public_core_lower_C3,hier)` when those admitted source owners are already restored.

The real dependency was one order upstream of p2: a fifth derivative of Q or the full axial inertial expression consumes the sixth derivative of the mean, energy and complete pressure. Earlier downstream arrays named axial4 had only entries0..4. Their internal Qpad zero was product bookkeeping. The new inlet source constructs the missing rows from actual analytic kernels; it never turns that zero into an extra derivative.

The six complete moment histories are retained on the original fixed integral domain:

```text
H=2*integral_0^1 t*Phi(rho*t,Z) dt
M=integral_0^1 U(rho*t,Z) dt
K=2*integral_0^1 t*Phi(rho*t,Z)*U(rho*t,Z) dt
A=integral_0^1 U(rho*t,Z)^2 dt
B=integral_0^1 t*Phi(rho*t,Z)^2 dt
C=integral_0^1 Phi(rho*t,Z)^2 dt
U=U0+epsilon*Psi, rho=4, R=4*epsilon
```

With `W=sqrt(2R)*F0/Pstar`, the actual generic normalized fields are:

```text
E=W*Phi, V=U/Pstar, m=M/Pstar
h=W*H/2, k=W*K/(2*Pstar)
e=(A-R*F0^2*B)/Pstar^2
p=R*F0^2*C/Pstar^2
pressure=P0+p
```

The inverse Pstar is applied once to physical axial velocity, mean and raw mixed history. The accepted physical primitive identities prove the h/k units. Independent analytic P0 is preserved separately from the pressure increment. Its original C3 handles remain, and new C4/C5/C6 nodes explicitly denote derivatives of that original function. An executable projection binds them to the same seven-row normalized pressure provider, the same fourteen atoms, analytic flatten contribution and source/datum hashes. They are not cloned value nodes, fitted pressure constants or selected enclosure endpoints.

The complete original generic equations retain every signed linear, quadratic, meridional and pressure term. The final recovery is:

```text
p2=R*(axial_linear+Pstar*axial_quadratic)/E
```

There is no extra scalar half on Pstar. The positive physical radius is applied once, after the unscaled quotient recurrence. Source-function jets carry ordinary derivatives and exact binomial factors. The numerical route converts ordinary rows to Taylor coefficients once, then executes the unchanged `recover_inputs` source equations with seven input rows. Its output is relabeled explicitly as axial5; old downstream producers are preserved.

Directed numeric views cover the whole continuous Z interval[-1,1] and the exact off-center point3/8. They use the actual degree24/depth6 core reconstruction and all nonlinear density-product tails, within one live context, five formal log bases and directed ledger. These are enclosures of the original inlet functions, not resolved point values or global mixed-norm accuracy certificates.

Focused validation checks actual source bindings, physical units, fixed-domain moment derivatives, nonzero sixth mean/energy/pressure inputs and all generic output derivatives against independent exact polynomial calculus. It also executes actual whole-Z/off-center source production and all higher P0 projections. No historical pressure grid or ancestor suite is regenerated.

Artifacts:

- `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source.py`
- Same stem `.json.gz`, `_check.py`, `_check.json`.
- `scripts/run_original_core_inlet_axial6_generic_source.py`.

Completed within this scope:

- [x] Preserve the admitted nonlinear common core, original C5 handles and complete sixth source derivatives.
- [x] Construct all six genuine C6 core moment histories and their rho4 inlet functions.
- [x] Bind P0 derivatives4..6 to the original independent function and executable order6 provider.
- [x] Prove physical E/V/m/h/k/e/p normalizations and once-applied inverse Pstar.
- [x] Construct genuine Q, all four signed inertial sectors and p2 through ordinary Z5.
- [x] Execute seven-input/six-output actual numerical production at whole Z and point3/8 with every density tail retained.
- [x] Preserve existing downstream arrays and keep all17/global/recursion gates open.

Next tasks, in dependency order:

1. [x] Propagate prescribed background Phi/U/mean and six moment C6 source functions from rho4 through both real micro windows and macro bridge to R100. The linked continuation differentiates original equations and retains histories; its quantitative containment and modified branch remain open.
2. [ ] Implement an actual versioned axial6 source on the Rm patch. Preserve original pressure datum, family, coordinate map, context, five log bases and ledger. Require seven rows for every E/V/m/h/k/e/p/P0 input before differentiating.
3. [ ] Recover Rm Q and all signed inertial/p2 axial5 outputs using the unchanged generic equations. Export six real entries under new keys. Preserve axial4 exports and the internal Qpad convention separately.
4. [ ] Carry the new actual sources into switch, reference-restore and long-reshape finite-N producers. Retain both linear and quadratic sectors, actual Pstar and one physical radius. Apply fifth-row factorial120 exactly once.
5. [ ] Extend the actual active loop/inverse/A/B to C5 with original flat support predicates and positive inverse Jacobian. Extend finite-N density and all five-rate Duhamel histories with every incoming memory.
6. [ ] Construct normalized complete n0 C4/C5 target functions with actual amplitude derivatives and correlated divided rows. Supply source-derived D4/D5 bounds; scalar matrix caps alone do not define target functions.
7. [ ] Extend the original compact implicit-control branch with the same Jacobian and selected N:

```text
J=B+DQ(h*)/N
J*h4=-d4-[4*D2Q(h1,h3)+3*D2Q(h2,h2)]/N
J*h5=-d5-[5*D2Q(h1,h4)+10*D2Q(h2,h3)]/N
```

8. [ ] Build actual active/compact leading U/mean C5, Q/radial Q C4 and full pre-Rp signed pressure C3. Keep leading weight1-delta separate from n1 recovery weight1+delta.
9. [ ] Compose complete five C3 defect functions with the existing core/lower/pulse layers. Differentiate the same physical reserved inverse, including actual Astar derivatives and finite-mu divided differences.
10. [ ] Recover repaired own mean_ZZ and Q/V Z1, required mixed/radial C3/C5 containment and all global joins. Preserve the already proved five terminal C1 zeros.
11. [ ] Finish compatible preheat/finite-energy/exact heat matching and admissible stress/flat remainder. Then implement true n>=2 coefficient equations, independent repair, finite remainders and smooth sum.
12. [ ] Construct oscillatory averaged-stress cancellation and independently validate corrected Cartesian NS with restricted forcing and both Linf/volume-L2 residual below1e-3.

The [preceding analytic core/lower C3 milestone](CURRENT_ORIGINAL_N1_ANALYTIC_CORE_LOWER_PRESSURE_C3_2026_10_13.md) and pulse C3/terminal C1 results remain valid. Genuine global coefficient recursion has not yet been completed.
