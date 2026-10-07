# Original q ordinary slow derivatives and source subdivisions

> Successor: [CURRENT_NATIVE_PHASE_FIRST_JETS_2026_10_07.md](CURRENT_NATIVE_PHASE_FIRST_JETS_2026_10_07.md) ([9be86964](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/9be869643f5379284db002c73a4330878f503bd8)) now implements original implicit phase/A/B first slow and fast derivatives, plus actual spatial first derivative chain on the declared domains and whole local integral cell. Higher primitive jets, density C1 integrals and inlet-to-Rc histories remain open.

Checked implementation: [8fc7c917](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/8fc7c9178d716025d231de36847bd3df12d6b049). **Original q_y, q_yy, q_Z, q_yZ and q_yyZ functions now execute alongside q on15 original native query boxes, the whole active integral radial/Z cell, and two new axial subdivisions.**20 source queries give18 enclosed queries and108 ordinary derivative enclosures. The two full-Z cutoff-crossing queries remain explicitly unresolved. Whole-chart coverage, implicit phase/A/B derivatives and C1 cumulative histories remain open.

## Original formula and derivative arithmetic

q=sigma(1-Delta/eta)*sqrt((2eta-Delta)/(2a)), with the same source-selected positive constant eta. The derivative input consists of the existing six ordinary original a,b,Delta root rows for y<=2,Z<=1. They are ordinary log-radius rows; hb/Pstar conversions are not applied again.

Delta>=eta uses the original lazy flat branch: q and every slow derivative are exactly zero, and the active square root is not evaluated. Delta<=0 uses original sigma=1, with exactly zero cutoff derivatives. Resolved positive transition boxes use the original sigma_jets Taylor coefficients multiplied by n! and the ordinary multivariate chain rule. Unresolved cutoff/normalization boxes return rows=None and a source-subdivision status.

The y2Z1 cutoff derivative is sigma_3*s_y^2*s_Z + sigma_2*(s_yy*s_Z+2*s_y*s_yZ) + sigma_1*s_yyZ, s=1-Delta/eta; sigma_k denotes the kth ordinary sigma derivative. The square-root jet is obtained by differentiating R^2=(2eta-Delta)/(2a) and dividing by the same positive2R. On the active source box, 2eta-Delta>eta; a<=a_max gives the rigorous positive root lower log(eta/(2a_max))/2. Endpoint bounds used here are proof coordinates, not field values.

## A real source-correlation fix

The new adapter forms Delta=(a-2)+b^2/a, cancelling the same source constant before adding the tiny quotient. The previous native C0 arithmetic formed(a+b^2/a)-2. In the original O2 axial/buffer charts a=2 exactly, so that arithmetic could lose the genuine tiny positive b^2/a lower bound. The formula and field are unchanged; the positive factor and cutoff branch are recovered more accurately.

The returned roots, loop, q and nested source record all use the refined expression. Pre-refinement branch/q metadata are explicitly labeled historical. The exact O3 excess from its original mu formula is retained rather than recomputed by subtraction.

Both added source subdivisions use original coordinate.1337 and Z in[.49,.51]: bridge_first and O2_axial now resolve as exact flat q with six zero derivative rows. Their original Z=[-1,1] queries remain unresolved because the broader source box crosses or fails to isolate the cutoff branch. Flat local subdivisions do not admit a full-Z flat field.

## Executed scope and API

- 15 of17 original native query boxes enclosed:13 active and2 O3 flat.
- The actual O2-slope integral cell coordinate[.13369999,.13370001], Z=[.49,.51], is enclosed and has all six derivative rows.
- Two additional Z=[.49,.51] source subdivisions are enclosed and flat.
- Total18 enclosed/20 queries,108 rows;4 flat queries;2 unresolved full-Z queries.

```python
from lei_ren_part1_paper_compliant_current_native_q_slow_jets import NativeQSlowJets

# correlated_owner is the accepted same live NativeCorrelatedShearQ.
backend = NativeQSlowJets(correlated_owner)
result = backend.query('O2_slope', ('.49', '.51'),
                       backend.ctx.mpf(('.13369999', '.13370001')))
q_Z = result['rows'][(0, 1)]
q_yyZ = result['rows'][(2, 1)]
source = result['source']  # refined same-source roots, loop and q
```

For downstream phase construction, use this returned source together with its q slow rows and the same selected d_star. Keep one shared source basis/ledger. Do not recompute a different source box and combine its derivatives with this box's C0 values. An unresolved result requires subdivision before derivative use.

Evidence: fresh same-source replay of20 actual queries and108 ordinary derivative enclosures;30 independent scalar multivariate derivative comparisons using the original scalar cutoff, plus six exact shared-source Delta=eta flat boundary rows; signed-b/transition/flat/Delta=0 cases; tiny positive q retained with genuine zero constant derivatives; regression retaining tiny positive b^2/a and its actual smaller-eta flat branch; nested refined branch/q provenance consistency. 1018 working/index hashes pass. Producer 34.125s; focused checker 34.906s on the already-live seed. Read-only review: **GPT-5.6 Luna / max**; math and provenance findings incorporated.

## Next production tasks

- [x] **LEFT4c2-q-slow-jets on admitted source boxes:** original six ordinary q rows, lazy flat cutoff, transition chain, positive root, tiny Delta=0 retention and explicit unresolved boxes.
- [x] **LEFT4c2-O2-excess-correlation:** exact a-2 cancellation before b^2/a; refined authoritative source record.
- [x] **LEFT4c2-local-open-box-subdivisions:** bridge_first/O2_axial flat source slices at Z=[.49,.51]; full-Z coverage remains open.
- [x] **LEFT4c2-implicit-phase-first-jets on declared domains (see successor; higher orders/whole coverage remain open):** use the same refined source/q jets to recover Phi_y,Phi_Z and implicit angle derivatives at fixed phase. Preserve positive rho, s and h^-1 in the signed Mobius chart, small-r series derivative tails and exact flat branches. Do not divide by a rounded1-r or by a zero r. Recover original A_y,A_Z,B_y,B_Z and higher needed mixed orders.
- [x] **LEFT4c2-actual-fast-phase-first-chain on declared domains (see successor; global common N remains open):** bind phi=N*log(R/r_minus) using the existing actual radius/Jacobians, then combine fast phase and slow profile derivatives exactly once. Radius phase is Z-independent. A free-phase derivative is not automatically an actual spatial derivative.
- [ ] **LEFT4c3-density-C1-and-integrals:** differentiate every original signed density kernel/expm1 factor, derive whole-cell C1 Z covers, and attach them to actual positive-weight integration. Preserve original cross terms and rates1,3/2,3/2,1,0.
- [ ] **LEFT4c2-full-source-subdivision:** cover the two unresolved full-Z boxes, including extremely narrow active/flat thresholds. Use retained log-factored b^2/a and eta; refine original Z/coordinate functions rather than selecting midpoint fields or replacing tiny positive scales by zero.
- [ ] **LEFT4c3-common-cell-basis-and-incoming:** rebase adjacent cells with exact source-factor cancellation, bind original five incoming histories and P0 at actual r_minus, and propagate quiet-region memory unchanged.
- [ ] **LEFT4c3-cumulative-Rc-functions:** all five actual signed histories and C1 Z jets throughRc, with controlled phase/oscillation/integration error and original chart lengths. Current local C0 contributions do not satisfy this requirement.
- [ ] **LEFT4d:** actual Rc target functions/A/A_Z/divided(J-M)/mu, reserved-band inverse, unique controls, five terminal Z-function identities and changed radial velocity/pressure/stress continuation.
- [ ] **HIGH / CONT / OUTER / ENERGY / LEFT4e:** higher derivatives, same-function seams, heat/pressure/energy compatibility, genuine common finite N and global modified admissible cone.
- [ ] **REC / WAVE / PHYS:** n-dependent coefficient recovery with independent repairs/smooth sum, oscillatory quadratic stress cancellation and corrected Cartesian NS/geometry/particle diagnostics.

Scoped gate: `current_original_native_q_ordinary_y2_Z1_slow_jet_functions_executed`. Every global completion gate remains false. Next work is actual implicit phase/A/B slow derivatives, not rerunning unchanged native C0/integral proofs. Full requirement/dependency detail remains in the [local integral handoff](CURRENT_NATIVE_LOCAL_SIGNED_INTEGRALS_2026_10_07.md).
