# Current variable O3 transition tensor (2026-10-06)

Implementation and scoped whole variable O3 transition tensor/power-join receipt: commit [b7fa407a](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/b7fa407a0acae14be5154fa5f052ad9ab29e87f2).

The actual O3 slope-mu transition now has a complete physical background stress tensor and momentum decomposition on the whole original offset [0,1]. Its offset 1 boundary is attached to the checked actual O3 power phase 0 by common source functions. The admitted tensor chain has **17 actual regions, 16 adjacent tensor joins and 4 internal end-support tensor traces**. The separate velocity/absolute-pressure atlas remains 14 adjacent / 8 internal.

Global cone/lift/NS, independent temporal flatness, prescribed-domain kinetic energy, resolved point values, axis regularity and actual n-dependent temporal recursion remain open. These regional counts are not an overall completion percentage. Full Gamma retains its regional exact NS identity.

## Entry points

- Producer: `experiments/root_st073/lei_ren_part1_paper_compliant_current_O3_transition_background_tensor.py`.
- Complete data: matching `.json.gz`, deterministic complete unpruned UTF-8 JSON; `read_producer()` reads it.
- Checker and receipt: matching `_check.py` and `_check.json`.
- Generic raw pre operator: `lei_ren_part1_paper_compliant_current_pre_pulse_stress_operator.py`.
- Focused controller stage: `currento3transitiontensor`.
- Owner: `CurrentO3TransitionBackgroundTensor(entrance_incoming_tensor=checked_current_entrance_incoming_tensor)`.

Reuse the same checked entrance/main/gap/end/flatten/heat graph, actual pre source, current selected C5 forcing, complete C5 future, full energy and current closed pressure. No historical pipeline restart is needed. Original kernels define functions; finite partitions only enclose them.

## Full variable source and raw units

The producer calls the actual `physical.pre.slope_mu(Z,offset)`. It preserves all five inherited raw histories, the original fourteen-atom analytic preheat pressure P0, full source ordinary derivatives and actual BASE spatial4/time1 source packet. The original log-amplitude jet is

```
D_y log u = -1/2 - mu*sigma(offset)
D_y^(j+1) log u = -mu*D_offset^j sigma(offset), j=1,2,3
u = Utheta/Pstar; V = Uz; y = log R
m = Mz/R
h = Mtheta/(sqrt(2)*R^(3/2)*Pstar)
k = Mtheta_z/(sqrt(2)*R^(3/2)*Pstar)
e = Mztheta/(R*Pstar^2)
P = raw_p + actual_P0
```

Bell derivatives reconstruct the full variable ordinary u jet through four. Only exact R powers and constant Pstar are factored from stress sectors. No constant pure-power velocity rate is substituted in the transition, and no source amplitude, absolute radius or microscopic cap is inverted in this stress adapter.

The raw five-history equations retain m_y=V-m, h_y=u-3h/2, k_y=uV-3k/2, e_y=V^2/Pstar^2-u^2/2-e and p_y=u^2/2. Local V is exactly zero in O3; inherited m/k/e and radial velocity are not deleted. Stress contains five theta and six axial sectors, including variable shear, retained moments, meridional transport and actual absolute pressure.

The radial jet is `shifted_rows(Q,.5,4)`: it represents D_y^j Ur divided by the **current-basepoint** sqrt(R/2), as in original `physical_mixed`. This does not multiply a numerical sqrt(R) value twice. Mode powers supply source values; ordinary coefficient jets already include the derivatives of those factors. For Q=q0*exp(-y), the first radial row is -Q/2, matching Ur proportional to exp(-y/2).

## Full physical transfer and source attachment

The raw stress operator is replayed against original equations (3.16)-(3.18) with arbitrary variable source functions and all moments: **20 stress mixed3 identities**. The velocity adapter is replayed against actual `physical_mixed`: **45 velocity mixed4 identities**. The actual remainder coefficients are compared with the original full remainder operator: **36 remainder mixed2 identities**, for a total of 81 velocity/remainder identities.

The full physical mapper keeps the original stress3, divergence2, completed diagonal2, Cartesian tensor/divergence, remainder2 and residual=-div(T)+E operators. Its single AST adaptation changes the remainder amplitude modes: radial has no Pstar, swirl has one Pstar, and local axial velocity has none. Radius powers, half normalizations, physical beta rates and coefficient formulas are unchanged. The checked original arbitrary-source full physical operator theorem supplies the NS decomposition and completed zero radial tensor divergence; it does not prove a global flat remainder.

At offset 1, original sigma is exactly one and its first four derivatives vanish. Original `pre.power` takes `pre.slope_mu(Z,1)` as its parent, and its phase 0 history kernels/decays are neutral. Thus u, all five histories, analytic P0 and the full source jets agree. Original BASE radius gives the same Rw on both sides. The checked actual power unit/formula theorem converts the right normalized layout to these same physical sources. The full original differential operators then give common tensor/decomposition traces.

Both sides contribute directed triangle bounds to **71 common full component functions**; 63 are nonzero in the saved join. Numerical interval overlap is not an identity proof. Nonzero boundary histories and pressure are preserved.

## Evidence and limits

- Six whole/endpoint/fresh views: **2088 physical contributions, 1512 nonzero source enclosures**; 348 contributions per view.
- Whole source offset [0,1], whole Z source box, exact endpoint flat jets and interior nonconstant log-amplitude derivatives pass.
- Actual spatial4/time1 source grid has 35 Cartesian multiindices, each with ux/uy/uz/p.
- Retained radial remainder is nonzero; the local V=0 axial viscosity is structurally zero.
- The 71-component power attachment passes both saved and fresh compact axial/time/viscosity sectors.
- Focused producer/checker/controller uses one checked parent; changed Python modules compile; staged whitespace and **691 working/index source hashes** pass.
- Read-only review uses **GPT-5.6 Luna / max**. The radial-normalization convention was explicitly resolved against original physical_mixed; root retains acceptance responsibility.

Scope is R>0, |Z|<1, tau>0, constant nu>0 and compact finite logtau sectors. Z=+/-1 boxes are source limits, not an axis theorem. Broad source bounds are not resolved physical point values, cone margins, total-energy certificates or temporal flatness. All corresponding acceptance gates remain false.

## Next coupled construction

- [x] F57C5a-O3-slope-mu: actual whole variable source tensor, full raw histories/absolute pressure and physical decomposition.
- [x] F57C5b-slope-mu-power: original offset1/phase0 full tensor function attachment before common bounds.
- [ ] **Next: F57C5a-O2-buffer + buffer/O3 join.** Consume this checked owner and actual `physical.pre.axial(Z,buffer_offset=offset)` on offset [0,11]. Reuse the raw operator with the original V ordinary jets, five histories and pressure. Preserve actual_y and BASE radius. Prove buffer offset11 / O3 offset0 common source and completed tensor traces. Export whole/endpoint/fresh sectors, directed bounds, receipt and focused controller; commit/push.
- [ ] F57C5a-O2-axial + turnoff/buffer join: construct the actual turnoff phase [0,1], using original ordinary logR derivatives of sigma(1-log(y)/Md), not selector-phase derivatives. Keep local V, radial recovery, all nonlinear transport, energy V^2 source, actual turnoff kernels and positive omitted tails. Prove phase1/offset0 attachment.
- [ ] F57C5a-O2-slope + slope/axial join: actual slope y[0,1], full sigma log-amplitude derivatives and unchanged slope integrals. Preserve five histories and actual pressure. Prove y1/turnoff phase0 full tensor attachment.
- [ ] F57C5a-reference + Rh/slope boundaries: actual repaired Rh reference extension; use the same five terminal moment functions and pressure datum. Enumerate original reference coordinate domain before construction. Connect Rh and reference/slope without replacing repaired moments by a template.
- [ ] F57C5a-core-retained: enumerate same fixed-point core, bridge, first/retained switch, long reshape, restore and actual moment-patch sectors. Construct full tensor/decomposition from live implicit source and analytic pressure. Prove every covered tensor attachment to the repaired Rh trace.
- [ ] F57C4e-axis: establish physical axis regularity and finite limiting tensor/remainder bounds for that same field. Source Z endpoint bounds alone are insufficient.
- [ ] F57C5b-angular-internal: upgrade four actual angular support endpoint/local-difference identities to full completed tensor traces, retaining inherited A/E/P and quadratic D/F contributions.
- [ ] F57C6a-global: enumerate complete prescribed chart/interface coverage; compose actual background tensors and full residual decomposition only on fully covered domains.
- [ ] F57C6b/F57E: independent physical temporal derivatives, high-order/flat remainder decay and tails. Spatial cutoff flatness and self-similar coordinates do not prove this.
- [ ] F57D/F57E/F57F: nonzero source-defined physical u,v,w,p, cone margins/admissible lift and prescribed-domain kinetic-energy integral from the same field.
- [ ] F58a-c: distinct n=1 / n>=2 recovery equations, shared inner domain, per-order five-moment repairs, finite-order residual and smooth summation.
- [ ] F59/F60/F61: both oscillatory families and mean corrections, averaged quadratic stress cancellation, independent corrected Cartesian NS, measured contraction/slenderness and material winding.

Build the next coupled region and boundary, reuse unchanged checked prerequisites, mark only finished scope done, update handoffs and commit/push. Keep the full long-term goal active.
