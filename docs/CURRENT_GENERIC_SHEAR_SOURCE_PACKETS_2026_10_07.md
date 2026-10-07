# Current original source packets and factored recovery, core through O2

Checked implementation: [07c3ca61](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/07c3ca615bd5d88a268604f1c4b7e6934c5ff731). **LEFT4c1-interface and factored recovery attachment are implemented.** Sixteen original chart covers now share a typed S=Pstar packet, including ordinary log-radius velocity and all five history rows, the separate original pressure datum, actual radius/selector provenance and unpruned factor bases. The generic recovery backend consumes these modes directly. Whole-current generic-loop input gates, actual changed field/history transport, terminal repair, new common N and scale recursion remain open. The long-term goal is active.

## Executable interface

`CurrentSourcePackets.saved(chart, view=None)` returns a `CurrentSourcePacket` for an explicitly saved box. Its default charts are core, bridge_first, bridge_second, bridge_macro, switch_first, switch_second, switch_power, reshape, inner_reference, axial_restore, restore_buffer, actual_patch, Rh_reference, O2_slope, O2_axial and O2_buffer.

`CurrentSourcePackets(owners={group: existing_checked_owner}).query(chart, Z, coordinate, ...)` is the arbitrary-coordinate route. It requires the exact admitted owner class, family/source/datum, checked receipt and live graph assertion. It never constructs an ancestor chain or falls back to a saved box. Successful live-owner queries were not exercised by this focused receipt; no live owner was constructed for the checks.

`packet.record()` exports the typed source cover. `packet.recover_original(delta)` recovers ordinary velocity/pressure/primitive rows through four and full signed stress rows through three from the original own histories. `FactoredRecoveryState(packet, original=..., defect=...).field_rows(...)` also supports supplied changed history covers in the same factor algebra. This is a recovery operator, not a newly admitted current changed-history transport or loop.

Normalized rows and recovered fields are in `experiments/root_st073/lei_ren_part1_paper_compliant_current_generic_shear_source_packets_views.json.gz`. The manifest and focused receipt bind its exact bytes and all source dependencies.

## Units and preserved source structure

    E=Utheta/S, V=Uz/S, Q=Ur/(S*sqrt(R/2)),
    m=Mz/(R*S), h=Mtheta/(sqrt(2)*R^1.5*S),
    k=Mtheta_z/(sqrt(2)*R^1.5*S^2), e=Mztheta/(R*S^2), p=Mp/S^2.

Native raw pre units already normalize theta, h, e, p and absolute pressure. Only native Uz, radial coefficient, m and k receive one S^-1 shift. The four original fixed factor logs are `(log hb, 2 log Pstar, 2 log F0base, 2 log swirl_unit_base)`, so division by S is the exact exponent shift `(0,-1/2,0,0)`. Neither S nor microscopic hb or absolute radius is materialized by production recovery.

The bridge, micro and core caches retain both the converted `current_unresolved_raw_source_rows` and `fixed_current_factored_source_log_bases`. The serialization removes the old algebra pointer, not the modes or coefficient data. The interface reconstructs a local FactoredAlgebra from the four logs and exact interval endpoint MP tuples. Half exponents remain exact. These objects are covers of the stored source box; they are not callable point source functions.

Checked raw adapters already applied microscopic hb^-j ordinary-y conversion and the original swirl factor. The interface rehydrates those converted rows directly, with neither conversion repeated. Bridge macro uses the original ordinary coordinate; its fraction only chooses coverage. The patch's D_y=x D_x and the core's Euler/Stirling conversion are likewise already present. Physical radial rows contain the one +1/2 derivative shift. Axial coefficients use Taylor factorial conventions, while radial rows use ordinary derivative conventions.

P0/S^2 is read from the original source's separately saved pressure coefficients. It is never reconstructed by subtracting two broad pressure intervals, divided by S twice, reset at a seam, or exterior-normalized before terminal repair.

## Recovery and receipt scope

The factored recovery extension compiles the already checked generic `field` and `field_rows` ASTs unchanged. Only finite-jet validation and the axial derivative dispatch accept FactoredJet. Its full signed linear/quadratic/shear modes, own-history recurrence, divergence recovery and absolute pressure formulas therefore retain the previous backend's mathematics. All factors stay symbolic through the products and derivatives.

Focused checks cover 16 actual saved charts, 9 exact unit identities, 9586 original velocity/history/pressure source coefficients and 543 coefficients for a nonzero changed-history fixture, plus 8 invalid input guards. Factored caches compare modes before resolution. Older plain energy caches already combine V^2/S^2 into their coefficient covers, so their range comparison uses total intervals after reciprocal unit normalization. These range comparisons are sanity checks; the bound source programs and receipts supply the function identities. No source factor is chosen as a point value. Working/index audit matched 955 bound hashes. Read-only review: GPT-5.6 Luna / max.

The default core cover has rho in [0.001,4], not the axis or a uniform rho->0 estimate. Saved boxes do not prove an arbitrary-coordinate function, a cross-view seam identity or complete original relaxed admissibility. The separate newer core-first interface/strict-collar receipts remain the source of that admitted join; the older three-bridge receipt intentionally has a smaller scope.

## Next production tasks

- [x] **LEFT4c1-interface:** typed actual source covers for core-to-O2 charts; common units, derivative coordinates, original radius tree, P0 and all five histories retained. Saved-box and injected-owner APIs have distinct contracts.
- [x] **LEFT4c3-factored-recovery:** attach unpruned factored rows to the exact own-field recovery backend, including nonzero supplied changed histories and full signed stress derivatives.
- [ ] **LEFT4c1-live-use:** when an existing checked owner graph is present, run selected query calls and bind the actual returned source/function domains. Do not rebuild the entire graph merely to test cache wrappers. Cache boxes are not a substitute for these function domains.
- [ ] **LEFT4c1-gates-a:** derive S/F=(-a,b) and I/F from the full signed current source sectors in a factored/log representation, retaining nonlinear meridional/energy/pressure terms. Do not use T/F for the generic loop input.
- [ ] **LEFT4c1-gates-b:** prove whole-source a>0, the original relaxed D and Q conditions and the equivalent H(t0)-2 input margin across every bridge, micro, power, reshape, restore, patch and outer chart. Source boxes may be subdivided for bounds; source functions cannot be replaced by endpoint fixtures.
- [ ] **LEFT4c1-gates-c:** collect actual lower a, interior margin, strict right-edge excess and p1/p2/t0 bounds into the generic conservative scales. Bind the same actual family/source/datum and the original exact width logs.
- [ ] **LEFT4c1-gates-d:** reserve the true power interval for the new terminal repair and derive a full-radius support/quiet-collar plan from the existing strict inner collar and strict right edge. Do not multiply the loop by an independent spatial taper where that destroys v>2.
- [ ] **LEFT4c2-current-phase:** install one shared phase N*log(R/r_minus) in the actual source functions across chart selectors, with the original radius tree and exact microscopic hb correlation retained.
- [ ] **LEFT4c2-current-derivatives:** construct current generic q/r/phase inversion and zero-mean primitives in a factored/log jet engine. Produce actual whole-cell E,V and increments through mixed4, with phase-held slow derivatives. Numerical point-loop examples are insufficient.
- [ ] **LEFT4c3-current-transport:** extend the Duhamel step to exact formal positive widths and actual changed source covers, then propagate all five defects from the unchanged inlet. Carry defects through quiet gaps and source seams, preserving P0. Per-chart algebra bases require explicit rebase identities, not pointer substitution.
- [ ] **LEFT4c3-current-tensor:** export the own recovered pressure/Ur/stress into the physical tensor mapper, including all nonzero signed sectors, Cartesian basis derivatives and time derivatives. Preserve the old source/function seam proofs until the changed source has replacements.
- [ ] **LEFT4d-repair-map:** reconstruct the new loop family's five-bump functional Jacobian, inverse, uniqueness and nonlinear error bounds; recover all five terminal moments as functions of Z with the same pressure datum.
- [ ] **LEFT4d-new-N:** derive loop/repair derivative and cone perturbation bounds from actual source factors, narrow widths, strict margins and repair powers. Select a new full finite N. Historical scoped N>=68,533,403 cannot be reused without these bounds.
- [ ] **LEFT4e / CONT / ENERGY:** close the whole modified stress cone and remaining physical interfaces; bound energy on the actual required spatial/time domain including the heat exterior.
- [ ] **REC / WAVE / PHYS:** implement the true n-dependent recovery equations and independent repairs, smooth coefficient sum, oscillatory/mean stress cancellation, corrected physical uvw/p/f and full Cartesian NS/dynamical diagnostics.

`current_original_upstream_common_unit_packet_interface_implemented=true` and `generic_recovery_accepts_factored_current_source_packets=true`. `whole_upstream_current_source_packet_assembled` remains false: original complete relaxed input and arbitrary-function assembly have not yet been admitted. The current generic loop, changed whole-family transport, repair/common N, global cone, coefficient recursion and corrected NS remain false. Registry counts and historical scoped N are unchanged.
