# Own generic shear histories, pressure and radial/full stress recovery

Checked implementation: [67b783db](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/67b783db67c8fa9e840474d9cf3643ba77d789ec). **LEFT4c3-backend is implemented:** all five changed histories propagate without resets, and their own pressure, divergence-free radial velocity and full signed stress can be recovered with ordinary log-radius/axial derivative rows. The backend is attached to saved original current O2 source covers. The new generic loop is not yet installed in the complete current source, so completed changed whole-family moments, independent repair and new common N remain open. The long-term goal is active.

## Common physical units and continuous history

GenericMomentRecovery in experiments/root_st073/lei_ren_part1_paper_compliant_current_generic_shear_moment_recovery.py keeps one formal velocity unit S=Pstar. Its variables are

    E=Utheta/S, V=Uz/S,
    m=Mz/(R*S), h=Mtheta/(sqrt(2)*R^1.5*S),
    k=Mtheta_z/(sqrt(2)*R^1.5*S^2), e=Mztheta/(R*S^2), p=Mp/S^2.

The actual ordinary y=log R equations are

    m_y=V-m, h_y=E-3h/2, k_y=E*V-3k/2,
    e_y=V^2-E^2/2-e, p_y=E^2/2.

The source AST binds these equations to the current pre_pulse_mixed_C4 physical_mixed program. Original histories and changed defects are kept separately; own() adds them. Given original E,V and increments dE,dV, the exact changed sources are

    dm source=dV,
    dh source=dE,
    dk source=V*dE+E*dV+dE*dV,
    de source=2V*dV+dV^2-E*dE-dE^2/2,
    dp source=E*dE+dE^2/2.

advance applies the continuous Duhamel formula with decay rates (1,3/2,3/2,1,0). The positive kernel mass is enclosed as width*integral_0^1 exp(-rate*width*t)dt, retaining extremely small positive widths without subtracting an exponential from 1. Inputs must cover the same actual functions and axial derivatives on the entire log-radius cell. Endpoint samples are insufficient.

On a quiet gap with dE=dV=0, the four radius-normalized defects decay with their rates; the pressure defect stays constant. In physical cumulative units those incoming defects persist. No pressure convention, moment vector or axis datum is reset when a local cutoff becomes zero.

## Recovered dependent fields and full derivative rows

The same unchanged P0/S^2 is added to own p. Before terminal repair, this pressure need not vanish at infinity; exterior normalization cannot be imposed prematurely. The radial coefficient is

    Q=Ur/(S*sqrt(R/2))
     =(2Z*V-(1-delta)Z*m-(1-Z^2)*m_Z)/(1-delta*Z^2).

The exact ordinary divergence identity follows from m_y=V-m. Pressure obeys d_y(P/S^2)=E^2/2. Both recoveries use the changed full history, rather than a local increment or copied original radial/pressure field.

field and field_rows recover the complete signed inertial and shear coefficients. Their physical modes are sqrt(R/2)*S for inertial linear terms, sqrt(R/2)*S^2 for nonlinear/energy/pressure terms, and S/sqrt(2R) for shear. Nonzero meridional transport, all mixed/energy histories and absolute pressure are retained. The current raw pre stress program is replayed symbolically to establish these units and signs.

field_rows requires actual ordinary E,V derivative rows through four. It builds own normalized moment derivatives, physical velocity/pressure and primitive derivatives through four, and signed stress rows through three, with axial Taylor jets. The physical radial and inertial +1/2 factors and shear -1/2 factor are differentiated exactly once. P0 contributes only to absolute-pressure row zero. The mixed4 triangle is available from these source rows; Cartesian lifting and full physical time derivatives remain separate downstream steps.

## Current source cache attachment

CurrentO2RecoveryCache loads checked current_O2_background_tensor.json.gz and its receipt, verifies dependency hashes and the common family/source/datum, and converts the four saved whole-chart source covers: Rh_reference, O2_slope, O2_axial and O2_buffer. It does not call predecessor constructors or return arbitrary point values.

Existing raw pre units require dividing only physical Uz, raw m and raw k by S. E,h,e,p and P0 are already in the correct units. The inverse S is represented in arbitrary-precision interval normalization; S itself, absolute current radius and Cstar are not materialized as physical field values. In particular P0 is not divided a second time. E derivative rows are recovered from the original E0 and exact log-Utheta derivative functions.

Full saved covers are in lei_ren_part1_paper_compliant_current_generic_shear_moment_recovery_views.json.gz, bound by the manifest/receipt hashes. Source covers remain marked as covers. This attachment is to the original field; it is not an assertion that the new Section 11 loop already supplies changed current whole-family inputs.

The source-entry scan confirms the common raw units through switch-power, reshape, restore, patch and O2. Bridge/microswitch caches require their checked phase-width/swirl-unit raw adapters before this normalization; their internal unresolved rows must not be relabeled. The newer current_core_first_interface_check and strict collar receipts have the core-first join/source attachment true. The older three-bridge functional receipt deliberately leaves that extra scope false; that older flag is not evidence that the newer core-first proof is absent.

## Evidence and scope

Focused receipt: 16 exact source/recovery identities, 60 independently integrated transport coefficients, 30 quiet memory coefficients, 1 positive thin-cell injection, 30 split continuous-transport coefficients and 10 independently recovered radial/pressure derivative values. The actual current original O2 comparison checks 840 saved velocity/pressure/primitive derivative coefficients, with source-program identities providing the unit/function equality. There are 7 invalid-input guards. Working/index audit matched 950 source/receipt hashes. Read-only review: GPT-5.6 Luna / max; no sign, normalization, pressure or prefactor error found.

The Duhamel interval operator is conditional on true whole-cell source covers. Its source identities and range checks do not prove an installed global generic-loop field, a new admissible N, terminal moment restoration, energy or corrected NS. Registry counts remain unchanged.

## Ordered next production tasks

- [x] **LEFT4a / LEFT4b:** original outer relaxed input and current strict inner support.
- [x] **LEFT4c2-kernel:** executable generic Section 11 curve, phase inversion and zero-mean primitives; see [CURRENT_GENERIC_SHEAR_LOOP_2026_10_07.md](CURRENT_GENERIC_SHEAR_LOOP_2026_10_07.md).
- [x] **LEFT4c3-backend:** own five-history transport, same P0, exact radial/divergence recovery, full signed stress and ordinary/axial derivative rows.
- [x] **LEFT4c1-O2-units:** actual saved original O2 common-unit conversion, all four charts and derivative source covers.
- [ ] **LEFT4c1-source:** implement one typed upstream packet interface. Apply checked bridge/microswitch adapters with their exact widths and swirl factors; normalize switch-power/reshape/restore/patch/O2 common rows. Keep coordinate selector, ordinary y derivatives, cumulative histories and P0 separate. Inject already existing live owners for arbitrary-coordinate queries; avoid default ancestor construction during focused checks.
- [ ] **LEFT4c1-gates:** assemble the complete original relaxed-input proof and whole-box a/H(t0)-2/p bounds, strict right edge and the reserved repair power interval. Use the current core-first attachment while respecting older scoped receipt flags.
- [ ] **LEFT4c2-current:** install the generic loop against actual source functions, one shared N*log(R/r_minus) phase and phase-held slow derivatives. Produce source-bound axial/mixed derivative covers for E,V and increments. The scalar engine's numerical examples are not those covers.
- [ ] **LEFT4c3-current:** feed those covers through this own-history backend from the actual unchanged inlet across every chart/support/gap, then use the recovered pressure/Ur/stress rows in the current tensor exporter. Carry outgoing defects into the reserved repair interval; no local-cutoff reset.
- [ ] **LEFT4d:** rebuild the new family's five-bump map, inverse/uniqueness and nonlinear error bounds; restore all five terminal moments with the same analytic P0. Bound the actual loop/repair norms and cone stability tolerance, including exact narrow source scales, then derive the full new finite-N requirement. Old N>=68,533,403 is not this bound.
- [ ] **LEFT4e / CONT / ENERGY:** global completed cone/common N, remaining physical interfaces and full required-volume energy/heat-tail contribution.
- [ ] **REC / WAVE / PHYS:** true n-dependent coefficient recovery and smooth sum, oscillatory/mean stress cancellation, corrected physical uvw/p/f, full Cartesian NS and measured contraction, relative elongation, scale recursion and material winding.

generic_shear_own_moment_transport_and_recovery_implemented=true and current_original_O2_common_velocity_unit_recovery_attached=true. Modified whole-current loop/moments/repair/N, global cone, coefficient recursion and full corrected NS remain false. The long-term goal stays active.
