# Successor: original Ra microscopic source functions implemented (2026-10-09)

[Current original micro functions](CURRENT_ORIGINAL_MICRO_FUNCTIONS_2026_10_09.md) provides the first conditional RM-W4a.2.2a/b source baseline at Z=0,.5. It retains actual core histories and both original coupled controls. The true finite-N micro-exit correction is still missing; typed source-join/Section11 input and both micro density integrations remain required. Full reconstruction **ACTIVE / INCOMPLETE**. Historical macro evidence below is unchanged.

---

# Current original frozen-macro finite-N source

Checked source [ee2d1a3b](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/ee2d1a3bdf7a2288ec9974ef259bb5751f713a2d). Full reconstruction **ACTIVE / INCOMPLETE**. Predecessor: [current original R100-to-R110 local drivers](CURRENT_ORIGINAL_SWITCH_FINITE_N_2026_10_09.md).

**RM-W4a.2.4 has its first complete conditional local implementation.** The original frozen-macro whole-cell background now supplies all five nonlinear finite-N density C0/Z drivers at Z=0,.5 / N257. Its local correction is composed with the accepted R100-to-Rm local source operator in one live algebra. The complete local-source operator now begins at the actual micro exit R0=Ra*exp(2hb), rather than R100.

The real finite-N correction at R0 is still **unsupplied**. Actual macro background inlets are retained, but they are not correction inlets. The two earlier bridge charts and their true Section 11 input remain required. Five-moment closure, whole-axis functions, stress cone/global N, heat/flat and actual n-dependent recursion remain OPEN.

## True whole-cell geometry

The original geometry is R0=Ra*exp(2hb), Y=log(100/Ra), L=Y-2hb>0. A whole exact rational rho cell uses

    R=Ra*exp(Y*rho)*exp(2hb*(1-rho))
    log(R/R0)=L*rho; dy=L*d_rho.

The accepted current MacroFlow source bases, ledger, coupled phi/V modes, micro angular/axial prefixes and six actual background inlet histories are reused. A shallow algebra view changes only the geometry provider so the original analytic field and Volterra expressions accept whole cells. This is not a hull of endpoints. The source radius/root radius and physical measure retain the nonzero formal hb. An ordinary outward S=L*rho cover is used only inside polynomial weights; it is never a selected hb or replacement radius. At rho=1, physical R=100 and sqrt(R)=10 exactly.

The original signed exponential-polynomial terms are integrated before enclosure. Full nonlinear exponential, field-product and finite Volterra remainder errors remain. Every six-history row has both its genuine incoming decay and its complete source integral. Histories are not reset to instantaneous forcing or comparison moments.

## Original macro source and common recovery

The frozen comparison modes are Dbar=sum_j d_j*R0^j*R^(1-j). Hydro/pressure drives have radial power1; swirl has power2. Full ordinary Z rows and signs remain. Source controls give

    a=hb*Dbar; phi_y=-a*phi/2
    V_y=-hb*(phi/barphi)*(G_hydro+Pstar^2*G_pressure+F0^2*G_swirl).

The same core phi and source quotient determine 1/barphi exactly once. Original source positivity of Dbar and E is checked from their current signed C0 rows. No earlier cone admission is borrowed. The macro hb factor and hydro/pressure/swirl scale assignment are AST-bound to the original coupled field formula and bridge chi=h control.

The accepted current physical-to-common converter retains sqrt(R), F0/F0^2 axial derivatives, separate P0 and all cumulative moment terms. It recovers the full signed pressure/meridional/inertial sectors before one actual R factor. Macro b=2V_common_y/E is generally nonzero; the original variable-a general q C0/Z source-domain union and all-signed-u primitive theorem are reused. All five exact nonlinear density kernels and Z terms remain. No N1024 report or old physical owner construction is used.

## Actual phase origin and retained unknown boundary

The source phase starts at the positive interior Section 11 point r_minus=Ra*exp(hb*s_c/2), not at Ra and not at an assumed s_c=1. The actual selected microscopic s_c, source family, implicit source, pressure datum and exact bridge-width tuple are bound to the checked current inner collar. Its positive log offset is retained separately rather than added to the enormous logRa:

    phase=frac(N*((Y-2hb)*rho+hb*(2-s_c/2))).

The original phase and r_minus assignments are AST-bound. Full-period phase covers are used, with exact phase_Z=0; numerical fast-phase cancellation is not claimed. The same selected s_c/collar hashes accompany every source packet.

Each macro cell uses its true positive Duhamel mass, own-rate incoming decay and downstream suffix, with rates (m,h,k,e,p)=(1,3/2,3/2,1,0). Bounded large-width masses are evaluated stably while their nonzero formal tail factors remain. Pressure has exact incoming memory1 and mass=L*cell_width.

    deltaH_j(R100)=exp(-lambda_j*L)*deltaH_j(R0)+macro_local_driver_j
    deltaH_j(Rm)=memory_j(R0,Rm)*deltaH_j(R0)+all_local_drivers_j(R0,Rm).

Accepted downstream rows are hash-bound and restored into the same live context/bases/ledger. Family, frame, candidate N and canonical P0 coefficient/formal-scale/exact-zero tuples must match. The centered magnitude rows enclose the signed integrals; they are not selected signed values. They do not provide a small defect, convergence, functional terminal repair or global frequency admission.

The original Section 11 zero initial condition is declared at r_minus in the earlier generic source contract. A zero background bridge increment at Ra does not prove a zero finite-N defect, and a zero local source does not remove incoming memory. The present result leaves deltaH(R0) unknown until the actual current r_minus input and both earlier bridge source windows are supplied and integrated.

## Scoped evidence

Artifacts: experiments/root_st073/lei_ren_part1_paper_compliant_current_original_macro_finite_N.py, .json.gz, _check.py, _check.json. Final producer 115.766s; checker 121.375s; both terminal exit0. **1084 staged exact dependency hashes PASS**. Only the new changed-scope source/checker and staged audit ran; no ancestor producer/checker was rerun. One reused **GPT-5.6 Luna / max** worker performed read-only source review; root owns edits/compute/Git.

Independent evidence: 183 exact nonlinear macro field/history/conversion coefficients, 30 true-measure weight comparisons and 3 nonzero axial-drive cases. The fixture has nonzero micro prefix and actual incoming histories, genuine F0 axial derivatives and separate P0. Interior source values and the exact R100 endpoint are checked against independently reconstructed full exponential solutions, not source polynomial truncations. Symbolic identities prove radius, dy, the selected-s_c phase origin, original phi/V_y equations and affine memory composition.

Live replay covers 8 whole source cells, 336 common-unit coefficients, 80 density rows and 10 retained macro memories. 8 invalid/mismatched inputs are rejected. Full goal remains active.

## Detailed next production tasks

- [x] **RM-W4a.2.4 FROZEN MACRO LOCAL DRIVER:** whole-cell actual phi/V and six histories, original Dbar/drives, full current recovery, nonlinear N257 C0/Z densities, true physical weights and composition through Rm. Conditional0,.5 only.
- [ ] **NEXT RM-W4a.2.2a FIRST BRIDGE FUNCTIONS:** actual R=Ra*exp(hb*s), s in[0,1], original chi=(1-sigma)+hb*sigma multiplying both angular and axial source equations. Reconstruct same-source actual phi/V and all cumulative six histories on whole cells; include the micro nonlinear integral errors and comparison derivatives. Preserve current pressure/amplitude data and dy=hb*ds. Do not use the later R100-switch owners here.
- [ ] **RM-W4a.2.2b SECOND BRIDGE FUNCTIONS:** R=Ra*exp(hb*s), s in[1,2], ending at R0. Use original shifted coordinate/cutoff and chi=hb, continuing actual first-exit fields and histories. Prove exact endpoint/source identity with the macro inlet; no independent cap-selected fields.
- [ ] **RM-W4a.2.1 TRUE r_minus INPUT:** bind the paper Section 11 declared zero initial condition at r_minus to this current source owner, candidate N, exact positive hb*s_c/2 offset, analytic P0 and defining modification support. Distinguish it from Ra seam-zero and old N1024 numerical vectors. If that source-owner identity is not established, retain an explicit unknown input.
- [ ] **RM-W4a.2.3a FIRST BRIDGE FINITE-N DRIVER:** begin at the actual s_c/2 input, recover complete variable a/b/full inertia, q C0/Z and original all-u primitives on every current whole source cell, then integrate all five nonlinear density terms with true microscopic own-rate masses. Retain fast phase origin and nonzero source errors.
- [ ] **RM-W4a.2.3b SECOND BRIDGE FINITE-N DRIVER:** propagate the true first-exit correction, integrate the second source window with actual shifted radius and hb measure, and export a typed complete R0 C0/Z correction packet.
- [ ] **RM-W4a.2.5 REAL R100 VECTOR:** apply the genuine source-bound R0 packet to this macro local operator. Bind family/source/datum/P0/bases/ledger/N and all preceding windows; reject arbitrary zero, missing-window and earlier-N substitutions.
- [ ] **RM-W4a.3 / W4b.4 / W4c.5 REAL R110 / Rm:** apply the true prefix to the complete current local R0-to-Rm operator and current Rm atlas. Local source completeness cannot admit the missing input.
- [ ] **RM-W5 SHARPENING:** preserve source a/q/u/radius correlations and original phase averaging through each complete window; compare improved covers with these genuine baselines. Bound changes alone are not measured defects/convergence.
- [ ] **RM-W6..W11:** complete Rc/all24 drivers and terminal functions; whole-Z/high jets; functional five moments, cone and one global N; pressure/heat/energy; admissible stress/flat and actual n-dependent recursion; original pulses/full corrected NS and quantitative dynamics.

Mark DONE with implementation, scoped report/receipt and commit. Prioritize the two bridge source windows and their true input over ancestor reruns.
