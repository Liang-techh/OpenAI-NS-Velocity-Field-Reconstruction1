# Actual Rm whole terminal five-density C0/Z Duhamel integrals

Checked source [a5057a30](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/a5057a30f06a0d31f8f71c135d37dacf3842e7ff). Full reconstruction **ACTIVE / INCOMPLETE**. Predecessor: [source-bound Rm phase and five Z densities](CURRENT_ORIGINAL_RM_PHASE_DENSITIES_2026_10_09.md). The new actual Rm source now has signed local five-density C0/Z integral enclosures over the full terminal window x in [71/40,e], at the two conditional axial frames0,.5 and candidate N=257. Four closed radial cells per frame are transported to the same exact Rh=e. This is a source-function enclosure and affine inlet operator, not an evaluated complete finite-N prefix, corrected Rh history, five-terminal closure, global N or true n-recursion.

## Original common units and defining five rates

The density interface already supplies E=Utheta/Pstar, V=Uz/Pstar and all signed increment densities in the original common units. No additional R, Pstar or N factor is applied here. The original own-rate rows are

    lambda_m=1; lambda_h=lambda_k=3/2; lambda_e=1; lambda_p=0
    m=deltaV; h=deltaE
    k=V*deltaE+E*deltaV+deltaE*deltaV
    e=2*V*deltaV+deltaV^2-E*deltaE-deltaE^2/2
    p=E*deltaE+deltaE^2/2.

AST bindings retain the defining signed returns in GenericMomentRecovery.increment_densities and the native signed density kernel, including both quadratic terms and every original nonzero V cross term. Z derivatives come directly from the accepted Rm density adapter. The separate patch repair-tail rates (2,1.6,1.2,.2 etc.) are not substituted for these five own rates.

Normalized histories have units m=Mz/(R*S), h=Mtheta/(sqrt(2)*R^1.5*S), k=Mtheta_z/(sqrt(2)*R^1.5*S^2), e=Mztheta/(R*S^2), p=Mp/S^2 with the same S=Pstar source definition. Physical patch rows require their original conversion before entering this density operator.

## Whole-source cell integral and downstream suffix

With y=logx, dy=dx/x. For a closed source cell [a,b] ending before target t=e,

    w=log(b/a); suffix=log(t/b)
    mass_lambda(w)=integral_0^w exp(-lambda*s) ds
    I_j_cell_to_target
       = exp(-lambda_j*suffix)
         * integral_a^b exp(-lambda_j*log(b/x))*f_j(x,Z)*dx/x.

The complete actual source/phase union encloses f_j and f_j_Z throughout the cell. Multiplying each signed range by its positive original mass gives a rigorous local integral enclosure; the same operation applies to Z derivatives. A phase union is hulled only in the shared original formal basis; there are no sampled source values or phase midpoints. The source radius, endpoints, weights and phase are Z independent, so differentiating under the finite integral introduces no boundary term.

Mass is computed stably as w*exp_average(-lambda*w), with exact mass=w for pressure rate0. Microscopic positive mass is retained rather than evaluating 1-exp(-lambda*w) by cancellation. Every local contribution receives its downstream suffix before summing. The current partition is [71/40,2,9/4,5/2,e], covering the entire terminal window; its source cells retain actual N-dependent phase and original derivatives.

The actual source family, closed-cell flag, exact cell geometry, live Rm radius factor and each phase cell's exact P0 object are bound. Each actual R equals the same factored Rm*x source in powers/offset/coefficient. Context, five bases and ledger are shared by all density, integral and incoming rows. The huge positive radius is never materialized.

## Incoming history is an explicit argument

For W=1-log(71/40), the correction history satisfies

    deltaH_j(Rh)=exp(-lambda_j*W)*deltaH_j(71/40)+I_j.

`transport_supplied_incoming` requires an explicitly supplied five-vector with C0 and ordinary Z rows in the same source algebra. An absent incoming vector is refused; no zero inlet is assumed. Its result is conditional on that supplied enclosure, not a proof that the preceding finite-N prefix was integrated. The producer deliberately leaves that finite-N inlet unsupplied.

The actual leading inlet and leading Rh histories are retained separately as source memory/comparison data. They are not the incoming finite-N correction. P0 stays the exact same live axial datum object and is not included in the radial p integral. Since lambda_p=0, incoming pressure memory is preserved with exact decay1.

## API and scoped evidence

`OriginalRmTerminalDensityIntegrals().contribution('0')` uses the declared full terminal partition at N=257. A custom strictly ordered rational terminal partition must end at exact Rh. `transport_supplied_incoming(label,incoming)` applies the explicit affine inlet operator. It neither invents incoming corrections nor admits global closure.

Files: experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rm_terminal_density_integrals.py, .json.gz, _check.py, _check.json. Producer 37.984s; checker 39.688s; terminal exit0. **1186 staged exact dependency hashes PASS**. Existing GPT-5.6 Luna/max reviewed rates, signs, units, suffixes and source scope read-only; root owns implementation/compute/Git.

Finite independent affine-in-y source fixtures at two Z values and both signs compare 120 cell integrals, downstream weighted sums and transport with genuinely nonzero incoming C0/Z values against exact high-precision closed formulas. Rate-zero pressure memory and tiny positive mass pass. Invalid partitions, negative transport width, foreign incoming algebra, foreign source family, non-cell cover and different P0 object are rejected. Finite fixture values are diagnostic only, not native field definitions.

Native live/report evaluation agrees on 20 whole-terminal C0/Z integral rows, 80 cell-to-Rh rows, 40 positive own-rate masses and 4 inlet/Rh memory packets. There are 8 closed source cells across the two conditional frames. All suffixes, units and P0 remain correct. Unspecified finite-N incoming histories are refused at both frames. Range-integral enclosures need not be sharp or sign-definite; no source/phase-mean cancellation or moment-closure accuracy is inferred from them.

## Detailed next production tasks

- [x] **RM-I1 WHOLE TERMINAL LOCAL C0/Z INTEGRALS.** Integrate all five actual source density functions with true dx/x measure and original own-rate masses; transport each closed cell to Rh before summing. Retain source phase, full signed terms and C0/Z rows.
- [x] **RM-I2 CONDITIONAL INLET TRANSPORT AND P0 MEMORY.** Implement the nonzero incoming affine operator; refuse missing inlet data, retain exact rate-zero pressure memory and preserve separate actual leading inlet/Rh and original P0 objects. Keep full-prefix and corrected-history flags false.
- [ ] **RM-I3 SHARPEN TERMINAL INTEGRAL BOUNDS.** Measure current enclosure widths in their exact formal units. Refine true phase/radius cells, use independently proved signed phase-mean cancellation where available, and retain both original N-dependent levels and peak/slow derivative terms. Full-period range covers are valid but do not prove a small defect. Never set an oscillatory integral to zero by assumption.
- [ ] **RM-I4 ACTUAL ACTIVE SOURCE ATLAS [1,71/40].** Extend the accepted narrow active cell to complete source coverage, partition all original support edges and signed-u/cutoff geometry changes, then apply this same dx/x operator. Use source-correlated denominators/branch proofs; do not extrapolate terminal-flat identities into active controls.
- [ ] **RM-I5 REAL FINITE-N INCOMING PREFIX.** Construct the upstream finite-N correction histories through core, micro, macro, switch, long-reshape, reference and restoration, preserving all prior tails and exact P0. Supply an independently bound correction vector at x=1 and then at71/40; actual leading histories are not this vector.
- [ ] **RM-I6 FULL PATCH TO Rc CUMULATIVE FUNCTIONS.** Combine actual active and terminal contributions with the proven prefix and outer windows, maintaining each own-rate suffix. Produce complete Rc_E/Rc_E_Z under the same family, pressure and explicit N, without reset at Rm or Rh. Assemble all24 required source windows and both original levels.
- [ ] **RM-I7 GENUINE MIXED y/Z DERIVATIVES.** Extend the source adapter beyond first Z with genuine fixed-phase y rows and N*A_phi/N*B_phi chains. Differentiate original cutoffs only on proven source branches; retain radial prefactor shifts once.
- [ ] **RM-I8 ACTUAL NEW-OWNER dstar/CONE.** Prove the full correlated signed H0,D,J and quadratic relaxed-cone inequalities with R factored. Selected eta/dstar remain definitions until this owner's margin theorem closes.
- [ ] **RM-I9 FINITE-N UNIQUE REPAIR, WHOLE-Z AND GLOBAL N.** Solve the actual B*h+N*r+Q/N=0 problem, not its leading surrogate. Prove five terminal conditions as functions of Z, include pole/midplane controls and choose one frequency compatible with every source, integral, cone, interface and tail bound.
- [ ] **RM-I10 HEAT/ENERGY/STRESS/FLAT AND REAL n-RECURSION.** Finish analytic preheat pressure/exact heat joins and finite energy; construct admissible stress and flat remainder; implement actual n=1/n>=2 recovery equations and independent moment repairs with divergence-free truncation/summation. Coordinate rescaling alone remains insufficient.
- [ ] **RM-I11 PULSES/CORRECTED NS AND DYNAMICS.** Add mean and oscillatory corrections, averaged stress cancellation and final forced Cartesian NS validation. Measure contraction, relative axial elongation, material winding and finite energy across physical times.

Mark DONE with defining source, scoped report/receipt and commit. Keep the full objective active. Proceed to source/integral production after necessary changed-scope checks.
