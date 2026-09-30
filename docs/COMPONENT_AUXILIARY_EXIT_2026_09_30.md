# Updated research objective and component auxiliary exit — 2026-09-30

The user-updated PROJECT_GOAL.md is authoritative. The project now distinguishes seven stages: geometry; leading core; five-moment/heat-exterior matched background; admissible stress; higher-order recursive background; oscillatory stress correction; independent full residual validation. Small background residual is not the leading-profile acceptance gate.

## Implemented in this lane

`ComponentExitComparison` carries the complete-preheat component core into the original Section 9.23 auxiliary comparison. It uses the same angular log-slope, axial ODE, pressure and five-moment equations and their Z jets. All divisions are evaluated in a declared finite Taylor ring in the pressure-tail parameter. Stress and frozen D/E driver coefficients use the same ring; no whole pressure/velocity atom is converted into a scalar before these equations.

The new pressure ring retains MP atoms, reciprocals, integer powers, exp/log and a declared order. This is finite-order pressure-parameter arithmetic, **not** the higher-order n-dependent temporal coefficient recursion. The finite radial core polynomial remains unchanged in meaning.

Independent resolved-tail replay at pressure parameter 0 and 1, in the active transition and frozen extension, agrees with the original scalar comparison for fields, pressure, five moments, stress and D/E drivers within a maximum scaled difference of 5.73825330761e-53 (order 9). This is a bounded implementation check; pressure-parameter and ODE remainders are not enclosed.

## Actual source replay

Command: `python experiments/root_st073/lei_ren_part1_paper_component_exit_comparison_check.py`.

The receipt uses the actual source delta=1e-200, logC=5e151, Lambda=1e36, degree=18, coherent waiting, complete analytic preheat datum and the source exit width exp(-100-1e154). It checks the core boundary, transition, endpoint and frozen auxiliary branch. The nonzero pressure tail and axial response survive all five probes. Zero-parameter scalar replay differs by at most 4.50111984858e-233. Read the generated JSON for retained pressure/velocity/moment atoms and replay details.

The tiny collar width itself lies far below a single atom's relative arithmetic precision. Width-dependent updates can therefore disappear inside a pressure-parameter coefficient. This limitation is flagged explicitly; preserving pressure powers alone does not certify the complete collar or its smoothness. A width-parameter representation or controlled bounds are still needed.

## Remaining paper dependencies

- Propagate the same atom data through the **prescribed** Section 9.25 exit bridge and its driver/normalization tangents; auxiliary comparison completion is not bridge completion.
- Retain tiny width dependence or enclose its effect; enclose parameter truncation and RK integration error.
- Carry component-valued fields through switching, long reshape and inner moment repair, then regenerate angular/axial/heat data and functional terminal identities together.
- Close pressure compatibility, radial finite energy and controlled exact heat exterior before declaring a matched background.
- Establish cone margins and flat remainder before genuine n-dependent recursion. Only after oscillatory averaged flux cancellation evaluate the full 1e-3 residual gate.

No prescribed bridge, full annulus, functional terminal closure, cone, finite energy or temporal recursion completion is claimed.
