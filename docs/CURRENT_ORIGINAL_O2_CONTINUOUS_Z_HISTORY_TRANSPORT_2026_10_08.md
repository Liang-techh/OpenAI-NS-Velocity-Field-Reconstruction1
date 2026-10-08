# Original O2 continuous axial history transport

Checked source: [01cb96ef](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/01cb96ef93663a0c36671dff1e3ac9a8625aeea1). Predecessor: [genuine fixed-Z density derivatives](CURRENT_ORIGINAL_O2_GENUINE_DENSITY_Z_INTEGRALS_2026_10_07.md). Existing read-only reviewer **GPT-5.6 Luna / max**; no new agent was spawned. The latest fetched documentation-main head f8e2b02e reviews earlier centered-phase/source work; it supplies no new original integral implementation to merge into this computation.

The original O2 five value/ordinary-Z contribution enclosures now cover **continuous rectangles**, y in[0,1], Z in[.36,.38] and[-.38,-.36], at candidate N7. Both 256/2048 refinements retain every source piece and genuine radius-phase union. No axial sample defines the function or substitutes for whole-Z coverage.

The source-defined original entrance histories are restored at nonzero Z. Set C=1/(1+Z²), S=Pstar:

| Channel | Original normalized inlet at y0 | Ordinary Z derivative |
|---|---|---|
| m | 4Z/Pstar | 4/Pstar |
| h | 5C/8 | -5ZC²/4 |
| k | 5ZC/(2Pstar) | 5(1-Z²)C²/(2Pstar) |
| e | 16Z²/Pstar²-5C²/12 | 32Z/Pstar²+5ZC³/3 |
| p | 5C²/2 | -10ZC³ |

These are **cumulative histories**, derived from the original slope/reference definition with J(0) and all three defining masses exactly zero. The local velocity densities E,V,EV,V²-E²/2,E²/2 are different functions and must not be substituted for inlet histories. The AST source guard and symbolic ordinary derivative identities bind the original formulas. Pstar inverse factors stay formal and strictly positive before normalized finite export; no microscopic source becomes exact zero.

The original unmodulated densities and their true ordinary derivatives are also integrated under rates (1,3/2,3/2,1,0). For each channel j,

```text
original_j(1,Z) = exp(-rate_j)*original_j(0,Z) + integral_original_j(Z)
own_j(1,Z) = original_j(1,Z) + integral_modulation_j(Z)
             + exp(-rate_j)*incoming_defect_j(0,Z).
```

The ordinary-Z relation uses the same Z-independent masses and decay, with every value replaced by its actual ordinary derivative. Incoming defects and their derivatives are required explicitly. This arithmetic API cannot certify a caller's arbitrary function covers; real upstream source/incoming functions still must be attached. Absolute P0 remains separate, and pressure's rate-zero memory is never reset.

## Executable interface and evidence

Producer: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_O2_continuous_Z_transport.py`; checker: `lei_ren_part1_paper_compliant_current_original_O2_continuous_Z_transport_check.py`; manifest: `lei_ren_part1_paper_compliant_current_original_O2_continuous_Z_transport.json`; accepted receipt: `lei_ren_part1_paper_compliant_current_original_O2_continuous_Z_transport_check.json`. Four separate lossless gzip archives remain below the GitHub 100MiB file limit. Together they are 192,473,558 bytes. The exact directed endpoints/MPF tuples and full source/phase/derivative records remain in those archives.

```python
import lei_ren_part1_paper_compliant_current_original_O2_continuous_Z_transport as source

owner = source.OriginalO2ContinuousZTransport()
record = owner.integrate(2048, Z_lower='.36', Z_upper='.38', N=7, bits=24)
# Supply all five real incoming defect functions and ordinary-Z function covers:
with source.mp.workdps(owner.c.dps + 40):
    own = source.apply_actual_incoming(owner.c, record,
        incoming=actual_incoming, incoming_Z=actual_incoming_Z,
        source_family=owner.family,
        original_P0_datum_sha256=owner.family['datum_enclosure_sha256'])
```

Focused checks cover 4,608 complete source cells, 4,640 actual original inverse/derivative phase-source records and 23,040 independent closed-form positive masses. Original endpoint values/ordinary derivatives independently overlap the accepted defining-mass H/D/P formulas; refinements have common intersections. Positive rectangles also intersect the previous genuine fixed-Z=.37 C0/Z cell results. Explicit nonzero synthetic incoming data test only affine wiring; they are not physical incoming data. Missing channels, changed datum, incorrect rates/window/source identity fail closed. Git-index dependency audit: 1094 hashes.

## Quantitative contribution intervals

Z in ['9/25', '19/50'], 2048 source cells:

| Channel | C0 modulation contribution | Ordinary-Z modulation contribution |
|---|---|---|
| m | [-5.27967744e-523, 5.66200849e-523] | [-2.22017861e-6, 2.22017861e-6] |
| h | [-5.34645969e-5, 6.22935505e-5] | [-7.33169087e-5, 6.85695614e-5] |
| k | [-1.85996565e-522, 1.76477647e-522] | [-1.52218697e-6, 1.52218697e-6] |
| e | [-9.22740728e-5, 7.91033349e-5] | [-0.000190042986, 0.000206438769] |
| p | [-0.000201031926, 0.000137768828] | [-0.000354095964, 0.000437952946] |

Z in ['-19/50', '-9/25'], 2048 source cells:

| Channel | C0 modulation contribution | Ordinary-Z modulation contribution |
|---|---|---|
| m | [-5.83395875e-523, 6.01971907e-523] | [-2.13157336e-6, 2.13157336e-6] |
| h | [-4.05599542e-5, 7.29830323e-5] | [-5.99049274e-5, 8.14563191e-5] |
| k | [-2.21616541e-522, 1.9536003e-522] | [-1.46118689e-6, 1.46118689e-6] |
| e | [-0.000125593669, 4.30595363e-5] | [-0.000251906673, 0.0001440292] |
| p | [-3.27405926e-5, 0.000301474214] | [-0.000220345728, 0.000571422338] |

These are normalized O2 contribution enclosures, not selected amplitudes, terminal closure defects, cone margins or full Cartesian residuals. Widths include variation throughout each Z tile and conservative source/inverse bounds. This advances source-functional transport; it does not solve the five terminal identities or establish coefficient/scale recursion.

## Executable next tasks

- [x] **CONTINUOUS-STRICT-SIGN-O2-C1:** entire y[0,1] / Z[.36,.38] and Z[-.38,-.36] rectangles, N7, 256/2048 refinements; original C0 and genuine density-Z graph, all original source/phase pieces and own-rate positive masses.
- [x] **SOURCE-DEFINED-NONMIDPLANE-INLET:** source cumulative history formulas, common Pstar unit, exact ordinary-Z jets, positive formal inverse Pstar/Pstar², no midplane constant substitution.
- [x] **ORIGINAL-HISTORY-C1-TRANSPORT:** original unmodulated value/Z densities, source inlet and independent defining-mass endpoint comparison. Separate P0 and cumulative pressure.
- [x] **EXPLICIT-INCOMING-C1-API:** all five incoming values and ordinary derivatives required; original + modulation + own-rate decayed incoming histories. Synthetic wiring test is not physical inlet data.
- [ ] **REAL-INCOMING-FUNCTION-ATTACHMENT:** attach the current inner/O1-to-O2 modulation inlet value and ordinary-Z functions, with same family, datum, full Z domain and source DAG/interval provenance. Pass them to apply_actual_incoming; no assumed zero incoming, and no source-bound original history reused as a modulation defect. Acceptance: real source functions enter all five channels and pressure memory at both interval boundaries and through quiet charts.
- [ ] **MICROSCOPIC-CROSSING-Z-ATLAS:** factor original odd p2=Z*Q with even Q and retain exact nonzero p2_Z. Build overlapping conditional signed-negative, regular |u|<=1/4 and signed-positive source domains around Z0, where u=p2*q/dstar. Keep original ordinary jets and every nonempty region. Acceptance: no unresolved y/Z state omitted, no numeric cutoff that removes native scales, no derivative zeroing at Z0.
- [ ] **CROSSING-Z-DERIVATIVE-UNITS:** ordinary derivatives at exact Z0 can be enormous in the current units. Keep signed/formal enclosures until source/density correlations justify finite transport. Derive joint expressions or rescale a reported bound with an explicit physical unit; do not clamp or reinterpret interval endpoints as source values.
- [ ] **AXIAL-TILE-REFINEMENT:** extend beyond these two tiles toward the full original Z domain with adaptive whole source rectangles. Join overlapping tile covers of the same function without summing duplicate contributions. Acceptance: coverage map, true source/inverse phase unions, common-family value/Z results and unresolved states listed explicitly.
- [ ] **DERIVATIVE-WIDTH-BUDGET:** isolate contribution widths from coefficient boxes, logq endpoints, inverse precision, source-u splits and rational direction intersections. Improve the limiting term using exact original identities/adaptive y-Z-phase refinement. Acceptance: reduced downstream target/control uncertainty on the same genuine source, not merely more checks.
- [ ] **ALL-CHART-ORIGINAL-INTEGRALS:** evaluate original and modulation functions over inner/O1, O2, O3, pulse/end, flatten and heat route under their own actual source units. Retain histories after local velocity modulation ends. Acceptance: all source windows and seams present; the shipped one-unit O2 integral is not a full Rc terminal integral.
- [ ] **CORRELATED-RC-TARGETS:** feed actual full-route original/incoming functions to the existing exact terminal graph. Retain joint C=D_k-A*D_m and its ordinary derivative before dividing by mu*A², and keep independent P0/P0_Z. Acceptance: real evaluated target functions with quantitative errors across the axial atlas.
- [ ] **FIVE-FUNCTIONAL-CONTROLS:** solve the original independent five bump/Gram response system with actual support, regularity and source datum. Acceptance: five terminal identities as Z functions, bounds for nonlinear repair terms, and an independent response/conditioning result. A handful of point fits is insufficient.
- [ ] **ONE-ADMITTED-GLOBAL-N:** N7 is only the genuine phase candidate used by this O2 receipt. For different N, regenerate actual radius-phase unions and downstream contribution covers; combine complete route amplitude, repair, cone and regularity constraints. Acceptance: one finite integer valid for the full original construction, with inherited local diagnostic frequencies excluded from that claim.
- [ ] **RADIAL-HIGHER-Z-JETS:** recover ordinary y and higher Z source/implicit primitive/history derivatives needed by pressure/stress and matching, preserving cross terms and shared local log bases. Acceptance: required mixed order at every chart and seam, not a derivative of a bound selector.
- [ ] **MATCHED-HEAT-BACKGROUND:** inner/annulus/pulse/flatten/heat joins, axis regularity, compatible analytic preheat datum, exact heat exterior and finite-energy radial tail. Acceptance: functional moment/pressure matching and correct derivatives through all interfaces.
- [ ] **GLOBAL-ADMISSIBLE-STRESS:** recover divergence-form T_B and separate E_B from the matched background; check the actual cone margin on the complete physical/similarity source domain, and high-order/flat remainder decay. A normalized O2 moment interval is not a Cartesian momentum residual.
- [ ] **TRUE-N-DEPENDENT-RECURSION:** implement original coefficient recovery separately for n=1 and n>=2, common inner domain, independent moment repair, finite-order remainder and smooth summation. Coordinate rescaling alone does not satisfy this task.
- [ ] **OSCILLATORY-STRESS-CANCELLATION:** mean corrections and both original pulse families; averaged quadratic stress realization/cancellation with actual parameters and remainder budget.
- [ ] **PHYSICAL-UVW-AND-DYNAMICS:** expose Cartesian [u,v,w](x,y,z,t) from the accepted complete corrected construction; independently evaluate forced NS residual, energy, core contraction/aspect, swirl/vorticity scaling and material winding. Full Linfinity/volume-L2 residual<1e-3 remains a later corrected-field target.

The long-term goal remains active and incomplete. The next practical step is real incoming function attachment plus microscopic crossing-Z coverage, followed by full-route terminal targets and independent repair.
