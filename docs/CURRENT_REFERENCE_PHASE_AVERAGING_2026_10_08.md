Current successor: [CURRENT_REFERENCE_MIXED_C1_2026_10_08.md](CURRENT_REFERENCE_MIXED_C1_2026_10_08.md) installs fixed nonzero-Z mixed source and C1 phase integration. The record below is historical; whole-Z/all-route closure remains open.

# Original reference phase averaging and actual slow radial jets

Checked source [56f3784a](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/56f3784ae318b49f23d638ff9ff9cc88f8a7339e). Predecessor: [CURRENT_REFERENCE_WHOLE_CELL_INTEGRALS_2026_10_08.md](CURRENT_REFERENCE_WHOLE_CELL_INTEGRALS_2026_10_08.md). Reused the existing GPT-5.6 Luna/max worker for read-only formula/source review; root implemented, computed, checked and published. No child was spawned and no ancestor producer was rebuilt.

## Executable change

**OriginalReferencePhaseAveraging** adds actual slow native log-radius derivatives to the original closed Rh_reference source. For each original term R^r*c, it computes

```
partial_y(R^r*c) = R^r*(r*c + f*c_f/10 + H*c_H/10 + D*c_D/5 + P*c_P/5)
```

Original P0 and its Z jets, Pstar,delta and L are fixed under this radial derivative. The same pressure-affine late error retains its Pstar^-1 factor in the differentiated term. The canonical atlas preserves exact source powers.

The accepted first-direction implicit chain rule consumes these actual radial rows at fixed original phase. Its derivative slots are reused with an explicit slow-y contract; they do not receive ordinary Z values or the total fast derivative N*f_phi. The bounded t*psi_y product is used just as in the previous original Z path, with the same exact rational-factor bound.

The original reflection identities A(1-phi)=-A(phi), B(1-phi)=-B(phi) give zero mean for the **linear leading densities only**. The full finite-N density remains

```
F_N = E*A*exprel(A/N) = E*A + Q_N/N
Q_N = E*A^2*R2(A/N)
R2(x) = integral_0^1 (1-t)*exp(t*x) dt
delta_j = L_j/N + S_j(N)/N^2

L = (B, E*A, V*E*A+E*B, 2*V*B-E^2*A, E^2*A)
S = (0, Q_N, V*Q_N+F_N*B, -E*Q_N+B^2-F_N^2/2, E*Q_N+F_N^2/2)
```

The row order is m,h,k,e,p. The nonlinear terms are evaluated from the original whole-cell source, with actual N in both analytic exponential factors. Their mean is not discarded. R2 uses a directed order64 integrated Taylor tail, including its exact zero-argument value1/2.

For each leading density, G(y,phi)=integral_0^phi L(y,s)ds is periodic with zero traces. The source-function bounds |G|<=sup|L|/2 and |G_y|<=sup|L_y|/2 enter

```
integral K(y)*L(y,frac(N*y+phi0))/N dy
  = ([K*G] - integral K*(G_y+rate*G)dy)/N^2
```

K(y)=exp(rate*y), native radius Jacobian1, rates m/e=1,h/k=3/2,p=0. The actual original phase origin, including its positive microscopic budget, is retained. The estimate applies to that phase without enumerating periods.

Only internal boundaries within this **one original closed reference source** telescope. The physical interval's two endpoint terms remain. G's periodic traces also handle modulo wraps. No cross-chart matching or terminal seam is assumed. Pressure rate-zero mass still totals5; P0 and unknown incoming histories are kept separate.

Direct source-range integration runs alongside averaging. The tighter valid source-factor bound is retained. This preserves the microscopic m/k C0 ranges when a derivative-based averaging bound would be weaker. Supremum bounds are used solely in error estimates, never as field values.

Four files with stem **lei_ren_part1_paper_compliant_current_original_reference_phase_averaged_integrals** are committed under experiments/root_st073: producer, compressed manifest, checker and receipt.

## Actual integral bounds

Complete reference-window[-5,0] contributions at Z=37/100 were computed with4/16 cells at N160 and16 cells at N320. Rounded outward absolute C0 bounds:

| Row | Direct,16 cells,N160 | Averaged,4 cells,N160 | Averaged,16 cells,N160 | Averaged,16 cells,N320 |
| --- | ---: | ---: | ---: | ---: |
| h |0.000971|0.000032489|0.000014791|0.000003698|
| e |0.001154|0.000037562|0.000015752|0.000003938|
| p |0.004394|0.000111770|0.000037532|0.000009382|

At N160/16 cells, the actual recorded bound reductions are approximately66,73 and117 times. N320 independently reevaluates Q_N and F_N; it is not a relabelled N160 result. These are source-window integral envelopes, not a terminal five-moment error or Cartesian NS residual. Their signs generally remain unresolved.

Ordinary Z contribution records from the previous N160 integration are retained exactly and labelled separately. **The current averaging tightening is C0.** Mixed-yZ source evaluation is still needed to tighten the Z integration by parts. No N320 Z result is claimed from the N160 data.

Producer **26.609s**, focused checker **26.610s**. Evidence:

- 18 independent defining weighted-exponential integral comparisons;
- 8 exact original radial factor/profile chain-rule identities and10 independent finite-unit original source derivative comparisons;
- 5 exact density split identities,5 linear reflection identities, weighted IBP and same-source endpoint telescoping;
- 12 independent weighted sinusoidal antiderivative comparisons with arbitrary phase origins, both frequencies and all three rates;
- actual source/phase/partition/zero-rate memory checks, strict h/e/p tightening, preserved microscopic m/k bounds and exact predecessor Z records;
- **1208 Git-index dependency hashes PASS**.

Finite diagnostic units are isolated from the original frame, source values and frequency choice. The existing scalar evaluator guards remain unchanged.

## Updated next actions

Read this handoff before historical queues. Mark DONE with implementation, scoped receipt and commit. Retain the full goal and remaining route in the predecessor queue; this milestone does not redefine completion.

- [x] **REFERENCE-ACTUAL-SLOW-Y-SOURCE:** differentiated full original reference terms, including radius powers and source pressure error factors.
- [x] **REFERENCE-LEADING-PHASE-CENTERING:** original A/B reflection and all five signed leading densities; no nonlinear zero-mean shortcut.
- [x] **EXACT-N-DEPENDENT-SECOND-REMAINDER:** actual Q_N/F_N and all cross/quadratic terms.
- [x] **REFERENCE-C0-ENDPOINT-IBP:** actual global endpoints, same-source internal telescoping, own-rate kernel and pressure memory.
- [x] **DUAL-DIRECT/AVERAGED-C0-BOUNDS:** retain whichever valid bound is tighter; do not degrade microscopic m/k results.
- [x] **INDEPENDENT-REFLECTION/RADIAL/REMAINDER/WEIGHTED-ENDPOINT-EVIDENCE:** committed with lossless source records.
- [ ] **REFERENCE-ACTUAL-MIXED-yZ-ROOTS:** differentiate the original ordinary-Z templates radially, including the radius factor. Preserve original L_Z terms already present in the Z row without differentiating them twice. Keep the shared P0 function and pressure-affine late factors. Acceptance: issued same-family C0,y,Z,yZ source rows over whole radial/Z cells, with independent mixed derivative references.
- [ ] **CORRELATED-PHASE-MIXED-yZ-EVALUATOR:** execute the exact original implicit mixed derivative graph, retaining T2_yZ and both first cross terms. Preserve positive rho/s/hinv and source factors across signed E/psi and small-u transitions. Do not multiply independently enormous direction/inverse ranges or differentiate saved bounds as functions. Acceptance: actual A_yZ/B_yZ enclosures, source branches/phase held fixed and independent finite-unit mixed references. If a genuine large factor remains, retain it formally rather than assert a finite bound.
- [ ] **REFERENCE-C1-ENDPOINT-IBP:** use true G_Z/G_yZ source bounds and the actual nonlinear remainder Z derivative. Keep physical endpoint terms, own rates and separate oracle/rounding errors. Acceptance: complete reference C0/Z averaging envelopes, not the current C0 result plus relabelled Z data.
- [ ] **WHOLE-Z/MIDPLANE-SIGNED-TRANSITION:** cover[-1,1] as Z functions, preserving p2=0/nonzero-p2_Z and the small-u branch. Acceptance: actual entire-window C1 function enclosures; fixed-Z sampling is insufficient.
- [ ] **SHARP-SIGNED-PERIOD-MEANS:** when needed, integrate the original finite-N nonlinear phase mean with controlled radial drift rather than only bounding it. Keep the original inverse/phase Jacobian and signed leading identities. Acceptance: genuine signed integral and separate mean/drift/endpoint/error records at the common candidate N.
- [ ] **O2-C1/AXIAL/11-UNIT-BUFFER:** apply the same actual-source method to accepted O2 kernels and inherited histories. Use exact native offsets and Jacobians; general a_y/b_y/q_y terms cannot inherit the reference's zero identities.
- [ ] **REMAINING-SOURCE-UNITS/O3/RESTORE/ALL17-24-ORACLE:** preserve the predecessor's detailed tasks. Integrate genuine source functions and incoming histories, with separately typed factored integral outputs and rejection of unsupported charts.
- [ ] **ACTUAL-FIVE-CONTROLS/TERMINAL-CLOSURE/EXECUTABLE-GLOBAL-N:** attach the complete C1 oracle to the accepted centered all-N/unit-ball bridge and numerical Picard replay. Candidate160/320 is not a global frequency selection. Keep numerical oracle errors separate from the Picard tail.
- [ ] **BACKGROUND-JOINS/EXACT-HEAT/STRESS/FLAT-REMAINDER/n-RECURSION/TWO-PULSES/CORRECTED-UVW:** full original objective and predecessor queue remain active. Actual n-dependent recursion and the final corrected Cartesian field are still open.

Mixed-yZ averaging, whole-Z terminal closure, all17/24 oracle, actual five controls, global N and recursive corrected-field flags remain **false**.
