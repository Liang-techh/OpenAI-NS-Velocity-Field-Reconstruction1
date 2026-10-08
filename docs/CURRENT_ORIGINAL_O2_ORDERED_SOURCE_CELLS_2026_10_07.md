# Original O2 ordered defining masses and nonzero-Z source coefficient cells

Checked source: [3eaaf5af](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/3eaaf5afa31224bce2c4bcb673a0ccf3169d7d36). Predecessor: [five midplane own integrals and incoming memory](CURRENT_ORIGINAL_O2_FIVE_OWN_INTEGRALS_2026_10_07.md). Focused checker PASS; Git-index dependency audit covers1078 hashes. The existing **GPT-5.6 Luna / max** read-only reviewer checked the source normalization, pressure correlation and factor bridge; no new worker, ancestor producer or repeated nested source quadrature was used.

Three original defining mass integrals now execute on ordered continuous source cells over the complete O2 y in[0,1] window. Their prefixes/suffixes produce whole-cell f/H/D/P ranges. Full original E/V/a/b/p1/p2 coefficient ranges and ordinary Z derivatives consume these ranges for fixed nonzero Z and entire Z intervals. A pressure suffix identity cancels the large baseline terms symbolically before interval arithmetic, retaining the original datum and positive late-source error.

This completes the defining-mass/source-coefficient step required before nonmidplane integrals. It does not complete those integrals, Z-functional terminal matching, all-chart histories, global N, five repair controls, matched stress or coefficient recursion. Some joint source cells still require signed-source refinement; their status is preserved explicitly.

## Executable source interfaces

Production/check: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_O2_ordered_source_cells.py` and `_check.py`. The complete producer report is `.json.gz`, about50.8 MB, losslessly compressed from210.6 MB. `_check.json` binds both compressed/decompressed digests and every dependency. No cell evidence is discarded. Source integration uses90 directed digits; native factor/conditioning uses260 digits to retain the original selected auxiliary log bits. The mass ranges are rounded outward into the native context, rather than mixing contexts in coefficient arithmetic.

```python
import lei_ren_part1_paper_compliant_current_original_O2_ordered_source_cells as source

owner = source.OriginalO2OrderedSourceCells()
saved = owner.parent.saved['actual_original_pressure_integral_refinements'][-1]
level = owner.integrate(saved)
cell = owner.coefficient_cell(8192, 4341, Z_lower='.36', Z_upper='.38')
```

`cell['record']` is serializable continuous-range evidence. `roots`, `q`, `kernel` and `ledger` expose the live conditioned-source bridge for the next inverse/integral work. A returned geometry requiring refinement is not an installed inverse. No source coefficient or enclosure midpoint is selected as a field point.

## Original masses and continuous profiles

The same accepted original monotone J source cells and exact J(1)=1/2 are reused. For each defining pair `(rate,power)=(8/5,1),(1/5,2),(6/5,2)`, integrate

```text
Mj(y) = integral_0^y exp(rate*s-(3/5)*power*J(s)) ds.
```

On [left,right], the positive exponential weight is enclosed directly from `(exp(rate*right)-exp(rate*left))/rate`. Multiplying it by the full cell's `exp(-(3/5)*power*J_cell)` encloses the continuous mass contribution. The finite shipped grid width is at least1/8192, so the directed subtraction at90 digits is well conditioned. This is an enclosure of the exact exponential integral, not a claim that finite endpoints are exact real values. Earlier own-rate mean-exponential masses were also enclosures of their exact Duhamel integrals.

Forward prefixes and backward suffixes retain all three positive original integrals. Monotonicity gives whole-cell mass ranges from adjacent node bounds. Full y/mass rectangles then enclose

```text
f = exp(y/10-(3/5)*J(y))
H = (5/8+M0(y))*exp(-3*y/2)
D = (5/12+M2(y)/2)*exp(-y)
P = 5/2+M1(y)/2
a = 4/5+(6/5)*sigma(y).
```

Endpoint H/D/f values alone are insufficient because these profiles need not be monotone. The nonzero reference inlets remain. H, -D and P at y=1 agree with the accepted original midplane own-history enclosures. These defining rates differ from the normalized own-history decay rates1,3/2,3/2,1,0.

Independent original defining J/M0/M1/M2 ODE references give terminal values approximately

| J(1) | M0(1) | M1(1) | M2(1) |
|---:|---:|---:|---:|
| 0.5 | 2.20790293024068 | 0.958925413993024 | 1.59322317005450 |

Two numerical tolerances differ by less than7.8e-15. Directed64/256/2048/8192 refinements contract for every defining mass. The checker verifies31680 continuous mass contributions,63360 profile covers and96 defining-ODE terminal/prefix comparisons. The floating reference corroborates directed source ranges and is not itself a functional proof.

## Correlated pressure datum and complete source factors

The actual baseline alpha source definition is explicitly bound to `current_original_pressure_point_jets.py`:

```text
alpha = 5/2+M1(1)/2+exp(-2/5)/2.
W(y) = integral_y^1 original_M1_density(s) ds/2+exp(-2/5)/2.
alpha = P(y)+W(y).
```

W is enclosed by backward suffix node integrals over the whole cell. Original P0 baseline rows use this exact identity:

```text
P0/S² = -(P+W)/(1+Z²)²
(P0/S²)_Z = 4*Z*(P+W)/(1+Z²)³
(P0/S²)_ZZ = (4-20*Z²)*(P+W)/(1+Z²)^4.
```

The coefficient substitution is verified symbolically against the original unmodified templates. Cumulative P and P0 remain distinct; the identity cancels correlated baseline terms in the stress coefficients before finite arithmetic. The true original pressure differs from this baseline by the previously accepted finite-end/all-late-source remainder. Its positive error logs for orders0/1/2 retain factors5/10/44 and are propagated through the original template sensitivities. Only exact Z0 odd pressure symmetry removes the order-one tail. Nonzero or crossing Z ranges retain it. Physical S² amplification is carried through the original factor sectors, not discarded.

At8192 cells, the baseline alpha source enclosure is approximately[3.31460519768,3.31464026009]. It overlaps the accepted independent alpha enclosure and contains its independent defining quadrature approximation. This narrows an enclosure of the same source function; it does not choose a new pressure datum.

Complete original p1/p2 templates keep the factor layout `(R,Pstar,delta,L)`, L=1-delta*Z², including p2's two Pstar sectors -1 and +1. Native arithmetic maps powers `(r,p,d,ell)` to `(p,d,ell,0,r)` in the basis `(logPstar,logdelta,logL,zero,logR)`. An AST binding checks that this is also the accepted original conditioned consumer's semantic order; the third slot is logL. The original q source has zero L/radius powers before being restored into this basis. This guard prevents silently reading the logL slot as an old zero slot.

E=f/(1+Z²), V=4Z/Pstar, a_Z=b=b_Z=t0=t0_Z=0; their ordinary Z rows remain. p2=V=0 only at exact Z0. The nonzero p2_Z and V_Z at Z0 are preserved. No first-derivative factorial division is applied. The independent checker uses the unmodified original coefficient templates and defining ODE values for1539 finite-coefficient diagnostics with a declared2e-12 floating comparison tolerance. Continuous-range certification comes from the directed covers and symbolic identities; these point diagnostics are not offered as its proof.

## Current conditioned geometry and remaining tasks

Four source queries, including the interval Z in[.36,.38] and its negative counterpart, admit the signed Mobius geometry on their continuous y cells. The exact midplane uses the regular small-r branch. Nine total source queries include all of Z in[-1,1], crossing-zero and near-end cells; four correctly remain `requires_signed_source_refinement`. A q enclosure whose lower endpoint is zero is not declared an exact flat branch. There is no nonmidplane inverse/integral admission in this checkpoint.

- [x] **ORDERED-DEFINING-MASSES:** all three actual original mass node prefixes/suffixes, whole y/mass profile covers, complete window/refinements and defining ODE reference; no per-node nested quadratures.
- [x] **PRESSURE-SUFFIX-CORRELATION:** same baseline alpha=P+W, symbolic pressure/stress coefficient cancellation, original datum/factor sectors and nonzero late errors.
- [x] **NONZERO/Z-RANGE-SOURCE-COEFFICIENTS:** complete original C0/ordinary-Z E/V/a/b/p1/p2 range bridge, positive L, both p2 Pstar sectors, exact Z0 parity/nonzero p2_Z and explicit semantic basis guard.
- [ ] **NEXT LOG-Q-WHOLE-CELL-ENDPOINT:** recover log(1-sigma) directly from the original flat cutoff. For y>=1/2, set lambda=1/y²-1/(1-y)² and log(1-sigma)=lambda-log(1+exp(lambda)), with exact endpoint extension. Enclose q²=(.6*(1-sigma)+eta)/a using the same native positive eta in log form. At a cell ending at1, the eta term supplies a finite positive logarithmic lower bound. Keep the full source range, possibly as a positive formal log interval with unit coefficient; never choose a q midpoint/cap or treat a zero-touching finite box as flat. Demonstrate that the signed-source conditioner consumes this complete positive source cover on the near-end nonzero-Z cells.
- [ ] **SIGNED-Z-CROSSING/REGULAR-COVER:** preserve the original odd p2 factor and its actual sign away from Z0. Use compatible regular/signed coordinates or justified source splits for joint boxes crossing p2=0; retain p2_Z. Do not choose a convenient signed coefficient value or discard the microscopic regular region. Publish any unresolved boxes explicitly.
- [ ] **NONZERO-Z-ACTUAL-PHASE-INTEGRALS:** use the complete continuous source bridge, actual common-N phase unions and original inverse/A-B formulas. Enclose all five signed Duhamel contributions over the complete window at fixed nonzero Z, then Z ranges. Preserve native q/rho/s positivity, original B/cross terms, source pressure correlations and incoming five histories/P0. A successful source geometry query alone is not a certified phase integral.
- [ ] **Z-FUNCTIONAL-OWN-HISTORIES:** integrate genuine density Z rows with original own rates, shared sources and nonzero p2_Z. Build correlated terminal histories as functions of Z and exact five functional controls. A finite axial sample grid cannot replace terminal identities.
- [ ] **RADIAL/HIGHER-JETS/ALL-CHART-ORACLE:** supply slow-y and fast-phase terms through the required ordinary recovery/stress orders; propagate actual incoming functions through all charts/gaps, flatten, collar and exterior. Keep all cumulative memory and same analytic datum. Install the full source/integral oracle only when the entire required domain and jets are available.
- [ ] **COMMON-N/FIVE-REPAIR/MATCHING:** admit one compatible N under simultaneous source/tail/stability constraints, solve independent repair functions and match all five terminal functions, analytic pressure, smooth annular/exterior joins and finite energy. N7 remains diagnostic.
- [ ] **STRESS/TRUE-RECURSION/OSCILLATORY/UVW:** actual global cone/flat remainder, correct n=1 and n>=2 coefficient equations with independent repair, divergence-preserving truncation and flat summation, two pulse families/mean corrections and averaged stress cancellation, then physical Cartesian uvw/full residual and core/winding diagnostics.

The full goal remains active. This is real advancement of continuous original sources and pressure correlation; actual coefficient/scale recursion and the complete forced-NS reconstruction remain open.
