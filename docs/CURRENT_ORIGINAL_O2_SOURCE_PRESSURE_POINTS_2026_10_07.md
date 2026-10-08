# Original O2 correlated parameters and normalized pressure point jets

Successor: [CURRENT_ORIGINAL_O2_FACTORED_INPUTS_PHASE_2026_10_07.md](CURRENT_ORIGINAL_O2_FACTORED_INPUTS_PHASE_2026_10_07.md) completes original O2 factored finite inputs with directed errors and conditional true radius phase points/covers. Next production is original conditioned loop scales/inverse/A-B primitives and one signed density integral; full oracle, controls and recursion remain open.

Checked source: [4ea85e94](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/4ea85e94699ea4cc923d2025ceb5901a824aef0f). Predecessor: [full original O2 inertial/pressure expressions](CURRENT_ORIGINAL_O2_INERTIAL_FUNCTIONS_2026_10_07.md). Parameter-frame and pressure-point checkers PASS. Git-index dependency audits cover 611 and 615 direct hashes. Existing read-only reviewer: **GPT-5.6 Luna / max**; no new workers.

Two production interfaces now bridge the expression layer toward actual original source evaluation. A common exact parameter frame binds the selected radius/scales and the prescribed raw waiting root. A normalized original pressure service returns approximate P0/Pstar² and its first two ordinary Z derivatives, with independently directed coefficient/arithmetic error and a source-bound all-late-stage error budget. It does not install the complete conditioned numerical p1/p2, phase, all-chart oracle, controls, common N or recursive field.

## One original source parameter frame

`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_O2_source_parameter_frame.py/.json` and `_check.py/.json` implement `OriginalO2SourceParameterFrame`. The same checked source family pins Md=40 and

```text
logPstar=exp(40)+11; yd=logPstar; Td=exp(40)+10
logmu=-log(1000)-4logPstar
logdelta=-4logPstar-30; logepsilon=logdelta-log(1000)
Tw=-60logmu; Ts=4(log(2)-logdelta); L=-30logmu; Tf=100
logRref=log(110)+10(logCstar+logPstar)
R=Rref*exp(y), y=log(R/Rref).
```

The min branch for delta is proved with integer inequalities: exp(40)>=841, so 4logPstar+30>=3438, whereas 200log(10)<600. No rounded interval branch is selected. The selected logCstar is the physical norm builder's explicit singleton dyadic choice, whose exact MPF tuple and canonical definition/family hash are verified. This is a permissible chosen source parameter; pressure, waiting and radius enclosure endpoints are not selected as field values. The exact radius formula is reconstructed from logCstar and exact logPstar rather than copying a rounded `selected_logRref` cap.

`ExactPositiveDyadic(mantissa,exponent)` denotes exactly mantissa*2^exponent without expanding an astronomical integer. `ExactSourceExponential(argument)` denotes the exact positive exp(argument), protected from automatic evalf. Ordinary SymPy assumptions otherwise tried to materialize a nested exponential and raised MemoryError. These protected operations are source definitions, not a floating-point evaluation backend.

`.functions(y)` ties all p1/p2 expression inputs, pressure parameters and radius to this frame. `.pressure_at(Z,order)` retains exact original integrals and their prescribed raw waiting root. All parameter symbols are eliminated; at a fixed radial coordinate the full p1/p2 value/Z expressions have only Z free. Approximate radial coefficients remain explicitly marked. The actual selected-radius core/outer field is still not numerically materialized or globally certified. Nonfinite MPF coordinates are rejected before tuple conversion, rather than being confused with zero.

## Original normalized pressure point service

`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_pressure_point_jets.py/.json` and `_check.py/.json` implement `OriginalNormalizedPressurePointJets(...).evaluate(Z)`. It computes the actual defining quadrature

```text
alpha = 5/2 + (1/2)*integral_0^1 exp(s/5-(6/5)*J_sigma(s)) ds
              + exp(-2/5)/2
approximate ordinary Z jets = -alpha * D_Z^j[(1+Z²)^-2], j=0,1,2.
```

The point coefficient is freshly evaluated from the source integrand. A separate 4096-cell directed calculation of the TRUE coefficient, using the existing original transition integrals, bounds its quadrature and arithmetic discrepancy. The point is never chosen as an interval midpoint/cap and is never asserted exact. The early reference mass is exactly5/2. Extending the finite axial segment to infinity creates a known removed tail; all eleven later original stages are retained in the error budget.

For real |Z|<=1, q=1+Z², the exact post-yd densities are algebraically matched to an integrable envelope

```text
envelope = exp(3/5-yd)*exp(-(y-yd))/(2q²), y>=yd
late density <= envelope/(1-epsilon)² <4*envelope.
```

The proof explicitly checks all eleven stage ratios and the removed finite axial-tail identity. It binds the original sigma/J definitions, phi exponent and unique positive raw waiting root to accepted source receipts. With 0<=J(t)<=t, 0<=sigma,phi<=1 and the fixed positive stage lengths, all remaining exponents are nonpositive. The source domain follows the Md=40 parameter contract.

The exact density derivative factors are checked through order2 on all fourteen stages. Beta-2 factors have bounds(1,2,4), beta-0 factors(1,0,0), and flatten factors(1,2,10). Including the removed axial extension, absolute late-source errors satisfy

```text
error_j <= (5,10,44)_j * exp(3/5-yd)/(2q²).
```

The returned budgets keep coefficient/roundoff error and this logarithmic source-tail error separate. Exact odd symmetry at Z=0 has an explicit `exact_zero` tail record. No positive tail is replaced by zero elsewhere. The exact fourteen-stage pressure expression remains available in the frame.

At Z=0 the approximate normalized datum is **-3.31462273001433254645298506290413385879609451632699097337247783174516159631**. The coefficient/arithmetic absolute error upper is **0.0000753670955476493668265612450611670892970945547325723890791551123793761365055**, plus a late-source error with log upper **-235385266837019994.891609178874879739325344405486544395786550031281026695734**. These are normalized pressure errors; physical pressure error is Pstar² times their sum. The normalized numbers alone do not bound physical pressure or the NS momentum residual, and cannot establish admissibility or terminal matching. Keep the original analytic pressure identity and correlated cumulative p when installing pressure in the field.

## Focused evidence and next production

The frame checker verifies parameter/radius/length correlations, the rigorous delta branch, safe opaque constants, complete pressure parameter binding at three Z orders, Z-independent substitution of the previously checked full I/F formulas, and seven source/domain guards. The pressure checker verifies four actual native normalized queries (twelve value/Z rows), signed symmetry and explicit zero budgets, five guards and one independent full fourteen-density integration in finite diagnostic units. The latter includes 28 nonzero density integrals through Z order2 and checks the general remainder lemma. Those diagnostic units do not select native parameters. The original full numerical oracle and all global field gates remain false.

- [x] **O2-SOURCE-PARAMETER-CONTRACT:** common exact/factored parameter definitions, selected logCstar/radius prescription, source-family/hash bindings, exact raw W and all stage lengths.
- [x] **O2-PRESSURE-NORMALIZED-POINTS:** approximate actual P0/Pstar² and first/second Z jets, independent directed coefficient/arithmetic bounds, source-bound finite-end/all-late errors and physical amplification accounting.
- [x] **O2-INERTIAL-EXPRESSION/PROFILES:** prior full p1/p2 expressions, original E/V and five-history coefficients remain available.
- [ ] **NEXT O2-FULL-FACTORED-POINT-ADAPTER:** produce complete original E,V,a,b,p1,p2 value/Z point objects with correlated R/Pstar powers and certified finite-coefficient errors. Extend the independent directed slope calculations to all radial f,H,D,P inputs. Combine pressure-point budgets with full original I/F recovery. Do not cast the astronomical factors to floats or install a full-oracle mode prematurely.
- [ ] **O2-CORRELATED-ABSOLUTE-PRESSURE:** keep P0 and cumulative p tied to their common defining integrals. Evaluate cancellation-sensitive sums by joint source integrals; independent forward enclosures or pressure approximations do not prove exact terminal closure. Propagate physical Pstar² and R/Pstar factors explicitly.
- [ ] **O2-CONDITIONED-LOOP-SCALES:** attach the original whole-input a/b/cone/margin contract and eta/d_star definitions. Keep microscopic excesses and Poisson widths correlated with p1/p2 and source scales.
- [ ] **O2-GLOBAL-PHASE-ORIGIN:** bind the actual original r_minus and integer-N fractional phase before conversion. Investigate exact removal of the selected logCstar's integer dyadic contribution to the phase; preserve the origin and any other correlated terms. This is not a new phase or selected common N.
- [ ] **O2-CONDITIONED-PHASE/ONE-INTEGRAL:** implement native factored angle inversion and A/B slow-Z primitive evaluation, then one genuine signed density own-rate integral with coefficient, phase, integration and roundoff errors. Reuse the existing scalar backend and exact N-dependent density formulas.
- [ ] **ORACLE-ALL-CHARTS/CONTROL:** extend source recipes to remaining fixed original cells, complete error/stability/tail contracts and whole-frequency constraints, choose one compatible global N, then install all five original controls and functional terminal bounds.
- [ ] **FIELD/OUTER/REC/WAVE/PHYS:** corrected field/pressure and joins, global stress cone/flat remainder/energy, genuine n-dependent coefficient recursion with repair and flat summation, oscillatory quadratic cancellation and physical Cartesian uvw/scale diagnostics.

Long-term goal remains active. These interfaces advance source computation; they do not establish a blow-up solution or completed scale recursion.
