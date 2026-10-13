# True original n1 common-inner numerical enclosures

Follow-on cutoff/own-history and five inner numerical ranges are now checked. Read [CURRENT_ORIGINAL_N1_CUTOFF_NUMERIC_BRIDGE_2026_10_13.md](CURRENT_ORIGINAL_N1_CUTOFF_NUMERIC_BRIDGE_2026_10_13.md) for the latest supported scope and next tasks.

The accepted true Volterra functions now have a directed numerical
enclosure interface. Three queries of the original source and parameters
have executed: an axis point, an interior core point and an actual
prescribed first-collar point. The full function graph remains the
definition; numerical centers, source bounds and tail bounds do not become
selected function values or coefficients.

The long-term goal remains **ACTIVE / INCOMPLETE**. This result covers
the canonical **unrepaired common-inner n1**. It does not evaluate the
complete cut-and-repaired family, close five terminal moments, construct
n>=2 recursion, or certify a corrected Navier–Stokes residual.

## Results and supported interface

| Chart | Z | Coordinate | Value/Z1 budget rows | Largest relative natural-component interval width |
| --- | --- | --- | --- | --- |
| Axis | .125 | rho=0 | Axis exact-zero branch; relative test inapplicable | Inapplicable |
| Core interior | .371 | rho=2 | 12/12 | 9.53021667222e-29 |
| Actual first collar | .371 | q=.625 | 12/12 | 4.51783342117e-10 |

Widths include the nonzero source errors and all remaining Volterra words.
They measure the pointwise value/first-Z enclosure in the accepted natural
component frame. They are not a momentum residual, a moment defect, or a
uniform error bound on the whole physical field. Ordinary Z2/Z3 intervals
are supplied by the accepted analytic bounds; comparable relative
precision is not claimed for those higher rows.

The checked quartet is
`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_numeric_bridge`
(`.py`, `.json.gz`, `_check.py`, `_check.json`). Its new gate is
`current_original_true_unrepaired_n1_common_inner_numeric_enclosure_bridge_installed`.
The receipt records `original_numeric_n1_solution_oracle_executed=True`
with that explicit local scope.

Archive-only restoration does not rebuild ancestors or run numerical
queries:

```powershell
python scripts/run_original_n1_inner_numeric_bridge.py
```

A fresh actual query restores the accepted numerical/source owners once:

```powershell
python scripts/run_original_n1_inner_numeric_bridge.py --evaluate --chart core --coordinate 2 --Z .371 --output tmp/n1_core_query.json
python scripts/run_original_n1_inner_numeric_bridge.py --evaluate --chart first --coordinate .625 --Z .371 --output tmp/n1_first_query.json
```

For multiple queries in one process, construct
`CurrentOriginalN1InnerNumericBridge(function_owner, numerical_owner, hierarchy_owner)`
once. `function_owner` is the checked
`CurrentOriginalN1FullDefectLinearFunctions`; `numerical_owner` is the
checked `CurrentOriginalN1RegularSeriesFunction`; `hierarchy_owner` is its
checked `CurrentOriginalCoreHierarchySource` with the identical live
physical/source objects. Use `.evaluate(Z,coordinate,chart)` and
`.report(packet)`. Optional `reuse_packet` accepts a live packet from that
same regular owner only after its exact Z/chart/coordinate guards pass.
It does not accept serialized results or a packet for another point.

All six physical components have ordinary Z0–Z3 enclosures. Supported
first radial mixed Z0–Z2 rows and own mean/Q/V are retained. At the axis,
rho derivatives of the last two odd-in-x components are not reported.
Their regular x derivatives and the accepted even-component axis limits
remain available in the numerical packet.

## Exact coordinate and physical identities

The numerical first coordinate is a fraction q in [0,1]. The true
function graph uses physical phase s in [0,sc/4]:

```text
s = q sc/4
rho = 4 exp(hb s)
R = epsilon rho
```

The original exact binary sc and positive source hb are retained. A
rounded rho enclosure near 4 does not redefine the first collar as the
core endpoint. Physical phase s must not be passed as the numerical q.
The inner cutoff is a separate operation, with argument 4q-2.

Exact execution of the accepted production mathematical formulas proves
all 72 B0/B1 entries and all three nonzero forcing rows under

```text
S = diag(F0 Lambda, epsilon, epsilon, Pstar^2 epsilon, F0 Lambda, epsilon)
B0hat = S^-1 B0 S + S^-1 B1 S_Z
B1hat = S^-1 B1 S
ghat = S^-1 g
F0_Z/F0 = -Lambda gradient
```

The -2chi principal contribution and the remaining Psi gauge add to the
full original gauge once. Sixteen exact ordinary amplitude derivative
rows through Z3 verify the real production Taylor recurrence. The n1
physical recovery keeps the weight **1+delta**:

```text
mean1 = Uz1 + K1
Q1 = [(1-delta) Z Uz1 - (1+delta) Z K1 - (1-Z^2) mean1_Z]/L
V1 = R Q1
```

The normalized principal-kernel series and the raw series solve this same
regular equation with the same zero-axis data. Their complete tail
admissions remain distinct; the normalized Kchi is not inserted into the
raw operator without conjugation. Unmaterializable F0/Pstar/Lambda/hb
sectors stay in signed factored/log frames, including their combined logs.

## Same-source proof and immutable ownership

The bridge retains the full accepted true graph and its FunctionRefs.
Its regular numerical owner supplies intervals only. Guards bind the
exact source family and source pins, actual twenty-term common fixed
point, original j/logC and anchored H root, exact parameter recipes,
same independent fourteen-atom P0 definition and callable AST, and the
original prescribed-shear collar and own inlet/history integrals.

The old numerical core is a cache context using the admitted Bessel/model
part plus correction bounds. The actual core oracle returns nonlinear
finite rows plus complete tails. Equality of source boxes is not the
argument connecting those two consumers: the exact producer records,
same operator/model/domain definitions, actual-map correction recipe and
accepted unique fixed-point theorem supply the connection.

The two interval contexts have different precisions. Sigma and correction
therefore need not have identical rounded endpoints. Each consumer must
reproduce `sigma=j/500` and
`correction=scaled_map_size_upper/(1-scaled_map_Lipschitz_upper)` in its own
directed context, using the same saved records. P0 cache masses are also
enclosures, never selected analytic masses.

Both amplitude consumers use the same H polynomial
`j+(9-delta)Z/2-j Z^2-4Z^3`. Their certified unique root in [-j,0] is the
same mathematical anchor. The complete legacy root certificate is bound
to its saved producer, and the actual core anchor is separately checked.
No overlap of root disks or selection of provisional centers establishes
function equality.

Changed owner/function/source metadata, conflicting pins, unsupported
domains, foreign packets and reuse at the wrong physical point are
rejected. Saved archive and receipt source identities and pins must agree.

## Completed tasks

- [x] Retain the 53,398-node true graph and full five C1 defect/selection owner.
- [x] Restore the accepted regular numerical owner and actual common nonlinear oracle without rerunning ancestor producers/checkers.
- [x] Bind original source/model/norm records and independent P0 callable identity.
- [x] Bind both unique H-root anchors without equating interval centers.
- [x] Prove the complete physical matrix/forcing conjugation and ordinary amplitude jets.
- [x] Install exact core-rho / first-fraction query guards and phase metadata.
- [x] Execute and save original axis, core interior and actual first-collar queries.
- [x] Retain all six-row feedback, nonzero source errors and complete word tails.
- [x] Export ordinary Z3 and supported first-radial Z2 intervals with the correct axis exclusions.
- [x] Check the public accepted constructor and archive-only CLI.
- [x] Preserve explicit false repair/global/recursion/NS gates.

## Next implementation queue

- [x] Bind numerical even rho2 F1/Uz1/mean1 rows to the accepted true radial functions. Away from the axis, use the exact equation and d_x components; do not differentiate numerical seed centers.
- [x] Bind second-even axis limits from the genuine forcing/equation jets, preserving ordinary factorials and the distinction between odd x and even rho derivatives.
- [x] Wrap the original F1/Uz1 cutoff numerically with q=s/(sc/4), retaining hb^-1/hb^-2 sectors in radial derivatives and the exact flat keep/zero branches.
- [x] Recover cutoff own mean/M/K/Q/V from its complete core-plus-collar history. Zero velocities outside support do not make the cumulative mean zero.
- [x] Recover cutoff own pressure from original P0/core inlet and the full N1p source. Do not discard pressure merely because cut velocities vanish.
- [x] Obtain interval endpoint and signed cumulative ranges for all five true inner moment functions at Rin. Preserve physical drho/dr conversion and all epsilon powers.
- [ ] Transfer complete lower/preRp/pulse signed pressure bounds to the new genuine full defect functions, retaining all 17 active charts and all finite endpoint terms.
- [ ] Transfer value/Z1 quantitative moment bounds through the actual stable five-by-five response inverse, keeping nonzero Astar, mu and the original physical normalizations.
- [ ] Extend genuine full m13/m15 derivative functions from C1 to the needed C3 order by differentiating original field, integral and implicit definitions.
- [ ] Differentiate the selected coefficient system, including matrix derivatives and Astar_Z, instead of using grid slopes or derivatives of error bounds.
- [ ] Install alpha/beta derivative functions needed for repaired Q/V derivatives and later temporal source terms.
- [ ] Construct repaired own M1 from cutoff plus compact axial bump primitives, retaining each primitive's incoming history and tail.
- [ ] Construct repaired own K1=mean1-Uz1 using the same physical cumulative M1.
- [ ] Construct repaired Q1 and V1 with the n1 weight 1+delta; do not reuse the leading n0 recovery weight.
- [ ] Construct repaired own P1 from the full original signed pressure integrand and angular repair primitives, preserving P0 and all exterior histories.
- [ ] State the regularity of each repaired component honestly: selected C1 bumps alone do not provide C3 Q/V or C5 coefficient bounds.
- [ ] Prove all five repaired terminal identities as functions of Z from the accepted response theorem and actual primitives, including value/derivative scope actually installed.
- [ ] Transfer axis, cutoff, compact-support and annular join identities to those repaired own fields.
- [ ] Transfer coefficient and field size/derivative bounds on original cells and holomorphic tubes, including all nonzero numerical and tail errors.
- [ ] Close positive-order five-moment repair only after the actual own-field construction, terminal identities and required quantitative bounds exist.
- [ ] Resume complete heat/preheat and global pressure/energy matching from the repaired family's own terminal data.
- [ ] Derive and evaluate admissible divergence-form stress and flat remainder separately on every required region.
- [ ] Implement genuine n>=2 source recursion with the correct n-dependent recovery and its own five-moment repair, using the completed prior coefficients.
- [ ] Bound finite-order remainders and implement the required smooth summation; coordinate-scaled copies are not temporal recursion.
- [ ] Construct the two-family oscillatory corrections and their mean terms from the admissible stress data.
- [ ] Evaluate averaged quadratic stress cancellation with retained oscillatory errors and flat remainder.
- [ ] Export the resulting 3D time-dependent Cartesian velocity/pressure with a fixed or restricted forcing definition independent of its residual.
- [ ] Validate the final corrected Cartesian momentum and divergence equations independently before applying the full residual targets.
- [ ] Measure time-dependent radial/axial core widths, aspect ratio, swirl/vorticity growth, material-line winding, energy and actual scale-recursion exponents from the resulting fields.

Each task should update the corresponding receipt, mark its own completed
checkbox and leave the next concrete action. Keep the full original goal
intact. Reuse accepted source owners and saved results; do not begin with
another historical pressure-grid run.
