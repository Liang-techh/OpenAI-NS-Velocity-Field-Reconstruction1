# True cut n1 own fields and five inner numerical moments

The numerical interface now follows the accepted true cutoff graph. It
cuts F1/Uz1, then encloses their own M1, mean1, K1, Q1, V1 and P1 with the
complete incoming history. All five signed inner cumulative functions,
including their Rin endpoint ranges, are exposed in physical factored/log
units. Numerical centers are not selected graph functions.

The full original goal remains **ACTIVE / INCOMPLETE**. This stage is the
cut, unrepaired n1 family. The reserved alpha/beta repair and two exterior
pressure contributions are separate remaining operations.

## Restore and query

The checked quartet under `experiments/root_st073` is
`lei_ren_part1_paper_compliant_current_original_n1_cutoff_numeric_bridge`
(`.py`, `.json.gz`, `_check.py`, `_check.json`). Its local gate is
`current_original_true_cutoff_n1_own_fields_and_five_inner_moment_enclosures_installed`.

```powershell
python scripts/run_original_n1_cutoff_numeric_bridge.py
python scripts/run_original_n1_cutoff_numeric_bridge.py --evaluate --chart first --coordinate 1 --Z .371 --output tmp/n1_cut_Rin.json
```

The first command restores the saved archive and checks its source pins.
The second restores accepted numerical owners once and executes a fresh
original query. For multiple queries, keep those owners alive and call
`CurrentOriginalN1CutoffNumericBridge.evaluate()/report()`.

Saved original cases at Z=.371 cover the axis, first inlet, cutoff
transition q=.625, velocity-zero boundary q=.75 and inner exit q=1.

| Case | Coordinate | Passing inner value/Z1 rows | Largest relative inner-moment interval width |
| --- | --- | --- | --- |
| Axis | core rho=0 | Exact-zero branch | Inapplicable |
| First inlet | first q=0, rho=4 | 10/10 | 2.68475822597e-28 |
| Cut transition | first q=.625 | 10/10 | 9.03566684437e-10 |
| Velocity-zero boundary | first q=.75 | 10/10 | 9.03566684437e-10 |
| Inner exit Rin | first q=1 | 10/10 | 9.03566684437e-10 |

These are pointwise relative widths of the accepted signed inner
moment value/first-Z intervals. They include source errors, range-cell
integration errors and complete feedback tails. They are not terminal
zero errors, complete pressure-moment errors or an NS residual.

Core coordinates are rho in [0,4]; first coordinates are q in [0,1], with
s=q sc/4 and rho=4 exp(hb s). In particular, **first q=0 is rho=4, not the
axis**. Exact sc is compared as an integer mantissa and binary exponent;
the code never allocates an integer of size 2 raised to that microscopic
scale's astronomical exponent.

The `cells` parameter partitions the cutoff-change integral. The accepted
base cumulative interface retains its default core/collar/source range
partitions. These are directed cell ranges, not point quadrature with
unaccounted errors. Increasing `cells` does not by itself certify whole-Z
or higher-derivative accuracy.

## Physical recovery and pressure

Only the two velocities are cut directly. The primitive and radial
velocity use

```text
M1(R) = integral_0^R Uz1_cut dR
mean1 = M1/R
K1 = mean1-Uz1_cut
Q1 = [2Z Uz1_cut-(1+delta)Z mean1-(1-Z^2)mean1_Z]/L
V1 = R Q1
```

Q1 is evaluated from the analytic recovery formula, including at the
axis. It is never obtained by dividing a boxed V1 by R=0. The true regular
axis values of F1/Uz1/mean1/K1 vanish; the n1 pressure also has the accepted
zero-axis O(x^2) condition. The independent leading P0 remains in the
leading source, rather than being added again to the n1 pressure.

After q=.75, the cut velocity and its supported radial jets are zero.
M1 keeps its nonzero terminal history, mean1=M1/R and the radial recovery
remain present. Pressure continues to use the unchanged leading N1p.
The exact same terminal M packet is used from the zero boundary to Rin.

The interface supplies F1/Uz1 rho0/rho1/rho2 with ordinary Z0/Z1,
M1/mean1/K1 through ordinary Z2, Q1/V1 through Z1, supported even mean/K
rho1/rho2 with Z1, and P1/P1_rho with Z1. It does not add pressure rho2,
unsupported odd axis derivatives, or a uniform C3/C5 numerical accuracy
claim. All source errors and the complete six-row Volterra tail remain.

## Five signed inner moments

The true signed densities per drho are

```text
m11: epsilon U1_cut
m12: 2 epsilon^2 rho F1_cut
m13: epsilon (2F0 F1_cut+N1p)
m14: 2 epsilon^2 rho (U0 F1_cut+U1_cut F0)
m15: 2 epsilon U0 U1_cut-epsilon^2 rho (2F0 F1_cut+N1p)
```

The old numerical cutoff increments do not contain N1p because it is
unchanged by cutting the velocities. The base cumulatives already contain
that pressure source. Exact production-formula checks verify base plus
increment equals all five true densities, with the m15 minus applied
once. Integration by parts for m15 preserves its signed pressure boundary
term. The two FTC cumulants use the actual own fields directly: m11=M1
and m13=P1.

At Rin, m11/m12/m14 cover the full velocity support. The pressure moments
still need

```text
m13_full = m13_cut(Rin)+integral_Rin^infinity N1p dR
m15_full = m15_cut(Rin)-integral_Rin^infinity R N1p dR
```

Their lower/preRp/pulse true function bodies already exist in the full
defect graph. This numerical bridge does not add or zero those exterior
integrals and does not establish terminal five-moment closure.

## Source and boundary guards

The accepted inner numerical bridge binds the true twenty-term core,
original anchored F0, independent P0, exact parameters and complete raw
equation. This cutoff bridge additionally binds the same live regular/
radial/cumulative owners, exact selected first phase, original sigma
value/derivative providers and the ordinary sigma'' factor
2*Taylor[2]. The original hb^-1/hb^-2 sectors and physical epsilon powers
are retained.

At first phase zero, the same regular core solution reaches rho=4. The
own mean/pressure histories have equal zero integral endpoints and preserve
their actual core inlet. Cumulative cutoff increments there are exactly
zero. Bounds may have different interval representations; interval overlap
does not prove this mathematical seam identity.

Wrong domains, interval cumulative endpoints, foreign packets, reuse at
another point or with a different cutoff partition, changed source
metadata and conflicting pins are rejected. Saved archive/receipt source
identities and exact byte pins must agree.

## Completed tasks

- [x] Bind supported even rho2/axis enclosures to the true regular n1 definitions.
- [x] Bind the original cutoff and hb^-1/hb^-2 ordinary radial derivatives.
- [x] Recover numerical own M/mean/K/Q/V/P without zeroing cumulative histories.
- [x] Expose all five genuine signed inner cumulative and Rin ranges with their retained errors.
- [x] Verify first q=0 uses the core rho4 inlet and zero first-history increments.
- [x] Execute original axis/inlet/transition/zero/Rin queries; all ten nonaxis inner value/Z1 budget rows pass.
- [x] Keep complete exterior pressure, reserved repair, C3/C5 accuracy, n>=2 and NS gates open.

## Next tasks

- [x] Bind true lower-collar m13/m15 C1 numerical enclosures and compose the actual Rin history through rminus. See [the lower numerical bridge](CURRENT_ORIGINAL_N1_LOWER_NUMERIC_BRIDGE_2026_10_13.md); the remaining preRp/pulse pressure transfer stays open.
- [ ] Attach the complete exterior pressure enclosers to the true lower/preRp/pulse C1 functions using the same source, original radii, native histories and actual finite-N fields.
- [ ] Compose numerical ranges of all five complete true C1 defect functions, retaining both nonzero pressure tails and all errors.
- [ ] Transfer these ranges through the physical reserved response inverse to the actual selected alpha/beta functions; do not choose interval centers as coefficients.
- [ ] Extend genuine pressure-defect and coefficient derivative functions to the orders needed for repaired Q/V and later temporal sources.
- [ ] Construct repaired own M/K/Q/V/P from cutoff plus compact axial/angular bump primitives with complete incoming histories.
- [ ] Prove actual repaired terminal five zeros, axis/FTC/join identities and required coefficient/field bounds before closing the n1 repair gate.
- [ ] Complete whole-region energy, heat/preheat pressure matching, admissible stress and flat remainder.
- [ ] Implement genuine n>=2 temporal recursion, its own moment repairs, remainder control and smooth summation.
- [ ] Build oscillatory stress cancellation, export 3D time-dependent corrected fields and independently validate the final restricted-forcing NS residual.

Reuse live numerical owners and saved archives. Do not start by rerunning
historical pressure grids or by dropping nonzero history/error sectors.
