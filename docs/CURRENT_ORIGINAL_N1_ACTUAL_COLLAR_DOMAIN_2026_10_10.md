> Subsequent progress: [the first coupled n=1 step](CURRENT_ORIGINAL_N1_COUPLED_FIRST_STEP_2026_10_11.md) now computes B0 pressure feedback and B1 derivative feedback, with a new nested-radius tail enclosing the remaining canonical inner series. Higher-Z/radial-velocity recovery, completed-leading support binding and positive-order repair remain open.

> Subsequent progress: [actual operator and computed initial regular integrals](CURRENT_ORIGINAL_N1_INITIAL_INTEGRAL_2026_10_10.md) now supplies source-owned actual inputs and Wbase=(I+J)Gg enclosures. The older bounds/axis scope below remains accepted. A full coupled n=1 field is still open.

# Genuine n=1: actual inner collar, common analytic bounds and integral tail

The next prerequisite for genuine temporal recursion is installed: same-source complex bounds for the analytic core and the **actual** first prescribed-shear collar, including the width-integrated n=1 forcing and a regular Neumann-series tail majorant. The existing finite n=1 axis solution remains valid. This checkpoint does not compute the full n=1 iteration prefix or deliver a full-interval positive-order velocity field. The long-term goal remains **ACTIVE / INCOMPLETE**.

Source quartet: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_analytic_domain` with `.py`, `.json.gz`, `_check.py`, `_check.json`. Its accepted gate is `current_original_n1_core_and_actual_first_collar_analytic_bounds_installed`.

## Correct inner source

The old core callback is jointly analytic on rho in [0,4.1]. For rho > 4 this is an analytic extension of the core. It is the comparison field in the first collar; the actual prescribed-shear bridge differs from it for every positive collar phase. Therefore `n1.system('4.1', Z)` supplies the analytic-core operator, not the actual completed-leading operator at that radius. Its axis coefficients and core operator at rho <= 4 are unaffected.

The new bounds use the actual source on

```text
0 <= rho <= 4                       analytic core
rho = 4 exp(hb*s), 0 <= s <= 1/2    actual prescribed-shear bridge
Rin = Ra exp(hb/2)
Rkeep = Ra exp(hb/4)
```

The exact positive hb is retained via its original log enclosure. These formal radii satisfy Ra < Rkeep < Rin; they must not be rounded to Ra. The numerical width cap is an error bound, never a selected field parameter. The same original K1 ledger binds hb = epsilon_b. Comparison alpha equals one on this prefix, while actual chi = 1-(1-epsilon_b)*sigma(s).

The new source is scoped to the canonical core/first-collar prefix. Before using this as the completed leading input in Assumption 14.1, independently bind the actual later moment-repair and shear-modification supports outside Rin, or choose a smaller formal Rin below their verified support endpoints. That completed-leading support gate remains false. No old closed pressure 4C/first-interface task is reopened.

## Common complex domain

The inherited admitted fixed-point norms Bphi and Bpsi belong to the same actual source, pressure datum and selected Cstar family. Paper (8.30) gives joint normal convergence for |rho| < 5 and local Z disks of radius h/2; its positive majorant is at most 4 times the Xh norm. The code uses radial Cauchy radius 1/2 and Z radius h/4.

Choose

```text
eta = a positive lower enclosure of h*phi_floor / [128*(1+Bphi)]
```

This radius is chosen for an error domain, not as a construction parameter or representative source value. It is below h/8. The core's complex perturbation from the real rectangle is at most 16 Bphi eta/h, which is below phi_floor/8. Thus |Phi_core| >= 7 phi_floor/8, and L has its own strictly positive modulus lower bound. The axis amplitude is the actual exponential exp(-selected_logCstar-Lambda G(Z)); the inherited selected-Cstar guard supplies its complex upper bound without materializing the extreme source amplitude.

The tubes reserve separate Cauchy margins:

| Use | Radius around the real Z interval |
| --- | --- |
| Core and actual-collar radial bounds | eta |
| Regular matrix and Q radial bounds | eta/2 |
| Known source and initial Gg | eta/4 |
| Recursive solution/tail estimate | eta/8 |

Q uses the actual radial mean of Uz, recovered without singular division at the axis. In the bridge, rho mean_rho = Uz-mean gives its radial derivative bounds from the actual velocity. Pressure retains the independent original P0 and the true integral of F squared. All fourteen positive pressure atoms and the exact flatten function remain in the source graph. The flatten factor is the fractional-log source

```text
2^(-2*sigma) * (1+Z^2)^(-2*(1-sigma)), 0 <= sigma <= 1.
```

Re(1+Z^2) >= 1-|Im Z|^2 > 0 fixes an analytic logarithm and gives the pressure modulus bound. This uses the actual outer definition rather than an older helper's explanatory rational-kernel comment. Original flatten branch premises and callable/AST bindings are checked from the original source records.

## Actual-collar radial derivatives

Let A = Phi_core_rho/Phi_core and r = Phi_bridge/Phi_core. The exact reduced-core equations and actual bridge ODE give

```text
Phi_bridge_rho = chi * Phi_bridge * A
Phi_bridge_rhorho = Phi_bridge * (chi_rho*A + chi^2*A^2 + chi*A_rho)
Uz_bridge_rho = chi * r * Uz_core_rho
Uz_bridge_rhorho = r * [chi_rho*Uz_core_rho
                        + chi*(chi-1)*A*Uz_core_rho
                        + chi*Uz_core_rhorho]
chi_rho = -(1-epsilon_b)*sigma_prime(s)/(hb*rho)
```

The exponential bridge has no zeros on the common domain. Its upper bound follows directly from the defining ODE integral; it is not a bridge/core closeness claim. The positive derivative recurrence retains inverse hb powers through radial order four. The same smooth radial ODE proves holomorphy on the same Z neighborhood for every fixed radial derivative; explicit computed majorants currently cover orders zero through four.

Bounds are positive Laurent polynomials sum c_p hb^p. A negative width power cannot be numerically evaluated. Ordinary R derivatives multiply the rho derivative by Lambda^i, with the exact scale kept separate. Ordinary Z derivatives use the appropriate reserved Cauchy margin.

## Width cancellation and recursive tail

Every actual-collar forcing row in the n=1 regular system has at most one inverse hb. The real integration measure is

```text
dx = (sqrt(rho)/2) hb ds,   0 <= s <= 1/2.
```

Multiplying before taking an upper bound cancels the inverse width exactly. The initial Gg bound is the sum of a core source integral bound and the collar's width-weighted integral bound. Remaining positive-width products are enclosed only after comparing their exact p*log(hb)+log(coefficient) against a safe numerical product cap. The original positive-width terms remain in the exported record. This avoids letting the generic 1e-50000 width cap swamp a much smaller true cancelled source.

The regular matrix has no inverse width. The code supplies infinity row-norm bounds C0 and C1, and inherits the exact nonadjacent B1 support: adjacent derivative blocks vanish through the diagonal Volterra kernel. For k further operator applications to Gg,

```text
Tk <= M0 (2*a*C0)^k/k! * max(1, alpha*ceil(k/2))^ceil(k/2)
alpha = C1/(C0*DeltaZ), a = sqrt(4.1), DeltaZ = eta/8.
```

The source is integrated first; k remaining nested radial integrations give the k! simplex denominator. This weighted L1 estimate is valid even though the collar's pointwise source supremum contains an inverse hb. The majorant for the tail after indices zero through 64 is about **3.95e-52 times M0**. This is a bound on omitted operator terms, conditional on correctly computing that prefix on the stated tubes. It is not the error of a computed field: no iteration terms or quadrature have yet been delivered.

Focused checks cover four independently differentiated ODE/width-measure identities, the scalar tail-ratio switch and independent high-precision tail sums, containment of accepted core matrix norms, original nozero/tube guards, and source/width/family/context ownership. One existing read-only worker, **GPT-5.6 Luna / max**, reviewed the source distinction and mathematical bounds. No new worker or Astra child was spawned.

## Live use

```python
from lei_ren_part1_paper_compliant_current_original_n1_inner_analytic_domain import CurrentOriginalN1InnerAnalyticDomain

domain = CurrentOriginalN1InnerAnalyticDomain(accepted_n1_owner)
packet = domain.domain()
view = domain.report(packet)
assert view['current_original_n1_core_and_actual_first_collar_analytic_bounds_installed']
assert not view['n1_Volterra_iteration_prefix_computed']
assert not view['n1_full_inner_interval_solution_certified']
```

The report contains bounds for implicit original functions. JSON cannot replace the live owners or hydrate selected source values.

## Next implementation tasks

- [x] **CURRENT-RP-N1-ACTUAL-COLLAR-DOMAIN-BOUNDS** Bind the same admitted analytic core, actual first collar, original width/cutoff relation, exact P0, Cstar guard and nozero domains. Deliver explicit radial majorants through order four and reserved Z Cauchy margins. Completion concerns this prefix-bound adapter only.
- [x] **CURRENT-RP-N1-COLLAR-SOURCE-MEASURE-CANCELLATION** Retain the inverse hb source sectors, multiply by dx before enclosure, and prove remaining positive-width product caps from their true source logs. Do not replace hb by the width cap.
- [x] **CURRENT-RP-N1-WEIGHTED-VOLTERRA-TAIL-MAJORANT** Supply C0/C1, weighted initial Gg and a factorial/Cauchy tail majorant. Keep its acceptance distinct from a computed quadrature/iteration prefix.
- [ ] **CURRENT-RP-N1-COMPLETED-LEADING-SUPPORT-GUARD** Read the actual source definitions of Rm and r_minus and bind them to the original completed-leading family. Prove Rin below both and the analytic/comparison limit, retaining strict formal inequalities despite extreme widths. If necessary choose a smaller positive formal phase. Close Assumption 14.1 only after the actual completed-leading dispatch and supports are bound; a comparison-core extension is insufficient.
- [x] **CURRENT-RP-N1-ACTUAL-COLLAR-OPERATOR-FUNCTIONS** Delivered as typed actual operator input enclosures in [the subsequent initial-integral checkpoint](CURRENT_ORIGINAL_N1_INITIAL_INTEGRAL_2026_10_10.md). Canonical actual field, own mean, independent pressure, mixed rho rows, genuine forcing and B0/B1/g retain source sectors and ordinary amplitude derivatives. Matrix Z2 and forcing Z1 are explicit; unavailable higher derivatives are rejected. Exact point moment-history reconstruction remains false.
- [ ] **CURRENT-RP-N1-CONTROLLED-Gg-PREFIX** Bounded progress: Wbase=(I+J)Gg now computes all six initial components and first Z derivatives with entire-cell intervals, width cancellation, zero-axis traces and exact integral core/inlet join. See the subsequent checkpoint. Complete this larger task by separating source/integration error, refining useful widths and supplying higher-Z/common analytic function control; the combined cell/source enclosure is not yet that complete error budget.
- [ ] **CURRENT-RP-N1-CONTROLLED-ITERATION-PREFIX** Compute the coupled integral iterations, including B1 partial_Z W or an explicit controlled bound for omitted derivative-block words. Preserve pressure feedback. Never label an iteration as zero because its scale underflows. Provide the actual retained terms, their quadrature/source uncertainty and compatible analytic domains; then add the existing tail majorant.
- [ ] **CURRENT-RP-N1-SOLVED-INNER-FIELD** Issue a callable F1/Uz1/K1/P1/V1 source on one genuine common interval, with interval and derivative enclosures. Verify the original equations, n=1 average/radial recovery, pressure and zero-axis data independently. Only then set the n=1 coefficient-solved/full-inner-interval gates.
- [ ] **CURRENT-RP-N1-FIVE-MOMENT-REPAIR** Compute five actual Z-function defects of the solved n=1 field; repair using original bump controls and controlled conditioning. Truncate streamfunction/vector potential before curl. Verify whole-Z terminal moments, pressure, joins and derivatives. Existing leading repair does not close positive-order repair.
- [ ] **CURRENT-RP-N2-AFTER-REPAIRED-N1** Build true n=2 sums and weights from solved/repaired lower temporal orders on the same interval. Radial Taylor degree two is not hierarchy n=2. Repeat solve, bounds and per-order repair.

Regional cone/flatness, actual core-width/material-winding measurements, physical finite-energy bounds and oscillatory stress cancellation remain part of the long-term goal. See [the physical/task checkpoint](CURRENT_ORIGINAL_RP_MULTITIME_REMAINDER_2026_10_10.md); this analytic-bound step does not promote those separate gates.
