# Genuine n=1: computed pressure/derivative feedback and analytic series tail

The first coupled regular iteration is now computed on the canonical core and the actual first prescribed-shear collar. It includes the unknown-field pressure feedback 2F0F1 and genuine B1 partial_Z W terms. The previous initial integral has been extended internally to whole radial endpoint cells, including cells that touch the axis. The long-term goal remains **ACTIVE / INCOMPLETE**.

This supersedes the next-action portion of [the initial-integral checkpoint](CURRENT_ORIGINAL_N1_INITIAL_INTEGRAL_2026_10_10.md). The source family, actual mean, original pressure, formal hb and common analytic domain remain unchanged. Do not rebuild the old ancestor validation chain.

## Artifacts and scope

Source quartet in `experiments/root_st073/`:

```text
lei_ren_part1_paper_compliant_current_original_n1_coupled_first_step.py
lei_ren_part1_paper_compliant_current_original_n1_coupled_first_step.json.gz
lei_ren_part1_paper_compliant_current_original_n1_coupled_first_step_check.py
lei_ren_part1_paper_compliant_current_original_n1_coupled_first_step_check.json
```

Accepted gate: `current_original_n1_first_coupled_regular_iteration_computed`.

The API computes all six state values for a supported inner query: F1, Uz1, K1, P1, d_x F1 and d_x Uz1. Initial values have first Z derivatives; the coupled prefix currently has Z order zero. Queries use rho in [0,4] or actual collar phase s in [0,1/2], with rho=4 exp(hb*s). It retains the actual cumulative mean and the independent pressure datum. It never uses the analytic-core extension as the actual collar.

The public report separates:

- `computed_prefix`: two genuinely computed preconditioned terms;
- `B0_feedback`: the off-kinematic matrix feedback, including pressure;
- `B1_derivative_feedback`: actual derivatives of the initial function;
- `canonical_inner_series_value_enclosures`: the computed prefix plus a uniform absolute bound for every omitted coupled series term.

The last values enclose the solution series of the **canonical unrepaired inner problem**. They do not certify the completed-leading support relation, higher-Z coefficient derivatives, repaired n=1 background, or useful relative error in every component. Inherited full-inner/coefficient-solved/positive-order-repair gates remain false. The older controlled-prefix gate is not promoted: its full source/integration and derivative budget remains open, and no original unpreconditioned 64-term prefix has been computed.

## Actual coupled computation

Let Akin have rows (1,5)=1, (2,6)=1 and (3,6)=-1. With the regular diagonal inverse G and J=G Akin,

```text
J^2=0
Wbase=(I+J)Gg
H=(I+J)G[(B0-Akin)+B1 partial_Z]
W=Wbase+H W
Wprefix=Wbase+H Wbase.
```

The exact nilpotent resummation keeps the kinematic rows while exposing the small remaining operator. It does not discard nonadjacent derivative words. B0 feedback and B1 feedback are integrated separately before being combined.

For a core endpoint X that ranges over an entire interval, use s=X*t. The ratio s/X is exactly t, including at X=0 by continuous extension. Thus the initial kernels do not divide by an interval that contains zero. Sixteen nested entire-core cells enclose Wbase and its first Z derivative for each outer core cell.

Eight outer core cells integrate the feedback. Two actual outer phase cells each use four nested collar cells and the true measure

```text
dx = sqrt(4 exp(hb*s))/2 * hb * ds.
```

Width factors combine before enclosure. All source bases `(hb,F0_anchor,Pstar_squared,Lambda)` remain formal; no source midpoint or numerical width representative is selected. All integrations use entire-cell enclosures rather than endpoint interpolation.

Published samples at Z=371/1000 are rho=4 and actual phase s=1/2. They contain 21 and 105 combined prefix source sectors respectively. The B1 pressure component is exactly zero by the operator support; the B0 pressure component retains the nonzero 2F0F1 contribution. Other B1 terms are computed explicitly, never replaced by zero because their scales are small.

## New tail proof

The old unpreconditioned 64-term bound cannot certify this differently organized prefix. The new bound treats all k derivative blocks as possible.

The middle-tube off-kinematic norm Cp is computed from the explicit absolute bounds of rows four, five and six, omitting the three Akin entries. It is not inferred by subtracting one from the previous full B0 norm. C1 is bounded on the same eta/2 tube. G has nonnegative diagonal kernels bounded by one; Akin has absolute row norm one and its map annihilates through every diagonal kernel. Therefore ||J||<=a, ||I+J||<=1+a, where a=sqrt(4.1).

Wbase is holomorphic on radius r0=eta/4, with Mbase<=(1+a)M0. For a word of length k, divide the available loss L=eta/8 into k gaps. Cauchy applied to the complete preceding product costs k/L at each derivative step. Derivatives of earlier matrix factors are thereby included. Radial Volterra ordering gives at least k simplex integrations; every optional J insertion adds one more. Summing these choices yields

```text
||H^k Wbase||_(eta/8)
 <= Mbase * [a*(1+a)*(Cp+k*C1/L)]^k / k!.
```

Since k! >= (k/e)^k, for k>=N+1 use

```text
beta_N = e*a*(1+a)*(Cp/(N+1)+C1/L)
sum_(k>N) ||H^k Wbase|| <= Mbase*beta_N^(N+1)/(1-beta_N).
```

The positive guard beta_N<1 is required. For the computed N=1 prefix, beta_1 is about 1.90e-408906090034569530 and the uniform omitted W-norm is about 1.16e-1022265225086424019 under the original source scales. These are absolute bounds for this inner coefficient problem. They are **not full Navier–Stokes residuals or guarantees of relative accuracy** for F1/Uz1. The prefix's source and cell-enclosure uncertainty is still combined; the remainder bound is added separately.

## Focused checks and use

Independent checks cover 84 exact normalized-axis/range and signed coupled polynomial rows, genuine Z-derivative feedback, the pressure unknown-field contribution, diagonal-kernel nilpotency, and 120 factorial/geometric tail comparisons. Prefix and full inner series traces vanish exactly at the axis. The accepted point initial adapter and normalized range adapter agree in their source enclosures. Owner/packet/range/chart/partition guards retain the scope and detect live prefix mutation.

The single reused read-only worker remains **GPT-5.6 Luna / max**. It reviewed the nested-radius proof and confirmed its validity with the explicit row-norm, Volterra ordering and tube guards recorded above. No new worker or Astra child was spawned.

```python
from lei_ren_part1_paper_compliant_current_original_n1_coupled_first_step import CurrentOriginalN1CoupledFirstStep

coupled = CurrentOriginalN1CoupledFirstStep(accepted_initial_integral_owner)
packet = coupled.evaluate('.371', '.5', 'first')
view = coupled.report(packet)
assert view['current_original_n1_first_coupled_regular_iteration_computed']
assert view['canonical_inner_series_values_include_remaining_uniform_absolute_tail']
assert not view['n1_full_inner_interval_solution_certified']
```

## Next executable tasks

- [x] **CURRENT-RP-N1-WHOLE-CELL-INITIAL-FUNCTION** Enclose Wbase and its first Z derivative over entire core/collar endpoint cells. Use normalized t kernels at the axis and actual phase integration in the collar. Completion covers this range input for the first feedback only.
- [x] **CURRENT-RP-N1-FIRST-COUPLED-FEEDBACK** Compute H Wbase separately for B0 and B1, retaining actual pressure feedback, own means and source sectors. Publish genuine core and actual-collar samples with both contributions.
- [x] **CURRENT-RP-N1-PRECONDITIONED-ANALYTIC-TAIL** Derive Cp row-wise, reserve nested Cauchy radii, include optional J integrations, prove the geometric guard and enclose every omitted H-series term. Keep the absolute W-norm tail separate from source/integration uncertainty and component relative accuracy.
- [ ] **CURRENT-RP-N1-COMPLETED-LEADING-SUPPORT-GUARD** Bind the original completed-leading dispatch, Rm and r_minus to this source family and prove all later modifications lie outside Rin. Retain strict formal radial inequalities. The canonical prefix alone does not close Assumption 14.1.
- [ ] **CURRENT-RP-N1-Z-DERIVATIVE-BUDGET-EXTENSION** Extend genuine leading mixed rows from total order three to the available total order four, forcing/g from Z1 to Z2 and matrix rows from Z2 to Z3, with matching derivative guards. Propagate Wbase Z2 and the first coupled prefix Z1. Verify ordinary F0 derivatives and own-mean derivatives; never add zero padding.
- [ ] **CURRENT-RP-N1-RADIAL-VELOCITY-RECOVERY** With certified Uz1/K1 Z derivatives, recover Q1=[(1-delta)Z Uz1-(1+delta)Z K1-(1-Z^2)partial_Z(Uz1+K1)]/L and V1=R Q1. Derive a compatible derivative tail by an additional Cauchy margin. Verify the exact n=1 mean/axis relation and pressure, then provide the complete inner coefficient field API.
- [ ] **CURRENT-RP-N1-COMPONENT-ERROR-BUDGET** Separate source uncertainty from integration error and refine normalized cell widths. Preserve amplitude-sensitive sectors; a uniform unweighted norm tail may swamp a much smaller individual swirl component. Supply field/derivative bounds with useful widths and independently assembled equation diagnostics before promoting coefficient-solved/full-inner gates.
- [ ] **CURRENT-RP-N1-FIVE-MOMENT-REPAIR** Compute five whole-Z defects of the solved n=1 field, apply its own original bump controls and pressure-compatible repair, truncate streamfunction/vector potential before curl, and verify joins/derivatives and finite-energy tail. Leading-order repair does not close this task.
- [ ] **CURRENT-RP-N2-AFTER-REPAIRED-N1** Build the true n=2 products and n-dependent recovery from solved/repaired n=1 on the same interval, then solve and repair independently. Continue higher temporal orders only from admissible lower orders.

Regional cone/flatness, physical core-width/material-winding measurements, finite-energy bounds, global/axis physical coverage and oscillatory stress cancellation remain in the full objective.
