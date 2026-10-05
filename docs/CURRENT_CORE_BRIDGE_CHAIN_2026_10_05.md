# Current core/axis and finite-width bridge integration - 2026-10-05

The active reconstruction goal remains the coupled, paper-faithful background followed by admissible stress, real order-dependent temporal recovery and oscillatory correction. Source-region counts are coverage counts, not a percentage of that goal.

## Accepted current core and axis

`CurrentCorePhysicalAssembly` adds `core` to the checked current29 heat/annular physical chain. Its provider is the existing `dispatch.anchor.patch.core` context adapter, whose `original` is the very same finite-width bridge's analytic core. No independent legacy core is constructed. The original nonlinear analytic derivative bounds, compatible pressure datum and nonsingular Cartesian map are retained.

Axis is a mode: `evaluate('core', Z, 0, axis=True)`. It is not a new chart owner. The original core domain is rho in [0,4], Z in [-1,1]; rho4.1 is reserved for analytic continuation inside bridge smoothing. Nonzero axis requests and unsupported chart labels are rejected.

Four whole-Z views cover whole core, axis inlet, core exit and axis mode. The scoped check passes300 profile source rows and1,008 spatial/time contributions. Current physical ownership is30, with6,486 regular source contributions,252 supplemental axis contributions and the previously retained216 pulse-gap contributions. Requested log(tau) is propagated to physical bounds, including a fresh unsaved Z. The existing independent nonsingular core fixture is retained by its checked hash instead of rerun.

Reproduce only the changed stage:

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage currentcorephysical
```

Receipt: `experiments/root_st073/lei_ren_part1_paper_compliant_current_core_physical_assembly_check.json`. Gates: `current_core_axis_cartesian_spatial4_time1_certified`, `current_thirty_source_physical_ownership_certified`.

## Current bridge implementation

`CurrentActualBridgeMixedC4` takes the same current graph's `history.upstream` finite-width integrals. The three original domains are first s in [0,1], second s in [1,2], and frozen macro fraction in [0,1]. One provider serves all three. `ActualR100BridgeAdapter.actual()` remains terminal-only and is never used to recover a microscopic coordinate from rounded theta.

A local AST replay replaces exactly the original bridge method's three parent acquisition calls. Every derivative, Cartesian-ready ledger formula, domain guard and source-width cancellation remains unchanged. Actual F/V, Q, pressure and six own histories come directly from `upstream.packet(Z,value,chart)`. Known comparison histories are acquired separately through axial order6: point micro/macro calls, or the admitted uniform micro cover for an interval. Actual histories remain axial5, with Q axial4. No comparison history is substituted for an actual moment.

The original source controls are bound before reuse of their operator evidence. The bridge keeps the nonlinear angular exponential, quadratic Volterra feedback, formal hb, Pstar-squared and F0-squared factors. An interval cap bounds an exact positive width; it does not define a selected source width. Graded mixed derivatives use total order at most4.

The new scoped receipt passes1,215 mixed4 rows and their final logR ledgers over nine whole-Z domain/endpoint views. It binds actual micro feedback, macro quadratic/Volterra histories and pressure dressing, retains the unchanged operator fixtures by hash, and exercises a fresh unsaved Z on every chart. Checked runtime loading also passes. Actual mixed4 feedback is now installed in this adapter's scope; the underlying history provider's older labels remain unchanged. The stage is:

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage currentbridgemixed
```

The new receipt is `lei_ren_part1_paper_compliant_current_actual_bridge_mixed_C4_check.json`, with gate `current_actual_three_bridge_charts_mixed4_available`. The legacy bridge receipt supplies unchanged operator evidence only. Generator mode disables only this new gate; its already checked current30 parent is required. There are now33 current profile-source owners across separately admitted core, bridge and downstream adapters. Current Cartesian/time physical ownership remains30 until a composed bridge physical receipt passes.

## Next executable work and limits

1. Prove current core/first, first/second, second/macro and R100/switch functional interfaces. Derive equality from source ODEs, identical complete histories and original flat endpoints; overlap of intervals is only a diagnostic. Reuse the now admitted bridge derivative/ledger source coverage.
2. Feed the admitted bridge ledger to the unchanged physical mapper. Reuse current30 physical evidence by hash and add the three bridge owners. Axis remains one core mode. Full33 physical ownership needs its own new scoped receipt.
3. Rebuild current heat pressure/stress companions using the checked current heat object. Global cone, tensor, flatness and required-domain energy remain separate obligations.
4. Recover actual nonlinear point histories with uncertainty and nontriviality preserved. Directed whole-domain source enclosures are not chosen point values for production [u,v,w,p].
5. Implement distinct n=1 and n>=2 coefficient recovery and per-order repairs once the leading/moment/stress dependencies hold. Spatial chart coverage and radial analytic recurrence do not establish this temporal recursion.

The current core/bridge interfaces, complete point field, quantitative native full pulse C4, global admissible tensor/cone, independent flat remainder, required-domain finite energy, temporal recursion and oscillatory stress correction are still open. No complete Navier-Stokes residual or measured blow-up claim is made here.
