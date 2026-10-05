# Actual finite-width histories installed in switch mixed4

Date: 2026-10-05. This completes the mixed4 installation left open in ACTUAL_BRIDGE_SWITCH_2026_10_05.md. The current actual bridge histories now drive the original first switch, second switch and postpower derivative algorithms. It does not complete the downstream background or select production point histories.

## Run and interfaces

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage actualswitchmixed
```

Producer/checker: `experiments/root_st073/lei_ren_part1_paper_compliant_actual_switch_mixed_C4{,_check}.py`, with matching JSON receipts.

`CompliantActualSwitchMixedC4` retains the original evaluate, postpower and width callables. Its bridge inlet is `CompliantActualBridgeSwitch`, consuming the current finite-width bridge source, own H/M/K/A/B/C histories and canonical P0. The known comparison histories continue to prescribe Dbar/Ebar; actual own moments recover the actual field and pressure.

`ActualSwitchSourceDispatcher.evaluate(chart,Z,coordinate)` owns exactly `switch_first`, `switch_second` and `switch_power`. The two switch coordinates are phase in [0,1] and [1,2]; the power coordinate is the original logarithmic fraction in [0,1]. It loads the new accepted receipt and checks current hashes before constructing its provider. This is not replacement of the full source dispatcher or a Cartesian field API.

## Source and derivative connection

The existing microscopic physical width factors and factored logarithmic source ledgers remain unchanged. The full F0 amplitude range stays separate from the extremely small width until final physical bounds. Neither width caps nor interval endpoints define a production field.

Ten AST bindings identify the original mixed4 input and control path. Seven further bindings identify the current actual contributions, field changes, exp(ell) and actual phi/V definitions. The shared axial source graph now belongs to this current source family and datum; it preserves the signed Ibridge/Ifirst definitions and the correlated V100, V110 and E primitive expressions.

At R100, exact phase-zero branches transfer the current actual phi/raw V, all six own moments, Q and pressure. Current source traces and the frozen comparison assignments are bound independently of interval overlap. The shared original physical primitive and 17 symbolic width-scaling identities relate bridge y derivatives to switch phase derivatives through order four. Width is independent of Z, so this also covers the required mixed rows.

The new aggregator therefore sets `actual_R100_functional_mixed4_join_certified=True`. The inherited packet's native `bridge_switch_inlet_join_certified=False` metadata is retained rather than silently relabeled. The claim here is the specifically traced current R100 source join; it does not certify every inner interface.

## Acceptance and scope

The focused checker passes with:

- 311 current input hashes and the original accepted generic fixture receipt.
- 1,080 first/second mixed derivative rows and factored log-radius ledgers.
- 405 complete postpower rows and phase1/R2/R110 functional joins.
- 10 mixed4 input/control AST bindings, 7 current field bindings and 7 exact R100 field/moment/Q/pressure trace groups.
- All three current dispatcher charts, with 405 derivative rows.
- Original source controls, canonical pressure and deferred width/amplitude caps preserved.

The original row verifier is replayed with the current provider and shared source owner. Unchanged generic comparison/control/physical/factored-scale fixtures are reused only through their current hash-bound receipt; they are not described as newly rerun independent comparisons. A separate smoke call also loads the acceptance receipt through the uncached dispatcher path and evaluates an interior first-switch point.

A bounded read-only GPT-5.6 Luna/max review found no material gap in this installation after the current field, phase-zero and control bindings were added.

Available: source-functional mixed4 enclosures for the actual R100-to-R110 history. Still open: exact production point history recovery; downstream reshape/restoration installation; recomputed leading five-defect input/Jacobian/remainders; full global tensor admissibility, temporal flatness, physical-volume norms and required-domain energy; actual n-dependent recursion; oscillatory correction; corrected residual and measured dynamics.

## Next implementation gate

1. Expose an exact bridge source-definition interface for long reshape. Its current `formal_axial_source()` expects `bridge.actual(...)[source_integral_definitions]`, absent from `ActualR100BridgeAdapter.actual()`. Bind any adapter interface to the current actual F/V/chi and the switch's exact Uz source.
2. Propagate the current provider through long reshape profiles/mixed4, reference restoration profiles/mixed4, and actual moment patch/mixed4. Use bounded local adapters or a deliberate injected constructor route. Legacy receipts continue to certify only the legacy construction.
3. Preserve the centered E source expression from V110 minus its known 4Z baseline before evaluation. Do not subtract wide independent baseline boxes, substitute a cap, or retain the old epsilon*Psi+j+6*cap formula as the defining E. Reference restoration must retain V=4Z+E*(1-sigma) and its exact final 4Z trace.
4. Keep the same core/family/source/datum and original Cstar/Pstar/delta/T/A. Preserve T=400A, the B and C2 bounds, logarithmic reference length and the original restoration offsets. Carry all six moment histories, canonical P0 and the five-patch x in [1,e], including E-squared tails.
5. Admit the actual downstream mixed derivatives and functional joins, then recompute shared five-defect inputs, the same implicit Jacobian and controlled nonlinear remainders. Only then advance the corresponding leading-background gate.
