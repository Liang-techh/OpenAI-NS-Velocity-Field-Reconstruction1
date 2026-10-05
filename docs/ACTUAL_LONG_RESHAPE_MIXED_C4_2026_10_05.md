# Current actual R110 histories installed through long reshape

Date: 2026-10-05. Current finite-width bridge/switch histories now continue through the original long reshape to Rsh with physical logR/Z derivatives of total order at most four. This closes the long-reshape transfer gate after ACTUAL_SWITCH_MIXED_C4_2026_10_05.md.

## Run and provider

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage actualreshapemixed
```

The producer/checker are `experiments/root_st073/lei_ren_part1_paper_compliant_actual_long_reshape_mixed_C4{,_check}.py`, with matching JSON receipts.

`CompliantActualLongReshapeProfiles` retains the original inputs/evaluate/report callables and original source-general theorem gates. Its inlet is the current accepted `CompliantActualSwitchMixedC4.switch`. The outer current switch provider remains available with its accepted source graph and hashes. Taylor context, cache, A/T, reference length and kernel-rate bounds are reconciled with this same source.

`CompliantActualLongReshapeMixedC4` invokes the original physical evaluate algorithm and unchanged R109-to-R110 power method. Its small evaluation wrapper adds current source ownership and explicitly leaves the current Rsh-to-reference interface uncertified. Existing legacy files and their receipts are preserved.

## Preserved construction

The original defining F/V/chi controls, core/family/datum, Cstar/Pstar/delta, A, T=400A, logref=10(logC+logP), normalization B=log(Cstar*u110*(1+Z^2)), cutoff and derivative factors remain. All six actual moment histories, raw V110 and canonical P0 pass into inherited kernels and primitive equations.

The source-general B C2 and shear bounds remain valid because this is the same defining original field, with tighter finite-width enclosures. Three explicit inlet normalization AST bindings supplement six profile input bindings, three constant-V transport bindings and original callable identities. No source parameter or width cap becomes a defining value.

A dedicated `bridge_source_definitions(provider)` interface reads exact F/V/chi from the original AST and binds them to the current provider. It avoids adding fields to previously accepted R100 packets. The signed Ibridge/Ifirst graph and centered E expression remain exactly the accepted current graph, including all hydro/pressure/swirl scales. Raw axial V is not the pressure primitive 4C. Centered E is still metadata at this stage, not a substitute for raw V in long reshape.

## Evidence and remaining scope

Focused PASS:313 current hashes;480 velocity/pressure rows and600 physical primitive rows;6 R110 input,3 normalization and3 transport AST bindings;7 current R110 source trace groups;135 exact two-sided R110 physical mixed rows;480 final positive-source cap proofs. Whole Z[-1,1], whole reshape, actual inlet/exit, interior packets and the original R109-to-R110 neighborhood are covered.

Unchanged independent kernel/decay/current-amplitude and physical primitive fixtures, symbolic moment coordinate changes and flat endpoint chain rules are consumed through their hash-current legacy receipts. They are not claimed as newly rerun fixtures. The new checks admit the current inlet/source wiring and exact R110 rows, not interval overlap as a functional join proof.

The actual long reshape is now installed. Current reference restoration and moment patch are not installed; their corresponding join flags remain false. Exact production point history selection, recomputed five-defect implicit inputs/Jacobian/remainders, completed global tensor/flatness/physical-volume/energy, n-dependent recursion, oscillatory correction and corrected residual/dynamics remain open.

## Next

1. Build the operational centered E from the same source: j+epsilon*Psi(4,Z)+actual signed bridge increment+actual first-switch increment. Use current correlated increments and their controlled enclosures; do not subtract independent broad V110/4Z boxes or retain the old 6*cap expression as the defining E.
2. Inject current long-reshape histories into reference/restoration profiles and mixed4, preserving V=4Z+E*(1-sigma), all centered moment source equations, canonical pressure, E-squared tails and original -8/-7/-6/-5 offsets.
3. Propagate to the actual five-moment patch, then recompute current shared defects, the same implicit Jacobian and controlled nonlinear remainders. Global and recursion gates follow separately.
