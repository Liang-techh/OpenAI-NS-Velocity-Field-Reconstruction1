# Current correlated E and reference/restoration mixed4 through Rm

Date: 2026-10-05. Current actual long-reshape histories now drive the original reference continuation, axial restoration and postrestoration to the unpatched Rm. The operational centered E uses the same fresh core source and current actual increments. This builds on ACTUAL_LONG_RESHAPE_MIXED_C4_2026_10_05.md.

## Run and providers

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage actualrestoremixed
```

Producer/checker: `experiments/root_st073/lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4{,_check}.py`, with matching JSON receipts.

`CompliantActualReferenceRestoreProfiles` consumes the current long-reshape provider. Original reference transport, restoration kernels, source equations, pressure datum, centered moment coordinates, offsets -8/-7/-6/-5 and final V=4Z remain. `CompliantActualReferenceRestoreMixedC4` reuses the original packet/reference/restoration/postrestore methods and physical mixed4 equations. Legacy modules and receipts are unchanged.

## Correlated current core and E

E is the exact source view j+epsilon*Psi(4,Z)+Ibridge+Ifirst. Directed enclosures of its components remain separate from the defining integrals. Neither independent subtraction of broad V110/4Z boxes nor the old six-cap formula defines E.

The legacy normalized_jets route and the fresh bridge inlet use different finite/tail enclosures. Their interval endpoints cannot establish an exact source identity. The new `current_centered_core` instead uses the SAME fresh recurrence packet and current source_profile_jet as the bridge:

1. Rebuild the current core packet at degree24/axial depth6.
2. Bind the original seed u=4Z+j, fixed U0 coefficients and initial Uz radial row0. Check its six axial coefficients against the exact affine baseline.
3. Make a local copy of the radial rows and replace only Uz row0 by [j,0,...], removing exact 4Z before summation. Every radial row n>=1 remains unchanged.
4. Use the same current source_profile_jet at rho4 and its unchanged epsilon*Psi tail. This encloses the same defining core field minus4Z.

Seed, inlet, fixed-row and eight local-copy/reconstruction AST bindings establish this source route. A read-only GPT-5.6 Luna/max review accepted the resolution of the core-source identity gap.

Add the current bridge macro packet's actual_delta_V and the current first-switch velocity_increment cover. The latter is an enclosure of Ifirst, not its defining integral. The common signed source graph, known comparison direction and original pressure datum remain.

Before applying the original E C2 theorem envelope, the current outward sum of ordinary derivative increment bounds through order2 is checked against rho_bridge, whose source theorem explicitly includes both short switches. Here the increment sum is at most about9.005e-49932, below the admitted budget about1.748e-502. These are bounds for this construction, not measured momentum residuals or chosen field values.

## Acceptance and limits

Focused PASS:317 current hashes;720 velocity/pressure rows and900 physical primitive rows;the current Rsh parent and all six centered histories retained;canonical P0 at each packet's actual Z;270 exact two-sided Rz/restore-exit rows;exact terminal Uz=4Z and its mixed derivatives;30 current unpatched five-defect coefficients at Rm.

Hash-current generic restoration/physical fixtures and structural identities are reused for the unchanged algorithms, rather than rerun or described as new actual production-point tests. The current source, component and C2 budget checks bind the new data route.

Current correlated E and reference/restoration mixed4 are installed. The full current Rsh source-functional mixed4 join is still an independent open gate; current patch, full implicit leading inputs/Jacobian/remainders, exact production points, global tensor/flatness/physical-volume/energy, n-dependent recursion, oscillatory correction and corrected residual/dynamics remain open. Native Rsh interface flags and the new aggregator Rsh flag stay false.

## Next

1. Establish the current Rsh long-reshape/reference functional mixed4 join using the same current primitive sources, centered moment coordinate identities, raw V110/E source relation and flat cutoff derivatives. Do not certify it by overlap.
2. Inject current restored Rm histories into the original actual five-moment patch and its mixed4 provider. Retain original bump functions, source integrals, E-squared terms, pressure datum and physical radius factors on x in [1,e].
3. Recompute common five-defect leading inputs, the same implicit Jacobian and controlled nonlinear remainders. Existing legacy receipts certify the legacy provider, not the current data route.
4. Extend current dispatcher ownership and then the independent global background gates before actual n-dependent recursion and oscillatory correction.
