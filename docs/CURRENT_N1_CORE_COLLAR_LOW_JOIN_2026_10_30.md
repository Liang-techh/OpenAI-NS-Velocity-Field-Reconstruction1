# No-reset core/collar first-word Z0/Z1 join, 2026-10-30

The saved same-source core-exit `W1_Z0/Z1` and the previously calculated local first-collar contributions now form a finite first-word six-component enclosure for real `Z in [-1,1]` and the full first-collar phase `0<=s<=sc/4`. The calculation retains the actual core history and adds each saved local sector exactly once. It is conditional on the authenticated parent source enclosures and is not a point value, full n1 solution, physical m15 integral or terminal five-moment closure.

## Transport and source preservation

With `d=hb*s`, the homogeneous transport of each core history is

```text
F = F_c + (1-exp(-d))*Fx_c
U = U_c + d*Ux_c
K = exp(-d)*K_c - (1-exp(-d))*Ux_c
P = P_c
Fx = exp(-3*d/2)*Fx_c
Ux = exp(-d/2)*Ux_c.
```

The microscopic `hb` remains a formal positive source factor. Directed bounds on the exponentials and `(1-exp(-d))/d` use the previously admitted `0<hb<1/2` geometry. Core powers `(F0,Lambda)` are lifted to the common four-atom basis `(hb,F0_anchor,Pstar_squared,Lambda)`. The saved local sectors already contain their own Green kernel and measure and receive no second integration or factor. The result retains the transported core sectors and local sectors separately; the finite sum is their unmerged concatenation, not a cancellation or a midpoint computation.

The receipt indexes all 3,408 input sectors by exact parent receipt path. Core entries retain their original outer-cell and matrix/W0 references through that index; their source-uncertainty detail was not separately classified in the parent and is explicitly labelled as such. Local entries preserve their original source IDs and error classes in addition to the new input ID. For Z0 the finite six-field group counts in `(F1,Uz1,K1,P1,d_x_F1,d_x_Uz1)` order are `313,437,437,31,225,317`; for Z1 they are `489,677,677,39,313,437`. These are unmerged interval terms rather than independent samples.

The previous finite Z2 receipt remains a separate reference. The earlier Cauchy local bound is an alternative bound on the same local quantity, and the all-word tail remains separate and unadded. No previous m15 term or error was deleted, canceled or duplicated; `P1` remains distinct from m15.

## Verification and next step

The one-entry consumer checked full source-family and atom records, parent hashes, the exact `sc` and phase interval, coordinate and ordinary-Z conventions, and every input's provenance. It also checks that each local sector is unchanged and appended once. A synthetic fixture passed 50 checks of the six transport components for both orders, full-phase inclusion, source IDs, local-once behavior and omitted-tail flags. Exponential references in this fixture use high-precision numerical evaluation and are cross-checks, not certified point values.

The first attempt stopped in the synthetic fixture because decimal reference truncation was too coarse for an exact rational interval. It did not perform the actual join. Its frozen sources and result are retained in `experiments/root_st073/n1_core_collar_low_join_v1_failure.json.gz`. The successful receipt is `experiments/root_st073/n1_core_collar_low_join_v2_receipt.json.gz`, SHA-256 `a2e25ccec22f97dc99d1eb80609f0d806d4276e7a75e763ca882d90532e7ce67`; the one-entry result SHA-256 is `ba6d0f948d08f9cfac5e9906122973c99a98ea01e4172f5dd795fc6c6252083e`.

Next, construct same-source physical m15 radial-integrand and endpoint contributions from the needed field and mean/pressure histories at ordinary Z0/Z1; this finite first word alone is insufficient for the full solution. Establish the corrected-reserved/current-pulse owner, measure and endpoint bridge before attempting correlated cancellation. Retain the previous 212 cumulative terms, all inherited errors and separately recorded all-word remainder.
