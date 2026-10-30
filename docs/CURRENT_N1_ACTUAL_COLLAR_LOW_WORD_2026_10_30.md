# Actual first-collar local W1 ordinary Z0/Z1, 2026-10-30

The actual same-source core `W0` inlet and first-collar forcing/matrix enclosures now yield the local first-word `W1_Z0` and `W1_Z1` Green contributions uniformly over the first collar and real `Z in [-1,1]`. This is a directed-interval enclosure conditional on the parent source assumptions. It is only the **local contribution**: the independent core-exit `W1_Z0/Z1` inlet is still missing and is neither set to zero nor included. The result is not a full field, a physical m15 integral or terminal closure.

## New calculation

At fixed phase the new feedbacks are

```text
S0 = C0 W0 + B1 W0_Z
S1 = C0_Z W0 + C0 W0_Z + B1_Z W0_Z + B1 W0_ZZ.
```

The local `W0_Z0..Z2` needed inside these new feedback integrals combines the saved same-source core `W0` inlet with the collar `g` Green integral. These `W0` quantities were intermediates of the earlier Z2 calculation and were not separately saved there; this new calculation forms the needed intermediates for the two new derivative orders, without rerunning an ancestor program or recomputing the previous Z2 output. It applies the six signed kernels and physical `dx=hb*u*dv/2` measure, retaining all four formal source atoms and their powers. One inner and one outer interval cell enclose the full phase; grouping is outward interval addition and does not prove same-source cancellation.

The local output has `137,197,197,23,137,197` compact groups for each order, in `(F1,Uz1,K1,P1,d_x_F1,d_x_Uz1)` order. It references 465 distinct source IDs, with the full 624 forcing/matrix input-group index retained in the receipt. Core whole-source, Cauchy-phase, mean-phase, cutoff-value and cutoff-derivative uncertainty classes remain in the output. The earlier finite `W1_Z2` receipt, alternative Cauchy bound and all-word tail remain separate references. No m15 term was added, removed or identified with `P1`.

## Checks and remaining work

The one-entry consumer authenticated the core inlet, collar forcing, phase, previous Z2 result, Cauchy and tail receipt chains. An independent synthetic fixture passed 39 checks: all six direct numerical Green-integral comparisons for each order, full-phase inclusion, source provenance, state/row mapping, and explicit missing-inlet and omitted-tail flags. Its high-precision quadrature is a cross-check, not a certified quadrature error bound. The actual inputs and output remain directed interval enclosures, not selected point values.

The result receipt is `experiments/root_st073/n1_actual_collar_low_word_v1_receipt.json.gz`, SHA-256 `f4ed142fd3c8ec79b919af92ad12f3e160b06a6c8b13a932bb4fc17df27df234`; the one-entry result SHA-256 is `28237f79860c2fd3a0c838d4fb7026194b5021d1e030a97b5ceadced21b1ed90`.

The prior core Z2 receipt retains the Z2 feedback products but does not retain the independent core `C0/B1` matrix factors; interval products cannot safely be divided to reconstruct the missing Z0/Z1 inlet. Next, construct those two core first-word histories from same-source core matrix data and already saved core `W0` outer-cell sectors, preserving the eight cell histories and avoiding a rerun of the old Z2 calculation. Then join them to these local contributions and construct the separate physical m15 radial integral and endpoint terms. The corrected reserved and current physical pulse sources still require an owner/measure/endpoint bridge, and all 212 existing terms and inherited errors must remain intact.
