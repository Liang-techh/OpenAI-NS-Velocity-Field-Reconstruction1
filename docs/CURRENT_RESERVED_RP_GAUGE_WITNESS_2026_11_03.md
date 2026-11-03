# Reserved Rp support and core gauge graph witness — 2026-11-03

The long-term Navier–Stokes objective remains **ACTIVE / INCOMPLETE**. This checkpoint advances two prerequisites: it checks the gauge constructor against edges in the saved graph, and implements an explicitly conditional bridge from the reserved compact bump support to a physical Rp endpoint. Neither result establishes a live common owner, the actual m15 endpoint comparison, or five-moment closure.

## Saved gauge graph: structural witness, not owner identity

The frozen raw inner-Volterra archive has SHA-256 `1b44be8c6cbfffb692c7a0430ecbc8a590f46001a233c99c59ca8c5e2c44a6fe`. After raw-byte authentication, the new consumer read its 49,024 serialized graph rows and the pinned core constructor source. It found this direct graph chain:

```
35646: definite_integral, integrand -> 35645, lower -> 35561, upper -> 20
35645: function_substitution, value -> 35644, variable -> 20
35644: bound_variable original_true_core_gauge_v
35561: true_core_anchor source parameter
```

The constructor AST separately contains the corresponding `v = g.symbol(...)`, substitution with `value=v.node`, and integral over `g.nodes[v.node]['name']`. This uses the actual substitution edge, not just a graph-wide name search. Two newly introduced negative structural cases, substituting the wrong binder and changing the upper endpoint, were rejected.

The result is a **candidate structural witness**. Source text and archived graph matching do not establish runtime control flow, a common live graph owner, archived-to-live owner identity, or exact source identity. No binder mapping was installed and no gauge function or numerical m15 was evaluated in this checkpoint.

## Conditional Rp support and ordinary Z1

The original reserved constructor defines five centers `j/40 * log(2)` for `j=12,...,16`, width `log(2)/400`, and angular channels `(1,2,3)`. Thus the last of those three angular bumps ends at `(151/400) log(2)` in `log(R/Rc)`, leaving the exact positive gap `(249/400) log(2)` before `2 Rc`. The separately saved C3 quiet source has `Rp=Rc*exp(Tw-2)`, `quiet_length=Tw-2-log(2)`, and a previously checked positive quiet interval. These observations concern independently pinned source records; they do not themselves bind their owners.

`reserved_rp_support_bridge_v1.py` requires explicit premises for shared source owners, positive quiet length, Rc mapping, and a single invertible angular matrix whose third row is the complete m15 row **as a full-Z C1 identity**. Under those unverified premises the full-support inverse action is `(0,0,1)`, with ordinary response-row derivative `(0,0,0)`. The source functions and their ordinary Z1 derivatives remain inputs. The certified interfaces always reject; the conditional interface requires explicit acknowledgment and cannot certify a physical row.

Six new synthetic checks passed, including rejection of a mismatched owner, a false matrix premise and certified-row use. A fresh algebra probe supplied nonzero pulse and derivative data, used the conditional row in the existing pure variable-Z evaluator, and recovered `D=J=13` and `D_Z=J_Z=17`; it also retained the `J-P` coefficient of a possible response-row derivative. These synthetic numbers are not physical data. No ancestor integral, producer, or old fixture was rerun.

## Provenance and next action

The component `reserved_rp_gauge_witness_v1` was frozen before its one entry. Its manifest SHA-256 is `ef62516df3ddbcb32ba806be00ebf0235621882257939850d8314d9b4dfc5e7c`; entry SHA-256 is `db58b15166c645daf3b447b187960310d4167da6167366d092f0c65855c524ca`; result SHA-256 is `f20af375ec2f15bc0254a9b9701977bf5d734790615bf50681e9b8cbd93ff352`. Full frozen inputs, source, entry and result are retained in [the receipt](../experiments/root_st073/reserved_rp_gauge_witness_v1_receipt.json.gz), SHA-256 `b8704427d20d0015da83afb5562524f35c60cf741980120419b3c97ddb795e6c`. Never re-enter this frozen component.

Next, bind the reserved support and quiet endpoint to the **same actual source owner and physical Rc**, and obtain the full-Z C1 matrix-row identity and invertibility from that same source. In parallel, connect the native 2Rc m15 ordinary-Z1 residual to corrected active and old compact sources using their actual endpoint and physical measure; the quiet Z1 defect is already zero in the saved source reduction. Only then may the conditional bridge be promoted and the reserved contribution added exactly once against the independent negative full-pulse m15 target.

Preserve the 212 forward terms, six signed inner/lower terms, all inherited errors and the independent pulse target. P1 is not m15; overlap and magnitude caps are not function values, and no old error has been removed. The exact heat exterior, stress/remainder, higher-order recursion, oscillatory correction and independent Cartesian validation remain open.
