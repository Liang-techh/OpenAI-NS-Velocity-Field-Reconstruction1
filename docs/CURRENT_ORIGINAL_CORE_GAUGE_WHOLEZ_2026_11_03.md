# Original core gauge on the admitted Z interval — 2026-11-03

The active long-term Navier–Stokes goal is **not complete**. This checkpoint adds a separately evaluable original **n=0 core gauge**, its ordinary first Z derivative, and rational interval enclosures on the admitted domain Z ∈ [-1,1]. It does not install a graph owner or construct the n=1 response, physical m15, or a five-moment closure.

## Definition and proof scope

The frozen source recipe uses the exact rational `j = 1/5721363714816000000000`, `sigma = j/500`, and the actual transcendental `delta = exp(-4*(exp(40)+11)-30)`. The evaluation encloses delta in `[0,10^-200]`; neither endpoint is substituted for its value. Define

```
H(v) = (1-delta)*v/2 + (1-v^2)*(4*v+j)
L(v) = 1-delta*v^2
G(Z) = integral_a^Z L(v)*H(v)/(H(v)^2+sigma^2) dv
G_Z(Z) = L(Z)*H(Z)/(H(Z)^2+sigma^2),
```

where `a` is the unique zero of H on `[-j,0]`. Exact bisection encloses both endpoint-delta roots, and their union encloses the actual anchor. The code proves the relevant signs on `[-1,1]`, uses geometric panels and rational interval rectangles, and rounds every panel outward to a dyadic grid. The anchor shift, rectangle enclosure width, and outward-rounding budgets are recorded separately. The entry called `integration_truncation` is an enclosure-width budget, not a signed error or a function value. Bounds and uniqueness are scoped to the stated domain; no all-real-axis assertion is made.

Source provenance includes `CURRENT_ORIGINAL_CORE_LOG_AMPLITUDE_2026_10_15.md`, the exact `core_transfer.json` data, and the source texts frozen in the receipt. This is a standalone evaluation of that recipe; the source graph binder/live owner is not proven by a name match.

## One-entry result and limitations

The source manifest was frozen first (SHA-256 `afbe0f06e105e2059798c0d501ec2b685f254f9679587b3754579282e62d35ab`); a single entry and result were saved before review. The new synthetic fixture ran once and all its checks passed, including a separate exact rational derivative witness, positivity witness, and rejection of an out-of-domain cell. One owner then evaluated four cells once. The exact rational results and their full budgets are in the authenticated receipt:

`experiments/root_st073/original_core_gauge_wholez_v1_receipt.json.gz`
SHA-256 `6335d097f5250a69519246ff0562c6779590b74c462ba69562de891b564d590a`.

The following decimals are rounded for orientation only; the exact rational enclosures in the receipt are authoritative. They are bounds, not function samples or midpoint values.

| Z cell | G interval, approximately | G_Z interval, approximately |
| --- | --- | --- |
| [-1,-1/2] | [4.46912, 40.4592] | [-9.19720, -0.0938490] |
| [-1/2,0] | [0, 36.2961] | [-1.43034e24, 1.43034e24] |
| [0,1/2] | [0.386700, 36.1746] | [0.239844, 9.93289e21] |
| [1/2,1] | [4.49328, 40.3377] | [0.0938490, 9.19720] |

The enclosed anchor width is about `7.13e-78`; the recorded anchor-uncertainty budget is about `1.31e-105` and cell rounding budgets are below `8.7e-114`. The much larger rectangle enclosure-width budgets are approximately 34.7–67.8. These conservative enclosures are too broad to settle the signed source correlation or internal endpoint cancellation. They are not a certificate of the claimed terminal conditions.

**Receipt metadata discrepancy:** the outer result's `whole_Z_coverage` field remains `false` from its initial default although its `coverage.whole_Z` and `leaf_result.whole_Z_coverage` fields are `true`; four returned cell domains cover `[-1,1]`. The sealed result is preserved unchanged. Readers must inspect the four saved cell records and this discrepancy rather than treating the outer flag as a successful coverage assertion.

## Preserved state and next task

The 212 accumulated m15 terms, the six signed inner/lower contributions, inherited errors, and all existing user modifications remain untouched. No previously completed producer, fixture, or ancestor integral was replayed. P1 and m15 remain distinct. No cancellation was inferred from overlap, and no old error was removed.

Next, recover constructor-level binder and live-owner witnesses for the actual gauge integral (the saved name-only candidate is not such a witness), and establish the native/corrected preRp same-source ordinary-C1 relation. Then integrate the independently defined reserved n1 own-swirl contribution exactly once against the separately defined negative full-pulse m15 integral, carrying physical scales and inherited errors. If tighter gauge enclosures are needed, plan a new frozen method with a priori error bounds; never re-enter this frozen component.
