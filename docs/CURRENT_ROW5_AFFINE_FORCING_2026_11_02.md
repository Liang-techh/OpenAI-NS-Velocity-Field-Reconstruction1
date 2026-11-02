# Source-linked row-5 affine forcing — 2026-11-02

The long-term construction remains **active and incomplete**. This checkpoint builds a conditional affine source expression using saved actual first-collar `g5`, `C0`, `B1`, and axis `chi` enclosures. It does not provide the unknown full state `W`, physical point values of `h5`, a completed collar solution or an NS validation. The 212-term m15 ledger, the six signed inner/lower terms and all previous errors remain untouched; P1 and m15 are distinct.

## Row-5 definition

For the raw six-row system `W_x + D W/x = B0 W + B1 W_Z + g`, `W1_x=W5`, and the only kinematic subtractions in `C0=B0-Akin` occur at `(1,5)`, `(2,6)`, `(3,6)`. Subject to the raw-system/frame premise, row 5 of C0 is therefore row 5 of B0. With `L0 W1=W1_xx+3 W1_x/x`, the *uncancelled* right side for `L0 W1+2 chi W1=h5` is

```text
h5 = g5 + C0_51 W1 + C0_52 W2 + C0_53 W3 + C0_55 W5
        + B1_51 W1_Z + B1_52 W2_Z + B1_53 W3_Z + 2 chi W1.
```

At fixed `x`, ordinary Z1 is formed with the product rule:

```text
h5_Z = g5_Z + sum_j(C0_5j,Z Wj + C0_5j Wj,Z)
              + sum_j(B1_5j,Z Wj,Z + B1_5j Wj,ZZ)
              + 2 chi_Z W1 + 2 chi W1,Z.
```

Here `j={1,2,3,5}` for C0 and `j={1,2,3}` for B1. In particular, the latter requires `W1_ZZ`, `W2_ZZ`, and `W3_ZZ`. The explicit `+2 chi W1` and its Z derivative remain separate: no common-owner gauge cancellation is assumed. Fixed-x differentiation and `rho=x²` are caller premises; no `dx` integration or coordinate transformation is performed here.

## Actual source inputs and single-entry result

The saved actual forcing receipt contains raw `g5` and dynamic `C0/B1` coefficients produced from the same source-input jets, with ordinary Z0/Z1 and formal powers `(hb,F0_anchor,Pstar_squared,Lambda)`. It also records that the relative F0 derivative was applied once. The new runner authenticates those saved bytes, axis-jet binding and source-family/log records before constructing the expression. Matching records are historical source evidence, not proof of a common live owner for the unknown W or of a physical `h5` frame.

The frozen `row5_affine_forcing_v1` was entered once. All **28 synthetic checks** passed, including independently computed rational h5/h5_Z results, negative formal powers, both chi product terms and rejection of a missing required `W3_ZZ`. On the actual saved coefficient packets the unevaluated expression contains:

| Output | Source g5 | C0 feedback terms | B1 feedback terms | Axis terms | Total |
|---|---:|---:|---:|---:|---:|
| h5 | 16 | 20 | 14 | 1 | 51 |
| ordinary h5_Z | 16 | 40 | 28 | 2 | 86 |

These are source-linked additive expression records, not point evaluations. They retain individual directed coefficient intervals, four formal powers, original paths, source IDs, uncertainty classes and explicit state-jet references. The actual packet was **not** passed to the evaluator and no W or formal-atom values were supplied. The evaluator's synthetic check only establishes the algebra for caller-supplied intervals; a caller must authenticate all state, source, axis, frame and geometry premises before interpreting any evaluated enclosure physically. No interval midpoint or bound is selected as a function value.

The complete [receipt](../experiments/root_st073/row5_affine_forcing_v1_receipt.json.gz), containing the frozen sources, entry and full result, has raw SHA-256 `b64d14b59b178502e562fe593d79419cc2e21200f80c9bf0006905d1b6e72515`. Its manifest SHA-256 is `c85998ecac56fab927a8ba0bcdd55e7311a23c57676514603661e639711733f9`; entry `7774b5c5e56f9b07ff83a9a8c5ce81d5f9d3c4239d199da57084c804c22b013a`; result `05be76b86e3c61e91b854839016ae337bb3290be15aff13fb80b15985479a03c`. Raw gzip inputs were authenticated before decoding; no ancestor calculation was rerun.

The leaf and fixture retain their frozen trailing blank lines. The staged whitespace check disables only `blank-at-eof`, and Git-index hashes for leaf, fixture and runner match the frozen manifest.

## Remaining executable work

Bind actual collar state enclosures and ordinary Z0/Z1/Z2 to this expression with the same source, scale and gauge; retain the distinction between finite W0/low-word components and the still-missing complete feedback state. Establish `rho=x²`, physical-x cell ownership and exactly one physical `dx` factor at integration. Only a same-owner identity can justify cancelling the principal `-2 chi` contribution against the explicit `+2 chi W1`; until then keep both. Next compare the reserved own-swirl contribution with the independent negative full pulse as functions of Z, without dropping inherited errors. Terminal five-moment closure, exact heat exterior, stress/remainder, time-order recursion, corrections and independent Cartesian NS validation remain open.

Workers were requested with the fixed Luna/max role. The spawning tool did not expose runtime model metadata, so independent verification of that metadata was unavailable.
