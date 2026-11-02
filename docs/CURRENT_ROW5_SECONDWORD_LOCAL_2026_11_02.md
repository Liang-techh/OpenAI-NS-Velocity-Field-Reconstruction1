# Raw-K0 row-5 local second-feedback seed — 2026-11-02

The long-term Navier–Stokes reconstruction remains **active and incomplete**. This checkpoint computes conditional source-bound coefficient enclosures for the local row-5-driven second K0 feedback word. It is not the complete second state, a bound on all of `R_hi`, a physical velocity or m15, or a proof of terminal five-moment closure. The inherited 212-term m15 ledger, six signed inner/lower terms and prior errors remain unchanged. P1 is not m15.

## Construction and limitations

For the saved finite firstword state, the raw row-5 feedback and its ordinary Z1 are

```text
H0 = C0 W_firstword + B1 W_firstword,Z
H1 = C0,Z W_firstword + C0 W_firstword,Z
     + B1,Z W_firstword,Z + B1 W_firstword,ZZ.
```

The source records are the previously emitted firstword product ledger. The consumer selects its C0/B1 terms and excludes standalone g5 and the separately introduced `+2 chi` axis terms. This is the raw K0 operator split; it is not a cancellation or a claim that the potential branch is zero. It forms the local responses `K0 H0`, `K0 H1`, and their normalized radial derivatives. The independent second-feedback core inlet, Kchi potential branches, other state components and all later feedback words remain separate and unknown.

For normalized `x=2 exp(hb*s/2)`, `s=sc*q/4`, cell upper endpoint `qhi`, and `S=(sc/4) qhi exp(sc/16)`, the geometric premises imply `x-2<=hb S`. The nonnegative K0 and d_x K0 prefix masses are bounded by `hb² S²/2` and `hb S`, respectively. The new leaf retains the hb powers formally. Each signed coefficient interval is represented as a product of its saved directed source coefficient, its integer scalar applied once, and the corresponding nonnegative kernel-mass interval; source IDs, state IDs, error classes and the original product record remain traceable. A separate monotone sum of absolute coefficient bounds is grouped by target cell, field, ordinary Z order, formal powers and conservative finite/error flag. It is a radius, not a signed field value; error classes remain linked to their original nodes. No midpoint, atom value, signed cancellation or old error deletion is used.

The same uniform whole-collar source enclosure is reused for all four target prefixes; these are not independent source-cell measurements or exact triangular cell masses. The common source owner, uniform state enclosure, fixed-normalized-x derivative and Z-independent geometry are caller premises, not proven by this result. In particular, Z1 needs the kernel and endpoints to be Z-independent. The normalized x is not physical Cartesian X, and no physical-volume integral or m15 contribution is computed.

## Run and evidence

The first frozen entry failed in its synthetic fixture before decoding actual receipts because the zero bound was supplied as a Python integer where an MPF endpoint tuple was required. The full failure, exclusive entry, fixture result and frozen code are retained in [the v1 failure archive](../experiments/root_st073/row5_secondword_local_v1_failure.json.gz), raw SHA-256 `c8446208821b7ab6055e4538e7a344113ee834755b6a625ae5cdac2410818144`; it was not re-entered. A separately versioned v2 uses fresh synthetic values and the saved q-cell partition.

The v2 entry stopped during its synthetic fixture because an exact MPF endpoint tuple was compared with an integer. Its [failure archive](../experiments/root_st073/row5_secondword_local_v2_failure.json.gz), raw SHA-256 `083b6cbe55af387de99177087b99400894197c575aa7036e2e5a2598bd3b7db8`, contains its traceback and frozen code. V3 passed 48 of 49 synthetic checks and stopped because its fixture compared list and tuple representations of identical endpoints; its [failure archive](../experiments/root_st073/row5_secondword_local_v3_failure.json.gz), raw SHA-256 `70cca25faaea2b9ad6ff7cf3c30db8776d99d76542ef687542a0c83911fe6005`, retains that result. Neither entry decoded actual source receipts, and neither was re-entered.

The frozen v4 entry returned successfully after **49 synthetic checks**, including exact endpoint comparisons, sign and scalar handling, error-class preservation, malformed input rejection, nonuniform synthetic cells, and independent closed-form K0/Fx prefix-mass examples. It consumed the authenticated saved source only once and selected **16,498** C0/B1 products for h5 and **36,900** for ordinary h5_Z1. The 1,147 g5/axis products were retained as excluded references, not canceled. The selected products form **224** formal-power/error groups and **1,792** separate target-cell/field/order absolute coefficient bounds. Error tags on source, state and product nodes are retained; classifying any tagged uncertainty as error is conservative and does not establish an error-only physical decomposition. The actual cells are the previously saved four equal q intervals; nonuniform cells appeared only in synthetic checks.

The complete [compact receipt](../experiments/root_st073/row5_secondword_local_v4_compact_receipt.json.gz) embeds the exact source, entry and full result JSON along with base64 frozen code; raw SHA-256 `70ce750f17202fc59a8f39d61eab7da864fc447be6a8ea9c79a3e625ea34ea55`. Its source manifest SHA-256 is `7902007d1883f076b970030a5d60f825b8953e3a28b74cb735d9b057003185e6`, entry SHA-256 `90c08ccb59b5c3d52051a514feb41249442b69743b2b4dc9c36c4322bd50eb2d`, and result SHA-256 `88702a0c24c4e453f1aeed77cffb6c941ec73ac0d8174c102734c5571bf20794`. All raw gzip inputs were authenticated before decoding; no ancestor producer or helper was invoked. The source rows retain the input's unproven common-owner and full-collar-uniformity assumptions.

The frozen leaf and fixture retain trailing blank lines. The staged whitespace check disables only `blank-at-eof`; staged executable source hashes match the frozen manifest. The compact receipt was byte-checked against the saved source, entry and complete result before publication.

## Remaining work

Firstword state/source ownership, fixed-x/phase geometry and gauge must be established before interpreting the local seed as part of a single physical solution. Independently bind the nonzero second-feedback core inlet, the omitted components and Kchi branches, then certify the whole feedback tail; the ordinary Z2 remainder requires at least a firstword Z3 input or a valid analytic Cauchy-tube bound. The local absolute radii cannot be added to an overlapping all-word bound as if the two were disjoint, and no finite seed may replace an old error without a same-source identity. The own-swirl/negative full-pulse functional comparison and later pressure, exterior heat, stress, time recursion, correction and Cartesian NS validation remain open.

Workers were requested using the fixed Luna/max role; runtime model metadata was not exposed for independent verification.
