# Local collar g5 word 2 and formal tail correction — 2026-11-02

This checkpoint concerns a conditional local source enclosure, not global five-moment closure or a completed physical Navier–Stokes field. The existing 212-term m15 ledger, its six signed inner/lower contributions, and all previously retained errors remain untouched. P1 and m15 are distinct data.

## Correction to the preceding formal collar checkpoint

The `formal_kchi_collar_v1` receipt from commit `e944f27a52774e4f9c9407c5d29158d6872aaa5f` used a tail factor `exp(sqrt(A)/2)`, where `A=2 sup|chi| S^2`, without an authenticated small-`A` premise. That is not a general upper bound. For the actual `j>=3` starting order, take `A=10000`: the positive series `sum(n>=3) A^n/(2n)!` exceeds `A^3/6! exp(sqrt(A)/2)`. Therefore the previous receipt's `j>=3` tail is **not certified** by its stated proof. Its saved input and finite-word records remain evidence; this correction does not delete or cancel any previous error. A separate corrected inlet-tail calculation is still needed.

## Inputs, ownership and scope

The saved actual local row `g5` has ordinary `Z0` and `Z1` grouped source intervals on the first collar. The axis potential `chi` and `chi_Z` are drawn from the separately bound whole-`Z` axis source; the phase-shear cutoff is not substituted for that potential. Source IDs, formal `(hb,F0_anchor,Pstar_squared,Lambda)` powers, directed interval endpoints and uncertainty classes are retained. The intervals are enclosures, not point values. The first-potential response and the existing local `j>=2` error envelope are referenced but not recomputed, added to the new alternative envelope, or removed.

Source-family and historical axis-binding evidence do not by themselves prove that the saved `(rho,Z)` raw `g5` sectors are the physical `h5(x,Z)` in the radial Green equation. A same-owner, same-frame adapter through `rho=x^2`, physical `dx`, and all formal scale factors remains necessary. No physical m15 projection or ledger append is asserted here.

For conditional source enclosures, `S=(sc/4)exp(sc/16)` and the physical collar length obeys `L<=hb S`, assuming `0<hb<1/2`. The n=2 expression is `4 G0 chi G0 chi G0 g5`. Its whole-collar bounds add `hb^6` in F and `hb^5` in Fx. Ordinary `Z1` keeps the source derivative and both separate `chi_Z` insertion paths. The new alternative tail starts at `j=3`, is symmetric error-only, and uses `exp(A/4)`; its `chi_Z` insertion factor is `(1+A/12)exp(A/4)`. The factorial comparison is `(r+2k)!/r! >= 4^k k!` for `r=7,8`, with the insertion weight `(3+k)/3`. This is a bound, not an assigned tail value.

## Execution

The first new version was frozen and entered once. Its synthetic fixture stopped at an undefined tail-power constant before actual receipts were consumed. The complete failure, entry, manifest and source bytes are archived in `experiments/root_st073/collar_g5_word2_v1_failure.json.gz` (SHA-256 `c9a53044de6c8623b9ad2e63742c329a47a98403035cf8fe5dfb8791b4e32eaa`). It must not be re-entered.

The independently frozen v2 used fresh synthetic inputs and returned successfully in one entry. All 44 fixture checks passed. The actual saved input contains 16 sectors each for ordinary Z0 and Z1. For each of F and Fx, the n=2 output contains 16 Z0 and 48 Z1 records; the alternative `j>=3` tail contains 16 Z0 and 32 Z1 error-only records. Thus there are 128 finite enclosure records and 96 tail records. The source role was not inferred from `uncertainty`: the actual output does not assert an authenticated signed-field role for these grouped sectors. The old j>=2 envelope remains separate and untouched.

The complete v2 receipt, including frozen sources, entry and result, is [`collar_g5_word2_v2_receipt.json.gz`](../experiments/root_st073/collar_g5_word2_v2_receipt.json.gz), SHA-256 `6072b9cdf00216d74fb18211bc730d6d9ee360c9fea945614830541d543339e0`. Its source manifest SHA-256 is `acf56688fde5f05158b67f25c20dab84f9770c4931d0e374eed95e229b8df23e`, entry SHA-256 `b5cfb34b661e11ae3eb8ba158f3793d37c1592966d50524ac91758ef15362543`, and result SHA-256 `cea5d737a8b763e9212a00e478788f30b6a8e0f145f97f1caa8cfba52076406c`. The large-parameter rational fixture checks finite positive prefixes; the stated all-order tail bound rests on the factorial comparison above, not on those finite checks alone.

The frozen Python bytes are preserved exactly, including trailing blank lines; the staged whitespace check disables only the blank-at-EOF rule for those immutable sources. Git-index hashes for the v2 leaf, fixture and runner match their frozen manifest.

## Next executable work

Use only the saved results and fresh versioned consumers. Establish the raw-g5 to physical-h5 same-source and scale identity, repair the previous formal inlet `j>=3` majorant without deleting the old error, and preserve ordinary `Z1` and source correlation. The independent terminal m15 matching, functional five-moment closure, heat exterior, stress/remainder, true time-order recursion, corrections and Cartesian validation remain open.
