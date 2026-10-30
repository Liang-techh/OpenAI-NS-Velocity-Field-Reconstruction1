# Local first-word cutoff enclosure (2026-10-30)

This is a local, conditional enclosure, not global five-moment closure or a
completed Navier–Stokes field. The preceding accepted state is
`CURRENT_COLLAR_RADIAL_Q_NATIVE_LINK_2026_10_30.md` at commit `3f92f62`.
The cumulative 212-term m15 calculation, its six signed inner/lower terms,
the separately reserved contribution, and all inherited errors remain as in
that state; none is recalculated or added here.

## New construction and one-entry result

`n1_uniform_first_word_certificate.py` builds a sector-preserving enclosure
from the saved q=1 finite-word output and the separate N=1 raw all-word tail.
Its one-entry caller `run_n1_uniform_first_word_once.py` checks frozen input
hashes, the accepted source-family and tail pins, the W0+W1 producer and
transport formulas, and the retained interval and cutoff premises. It imports
no ancestor producer. The recorded invocation returned successfully with one
consumer call. Source, entry, complete result, hashes and exact interval
endpoints are in `experiments/root_st073/n1_uniform_first_word_v1_receipt.json`.

The saved Z interval is the very narrow directed enclosure around 0.371
(width about 6.10e-502); it is not an interval covering the full axial domain.
For each of ordinary Z orders 0 and 1, the result retains three signed formal
source sectors for F1 and four for Uz1. Their coefficients and source-log
intervals remain separate; the positive source exponentials are not
materialized. Local W0 and W1 error radii and the raw omitted-word tail are
retained separately, then added once to form the positive enclosure radius.
The tail is a vector-norm majorant, not a signed field value. The exact MPF
endpoint tuples in the receipt are authoritative: scientific displays use
floating-point logarithms and can be inaccurate when exponents are enormous.

For the preserved source premises, q=1 uses the transport box
`0 <= y <= sc/8`; positivity and monotonicity of the local radii make this
one box cover `0 <= q <= 1`. This is a uniform enclosure, not an evaluation
of the field at q=1 used as a value at other q. The cutoff
`eta(q)=1-sigma(4q-2)` is bounded by [0,1], with exact zero for q>=3/4;
ordinary Z differentiation does not differentiate eta. The cutoff enclosure
here concerns only F1 and Uz1. Mean and pressure histories, including N1p,
are not set to zero when eta vanishes. The output is not an exact q-dependent
cutoff oracle, an integral, or a signed terminal value.

This representation differs from the earlier q=0 inlet plus q=1 variation
bound: it retains the q=1 transport's formal source sectors. It must not be
appended as another m15 contribution. P1 and m15 remain distinct quantities.

## Remaining work

Connect the source-bound raw finite prefix and its all-word/local errors to
the same-source mean and pressure histories and physical signed m15 integral,
including ordinary Z1 and physical factors exactly once. Establish the
finite-prefix residual and compare it with the canonical Volterra source;
the old finite polynomial and a majorant cannot substitute for that equality.
Then prove the relevant same-source endpoint identities before cancelling any
large inherited error. Independent terminal matching as a function of Z,
the exact heat exterior, stress/remainder, higher-order recursion, corrections,
and independent Cartesian NS validation remain open.
