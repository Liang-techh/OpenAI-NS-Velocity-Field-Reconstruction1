# Kchi cell geometry and compensation kernel errors — 2026-11-02

The long-term construction remains **active and incomplete**. This checkpoint
verifies coordinates in two saved eight-cell packets and propagates previously
uncomputed `+2 chi E51` kernel-error terms over the target cells. It does not
prove common live source ownership, complete W1, m15, or terminal matching.

## Saved-cell geometry

For each index `i=0,...,7`, the authenticated Kchi and K0 records have exactly
`x=[i/4,(i+1)/4]`, `t=[i/8,(i+1)/8]`, and
`rho=[i^2/16,(i+1)^2/16]`. Thus the saved coordinate records satisfy
`x=2t`, `rho=x^2=4t^2` and the source cells retain `j<=i` with the diagonal
triangular. The bridge also exports exact rational target weights
`8*(t_hi^2-t_lo^2)` and `-16*(t_hi^4-t_lo^4)`, corresponding to the existing
axial and swirl prefactors; `d rho=8t dt=2x dx`. These are weights, not field
integrals. Matching coordinate records and family/log labels does not prove a
common callable, cutoff owner or source identity.

The first frozen geometry entry stopped at its synthetic fixture: one negative
test accidentally replaced cell 4's interval with the identical interval.
The full failure, including exact frozen sources, is archived; v1 must not be
re-entered. A separate v2 used fresh synthetic records, fixed that test, passed
23 checks and verified the saved records. Its source manifest SHA-256 is
`8354001d77371d9ab20da85c7dff0676fc1d8e97cccd1728d916f6f973d26230`;
entry SHA-256 `599ff123b37ba2c8b7fc4f44003c21b611d323674ad4bd83edde12c4dfb013e8`;
result SHA-256 `ec4bfe1fbb84148a4386ec06376670a3c613d4bf04820d86dedc4dd04f1f93e9`.
The complete v2 archive raw SHA-256 is
`88e5aeecce42f34e758a429862744c096b47cd00051876495ffe1b8b11d74d84`;
the v1 failure archive is
`6e7c00af611821ceba06d5115e3165925b1c40af00a1f05cbdfd0c59aba0a8b1`.

## Target-cell error propagation

The earlier raw-Q feedback and the four finite `P=+2 chi E51` compensation
branches were already computed and are not repeated. Their fixed `x=2`
potential-tail errors were also already computed. The new consumer instead
produces target-cell conditional error enclosures for

`G0 P Rg`, `D P Rg`, `R P G0g`, `R P Dg`, and `R P Rg`,

where `D=(-2 chi) T K0` and `R` is the saved j>=2 kernel remainder. Both `F1`
and output `d_x_F1` are retained through ordinary Z0/Z1. Since `E51` consumes
the input F1 component, both outputs use its `Rg` enclosure; their outer kernel
masses and remainder bounds are field-specific. The new calculation does not
compute another fixed exit or m15 integral.

The 21-check independent small-rational fixture passed, including nonzero
`chi_Z`, the ordinary Z1 product-rule routes, five explicit rational Z0
checks, and the column-1 distinction. The authenticated saved inputs produced
2,496 unaggregated terms. A post-save audit found every coefficient symmetric
and every term marked error-only. Per field, the Z0/Z1 counts are:

| Branch | Z0 | Z1 |
| --- | ---: | ---: |
| `G0 P Rg` | 36 | 108 |
| `D P Rg` | 36 | 144 |
| `R P G0g` | 36 | 108 |
| `R P Dg` | 120 | 480 |
| `R P Rg` | 36 | 144 |

The G0 estimate uses conservative whole-triangle upper masses `x_hi^2/8` for
F1 and `x_hi/4` for its x derivative, separately for each source cell; these
can be wide. D reuses saved target/source `T K0` masses. The outer R bound
uses the saved field-specific kernel majorant and the source-cell width at
most `1/4` once; the width already contained in the incoming Rg bound is
not added again. Ordinary Z1 applies the product rule to P, D and R without
turning an error bound into a signed value or removing an existing error.
The enclosures are conditional on the saved source intervals, kernel bounds,
and physical/gauge interpretation; the remaining same-owner gauge identity is
not proved here.

The target-tail source manifest SHA-256 is
`9c80ee4decb9c81695c5bab9354d0d088476490c1d7c56263aa73944fd676bb5`;
entry SHA-256 `b3e1dc2d5137fef107472af6c9cadc14a03e6301bed346ec5b8ccafcaf13d173`;
result SHA-256 `2f46d930b007a6640859cb4f227757478994fdb7f1f154e8cb17b4cc34a61c159`.
The complete receipt, including exact frozen sources, has raw SHA-256
`7e50ace01f0618a39a89528b487ee1bd959d817ad5eddc0ad18a4c6032d7c9fa`.
Authenticate raw gzip bytes before decompression.

## Remaining work

The target-tail terms cover only the listed potential-kernel truncations;
they do not cover the whole feedback-word remainder, collar transport, pressure
or later physical layers. They have not been integrated into m15, combined
with signed finite terms, or appended to the inherited 212-term ledger. All
six earlier signed inner/lower contributions and inherited errors remain.
Next, bind these distinct error paths and already-saved raw-Q/P finite paths
to a same-source Kchi-aware collar transport, establish the gauge and seam
owner identities, then audit one-time ledger correspondence and compare the
independent negative complete pulse as a function of Z. P1 remains distinct
from m15. No whole-domain five-moment closure or Cartesian NS validation is
claimed.

Workers were requested using the fixed Luna/max role; the spawn responses did
not expose model metadata, so independent runtime-model verification is not
available.
