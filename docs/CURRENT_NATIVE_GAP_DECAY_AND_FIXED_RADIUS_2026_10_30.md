# Native gap decay and fixed-radius transformation — 2026-10-30

Status: source-record and conditional mathematics only.  The fixed-radius v1
entry was attempted once and failed before the calculation function body
executed; an independent v2 now has a five-case synthetic result.  This note
does not promote a bound to a function value and does not delete, weaken, or
replace any historical error ledger.

## 1. Canonical native gap decay

The original `whole_Z` producer defines the canonical parameter

```text
mu = (1/1000) * exp(-4 * (exp(40) + 11)),
D  = 2/mu - 4.
```

The source location is the precise path
`whole_Z_all_N_control_connection.py:308-316`.  Since the exponential is
strictly positive and its displayed exponent is strictly negative,

```text
0 < mu < 1/1000,
```

and hence

```text
D = 2/mu - 4 > 1996.
```

The positive defining series for `e` gives

```text
e = 1 + 1 + 1/2! + 1/3! + ... > 2.
```

Therefore, for `D > 1996`,

```text
exp(-D) < 2^(-1996) <= 2^(-128),
```

so the independent gap factor `g(D) = 1 - exp(-D)` has the closed rational
directed enclosure

```text
1 - 2^(-128) <= g(D) <= 1.
```

The displayed inequality is a canonical-original-parameter rational
enclosure, not a function value, midpoint, or evaluated exponential.  It is
conditional on the caller establishing that the `mu` object really is the
source definition above.  It does not prove F ownership, an F/source
identity, a physical-scale identity, an endpoint identity, or removal of any
old error.

The source audit recorded one source-AST match followed by eight leaf calls:
four returned intervals and four rejected invalid inputs as expected.  The frozen provenance pins are:

```text
source checkpoint SHA256: 2c85de04aae4cdf606a7953f5c37eee5295c632c43826ec8a893e9d1e1b1b48a
entry SHA256:             83f53cfb2b734049deab9eb571e6c8e2e1154344cb20ee4b02cc926e0f03fb82
result SHA256:            7a8e2a66c4c7b9979a3c541a6760c5ff309aa7388b1964b6a3eb57d92e658efb
receipt SHA256:           fe275896ac60cc3302147ff16a0cb50a70991db768b0c67f1c0ca5bff5d79445
```

The independent implementation in
`src/openai_ns_reconstruction/gap_decay_directed_interval.py` keeps the two
directed routes separate:

* For small and moderate arguments, it reduces `x` to `t=x/2^k` in
  `[0,1]`, uses adjacent odd/even alternating Taylor partial sums as closed
  rational lower/upper bounds, and raises those bounds by the exact integer
  `2^k`.  No midpoint is formed.
* For large arguments, it uses `exp(-x) < 2^(-floor(x))` from `e>2`, with a
  caller-visible dyadic cap (`max_bits` and `tail_exponent_cap`).  A fixed
  conservative 65,536-bit intermediate-work guard estimates the input
  numerator/denominator lengths, term count, binary term-count length,
  reduction exponent, and final powering scale before any Taylor Fraction
  work.  An over-budget series request falls back to a valid dyadic tail (or the
  full interval when the cap is zero);
  it never allocates an unbounded rational intermediate.  `max_bits` controls
  the tail denominator budget only and is not a Taylor-accuracy promise.

## 2. Fixed-radius jet transformation

Let

```text
r   = (log Rp)_Z,
r_Z = (log Rp)_ZZ,
m_R = m1 + r*m0.
```

After reconstructing the fixed-radius quantity, the ordinary `Z` jets are

```text
Q_RZ  = q1 + r*q0,
Q_RZZ = q2 + 2*r*q1 + (r_Z + r^2)*q0.
```

Here `q0`, `q1`, and `q2` are the actual `Z` jets of the reconstructed
fixed-`y` `Q_R`.  They are not the original fixed-`y` `Q` jets.  The
coefficient roles, source function, F owner, and physical scale remain
unresolved; this is a coordinate/jet transformation only.

The retained `r=1` counterexample is essential.  In that fixture the correct
reconstructed value is

```text
Q_RZZ = -3.
```

Dropping the term `-2*s*r^2*m1` changes the reported result to `-1`.  Thus the
term is not an optional simplification: omitting it changes the jet by two.
The counterexample does not establish an F identity or a physical endpoint
identity.

### Fixed-radius v1 disposition and v2 synthetic result

The fixed-radius v1 entry failed once before any calculation function body
executed.  The runner incorrectly used the `from_rows` row-wise interface and
passed a tuple, while the fixed-radius calculation requires the separate
`r`, `r_Z`, and `r_ZZ` arguments.  This is an interface-entry failure only;
there is no v1 numerical result or accepted jet computation.

The frozen v1 provenance pins are:

```text
source SHA256: 1f39b3d83c85988af206d9333f22cc4e9f276234afff4a430f7a9c195fcac8c4
entry SHA256:  b47074e9986219d36619174ddaa311a0eca6b8ca9b9b69d88103d5edf7c9b887
result SHA256: 7b1f701bc99e531a47dd78ce67df638fbd400be7b1ff70a0b522777a8906156e
receipt SHA256: 8df1d6f26519730301b623d6f4c43b2722d87288f3e3aa371dc2204aec6be34a
```

An independent v2 entry used the corrected argument contract.  Component
`gap_fixed_radius_transport_v2` returned five synthetic cases.  The principal
counterexample returned

```text
Q_R   = -1,
Q_RZ  = -2,
Q_RZZ = -3.
```

The joint nonzero synthetic case with `r=1`, `r_Z=2`, and `r_ZZ=3` returned
the tuple `(-2, -7, -23)`.  The `r=0` cases provide the expected degenerate
checks.  These are transport-algebra fixtures only.

The frozen v2 provenance pins are:

```text
component:    gap_fixed_radius_transport_v2
source SHA256: cf8a8fc176b54d8f803d2278558e76e6340ff7cc575cf65910716aebae2f8f38
entry SHA256:  5ceef90c65431a4630ec2beb54e92915a4ec015ba4dca7f5dee6c3495137c3dc
result SHA256: 8c1d8120356b7c6bdd2f62ceb155c3b8dce64a1717353d97f0781c9ce9650a48
receipt SHA256: 76e598257a27048eb7362653c9132780345bff066e226676f730ac81697a4e86
```

The v2 result is synthetic and does not authenticate an actual F physical
source, F owner, physical scale, endpoint identity, or certified error.  The
next step is to supply authenticated actual inputs and prove their same-source
ownership before attaching this transport to any F or physical claim.

## 3. Frozen unresolved relations and preservation pointers

The following remain explicitly unproved and must stay separate from the
conditional native scalar bound:

* The F gap/source relation, F owner identity, physical-scale binding,
  moving-endpoint identity, signed transport, and final error closure remain
  open.  The current handoffs are
  `CURRENT_F_GAP_REDUCED_HISTORY_CLOSED_FORM_2026_10_30.md`,
  `CURRENT_F_GAP_ROW_REDUCTION_AND_SIGNED_Z_TRANSPORT_2026_10_30.md`,
  `CURRENT_F_GAP_SOURCE_BINDING_LOG_MEAN_2026_10_30.md`,
  `CURRENT_F_Q_JETS_AND_SIGNED_ENDPOINT_2026_10_30.md`,
  `CURRENT_F_GAP_EMITTED_DAG_2026_10_30.md`, and
  `CURRENT_F_PRESSURE_RADIUS_PATH_2026_10_30.md`.
* The `212` cumulative, `213` raw-signed, and `183` reduced ledgers remain
  frozen, together with all six signed channels, all old errors, and the
  single once-only reserved contribution.  See
  `CURRENT_ORIGINAL_NATIVE_ACTIVE_S1_WEIGHTED_IBP_2026_10_21.md`,
  `CURRENT_ORIGINAL_NATIVE_REMAINING_B_SEAM_SOURCES_2026_10_21.md`, and
  `CURRENT_ORIGINAL_N1_SIGNED_LOWER_M15_2026_10_17.md`.
* `P1` is not `m15`; the distinction and the reserved contribution remain
  explicit in `CURRENT_ORIGINAL_CORRECTED_RESERVED_M15_2026_10_18.md`.
  The corrected reserved source is a once-only registration, not a license to
  duplicate or remove an inherited term.
* Prior deviation, failed-entry, and replay/audit records remain preserved as
  do-not-reenter history.  Relevant paths are
  `CURRENT_ORIGINAL_MICRO_AFFINE_SOURCE_IDENTITY_2026_10_20.md` (worker brief
  deviation),
  `CURRENT_ORIGINAL_MICRO_MACRO_CORRELATED_M15_2026_10_17.md` (discarded
  replay recorded in the execution audit),
  `CURRENT_PULSE_CORRELATED_DENSITY_CHART_ATTEMPT_2026_10_28.md` (failed
  exact-chart fixture),
  `CURRENT_COMPACT_M_RECIPE_FAILURE_PULSE_GAP_2026_10_28.md` (entered-once
  failed recipe), and
  `CURRENT_ORIGINAL_NATIVE_LOGRP_SOURCE_BRIDGE_2026_10_28.md` (partial worker
  comparisons without the required checkpoint).

All of those records remain historical evidence and unresolved gates.  This
document adds no numerical F result, no physical completion claim, and no
permission to replay a frozen ancestor or read the prohibited large JSON
archives.

The v2 leaf and runner were frozen with an extra blank line at EOF. Their raw
bytes and hashes are preserved; `git diff --cached --check` reports only those
two frozen EOF lines, and the same check with `-blank-at-eof` passes.
