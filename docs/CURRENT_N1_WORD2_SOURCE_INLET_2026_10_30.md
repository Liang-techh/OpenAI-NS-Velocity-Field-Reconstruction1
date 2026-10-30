# Signed second-Z first-word source at the inlet

At the actual core/first-collar seam `rho=4`, `s=0` and the saved narrow
real Z interval around 0.371, the new source-jet adapter computes the signed
known part of the second ordinary Z derivative of the first-word feedback
source. It uses the saved six-component `W0` jets through Z2 and the actual
first-inlet raw `B0/B1` matrix jets through Z2. The same source family,
four formal source-log anchors, physical inlet and selected Z interval are
checked before combining them. Core-cell quadrature and source uncertainty
remain inside the original directed coefficient intervals.

Write `C0=B0-Akin`, `W1=(I+J)G(C0 W0+B1 W0_Z)`. At fixed radial coordinate,
the differentiated feedback source is

```
C0_ZZ W0 + 2 C0_Z W0_Z + C0 W0_ZZ
+ B1_ZZ W0_Z + 2 B1_Z W0_ZZ + B1 W0_ZZZ.
```

The new leaf `n1_word2_source_integrand.py` computes the first five terms
using the saved ordinary Z derivatives. It retains every signed product
sector and its own combined source-log interval, without materializing the
source exponentials or choosing a midpoint. It returns `B1 W0_ZZZ`
separately as an uncomputed, non-omittable term. The one-entry actual-data
run returned 3, 28, and 39 known product sectors for dynamic rows 4, 5, 6;
the missing term has respectively 0, 3, and 3 formal B1 factors. The row 4
feedback source has no B1 contribution.

The values are **integrand enclosures at the seam**, not the Green integral
for `W1_ZZ`. The accepted canonical source constructs ordinary Z0–Z3 function
graph handles for the finite words, including `W0_ZZZ`; those exact function
definitions are distinct from a computed signed numerical sector. The saved
coupled numeric input here stops at `W0_ZZ`, and its cell ledger does not
contain the higher-order sector data needed to integrate the full expression.
The previously certified Cauchy bound remains a separate positive error
bound and is not replaced by a midpoint or by the partial signed source.

The exact source manifest, entry record, directed endpoint tuples and
individual product factors are in
`experiments/root_st073/n1_word2_source_inlet_v1_receipt.json.gz`; the runner
is `scripts/run_n1_word2_source_inlet_once.py`. No m15 contribution was
recomputed or appended. The 212 cumulative terms, six signed inner/lower
contributions, reserved term and inherited errors remain unchanged; P1 and
m15 remain different quantities. Neither an all-Z terminal identity nor a
completed positive-order field follows from this local result.

Next, compute actual source-bound `g_Z3`/`W0_Z3` sector data on radial cells
and the corresponding W1_Z2 Green integral without replaying accepted
ancestors. Use the saved exact graph definitions and existing source data as
inputs, keep the core-to-collar history and fixed-coordinate derivative
semantics, and separately retain the Cauchy fallback and all-word tail. The
subsequent signed m15 propagation still requires same-source history and
endpoint identities.
