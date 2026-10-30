# Current F pressure integral roles — 2026-10-30

Status: **ACTIVE / INCOMPLETE**.  This is a source-aware handoff for the two
saved F pressure paths.  One bounded consumer run is recorded as
`returned`; this update does not rerun it.  The returned receipt has manifest
SHA-256
`8843a60f00b9e8417b5722ab7cd2244dd8c8ca5e03b817350c6636dfc28eb2b4`, result
SHA-256
`ecbaeda54bdfa82952e09ed1639498ea0bb057716ffc087ab078aa31673c780b`, and
`7` authenticated inputs.  The receipt is source/role evidence only: it does
not evaluate a source row, an integral, a numerical value, or a finite sample.

## Saved paths and what they establish

The saved path archive contains two roots:

| root | saved path shape | current structural result |
|---|---|---|
| `F.pressure.m15[Z0]` (`279426`) | 170 nodes, 304 edges | reaches protected rewritten target `33982` |
| `F.pressure.m15[Z1]` (`279542`) | 198 nodes, 357 edges | reaches protected rewritten target `33982` |

The target is not dereferenced.  The complete F archive JSON packet was
decoded, but there was no graph traversal, old-ancestor replay, source-value
selection, numerical evaluation, or integral evaluation.  The paths therefore
establish only bounded operand reachability.  They do not establish that the Z
1 graph is the ordinary-Z derivative of the Z 0 graph, that both roots use the
same integral definition, or that either root owns the physical F radius.

The bounded direct read set is exactly six graph rows:

- integrals `279383` and `279479`;
- integrands `279382` and `279478`;
- shared product `34721`; and
- its one non-protected sibling `34710`.

Only the sibling operation `positive_quotient` was checked.  That does not
prove the sibling is the physical Jacobian, nor does it prove a radius owner.

## Integral chart, measure, and original I 1 roles

The original pressure producer
(`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_functions.py`)
loops over `main`, `gap`, and `end`.  It takes the chart query `x` and the
chart field radius `R`, then uses

| chart | integration variable and bounds | radial Jacobian | original I 1 measure |
|---|---|---|---|
| `main` | `xi` from `0` to `11` | `R / mu` | `(R^2 / mu) dxi` |
| `gap` | `xi` from `11` to `13 - 4*mu` | `R / mu` | `(R^2 / mu) dxi` |
| `end` | `s` from `-4` to `0` | `R` | `R^2 ds` |

For an `I0_*` row the producer uses only the radial Jacobian.  For an
`I1_*` row it multiplies by one more `R`; this is why the `I1` rows carry the
`R^2 / mu` or `R^2` measure.  The original `I1_q`, `I1_q_Z`,
`I1_q_squared`, `I1_U_q_Z`, `I1_U_squared`, and `I1_U_U_Z` roles feed the
reduced `S1` expression.  The original pressure reduction is

```text
m15 = 1/2*S1 - (R_q_y at upper - R_q_y at lower).
```

Thus an integral-role certificate must preserve the I 1 row identity, its
chart variable, its physical measure, its signed orientation, and its
connection to the `S1` row.  A common target node or a common row name is not
enough.

For each value/Z 1 pair the old producer stores the same variable and bounds
with `ordinary_slow_Z_derivative_order = 0` or `1`, and declares
`endpoint_independent_of_Z = True`.  This declaration is a property asserted
by that producer.  It is not a source proof for the F rewrite, and it does
not by itself prove that the radius, endpoint fields, or source substitution
are Z-independent.

The returned F mapping checked all `36` source references.  The two selected
F rows are uniquely mapped as

```text
Z0 -> gap / I1_q_squared[0]
Z1 -> gap / I1_q_squared[1]
```

Their formal m 15 coefficient is `c = -1/4`, and their formal integrand role
is `Q^2`.  Both rows report the same metadata:

```text
chart  = gap
measure = (R^2/mu)*dxi
variable = original_own_pulse_gap_xi
ordinary order = 0 / 1
endpoint_independent_of_Z = True
physical_Pstar_not_reapplied = True
```

Their lower and upper values are opaque references `4` and `33889` in both
rows.  Equal opaque references establish only equal stored references; they do
not establish endpoint ownership, axial constancy, or a source-level
`xi_u_Z = 0` identity.

The shared product row `34721` occurs as `arguments[0]` in both selected
integrand rows.  Its `arguments[0]` is the opaque radius target `33982`; the
non-protected sibling `34710` is a `positive_quotient`.  These structural
facts do not elevate the sibling to a proved Jacobian or the target to a
proved F radius.

## Complete moving-endpoint derivative

For an F chart integral

```text
I_F(Z) = integral from xi_l(Z) to xi_u(Z) of f_F(Z, xi) dxi,
```

the required formal derivative is

```text
I_F_Z = integral of partial_Z f_F dxi
        + f_F(Z, xi_u) * xi_u_Z
        - f_F(Z, xi_l) * xi_l_Z.
```

For the saved gap chart,

```text
xi_l = 11,
xi_u = 13 - 4*mu_F(Z),
xi_l_Z = 0 only after the literal lower endpoint is source-bound,
xi_u_Z = -4*mu_F_Z.
```

The upper boundary term must remain explicit until the F source owner proves
`mu_F_Z = 0`.  A source-owned, structurally Z-free canonical `mu` expression
could discharge that term, but the canonical expression alone does not prove
that the protected F node has that owner or that the rewrite preserves the
identity.

If a pressure term is written as `c(Z) * I(Z)`, its derivative also has the
full product rule

```text
(c*I)_Z = c_Z*I + c*I_Z.
```

For the returned `gap / I1_q_squared` role the formal coefficient is the
constant `c = -1/4`, so its formal coefficient derivative is `c_Z = 0`.
That only removes the coefficient product-rule term.  It does not remove the
unproved integral boundary terms, the F source derivative, capture-free
transport, or endpoint transport.

The `c_Z*I` term cannot be dropped when `c` contains `mu`, `delta`, `L`, a
radius factor, or an endpoint normalization.  For example,
`L = 1 - delta*Z^2` requires
`L_Z = -(delta_Z*Z^2 + 2*delta*Z)`; `delta_Z = 0` must be source-proved,
not supplied as a caller convention.

At a moving endpoint, an endpoint field also needs the total chain rule:

```text
d/dZ f(Z, xi_u(Z)) = (partial_Z f)_xi + f_xi * xi_u_Z.
```

Consequently, the endpoint term `R*Q_y` requires both the fixed-coordinate
Z row and the source-owned `xi` derivative of `Q_y`, as well as the total
radius derivative.  The current candidate leaves these `partial_xi` rows and
the F endpoint transport unresolved.

## Radius and source-ownership requirements

The candidate chart relation is

```text
R_F = exp(logRp_F + xi/mu_F),

(partial_Z R_F)_xi
  = R_F * (logRp_F_Z - xi*mu_F_Z/mu_F^2),

d_Z R_F(Z, xi(Z))
  = R_F * (logRp_F_Z + xi_Z/mu_F - xi*mu_F_Z/mu_F^2).
```

These formulas become a certificate only when all of the following are
source-bound:

1. the radius used by the pressure integrand and by both endpoint tensors is
   the same F field owner;
2. `logRp_F`, `mu_F`, and `delta_F` have authenticated Z rows and the fixed
   `xi` derivative semantics;
3. the substitution is capture-free and preserves the derivative operation
   tree, including protected leaves;
4. lower and upper endpoint owners, chart orientation, measure, and `Pstar`
   are transported from that same source;
5. the Z 1 integral row is the exact ordinary-Z derivative of the Z 0 row,
   including all boundary and endpoint-chain terms.

The separately emitted F `pulse_coordinates.radius` is current-geometry
metadata.  It is not yet proved identical to the protected radius node in the
rewritten pressure graph.  Similarly, using a historical-template `logRp`
with F `mu` is not a same-owner identity.

At the `gap` / `end` seam, matching formal log-radius expressions is only
valid after matching the source owner, `mu`, `logRp`, coordinate, orientation,
and endpoint.  Converting the measures additionally requires
`dxi = mu*ds`; reachability at the seam does not prove this transport.

## Preserved differences, ledgers, and errors

No cancellation is admitted in this leaf.  Preserve all of the following
while the source-owner and derivative bridges are open:

- the `212` cumulative, `213` raw-signed, and `183` reduced ledgers;
- six signed channels, the once-only reserved contribution, the `17` K terms,
  the `2` outer-B terms, and the archived `16` J terms where those packets
  apply;
- the four signed same-source m15 seam corrections.  The local-minus-global
  endpoint form remains

  ```text
  [K(global.lower) - K(main.lower)]
  + [K(main.upper) - K(gap.lower)]
  + [K(gap.upper) - K(end.lower)]
  + [K(end.upper) - K(global.upper)].
  ```

- the independent reference `offset` and the physical pulse correction;
- the formal difference `(P-G) + (G-F)`, without identifying P, G, and F;
- `P1 != m15`; a P 1 row is not a substitute for the signed m15 source;
- every prior quantitative and analytic error packet.

The four corrections and the `(P-G) + (G-F)` decomposition remain signed,
unevaluated source terms.  Bounds, midpoint overlap, finite samples, or a
common protected node cannot remove them.

## Required certificate status

A new source-aware integral certificate should report, separately for Z 0 and
Z 1, the chart/measure fingerprint, ordinary-Z declaration, exact integral
provenance, bound derivatives, radius owner, capture-free substitution, and
endpoint transport.  It should fail closed with statuses such as
`moving_bound_formal_only`, `unresolved_F_radius_owner`, or
`ordinary_Z_transport_unproved` until those fields are source-authenticated.
A Z 1 integrand derivative certificate alone does not discharge the
`endpoint_independent_of_Z` declaration: bound ownership, bound derivatives,
and capture-free transport require separate evidence.

The eight synthetic guards in the saved consumer passed.  They are guard
tests only.  The source gates remain false for source identity,
capture-free substitution, actual ordinary-Z transport, measure ownership,
endpoint transport, and error transport.  The next substantive work is to
authenticate the Z 1 minus Z 0 integrand relation and the F source for
`mu_Z` at the upper bound, including its Leibniz boundary contribution.  A
real m 15 comparison remains incomplete.

This update has not rerun the project producer or checker.  No source value,
quadrature, finite-sampling derivative, or new prior mathematical result is
admitted here.  All previously retained ledgers, corrections, offsets, and
error packets remain unchanged.

Published receipt `experiments/root_st073/f_pressure_path_integral_roles_v1_receipt.json` preserves the new code, runner, checkpoint, entry, result and prior result bytes; SHA-256 `7be3f0d22defb8a4f067891566bec4348e59103c9d63a8cb7d45e3385f4843c7`. The tracked F archive and source anchors are pinned by raw SHA in the checkpoint. Do not re-enter this component.
