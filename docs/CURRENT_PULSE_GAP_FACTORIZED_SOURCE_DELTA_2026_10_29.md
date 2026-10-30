# Actual-source conditional factorization of the pulse gap — 2026-10-29

Goal ACTIVE / INCOMPLETE. A new one-entry source-expression builder authenticates the saved P, F and original-kernel inputs, then records their direct source roles without traversing graph edges or evaluating a source function. Its output is a conditional, unevaluated P-minus-F gap expression, not admissible for numerical use. It keeps the nonzero incoming history, radius difference and finite-endpoint term. No ancestor constructor, quadrature or source-value reducer was called.

The saved F selected `C1` Z0–Z3 references are `[275936,275938,275940,275942]` and C2 references `[275946,275948,275950,275952]`. P records C1 `[282525,282527,282529,282531]` and C2 `[282535,282537,282539,282541]`. Both saved current-frame `u` lists record `[265191,265192,265193,265194]`; this metadata alone is not a cross-owner semantic or derivative proof. Original kernel first-row weight references are `[31242,31255]`; they are two different bump weights, **not** Z0/Z1 rows or the energy quadratic's A0/A1. Each reference carries its source archive/owner namespace and is never used as a numerical value.

Under explicit common source, coordinate, ordinary-jet and first-row assumptions, let `dCj=Cj_P-Cj_F`, `cF=A0*C1_F+A1*C2_F`, `dc=A0*dC1+A1*dC2`. The constructed ordinary jets are

```
mF_j = -sum(k=0..j) binom(j,k) u_F,k * cF_(j-k)
dm_j = -sum(k=0..j) binom(j,k) u_k * dc_(j-k),  j=0..3.
```

The named expression DAG applies the normalized Q recovery to these m jets through Z2. It constructs the linear and quadratic signed reduced-S differences, including `delta_squad=-(qF+dq/2)*dq` and its ordinary Z1. With a common Z-independent field scale `T`, the physical density difference is `DeltaS=T*delta_slinear+T²*delta_squad`; the density retains `R_P²*DeltaS+(R_P-R_F)(R_P+R_F)*S_F`, and the endpoint retains `R_P*DeltaQ_y+(R_P-R_F)*Q_yF`. Physical m15 still requires the directed `1/2 dxi/mu` measure and **minus (upper minus lower)** endpoint combination. F's current-radius metadata is retained only as metadata inside an explicitly unresolved historical field-radius owner; it is not used as the actual F radius.

The new source-level expressions still define `dCj` as unevaluated P/F function differences. They do not yet provide a cancellation-stable numerical difference provider, and equal raw reference integers do not establish a common function. Common `u`, first-row weights, mu/delta, field scale, coordinate, ordinary derivatives and source substitution remain proof gates. P/F radii and endpoint ownership are not assumed equal. No source numerical value, certified error bound, pulse integral or terminal identity has been obtained; P-F cannot be removed.

An independent exact-polynomial model substituted synthetic functions into the emitted DAG and passed **20/20** value, ordinary-Z, physical density and endpoint checks, including nonzero radius derivatives. This verifies the constructed algebra at that synthetic case, not all-Z source identity or error bounds. Runtime model metadata was unavailable (`null`).

The initial v1 attempt failed **before entry** because its prior-result SHA literal omitted two characters. No entry marker, source consumer or old calculation ran. Its code, frozen source checkpoint and preflight failure metadata are embedded alongside the corrected, independently frozen v2 code, entry and result in `experiments/root_st073/lei_ren_part1_paper_compliant_current_pulse_gap_source_factorized_delta.json`. V1 was not reentered.

All 212/213/183 inherited ledger counts, six signed channels, once-only reserved contribution, four physical correction roots, independent reference offset and prior errors remain unchanged; P1 is not m15. Global five-moment closure, heat, admissible stress/remainder, actual temporal recursion and corrected Cartesian NS remain incomplete.

Next: establish the same-source ordinary Z0–Z3 `u/C` binding and common first-row/mu/delta/scale premises, obtain a direct source-owned correlated dC provider (for example through the original selector's future-energy difference), and bind historical F radius/endpoints. Transport the two radius-dependent terms and certified errors before using them in the independent terminal comparison. The selector's rationalized `d_ap=d_future/[A2*(sqrt(D_P)+sqrt(D_F))]` is valid only when all its other source inputs and positive-root branches agree; those conditions remain unproved here.

V2 source checkpoint SHA-256: `34ac697f5690d78f2c651abf9172ec67717176b5f39bf4cf9e84b5b0bee1086d`. V2 complete result SHA-256: `a2f1e5edc3ebe60f237b6aafc57fd7bd36536fca7ab5697a4ea3e37a797e56e6`.
