# Radius-sensitive m15 endpoint Z1 — 2026-10-29

Goal ACTIVE / INCOMPLETE. A new caller-owned point-arithmetic endpoint implementation retains the radius ordinary-Z derivative instead of silently fixing it. Its synthetic fixture was frozen and entered once; 12/12 controls passed. No original producer, saved graph, old numerical source or ancestor computation was run. The implementation is not a source identity or a physical evaluation.

For one fixed common endpoint, `src/openai_ns_reconstruction/m15_radius_sensitive_endpoint.py` requires caller-supplied `R,R_Z,U,U_Z,Q,Q_Z,Q_y,Q_yZ,Z,delta,delta_Z`. It forms `A=R^2 Q`, `B=R^2 Q^2`, `C=R^2 U^2`, `D=R Q_y`, `L=1-delta Z^2`, and the existing signed endpoint expression

```
K = A/(2L) + B/4 - Z^2 C/L^2 - D.
```

It uses `L_Z=-delta_Z Z^2-2 delta Z` and applies the ordinary product rule to all four endpoint tensors, including `2 R R_Z Q`, `2 R R_Z Q^2`, `2 R R_Z U^2`, and `R_Z Q_y`. No `R_Z` or `delta_Z` is defaulted to zero. The difference between the full derivative and the same derivative with only `R_Z` set to zero is

```
R_Z * (2 R * (Q/(2L) + Q^2/4 - Z^2 U^2/L^2) - Q_y).
```

This shows the exact term the historical `a.fixed(radius)` specialization would omit if the applicable radius were not Z-independent. It does **not** establish that the actual F radius has nonzero R_Z, that it is independent of Z, or that F/P/G endpoint functions agree. The callable rejects nonpositive R, nonfinite inputs, point L=0, and interval inputs; its rounded inputs and outputs are not certified. It assumes a fixed common endpoint or caller-provided total endpoint jets; a moving endpoint's chain rule must be established externally.

The fixture uses Fraction polynomials and a separate rational finite-difference reference, exercises nonzero R_Z, a nonzero delta_Z, signs, the R_Z=0 special case, singular input rejection and the full omission identity. Its nonzero omitted term is a synthetic counterexample to dropping the radius term without a premise, not an actual physical estimate. There is no actual source value, signed m15 integral, seam proof, source-provider equivalence, certified bound or error reduction.

Static source review traces historical Gamma's mu formula to constants, while the historical `logRp` expression includes an opaque `phase_rc` source handle. Absence of a literal Z in the producer does not prove that handle's semantic Z independence. Historical F radius/finite-endpoint owners and capture-free inlet substitution remain unbound; the scoped current geometry radius metadata cannot substitute for F's source owner. P1 and M are distinct from m15. Keep `(P-G)+(G-F)`, four physical pulse corrections, 212/213/183 inherited counts, six signed channels, once-only reserved term, independent reference offset and all old errors.

A read-only worker search accidentally included `--glob '*.json'` across `experiments/root_st073` while looking for Gamma schedule references. Its output showed no `Rm_loop_density_identity.json` filename or content, but the broad glob may have caused `rg` to enumerate or read that prohibited file; this cannot be excluded. The worker did not decompress archives, parse a graph or run a source consumer, and the present endpoint result does not rely on that search. The exact command was `rg -n -i --glob '*.md' --glob '*.txt' --glob '*.json' 'Gamma\.schedule|schedule\s*\[.?mu|schedule\s*\[.?logRp|Z[- ]?independence.*Gamma|Gamma.*Z[- ]?independence|historical.*Gamma' docs experiments/root_st073 2>$null`. Do not repeat it; explicitly exclude the prohibited file in future searches.

Next: bind the actual F radius, `mu`, `logRp`, endpoint substitution and ordinary Z1 on a common source chart, including capture freedom, moving-endpoint behavior, measure/orientation and Pstar; supply the new full-radius expression only with those source-owned jets. Prove the four same-source seam identities and certified error transport before any cancellation. Global five-moment closure, heat/stress, order-dependent time recursion and Cartesian NS validation remain incomplete. Runtime model metadata was unavailable (`null`).
