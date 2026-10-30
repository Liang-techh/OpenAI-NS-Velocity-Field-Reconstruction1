# Original gap scaled field — 2026-10-29

Goal ACTIVE / INCOMPLETE. A new point interface implements the original gap formulas with caller-owned ordinary Z0–Z3 jets, retaining the nonzero accumulated mean, auxiliary Q and Q_y. It returns normalized `m, Q, Q_Z, Q_ZZ, Q_y, Q_yZ` and signed m15 reduced-S linear/quadratic components with Z1. The original field, kernel and pressure producer bytes are pinned; static reading supplies the mathematical derivation, but this callable does not authenticate supplied jets as actual historical F or current physical P functions. It is not a certified bound or a physical-source evaluation.

For `lambda=1/2-mu`, `b=1/2+mu`, `xi=anchor+mu*offset`, and the original first moment row `(A0,A1)`, the normalized gap mean is `-u*(A0*C1+A1*C2)`. The single common logarithmic field scale is retained as separate terms

```
(11/2-anchor)/mu,  -offset,  -13,  -3*(4mu)^(-2/3),
-(4/3)log(mu),  (2/3)log(4)-(1/2)log(6),  logPstar.
```

The endpoint coordinate is stored as anchor plus offset rather than relying on the rounded sum. The input `xi` and parameter independence declarations are caller premises; interval values are rejected and the displayed rounded xi is diagnostic only. The output includes separate, unsummed physical scales for `R*Q_y`, `R²*Slinear` and `R²*Squadratic`, with `logRp` and `logPstar` explicitly supplied. The first inverse-mu terms are respectively `(11/2)/mu`, `(11/2+anchor)/mu`, and `11/mu`; the linear R² term additionally retains `offset`. No scale is exponentiated by the point interface. A single output factor must not be multiplied by its field scale again.

For `L=1-delta*Z²`, the gap has `U=U_y=0`, `m_y=-m`, and `Q_y=-Q`. Thus `Slinear=(1-delta)Z Q_Z/(2L)` and `Squadratic=-Q²/2`; their ordinary Z1 includes the derivative of L. The physical contribution is `1/2 integral R²*S dxi/mu - [R*Q_y]`, so the integration, endpoint signs and nonzero gap history are retained rather than claimed complete. Input mu, delta, coordinates, log scales and row weights must be Z-independent for the returned ordinary derivatives; these declarations are not source proofs.

A single frozen fixture returned **38/38 checks**, including a medium-mu comparison against the original unfactored expression, independent finite differences of Q_Z, Q_ZZ and linear S_Z, consistency of all three physical log scales, and a `mu=10^-400` upper-endpoint/outside-domain check. It uses synthetic u/C jets and a synthetic mu. The real source mu is recorded as `exp(-4*(exp(40)+11))/1000` but was not numerically instantiated. There is no actual P/F signed integral, quadrature error certificate, source-owner identity, seam cancellation or removal of old errors. Runtime model metadata was unavailable (`null`).

Preserve 212 cumulative / 213 raw / 183 reduced entries, six signed channels, once-only reserved contribution, four physical corrections, `(P-G)+(G-F)`, independent reference offset and every old error. P1 is distinct from m15. Global five-moment closure, exact heat, stress/remainder, order-dependent recursion and independent Cartesian NS validation remain incomplete.

Next: bind actual F and P gap inlet and selected coefficient ordinary-Z jets to the same-source formula, prove radius/mu/scale independence and common coordinate/measure with capture-free substitution, then feed each independently scaled linear, quadratic and endpoint term to a directed signed comparison with certified error transport. Extend to main/end before asserting terminal matching. Do not reenter this or any old leaf.

Source checkpoint SHA-256: `8415e90737ae14cb45ab1f2cf756fb2322b5eeb5d7197ed94269b8297322fead`. Complete result SHA-256: `ff9215f1f5935f17cf62394a55cfd29261c38ea258d8db699d42abfd065aaf9b`.

Independent read-only review confirmed the derivative and scale algebra. The runner docstring incorrectly says it performs a static fragment audit; the executed runner only pins the source bytes, as its saved result reports. `scale_coordinate_and_row_weights_Z_independent` is a required caller declaration, not a proved source property. The fixture does not test the lower-domain edge at tiny mu or mu=1/2. Entry creation is exclusive and flushed, but does not fsync; the runner does not claim protection from concurrent file changes between hash checking and import.
