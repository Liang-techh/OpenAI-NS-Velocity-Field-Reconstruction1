# Gamma weighted energy difference — 2026-10-29

Goal ACTIVE / INCOMPLETE. This source-shaped algebra retains the explicit S factor in the Gamma heat integrand and its ordinary Z1. It does not evaluate actual P/F source values or establish global five-moment closure. The prior 212/213/183 entries, six signed channels, once-only reserved contribution, four pulse corrections, `(P-G)+(G-F)`, independent reference offset and all inherited errors remain untouched; P1 is not m15.

For common a, pre and A, set `s=S_F`, `h=S_P-S_F`, `D=D_F`, `K=(D_P-D_F)/h=A*K0`, `X=s*D`, and `Y=D+s*K`. These common-source relations are gates, not proved P/F identities. Formally, with `square=D*(2*pre-a*S*D)`,

```
Delta(S*square) = h*(G0+h*G1+h^2*G2+h^3*G3)
G0=2*(pre-a*X)*Y
G1=2*(pre-a*X)*K-a*Y^2
G2=-2*a*Y*K
G3=-a*K^2
```

The ordinary Z derivative includes both `h_Z*G` and `h*G_Z`, including `a_Z`, `pre_Z`, `s_Z`, `D_Z`, `h_Z` and total `K_Z`. The old fixed-scale K1 is not a total K_Z when scale varies. Missing derivatives are rejected. For differing branches the separate exact zero-order identity retains `DeltaD`, `Delta pre`, and `Delta a`; a general non-common-branch Z1 implementation is still needed.

The source representation separates collar `[0,3]` and exterior `[3,infinity)` and preserves the time-weight difference and its Z derivative. The latter differentiates `expm1` using `exp`, not `exprel`. The integrals are formal, not numerically computed. For `M=Nrel*tail_multiplier*a`, retain `Delta heat=M_P*Delta(S E)+Delta M*(S E)_F` and its four Z1 product routes; `Nrel_Z/Nrel=4Z/(1+Z^2)` applies to the displayed source form. The weighted heat contribution carries `B=-mu*exp(-26)/2`, with the `Delta B` term retained. This is only the heat contribution to the future, not the whole future or a signed m15 integral. The caller must avoid applying the heat minus sign twice.

The frozen one-shot exact SymPy fixture passed all 18 checks, including non-common `a/pre/D`, nonzero Z jets, both heat and weighted four-route product derivatives, the time-weight derivative, and missing-derivative rejection. It used synthetic inputs only; no actual Gamma value, physical integral, certified bound, or ancestor computation was produced. Raw input hashes were checked before gzip decoding. Direct saved Gamma roles include `S=31794`, `Dhat=[31833,31840]`, collar energy `[31909,31910]`, exterior energy `[31913,31914]`; complete-future roles include energy `[32184,32185]`, weighted `[32188,32189]`, and heat difference `[32172,32175]`. These are role references, not function values or P/F equivalence proofs. Both pinned receipt archive hashes matched; all direct declared role pointers were present. The runner does not enforce receipt acceptance or canonical paths, and fixture success alone does not certify source binding: these gates are manually reported from this frozen run and must be tightened in a new leaf before reuse. Caller-provided common-source booleans are declarations, not evidence. Runtime model metadata was unavailable (`null`).

The previous Gamma leaf's exposed `ordinary_z_derivative_kernel` converts inputs before entering its requested mp work precision, so direct use can lose caller input precision. The previous synthetic fixture remains only a synthetic diagnostic; do not reenter that helper or claim a certified bound from it. Existing point arithmetic can round small contributions despite their structural retention.

Next: in a separately frozen new leaf, implement and independently check general non-common P/F Z1, bind both branch owners and their actual S, a, waiting, pre, A, delta and physical radius, and enforce canonical archive/receipt gates without traversing or rerunning old graphs. Keep separate `delta_P_Z` and `delta_F_Z` declarations. Carry all five future terms, the independent F endpoint/radius and all signed error channels before attempting `(P-G)+(G-F)`. No exact heat exterior, stress admissibility, temporal recursion or corrected Cartesian NS validation is claimed.

Frozen checkpoint SHA-256: `3d9788787d7f966d6a92cc1137f6c77dbbdb65aa577f99251b802163e7671106`. Complete one-shot result SHA-256: `949ec01cb201c1c5d1b43edba1af23e4a18510cadd42676e305b13cbd4fe3c45`. The publication receipt embeds exact code, runner, checkpoint, entry, result and prior ephemeral receipt bytes.
