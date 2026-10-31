# Fixed E00 collar Volterra tail — 2026-10-31

This handoff adds a conditional, source-indexed bound for the `j>=2` Kchi Volterra words applied to the fixed E00 inlet of `CURRENT_N1_E00_FIXED_COLLAR_SEAM_2026_10_31.md` (published predecessor `58c042ad5b50550de32c283ec8277dc366d055c0`). It is not a complete n=1 inlet or collar solution, a terminal five-moment match, or an NS validation. The 212-term m15 ledger and all previously retained errors remain unchanged.

## Bound and assumptions

Let `a=2`, `x(q)=2 exp(hb sc q/8)`, `0<=q<=1`, `0<hb<1/2`, and `L=x-a`. Write `B=(sc q/4) exp(sc q/16)`; then `L<=hb B`. For the fixed E00 inlet `H(t)=Fc+(1-4/t²) Fxc`, `|H|<=|Fc|+hb B|Fxc|`. With `M>=sup|chi|` and `Mz>=sup|chi_Z|` on the source-bound collar, the kernel and its physical x derivative satisfy

```
0 <= K0(x,t) = t/2*(1-t²/x²) <= x-t,
0 <= d_x K0(x,t) = (t/x)³ <= 1.
```

For `A f = -2 integral_a^x K0(x,t) chi(t,Z) f(t) dt`, ordered integration yields `|A^j H| <= |H| (2M)^j L^(2j)/(2j)!` and `|d_x A^j H| <= |H| (2M)^j L^(2j-1)/(2j-1)!`. Summing `j>=2` using factorial majorants gives coefficient bounds

```
F:  (2M)^2 B^4 exp(M B²/2) / 24
Fx: (2M)^2 B^3 exp(M B²/2) / 6
```

with formal `hb^4` and `hb^3` respectively. The Fxc source adds one formal `hb` and one factor `B`. For ordinary Z1, the incoming Z1 term uses the same bound; the terms obtained by differentiating a chi factor are retained separately with the coefficient factors

```
F:  8 M Mz B^4 exp(M B²) / 24
Fx: 8 M Mz B^3 exp(M B²) / 6.
```

These include the `j` possible placements of `chi_Z`: use `(j-2)+2 <= 2*2^(j-2)` and the factorial majorant. They do not divide by `M`, so `M=0` with `Mz!=0` is supported. Exponent bounds use the upper closure `hb<=1/2` only; hb remains a formal power and no midpoint is substituted. The source-bound chi and chi_Z must have zero formal powers. All outputs are symmetric error sectors, not signed evaluations of a tail function. The four formal hb increments for `(F_from_Fc, F_from_Fxc, Fx_from_Fc, Fx_from_Fxc)` are `(4,5,3,4)`.

## Execution and provenance

A new consumer reads the authenticated saved V4 q=0 identity routes, preserving the original E00 source index and paths, and the separately pinned axis-potential chi and chi_Z. It does not use the first-potential output as an inlet. Its 528 input terms comprise F1 and d_x_F1 at each ordinary Z order: 88 Z0 and 176 Z1 per field. Each tail term retains its input ID and its chi ID(s), source classes, ordinary Z route, and formal powers. Axis chi is not the phase shear cutoff.

The fresh synthetic fixture passed **28/28** checks. The one saved actual-source execution emitted **5,632** symmetric tail sectors over the exact seam and three q cells `[0,1/4]`, `[1/4,1/2]`, `[1/2,1]`; **1,408** seam sectors are exactly zero. The result indexes 530 source IDs, including chi and chi_Z. These are enclosures conditional on the saved input intervals and source bindings, not numerical values of the unknown source functions.

| Artifact | SHA-256 |
| --- | --- |
| `experiments/root_st073/n1_e00_collar_j2_volterra_v2_receipt.json.gz` (raw gzip) | `94ca27da585d7e3d252ba94bac1453f037f9601e5fd3d14a9b1b30d2de4f11c5` |
| Frozen source manifest | `82d4e11746bd4ea717baeceec410d5b424258b3e5147fdaafd8932081a3cde34` |
| Entry | `aca536a50d73ea3e3489a2f4685274b47aa98343ae3bed4e0915e7387f0039c6` |
| Complete result | `15e8dc8c5426a77082c9ec39254228630c9730bd994e764342a6932a6363123b` |
| V2 leaf | `af988c8f4032ece433cc901c7a66d926ad93a425cbc4c6c168d12a2dff3ee29e` |
| V2 fixture | `cb41d3be4b1743cccd9485a000de87e98e5809f52913a0eb1faebc0eb78f81f8` |
| V2 runner | `2eb251f1fa58054df6bb1709cfc0c18cbe23c12f52c84b0b2d01502ab9326577` |

V1 failed in its independent synthetic fixture before decoding actual receipts: integer endpoints were passed to the interval constructor without conversion. It was not reentered. Its frozen source, entry, full failure result and code are preserved in `experiments/root_st073/n1_e00_collar_j2_volterra_v1_failure.json.gz`, raw SHA-256 `623cc3617e8c342b6fdef5e4a0cfc8855521e6202e2634ff0bb2fb5896174903`; V2 uses fresh synthetic data and correct endpoint conversion.

## Remaining work

This tail applies only to homogeneous fixed-E00 inlet feedback. Positive-q local forcing, remaining inlet feedback and Q/P branches are not assembled by this result. The old W1/low and new fixed-E00 endpoint records share source ancestry but have no proved term-level function identity or disjointness, so interval containment cannot eliminate any old error. No sector has been appended to the 212-item m15 ledger. Reserved own_swirl, the independently defined negative full-pulse m15 comparison, terminal closure and all subsequent physical validation layers remain open.

Next, use the V4 and this tail receipt to construct a separately frozen, source-bound full fixed-E00 Kchi interface including positive-q local forcing, or prove the old/new endpoint operator identity from common source functions and exact ownership; neither task can be replaced by interval overlap. Preserve ordinary Z1, physical derivatives and all existing signed sectors and error bounds, then continue the reserved/pulse bridge.
