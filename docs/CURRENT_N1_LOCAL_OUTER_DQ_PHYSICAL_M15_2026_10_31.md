# Local raw-Q outer potential in physical m15 — 2026-10-31

The saved actual-collar local first-feedback word is now carried through one additional outer axis-potential kernel and projected to physical swirl m15. This is a new finite raw-Q local branch; it continues `CURRENT_N1_LOCAL_G5_POTENTIAL_PHYSICAL_M15_2026_10_31.md` (published predecessor `a603e5b02ed3eccffbf50255a4d0c3c3c38ca14b`). It is not the full Kchi feedback, complete m15, reserved/pulse comparison, terminal closure, or an NS validation. The 212-term ledger and all inherited errors remain unchanged.

## Source and operator

The input is the already saved whole-collar local word

```
L0 = [G0 (C0 W0 + B1 W0_Z)]_F,Z0
L1 = [G0 (C0 W0 + B1 W0_Z)]_F,Z1.
```

`n1_actual_collar_low_word_v1_receipt.json.gz` supplies its 137 Z0 and 137 Z1 directed sectors; it excludes the independent core W1 inlet and does not include an all-word tail. These are whole-phase enclosures, not point samples or a resolved q profile. Their original source IDs remain linked to the saved core-inlet and local-forcing parent indices. The local W1 Green integral is not recomputed.

With the separately bound axis potential `chi` and `D=-2 chi K0`, this consumer bounds the local outer word

```
D L0 = -2 K0(chi L0)
(D L)_Z1 = -2 K0(chi_Z L0 + chi L1).
```

The axis potential and its ordinary Z derivative are the saved native axis jets, not the phase shear cutoff. With `d=hb sc q/4` and `x=2 exp(d/2)`, the positive kernel obeys

```
integral_2^x K0(x,t) dt = cosh(d)-1
 <= hb^2 * (sc/4)^2 * exp(sc/8)/2,  0<=q<=1.
```

Each signed source sector is multiplied by a directed `[0, upper]` kernel-mass enclosure; the upper bound is not a function value. The new outer Green measure is applied once and contributes formal `hb^2`, in addition to the powers already carried by L. `q=0` is included as zero.

The physical swirl projection uses the original cutoff mass `5/8` and the radial factor

```
-2 Lambda^-2 integral eta rho Flead (D L) d rho,
rho=4 exp(hb sc q/4),  d rho=hb(sc/4)rho dq,
physical coefficient factor = -5 sc [1,exp(sc/4)].
```

The physical measure contributes formal `hb Lambda^-2` once. Thus the combined new outer-kernel and m15 shifts are `(hb,F0_anchor,Pstar_squared,Lambda) += (3,0,0,-2)`. Ordinary Z1 retains `Flead_Z*chi*L0`, `Flead*chi_Z*L0`, and `Flead*chi*L1` separately; the F0 connection in Flead is not reapplied. Neither P1 nor the P compensation is physical pressure m15.

## Saved execution

The one-entry consumer authenticated the saved low-word, core parent index, axis-potential binding, source family and logs, exact sc/geometry, original cutoff lineage, leading-F source paths and the F0 provenance correction. It did not enter any old producer. The independent synthetic fixture passed **25/25** checks, including nontrivial kernel and radial-integral numerical comparisons, ordinary Z1 paths, physical powers, source paths and error handling. These synthetic quadrature comparisons are cross-checks, not certified values of the actual source.

The actual-source result contains **137 field Z0** and **274 field Z1** sectors; the physical projection contains **274 m15 Z0** and **1,096 m15 Z1** sectors. Its three Z1 paths contain 548 leading-Z, 274 chi-Z and 274 input-Z sectors. All 282 direct source records are referenced. A saved-result audit found no unclassified uncertain terms: axis-potential, core-source, Cauchy-phase, cutoff-value/derivative, mean-phase and F0-connection classes remain represented. All are conditional interval enclosures, not selected signed values.

| Artifact | SHA-256 |
| --- | --- |
| `experiments/root_st073/n1_local_outer_dq_m15_v1_receipt.json.gz` (raw gzip) | `8c249ce31c9e0c6abfc9fd6d7a293a76d824004ee6776a1f6da440bc0a6de512` |
| Frozen source manifest | `ac393c2963c8ceea3d6da7117970fda61b239f35c3aaaac5283a263eb31e1f62` |
| Entry record | `5050a919eed29ed9675b82eaf2071ca9da260bafec199d1ce8f1e2b563ae989d` |
| Complete result | `95c18b251336aa68ec876f759565507ef3dd3cf8253b1c2d199cf2074d3283a6` |
| Leaf | `18b084b5b1c18d2e9adbccd3146339d5ebdd06de7b3eef8f8a8fc11c71a82a5e` |
| Fixture | `ddd584b79be3314835d90a8bd8a498c778567fc13aa6759519763bcebedc7216` |
| Runner | `66060094501b2f4fe33503009050dd8b8f5c766c4f0a9ce53a7bd42deaf1d0de` |

## Limits and next work

This is the outer D applied to the saved local raw-Q W1 only. It is distinct by operator and input path from direct local `Dg`, previously projected K0 W1, and core-exit E10, but no global owner disjointness or cancellation is claimed. The physical `+2 chi E51` compensation required when moving the extracted potential across the operator, inner-potential feedback corrections, `j>=2` outer tails, other Q/P terms and full feedback remain separate and uncomputed here. The output is not added to the 212-term ledger.

Next construct a separately frozen same-source local physical P-compensation contribution, or the inner-potential feedback branch with its Z1 and errors, and preserve it separately from this raw-Q term. In parallel, the corrected reserved own-swirl full-Z signed source evaluator and the independent negative pulse/endpoint ownership bridge remain required. Only a proved same-source operator/endpoint identity can justify a later assembly or remove any old error.
