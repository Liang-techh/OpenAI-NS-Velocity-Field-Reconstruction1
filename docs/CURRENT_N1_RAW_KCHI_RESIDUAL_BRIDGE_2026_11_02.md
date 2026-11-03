# Conditional raw--Kchi residual bridge — 2026-11-02

## Status and scope

The new `src/openai_ns_reconstruction/n1_raw_kchi_residual_bridge_v1.py` is a caller-owned conditional source-operator leaf. It applies caller-supplied operators; it does not itself supply actual source functions, a numerical physical result, an error bound, or a completed finite-physical-to-canonical identity. It does not alter the 212-term ledger, the six signed inner/lower terms, inherited errors, or the separate `P1` datum; `P1` is not `m15`.

The same live source owner, domain, units, initial data, full `K` resolvent, and fixed-`x` `G_Z=0` are caller premises. Labels, matching source IDs, interval overlap, or equal-looking values do not prove those premises. The implementation records `actual_source_verified=False`.

## Conditional operator identity

With the same raw source owner and caller-supplied operators,

```text
d_raw = W - G(g + QW)
d_chi = W - K(g + (Q-A)W)
A     = -2*chi*E51
K     = (I-GA)^-1 G
```

the retained resolvent relation is

```text
(I-GA)d_chi = d_raw
d_chi = (I+KA)d_raw
```

At fixed `x`, when the caller certifies `G_Z=0` and supplies the ordinary `Z1` residual,

```text
d_chi_Z = (I+KA)(d_raw_Z + G A_Z d_chi)
A_Z     = -2*chi_Z*E51
```

For `Q=C0+B1*d/dZ`, the ordinary derivative is kept as

```text
d_Z(QW) = C0_Z W + C0 W_Z + B1_Z W_Z + B1 W_ZZ
```

In particular, the `B1*W_ZZ` term is explicit and is not replaced by a lower-order or interval proxy. The bridge applies operators in the non-commuting order shown above; `K*A` cannot be replaced by `A*K`.

## What the synthetic witness establishes

The companion `src/openai_ns_reconstruction/n1_raw_kchi_residual_bridge_fixture_v1.py` is a fresh, non-commuting exact `Fraction` witness. The saved receipt records one 41/41 check set over 2 × 2 and 3 × 3 cases, with nonzero `chi_Z` and nonzero coefficient `Z` derivatives. It rejects a false owner, missing units or initial data, wrong `A`/`A_Z` signs, unknown or nonzero `G_Z`, missing Z1, and an incomplete `B1*W_ZZ` contract.

This is only an algebraic witness. It does not evaluate a real raw source, a finite Kchi packet, an old producer, a historical owner, or a mathematical reconstruction.

The source and receipt pins are:

- source checkpoint manifest: SHA-256 `e384d9491346faededcf4cc6af3866c76349af1de8c824e8d47804a9c9bb9324`;
- leaf: `src/openai_ns_reconstruction/n1_raw_kchi_residual_bridge_v1.py`, SHA-256 `616fdddb607fbfbff7fca942b68329017656578c2fa6d7b250cedd2dea33937c`;
- fixture: `src/openai_ns_reconstruction/n1_raw_kchi_residual_bridge_fixture_v1.py`;
- entry: `fdb8b1a0e5aeab03ad3cef5190d45c041c5b36aabfd25f1d227b5f42574bece7`;
- result: `9f7564785800a9731da34a1809bb856dd7e4e2a9092347f9476dbcc8f6086812`;
- full receipt: `experiments/root_st073/n1_raw_kchi_residual_bridge_v1_receipt.json.gz`, SHA-256 `18681f0532f41e502ddc70610624b99bf7aed21837019bfc9be1ac54c9fed7df`.

No real source function was evaluated in this source review, and no old production, import, or mathematical rerun is claimed here.

## Boundary with existing finite Kchi code

The bridge must not promote the finite `G0`/`D` pieces to the full resolvent `K`. In particular, the existing first correction `D=(-2*chi) T K0` is only a finite first word. It is not evidence for `K=(I-GA)^-1G`, nor does it certify the omitted Kchi words or their tail.

Likewise, the natural-component frame must not be identified with the raw physical frame. The raw-to-forward scale and its connection, including the `F0_Z/F0` contribution, must remain explicit. Dropping the `F0_Z` connection changes the source/operator identity even if formal powers or field labels coincide.

The bridge also does not by itself connect the existing finite `W0`, `DeltaW0`, `E00/E01/E10/E11`, or `P` packets to canonical raw `g`, `C0`, and `B1`. Their source ownership, units, support, and residual provenance remain separate until an actual adapter proves the relation. No residual numerical value, error upper bound, all-word tail contraction, or historical cancellation is established here.

## Next concrete implementation

Build one source-owner adapter that supplies the actual raw `g/C0/B1` and the Kchi `G`, `K`, `A`, and `A_Z` operators to this leaf. The adapter must certify, on the same domain and ordinary `Z0/Z1` jets:

1. one live owner and one initial datum for raw and Kchi packets;
2. raw physical units, including the `F0` amplitude connection and the `Pstar`/`Lambda` factors;
3. the full resolvent relation for `K`, rather than substituting finite `G0`/`D` pieces;
4. `G_Z=0` only in the frame where it is actually true;
5. `Q=C0+B1*d/dZ` together with `C0_Z W+C0 W_Z+B1_Z W_Z+B1 W_ZZ`;
6. explicit finite-prefix, omitted-word, and all-word-tail residual packets, without merging them or deleting inherited errors.

Do not join the ledger or claim cancellation until this source-owner and unit-level identity is proved from the actual source objects; source labels alone are insufficient.
