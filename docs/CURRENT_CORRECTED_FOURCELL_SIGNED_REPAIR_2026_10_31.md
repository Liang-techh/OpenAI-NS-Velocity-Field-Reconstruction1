# Conditional signed four-cell lower-endpoint repair — 2026-10-31

The saved physical lower `V_y` source enclosures now feed a directed, factored interval calculation on all four `Z` cells. The calculation returns four Z0 terms and eight separate ordinary Z1 terms. Its source-basis and parameter connections are **explicit caller hypotheses**, not a proved native/corrected endpoint identity. No point value, reserved response, complete m15, pulse match, or global five-moment closure is claimed. The historical 212 terms and all prior errors remain unchanged.

## Definition and directed enclosure

On a fixed cell and under one physical source, write `Y0=V_y`, `Y1=partial_Z V_y`, `L=1-delta*Z^2`, `F=2Z/L`, and `Fz=2(1+delta*Z^2)/L^2`. With the original positive radius `R` and `Pstar` independent of `Z`, the lower-endpoint excess removed by the corrected source is

```
repair0 = R*(Pstar-1)*F*Y0
repair1 = R*(Pstar-1)*(F*Y1 + Fz*Y0).
```

The new implementation retains the two Z1 terms separately, including their signed source coefficient intervals and formal positive scales. It composes powers and offsets only within each cell; source exponents are opaque coordinates whose basis names and cross-cell object identities were not serialized. No cross-cell sum or independently boxed cancellation is performed. The original physical `Y` already contains its `Pstar` factor. In the identity below, the additional `Pstar` is part of the **old excess** `Pstar-1`, not a new normalization of `Y`:

```
Pstar-1 = Pstar*beta,   beta = 1-1/Pstar.
```

For `p=log(Pstar)>0`, the directed calculation uses `beta in [p/(1+p),1]`, because `exp(p)>=1+p`. The input `p` interval is checked to have strictly positive lower endpoint; the ratio is rounded outward. This is an enclosure of the same shared factor, not a selected beta value, and it avoids materializing an enormous exponential. The output is a factored enclosure for physical repair terms; it has not been divided by the common m15 unit or inserted into a ledger.

## One-entry run and limits

The frozen consumer authenticated the previous four-cell source receipt and corrected source report before decoding. It checked the saved `[0,.5]` lower rows and proof against the corrected input, retained the other cell records, used the saved delta and log-Pstar intervals, and verified positive coefficient enclosures for `R` and `Pstar`. In this arithmetic run all four `L` intervals have strictly positive lower bounds. This is a statement about the supplied delta interval; its identity with the physical parameter on every corrected cell remains a premise.

The new independent synthetic fixture passed **23/23** checks. It exercises whole-cell geometry and sample containment, nonzero ordinary Z1, source-specific signs, half-integer formal powers, a shared beta factor, and rejection of an unproved Pstar or a denominator reaching zero. The actual output has four Z0 terms and eight Z1 terms, each retaining a signed coefficient interval and a separate positive formal scale. No extra force is introduced and no error is removed.

The saved records do not themselves serialize the common live basis/ledger object for `R`, `Pstar`, `Y0`, and `Y1`. The diagnostic also did not find explicit saved flags proving `delta_Z_independent`, `R_Z_independent`, and `Pstar_Z_independent` across this corrected source. The runner therefore keeps the genuine diagnostic separately and invokes the arithmetic leaf under an explicitly labeled **assumptions-only certificate**. Its `passed` flag means the computation succeeded under that certificate, not that those physical owner hypotheses were proved. In particular, `[0,.5]` has the exact corrected-row/proof binding; the three other cells are not independently rebound to the complete corrected reserved source.

## Remaining endpoint identity

Static source formulas give, with `W=U-m`,

```
N = 2Z*Y - (1-delta)*Z*W - (1-Z^2)*W_Z
Q_y = N/L.
```

The native `Q_y` formula and the corrected `qplus-q` formula agree algebraically under common jets. If both use the same fixed physical radius and the same delta, their difference is `R*DeltaN/L`, with ordinary derivative

```
DeltaN_Z = 2*DeltaY + 2Z*DeltaY_Z
           -(1-delta)*(DeltaW+Z*DeltaW_Z)
           +2Z*DeltaW_Z -(1-Z^2)*DeltaW_ZZ
partial_Z(R*DeltaN/L) = R*(DeltaN_Z/L + 2*delta*Z*DeltaN/L^2).
```

This becomes zero only after proving the required native/corrected `U,m,V_y` C1/C2 source jets, common radius and incoming-history identity. Native initial `R Q_y` `[129567,129568]` and corrected lower `[53450,53462]` remain distinct ports. Different integer handles, coincident boxes, or shared source family do not prove the identity; no endpoint cancellation or error deletion is performed here. The separate Rm-to-Rp transport, reserved own-swirl binding, and independent negative-pulse comparison remain open.

## Reproducibility

| Artifact | SHA-256 |
| --- | --- |
| `experiments/root_st073/original_corrected_fourcell_signed_repair_v1_receipt.json.gz` (raw gzip) | `5876bf565709eab3e280b3e49b68cff57c30e5c772cb22a0c6c9f8b3f3b65b78` |
| Frozen source manifest | `d2bfd728eafd1b2a0552335c423be6bb899bd076181bf213048f6964e19c5cc2` |
| Entry | `266c054928701653a42cefdb77a2dfa9c27644bb138ec1280c923f718634dac4` |
| Complete result | `94d26ba64112c01d58f50144b8d19ea25d683d9ed1f83cd3e8872addfae497c8` |
| Leaf | `330c15e10a4fe7fe1ca114c7b7cb3d5e1303f70c75e05a40038150bd54bda537` |
| Fixture | `ff82b2b0652833b2a26388af31f640d3fab7acd7b3efe9c3f4c66642b8917170` |
| Runner | `17117d37d81e3a794d6579acbcc88c7d46e9f6eedf5c88e7edbfa9b728eb55aa` |
