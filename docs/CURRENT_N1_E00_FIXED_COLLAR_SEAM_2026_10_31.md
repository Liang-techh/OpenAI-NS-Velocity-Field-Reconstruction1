# Fixed E00 core-to-collar seam — 2026-10-31

This is a conditional, finite-sector handoff, not a completed n=1 inlet, a five-moment terminal match, or a Navier–Stokes validation. The published predecessor is `11090af3936f78ddd1772d703a8469f9602a0161` and `CURRENT_N1_CORE_E00_FIXED_FOUR_FIELDS_2026_10_31.md`. The 212-term m15 accumulation and its existing errors are unchanged.

## What was implemented

`n1_e00_collar_fixed_seam_v4.py` propagates the fixed-target K0 first-feedback E00 inlet at `x=2`, `rho=x^2=4`, separately from the old coarse W1 and Q/P records. It retains ordinary `Z0` and `Z1`, original sector intervals, source IDs and indexed parent paths, and four formal powers `(hb, F0_anchor, Pstar_squared, Lambda)`. The 624 saved F1/Uz1 sectors are references to the earlier fixed-target profile, and the 1,008 sectors for K1/P1/d_x_F1/d_x_Uz1 are the fixed E00 supplement. They are input records of this separate calculation, not new additions to the 212-term ledger. P1 here is a state component, not physical pressure or m15.

With `x(q)=2 exp(hb sc q/8)` and `d=hb sc q/4`, the homogeneous transfer uses `F=Fc+(1-exp(-d)) Fxc`, `U=Uc+d Uxc`, `K=exp(-d) Kc-(1-exp(-d)) Uxc`, `P=Pc`, `Fx=exp(-3d/2) Fxc`, and `Ux=exp(-d/2) Uxc`. At `q=0` this is the exact algebraic identity on all six components, and, for the authenticated Z-independent geometry, on their ordinary first Z derivatives. The zero-length collar integrals vanish there; this does not assert that their local source functions vanish. Radial derivatives are physical `d_x`: at the seam `d_rho=d_x/4`.

For `q` in `[0,1/4]`, `[1/4,1/2]`, and `[1/2,1]`, the new leaf separately encloses this homogeneous transfer and its first `-2 chi K0 H0` potential word. `chi` is the source-bound axis potential, not the phase shear cutoff. The ordinary Z1 word keeps both `chi_Z * inlet_Z0` and `chi * inlet_Z1`. The four first-word branches add formal hb powers `2, 3, 1, 2`, respectively, for `F_from_Fc`, `F_from_Fxc`, `Fx_from_Fc`, `Fx_from_Fxc`. The `[0,1/2]` hb interval is used only to bound coefficients for the original strict range; hb powers remain formal and no midpoint is selected. These are conditional interval enclosures of the stated sectors, not evaluated source functions.

## Saved execution and evidence

The source manifest and code were frozen before the one permitted V4 entry. The independent synthetic fixture passed **47/47** checks. The saved actual-source result then checked **1,632** identity source routes and **984** zero cross routes at `q=0`, with **1,634** indexed IDs including two chi jets. Three q cells were emitted. Each of the six fields retains Z0 and Z1 and source-parent indexing. The full result was saved before inspection.

| Artifact | SHA-256 |
| --- | --- |
| `experiments/root_st073/n1_e00_collar_fixed_seam_v4_receipt.json.gz` (raw gzip) | `6f8cc4c4fb4baedb9d14274c2f8ddc92f4ae536766d3d911a807d6f626539766` |
| Frozen source manifest | `a4295473c53a43ab183bf10a23ce3be9581477306377ddeb206968342d74853f` |
| Entry record | `9a6d00e6a2817a9c7db1e79f79dc3a0b75aa63b9137857f463da4ada15cde460` |
| Full result | `c15dff1be792915ced695de04f06f445a5161d8082a9ee417cda737e362bdf2e` |
| V4 leaf | `19379232a77926759ee3f955502aa0dad392494ef9d7a280203d64f124039ca9` |
| V4 fixture | `45656be7bde198a8f5ce3553bf94fed4d03ec6ce995070519d02176b25bd7555` |
| V4 runner | `1a871a7e5bfa14e033c9858c13b6bf1171aa3f809dd02faae8ac786db07efa8d` |

Earlier one-entry versions V1–V3 failed inside synthetic fixtures before actual receipt decoding and are not rerun. V1 checked an unprefixed error label against a prefixed label; V2 mixed integer and public Z keys; V3 mixed integer and branch Z keys. Their exact manifests, entry records, full failure results, and frozen code are saved in `n1_e00_collar_fixed_seam_v{1,2,3}_failure.json.gz`, with raw hashes `44de3a7ce09e1aadd4f9b01c8cdbc589a7fd7ff351f04e0b87801c97cd6e19b0`, `57cd598d34b227d368c8981d9791766ed0ffa19806b168f04cbef2bd549ec5d4`, and `e1880de571e505a9a4f92be016f5b84c605f1d0f499c6065d53730838a4f8af2` respectively. V4 also supplies explicit top-level source IDs required by its adapter.

## Limits and next task

The seam identity acts on the supplied fixed E00 packet; it does **not** prove that the old W1/low/active/compact/quiet records represent the same endpoint functions, nor that old and new records are disjoint. No error interval was deleted by overlap or midpoint. The positive-q local forcing and `j>=2` potential remainder are not computed here; the first word alone is not a complete collar solution. The reserved own_swirl contribution, independent negative full-pulse m15 comparison, global five-moment closure, later physical reconstruction layers and independent Cartesian NS validation remain open.

Next, use the frozen V4 source paths to construct a separately versioned same-source endpoint identity or a justified `j>=2` remainder for this fixed inlet, without re-entering any saved ancestor. Preserve the existing six signed core/lower sectors, all accumulated errors, ordinary Z1 and physical scales; do not append sectors to the 212-term ledger until exact ownership and no-duplication are proved. Then continue the independently defined reserved/pulse comparison.
