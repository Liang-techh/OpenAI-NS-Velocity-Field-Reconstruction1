# Fixed E00 physical velocity m15 cells — 2026-10-31

This handoff projects the saved fixed-E00 homogeneous inlet contribution to physical first-collar velocity m15. It consumes the K0 and J1 field sectors from `CURRENT_N1_E00_FIXED_COLLAR_SEAM_2026_10_31.md` and the independent J2 symmetric field tail from `CURRENT_N1_E00_FIXED_J2_VOLTERRA_TAIL_2026_10_31.md` (published predecessor `6a80ef10c8a453f63b56bfce3174e7b3a386fad1`). This is a conditional partial integral enclosure, not the complete m15, a signed evaluated target, five-moment closure, or an NS validation. No output has been appended to the existing 212-term ledger and no old uncertainty has been removed.

## Physical construction

For `rho(q)=4 exp(hb sc q/4)` and `d rho=hb(sc/4)rho dq`, the source-bound leading fields and fixed E00 corrections use

```
axial:  +2 Lambda^-1 integral eta Ulead U1 d rho
swirl:  -2 Lambda^-2 integral eta rho Flead F1 d rho.
```

The K0 part contributes separate axial U and swirl F channels. The J1 finite potential word and J2 symmetric tail contribute only swirl F; `d_x F` is a derivative field, not another swirl m15 term. The physical radial measure is applied once. Each sector retains its original four formal powers and the measure adds `hb^1 Lambda^-1` in axial or `hb^1 Lambda^-2` in swirl. In particular, the J2 Fc and Fxc swirl sectors acquire `hb^5` and `hb^6` in total relative to their fixed inlet.

The original velocity cutoff `eta=1-sigma(4q-2)` has exact masses `1/4, 1/4, 1/8, 0` on four quarter-cells. Input field and leading enclosures cover `[0,1/4]`, `[1/4,1/2]`, `[1/2,1]`; the last enclosure safely covers the active `[1/2,3/4]` region and can be wider than needed. The cutoff mass is used as a positive integral weight, not as a selected q value. The `[3/4,1]` zero applies to velocity only and makes no pressure claim.

Because geometry and eta are Z-independent, ordinary Z1 is assembled as `leading_Z * field_Z0 + leading * field_Z1`. The J1/J2 field Z1 input already contains its separate `chi_Z * inlet_Z0` and `chi * inlet_Z1` routes; neither chi nor a potential operator is applied again. J2 products remain symmetric `error_only` enclosures. Fixed E00 P1 is a state component and is not physical pressure or m15.

## Saved execution

The new one-entry consumer authenticated raw input hashes, common source family/logs, identical fixed-E00 tail source records, exact q-cell agreement, original cutoff/leading-source lineage, and the F0 connection correction before projecting. It preserves a 1,644-entry source index, with 1,260 IDs used by these products. The actual source boxes and results are directed intervals, not point values.

The new synthetic fixture passed **83/83** checks, including signed channels, distinct U/F ownership, ordinary Z1 routes, formal physical powers, symmetric J2 handling, zero velocity cutoff tail, and independent nonconstant-rho constant-cell integral comparisons. A saved-result audit found all **6,336** J2 sectors symmetric and error-only, including **1,056** terms carrying the chi_Z source ID. The actual-source result contains **23,712** ungrouped sectors over the three active cells:

| Component | Z0 | Z1 | Total |
| --- | ---: | ---: | ---: |
| K0 axial | 1,440 | 4,320 | 5,760 |
| K0 swirl | 1,056 | 4,224 | 5,280 |
| J1 swirl | 1,056 | 5,280 | 6,336 |
| J2 swirl (symmetric error) | 1,056 | 5,280 | 6,336 |

| Artifact | SHA-256 |
| --- | --- |
| `experiments/root_st073/n1_e00_physical_m15_cells_v2_receipt.json.gz` (raw gzip) | `209218dbe33f3f17c4e78d8e31c54d6398e0e13841a5993c48b9c099d64bd2c1` |
| Frozen source manifest | `11fa2cab21e8959d3ca405d5802732115129a1b5a08f87752340a6cd3e4499d6` |
| Entry record | `9401b6e13ce9072aa195f5c71169cbcfd40646226b6e3415dbfb7575dc3ea0fe` |
| Complete result | `34c0456c169e230672fe5bf2cc2d36f3fd553d46557210af404459e7b86c218c` |
| V2 leaf | `11c195d018d0e883149b8642fa4d8ce303acb7cdf4bc27069c45f91db6d9ff17` |
| V2 fixture | `3b14151d71049bbad1c7cf89dde9fb122f9ba79ee5cd508fa81b3be68c9b2fb0` |
| V2 runner | `47c0227d0651e6cdc632d28933ba0c5edd9d67329310a3d663a365d8e8a9fe19` |

V1 passed its fixture and reached actual-source assembly but failed while saving because its adapter expected an `axial` key in the J1/J2 channel map, where those channels are structurally absent. It was not reentered. Its frozen code, input manifest, entry and full failure result are preserved in `experiments/root_st073/n1_e00_physical_m15_cells_v1_failure.json.gz` (raw SHA-256 `48f0382f6a259b5756d325546a2c2674f82a3fc09934e2aa453a258ea9816659`); V2 uses fresh synthetic inputs and explicitly treats these absent channels as empty.

## Limits and next executable work

This does not assemble positive-q local `g/C0/B1` forcing, its omitted words, pressure, higher-feedback errors or a full Kchi solution. It also does not prove that fixed E00 and old W1/low/active/compact/quiet sectors are disjoint or identical. Therefore its 23,712 sectors are not added to the earlier ledger and interval overlap is not used to cancel errors. The separately defined reserved n1 own_swirl and negative full-pulse m15 comparison, function-level terminal closure and subsequent physical validation remain open.

Next construct a separately frozen source-bound positive-q local-forcing interface, preserving the local G0g/C0-B1 word, its ordinary Z1, physical measure and missing higher-feedback tails separately from the fixed E00 inlet. In parallel, pursue a term-level old/new endpoint identity and the full-Z reserved/pulse source evaluator. Only a proved owner and no-duplication relation permits any ledger join or error cancellation.
