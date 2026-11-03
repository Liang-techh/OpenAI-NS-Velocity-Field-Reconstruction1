# Reserved Rp source geometry and general C1 matrix action — 2026-11-03

The long-term Navier–Stokes objective remains **ACTIVE / INCOMPLETE**. This checkpoint establishes the reserved/quiet radius relationship at the authenticated source-recipe level and adds an ordinary-Z derivative implementation for a varying three-channel matrix. It does not establish a common live owner, the complete m15 C1 row identity, or five-moment closure.

## Actual reserved and quiet source geometry

The new source consumer read three previously saved graph archives only after their raw SHA-256 values were checked. The reserved and quiet graphs have the same first 24,952 serialized rows, including the shared `original_Rc_offset` node 4893. **That prefix comparison was only between reserved and quiet**; the native-frame graph was checked through its source family and specific expression rows, not compared row by row with that prefix.

Write `O` for the shared offset and `Lm = log(4)-4 logP-1000+hbB*sc/2`. The checked graph recipes give

```
Rc_quiet = exp(O)
Rc_physical = exp(O + Lm)
Rp_quiet = 2 exp(O) exp(Tw - 2 - log(2))
log(Rp_physical) = log(Rc_physical) + Tw - 2
Rp_physical / Rc_physical = exp(Tw - 2).
```

Thus the quiet radius and the physical reserved radius differ by the common factor `exp(Lm)`; the factor cancels in the ratio. These are graph/source expression identities, not equalities inferred from interval overlap. The previously checked quiet archive contains a strictly positive lower bound for `Tw-2-log(2)` on each of four axial cells. This is reuse of the pinned directed-range evidence, not selection of an endpoint as a function value.

In units of `log(2)`, the five source-support intervals end at `121/400`, `131/400`, `141/400`, `151/400`, and `161/400`. All five lie before `2Rc`, with exact margin `239/400 log(2)`; the three angular bumps have margin `249/400 log(2)`. Under the authenticated source recipes and saved quiet-length evidence, their complete support lies before `Rp`. Two modified-graph negative cases, changing the quiet Rc or Tw edge, were rejected. This is source-recipe geometry: the final graph roots have different suffixes and the physical live owners are not identified by this result.

The one-entry geometry run passed. Source SHA-256 `1665259fad1025f8e01aad4e6e518ac8f72cb5bf31b20b67317e831823f44eaf`; entry `6ab80a7da557b6a8b8fd5c0f9a8d1b80b46865b18802d8e280dac6eb57853f9c`; result `7fcaa2ed6a0ccbf6c2c8c4b4a2d2687ed19ac195ce4034ca6afd932ac4fcebe6`. Full input bytes and output are preserved in [the geometry receipt](../experiments/root_st073/reserved_rp_source_geometry_v1_receipt.json.gz), SHA-256 `e57cabcf65fd56baafd24b7e0d6202a5a5b832289c176720747c693bcbe152fd`.

## General ordinary-Z matrix derivative

The existing reserved coefficient action applies a fixed inverse to the derivative of its right-hand side. That is valid only when the actual matrix is Z-independent. The new `reserved_matrix_c1_action_v1.py` also handles a general differentiable matrix using

```
beta   = M_inverse y
beta_Z = M_inverse (y_Z - M_Z beta).
```

For a complete third row it constructs both `M[2,:] beta` and `M_Z[2,:] beta + M[2,:] beta_Z` from the same graph handles. The physical m15 minus sign remains a caller responsibility. The inverse, the derivative pairing and the actual source owner remain caller premises; this module does not claim to have proved `M M_inverse = I` or bound the new action to the real m15 owner. `A_Z` is retained through `y_Z` and `mu_Z` through `M_Z`.

A new exact-rational synthetic test with nonzero `M_Z` passed ten checks, including foreign-owner and incorrect-row rejection. It obtained `beta=(2,-1,3)`, `beta_Z=(9/2,-7/2,3/2)`, and third-row value/derivative `(3,2)`. Omitting `M_Z` instead gives `(7/2,-3,2)` for `beta_Z`, confirming that the extra term matters. These are synthetic values, not physical m15 results.

The first frozen matrix runner failed before any mathematical check because its component name did not match the module schema. Its failure and inputs remain in [the v1 receipt](../experiments/root_st073/reserved_matrix_c1_source_v1_receipt.json.gz), SHA-256 `3b47c799b96ec624d065ff58c7573e2c75888f5f4aaf6ee17605d6446e5331c0`. The corrected distinct v2 run passed and is preserved in [the v2 receipt](../experiments/root_st073/reserved_matrix_c1_source_v2_receipt.json.gz), SHA-256 `e0dd5063f28f1df94b8f6aab2c4ce392eff9e67be85305ac5ece854ea17e067e`. No old producer, integral, or fixture was rerun.

## Next action and preserved state

Bind the checked source-recipe geometry to the actual reserved and quiet live source owners; independently establish that the complete physical m15 row and the third angular row are the same full-Z C1 functions, with their common inverse and sign/scale. Then connect the native-versus-corrected preRp ordinary Z1 endpoint and measure differences to the same source, retain the independently defined negative full pulse, and add reserved own-swirl exactly once.

The 212 forward terms, six signed inner/lower terms, all inherited errors and independent target remain intact. P1 remains distinct from m15. This checkpoint did not compute physical m15 values, contract old errors, or close any terminal moment, exact heat exterior, stress/remainder, higher-order recursion, oscillatory correction or Cartesian NS validation.
