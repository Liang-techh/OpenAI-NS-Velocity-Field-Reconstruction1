# Core K0 nested feedback m15 — 2026-10-30

## Result and scope

The frozen one-entry `n1_core_k0_feedback_m15_v2` returned directed interval sectors for the target-radius profile and physical axial/swirl m15 contributions of the **nested K0 finite subterm**. The receipt is [n1_core_k0_feedback_m15_v2_receipt.json.gz](../experiments/root_st073/n1_core_k0_feedback_m15_v2_receipt.json.gz), SHA-256 `ca947df9fdf7765d9211f428375d81956af783d89074f9e782ba56aeff91dad0`; its result SHA-256 is `3f7f3398bcad74e43acfeebe6eef64371e5156e77240c2cbe34742b776450952`. This is not a complete original W1 enclosure or a completed m15.

The saved source defines the regular principal inverse using `Kchi = sum_j (-2 chi)^j T^j K0` for F1 and its radial derivative, while the saved core W0/F1 and first-feedback constructions used only K0. The new leaf deliberately uses those saved finite sectors; it does not add `Kchi-K0`, its ordinary Z derivatives, cross terms from the omitted W0 sectors, or an all-word remainder. Consequently even its Uz1 feedback is only a nested finite subterm: the omitted F1 input can feed the Uz1 row. The prior core W0 m15 handoff now carries this correction; its axial W0 and direct N1p sectors are separately identified there. The correction does not remove an error or prove any cancellation.

## Computation

The leaf reuses the saved W0 ordinary Z0/Z1/Z2 cells and C0/B1 matrix cells. On target cell i it integrates source cells j<i and the triangular part of j=i; j>i is absent. Each Green mass is a directed source-cell integral, with the removable axis limit and no division by zero. In ordinary derivatives the finite feedback source is

`Z0: C0 W0 + B1 W0_Z`

`Z1: C0_Z W0 + C0 W0_Z + B1_Z W0_Z + B1 W0_ZZ`.

The source cells and target cells are retained, along with matrix and W0 indices, formal powers, and signed interval coefficients. Across the eight target cells the F1/Uz1 unmerged term counts are 396/540 for Z0 and 792/1080 for Z1. The physical target-radius integration produces axial Z0/Z1 sector counts 32/32 and swirl counts 16/32, with 461 source references. The F0 connection is reused from the same authenticated leading source. This velocity adapter does not consume g4 and adds no pressure contribution.

A separate x=2 computation uses fixed-target source-cell masses. Relative to the prior saved K0 exit sectors, all term keys match and all 88 F1 and 120 Uz1 Z0 terms, and all 176 F1 and 240 Uz1 Z1 terms, are contained in their corresponding older interval terms. This is an enclosure consistency check, not a point value, a same-source cancellation, or proof of equality of the full functions. The synthetic profile fixture passed 25 checks and the independent variable-radius m15 fixture passed 102 checks. The first frozen attempt failed all 16 synthetic swirl reference checks because its reference omitted the formal `F0_anchor^2` factor; it stopped before actual profile/m15 computation. Its exact source, entry and failed result are preserved in [the v1 failure archive](../experiments/root_st073/n1_core_k0_feedback_m15_v1_failure.json.gz), SHA-256 `974d9a88e5671113e42eb6192807fe9d83043f1279fe5751e2821041b44f34c5`. The corrected v2 was frozen separately.

## Execution incident

A parallel read-only audit decoded and traversed prior v2/v3 and collar receipts before freezing a new consumer and before the raw-hash preflight required by the local continuation rules. The agent later reported 12 temporary Python receipt reads/traversals, including an eight-cell comparison, and one raw file-hash command. It reported no writes, old producer imports, or ancestor calculation. Full command transcripts for all traversals were not retained; the comparison result is not used as a certification. The new v2 consumer independently authenticated its pinned raw inputs before decoding and persisted its own entry before its calculations. The incident does not authorize replaying the uncheckpointed traversals.

## Remaining work

The 212-term ledger, its six earlier signed inner/lower contributions, and old error terms remain unchanged. The finite terms here have not been appended. Before treating the response as the original regular solution, construct and source-bind Kchi-K0 and its required ordinary Z derivatives for both the inner W0 and outer feedback kernel, including cross terms and a separately certified remainder; do not substitute a tail bound or midpoint for a signed value. A whole-Z omitted-word bound, a no-duplication ledger correspondence and same-source terminal identity are also still required. Five-moment closure and the later physical construction layers remain open.
