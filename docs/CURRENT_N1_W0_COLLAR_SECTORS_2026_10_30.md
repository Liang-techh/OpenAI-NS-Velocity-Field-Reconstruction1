# Same-source W0 first-collar sectors in four phase cells, 2026-10-30

The saved core `W0` transport and the actual collar `g` enclosures now give separately stored local `W0_Z0/Z1` sectors for all six fields on four phase cells spanning the first collar. Their sum is a finite zeroth-word interval enclosure under the parent analytic/source assumptions. It is not an exact point value, all-word n1 field, physical m15 calculation or terminal matching.

## Construction

The source `g_Z0/Z1` is an interval enclosure on the entire original collar, not a set of four independent samples. The new calculation bounds its Green integral for each phase cell `q in [i/4,(i+1)/4]`, `s=(sc/4)q`. The previously saved core `W0` homogeneous transport already covers the full phase and is reused unchanged; a cell's finite W0 is defined as that core sector list plus the corresponding local list once. The local Green terms were intermediate quantities in the earlier Z2 computation but were not exported there. They are now calculated directly as a new cell enclosure; the earlier program and Z2 output were not rerun.

With `u=2 exp(hb*v/2)` and `d=hb*(s-v)`, the six kernels, in `(F1,Uz1,K1,P1,d_x_F1,d_x_Uz1)` order, are `(u/2)*hb*(s-v)*exprel(-d)`, `(u/2)*hb*(s-v)`, the negative of the first, `1`, `exp(-3d/2)` and `exp(-d/2)`. The physical Green measure is `dx=hb*(u/2)dv`. Positive kernel-mass bounds use `0<hb<1/2`, the exact cell phase bounds and outward interval arithmetic; K keeps the negative kernel. Consequently the first three local groups acquire formal `hb^2` and the last three `hb^1`. `hb`, `F0`, `Pstar_squared` and `Lambda` are not materialized.

Each cell's six local fields contain `16,43,43,43,16,43` source sectors for each of Z0 and Z1, respectively. The saved core branch remains separate; all 234 new input IDs map to exact parent receipt entries and retain the original parent source IDs and uncertainty classes. A single whole-phase source enclosure is reused on every cell; the result does not claim new source samples, same-source endpoint cancellation or independent cell errors. The existing first-word axial/swirl m15 contributions, Z2, Cauchy alternative, all-word tail and 212-term ledger remain separate and unchanged.

## Verification and next step

The one-entry consumer authenticated the saved core transport, actual forcing, prior physical contributions, same source family and four formal atom records, and the exact geometry. The synthetic fixture passed 246 checks, including independent high-precision direct integration of six kernels at a representative point within each cell for both ordinary orders, formal powers, signed K, cell geometry, input provenance and scope guards. These numerical cross-checks are not a proof by sampling of a whole-Z identity; the actual whole-cell bounds come from the analytic positive-kernel inequalities and directed input intervals.

The receipt is `experiments/root_st073/n1_w0_collar_sectors_v1_receipt.json.gz`, SHA-256 `4229c29a2a21207f8b8da3ba53d70bbf6f50e09eb67110bd3b7cb30fa6577c06`; the one-entry result SHA-256 is `163799c9fb11c10d95ded8caab5e6118de7b23dbaab01b46af77f327b7568d33`.

Next, integrate the W0 axial and swirl m15 contributions using these cell-specific sectors and the same physical leading sources, with the exact own-cutoff and ordinary Z1. Keep the W1 contributions separate until a nonduplicating finite-prefix assembly is justified. Pressure history/endpoint, all-word remainder and corrected-reserved/current-pulse owner and measure bridges remain open; `P1` is not m15.
