# Q/P inlet first-potential physical m15 cells — 2026-10-31

This one-entry calculation supplies the first collar potential word for the same 17 Q/P fixed-core-exit inlet groups used by [the K0 physical-cell result](CURRENT_N1_COLLAR_QP_K0_PHYSICAL_M15_CELLS_2026_10_30.md). It applies the variable-target potential kernel to the original `x=2` inlet sectors; the earlier `q=1` result is referenced, not used as a phase profile. It does not recompute K0, E00, the preceding endpoint word, the local forcing or pressure. The 212-term ledger is unchanged.

With `r=hb sc q/4`, `X=2 exp(r/2)`, `K0(X,t)=t(1-t²/X²)/2`, the two positive source-kernel masses are

```text
Ic(q)  = integral_2^X K0(X,t) dt = cosh(r)-1
Ixc(q) = integral_2^X K0(X,t)(1-4/t²) dt = sinh(r)-r
F_j1(q) = -2 chi [Fc Ic(q) + Fxc Ixc(q)].
```

The physical swirl contribution is `-2 Lambda^-2 integral eta rho Flead F_j1 d-rho`, `rho=4 exp(r)`, `d-rho=hb(sc/4)rho dq`. The two branch powers gain formal `hb^3` and `hb^4`, respectively, including the physical radial measure. Whole-cell directed bounds on positive kernel masses and radial factors use the exact velocity cutoff masses `(1/4,1/4,1/8,0)`; the final zero is for cutoff velocity only and does not suppress pressure. The leading F0 connection is reused once. Ordinary Z1 retains the three products `Flead_Z*chi*inlet`, `Flead*chi_Z*inlet`, and `Flead*chi*inlet_Z`. The axis chi is the separately authenticated potential, not the phase shear cutoff.

The result has 684 grouped interval sectors, 228 Z0 and 456 Z1, separated by source group, phase cell, and Fc/Fxc branch. All 3,984 inlet source paths and their parent records match the saved Q/P receipt; the result index contains 3,992 source paths and 3,736 are used. All inherited error-only outputs passed an independent saved-result audit for symmetric endpoints and error labels. The synthetic fixture passed 21 checks, including high-precision cellwise quadrature comparisons; those comparisons are numerical cross-checks, not certified quadrature bounds. Source and result boxes are enclosures, not actual function values.

The Q/P first word is a finite partial contribution. It does not include `j>=2`, complete feedback, local source history, or a full m15 function. `P1` is not the physical pressure density. No cancellation with errors or old terms is claimed, and disjointness from the old W0/W1 prefix has not been established by a source identity.

## Evidence

- Compressed receipt `experiments/root_st073/n1_collar_qp_j1_m15_cells_v1_receipt.json.gz`, raw SHA-256 `d03bb6c2ff6b0d26a74637ba9fb3428f55a1192498972ad3b65c823b6fd18144`.
- Source `fb228fc5edcebb9519731a1fce6828e2af4e1cb2d8037fec73e5b4b9e5f6e2e1`; entry `e24b375585ca9bc6691334a59fa72e563dfc9eee08e4f606a6e8ccfa1ae74a1c`; complete result `cbae3a6144180aca7349b6523c5bd005c9f643885430f5afe8057cbf1da90887`.
- Frozen leaf `aa70ed02c12fa381ed1e151b979cec77a3b24a858ac1eed784b8fe29a8523a9b`; fixture `626fd5bde096357b3cc74d4515f0a07e5a91c818b0788e7f8678022bc678e18a`; runner `e8393c81ba462b1a69e8448f0246fd3f046c4643fcb6819d40238d1e844a1e71`.

The receipt retains full source, entry, result and exact code. Authenticate raw gzip bytes before decoding; never reenter this or earlier mathematical families. Next compute a separate cellwise symmetric j>=2 error for the same Q/P inlet, then establish old-prefix ownership, the complete feedback source and its higher-Z control, and the independent reserved/pulse endpoint identity before any ledger join.
