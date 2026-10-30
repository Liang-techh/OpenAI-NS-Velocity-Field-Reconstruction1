# Q/P inlet K0 physical m15 by collar cell — 2026-10-30

This one-entry calculation carries the 17 saved Q/P core-exit inlet groups into physical first-collar axial and swirl m15 contributions. It uses the original fixed `x=2` inlet interval sectors and the no-reset homogeneous `K0` map at each collar phase, not the earlier `q=1` endpoint output as a profile. The sources are authenticated; all 3,984 saved inlet paths and their parent records match the earlier Q/P endpoint receipt. Shared leading U/F sources are reused as multiplicative inputs. The input interval boxes enclose the indicated source data; neither they nor the output are point values or an exact signed function.

For `q in [0,1]`, `rho=4 exp(hb sc q/4)`, `epsilon=Lambda^-1`, and `d=hb sc q/4`, the inlet map uses

```text
H0U = Uc + d Uxc
H0F = Fc + (1-exp(-d)) Fxc.
```

The physical channels are `+2 epsilon Ulead H0U d-rho` and `-2 epsilon^2 rho Flead H0F d-rho`, with `d-rho=hb(sc/4)rho dq`. Ordinary `Z1` keeps `leading_Z*inlet_Z0 + leading*inlet_Z1`. The leading F0 connection is reused once. `P1` and `K1` are not treated as m15 or pressure; the separately saved direct `g4/N1p` pressure contribution is referenced and is not recomputed or added.

The original velocity cutoff has exact masses `(1/4,1/4,1/8,0)` on the four quarter-cells. The final zero applies only to the cutoff velocity terms, not to pressure. Directed phase-cell factors preserve the physical measure and source powers `(hb,F0_anchor,Pstar_squared,Lambda)`. In the Fxc branch `phi(d)=(1-exp(-d))/d` is bounded separately, using `exp(-d_max)<=phi(d)<=1`; hb remains formal. Finite signed sectors and the ten inherited error-only groups remain separate, and no correlated endpoint cancellation is asserted.

The result contains 732 grouped output sectors: 252 ordinary Z0 and 480 Z1, kept by input group, phase cell, channel, and inlet branch. The source index contains 3,994 entries, of which 3,866 are used in velocity products; the remainder retain unused nonvelocity or structural-source paths. All error-only output sectors passed the saved-result symmetric-endpoint and error-label audit. An independent synthetic fixture passed 21 checks, including a nonzero-phase derivative-branch quadrature cross-check. Those numerical checks do not certify the actual source boxes as point evaluations.

The old E00 and first-potential collar word are not recomputed. The new groups are not appended to the old W0/W1 prefix or the 212-term cumulative ledger: disjointness from that prefix has not been proved at the required source-identity level. The calculation supplies neither the Q/P first-potential cell contributions nor complete feedback and all-word errors. It therefore does not establish the full reserved own-swirl m15, its pulse comparison, or global five-moment closure.

## Evidence

- Receipt `experiments/root_st073/n1_collar_qp_k0_m15_cells_v1_receipt.json.gz`: raw SHA-256 `48409ed314e95fe0f3117f6c27ed74b5aa8b56822c810ba5050d2b000b08e5c2`.
- Source checkpoint `fcc3f39757cd7f923105b5388d804b24eb13ee226ce227c2965fa3742e269585`; entered record `9a9f76d8f09d66a3df10fcb7099fb5f74d796c52f87b75722300eb7f007281f4`; complete result `727718a05ea68fb7e35ab978aa0654bb2ee2c201a2a4c4d336547909d4c374b3`.
- Frozen leaf `406bd7f3cf3cb55c1f793eca1db9f5f13c21118f3f3d2278299f80d89734f4cb`; fixture `bcd09b9ac99a4b57f38e6c1baa3ec0778cd3400279cbacaae3a7bf0171ade54d`; runner `4c4358bd3ab7669d5cc4793999950a3420784337fa003bd99ee08bb00113f16b`.

The receipt contains the complete result, checkpoint, entry and exact frozen code bytes. Authenticate raw gzip bytes before decoding; do not reenter this consumer or prior mathematical families.

## Next work

Construct cellwise first-potential Q/P propagation with both `chi_Z*inlet_Z0` and `chi*inlet_Z1`, and its separate j>=2 uncertainty; then bind the local source and existing W0/W1/pressure ownership without assuming disjointness from interval overlap. Establish common physical phase-cell sources and a collar-wide feedback bound requiring higher ordinary Z jets or a certified analytic tube. Only a same-source endpoint/measure identity permits a later reserved own-swirl versus independent pulse comparison and ledger join.
