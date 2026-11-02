# Safe first-collar inlet tail — 2026-11-02

The long-term reconstruction remains **active and incomplete**. The previously uncertified `j>=3` tail in `formal_kchi_collar_v1` has a new, separately computed conditional upper bound for the same saved baseline inlet. This does not establish complete collar transport, functional five-moment closure, physical m15, or a Navier–Stokes solution. The existing 212-term ledger, six signed inner/lower terms and all prior errors remain unchanged. P1 is not m15.

## Safe upper bound

On the saved first-collar geometry `x(q)=2 exp(hb sc q/8)`, assume `0<hb<1/2`, Z-independent `hb,sc`, `S=(sc/4)exp(sc/16)` and `L<=hb*S`. The causal Green kernels satisfy `0<=K0(x,t)<=x-t` and `0<=K0_x(x,t)<=1`. Put `M=sup|chi|`, `Mz=sup|chi_Z|`, and `A=2 M S²`. For the omitted `j>=3` inlet words, the value radii before the recorded formal `hb` factor are

```text
F:  (2M)^3 S^6 / 6! * exp(A/4)
Fx: (2M)^3 S^5 / 5! * exp(A/4)
```

The separate ordinary-Z1 insertion radii replace `(2M)^3` by `3(2Mz)(2M)^2` and multiply by `(1+A/12)`. The proof compares `(r+2k)!/r! >= 4^k k!` for `r=6,5` and sums the insertion weight `(3+k)/3`; no small-A assumption is used. A separate geometric factor `[0,sc/4]` and one additional formal `hb` power apply to the saved Fx inlet. Consequently the respective formal power increments for `F<-F`, `F<-Fx`, `Fx<-F`, `Fx<-Fx` are `6,7,5,6`. No `hb` or other formal source atom is materialized.

Each source sector and Z derivative route is retained separately. Outputs are symmetric error-only enclosures, never actual tail values. The earlier invalid `exp(sqrt(A)/2)` proof is not reinstated; the earlier `j>=2` envelope and the new `n=2 + j>=3` decomposition are alternatives, not additive independent terms. No earlier error is subtracted or erased.

## Execution and limitations

The independently frozen `inlet_kchi_tail_safe_v1` ran once. Its 39 synthetic checks passed, including fresh signed and error-only inputs, zero and nonzero axis potential, both ordinary-Z1 routes and large-A rational prefix checks. The saved source input contains 32 F/Fx, Z0/Z1 parent records. For each output field F and Fx there are 16 Z0 and 32 Z1 records: 96 in total, all symmetric/error-only. The fresh rational checks exercise finite prefixes; the all-order claim uses the factorial comparison above.

The complete [receipt](../experiments/root_st073/inlet_kchi_tail_safe_v1_receipt.json.gz), with manifest, entry, result and exact frozen Python bytes, has raw SHA-256 `3720e242c83372be22f5924660665d8b6e6b2a2d617d54c2269896113e9b715d`. Manifest SHA-256 is `92dc34104931deec411ca508e8b037d5971038a0fabc73f25679c2c2b3bdde51`; entry `0702df794fe4a8c23c40483970826d100930e466ce72db13b85faf9395996609`; result `ae7fb40d3b0fbea7c4f68df46f23e59a9ce6946957ead3ae9ee22b3035d6dcf7`. Raw receipt hashes were authenticated before gzip decode, and no ancestor function was invoked.

This is a conditional, source-linked bound. The frozen runner checks the actual saved source family, axis binding and geometry certificate; it does not prove a common live owner. Direct callers of the leaf must independently supply the same source, axis and geometry premises: its API does not by itself authenticate the physical owner, a whole-Z source ID, or the support-measure certificate. Do not treat the leaf's metadata alone as such a certificate.

## Source equation and next task

The raw six-row system has `W_x + D W/x = B0 W + B1 D_Z W + g`, with `W1_x=W5`. Hence the row-5 principal left side is `L0 W1=W1_xx+3 W1_x/x`. The raw source-only Green response directly integrates `g5` with `K0` and `K0_x`, but it does not apply the B0/B1 feedback. The full Kchi equation `L0 W1+2 chi W1=h5` consequently needs a same-owner decomposition of the complete row-5 right side. Algebraically, before any cancellation is proved, its right side contains `g5`, the row-5 B0/B1 couplings, and the added `2 chi W1` term. Neither `h5=g5` nor cancellation of that term follows from matching receipt labels.

The next executable work is to build and authenticate an explicit same-owner row-5 coupling/potential adapter, retaining ordinary Z1, signed source identities, `rho=x²`, physical `dx`, and formal scales. It must distinguish the standalone g5 Green contribution from the full `h5` forcing and preserve all existing errors until a same-source identity justifies any cancellation. The exact heat exterior, admissible stress and flat remainder, time-order recursion, corrections and independent Cartesian validation remain open.

Workers were requested with the fixed Luna/max role. The spawning tool did not expose runtime model metadata, so independent verification of that metadata was unavailable.
