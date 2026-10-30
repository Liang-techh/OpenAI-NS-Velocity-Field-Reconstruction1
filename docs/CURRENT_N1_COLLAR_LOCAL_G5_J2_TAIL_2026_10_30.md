# Local collar g5 potential tail — 2026-10-30

The new one-entry consumer bounds the omitted `j>=2` potential words of the local collar row-5 source, separately from the previously saved core-inlet and Q/P groups. It authenticates the same saved raw `g5`, axis `chi`/ordinary `chi_Z`, geometry and first-local-word receipts. The previous local `K0 g5` and first potential word are referenced, not recomputed; the 212-term m15 ledger remains untouched. This is an enclosure of a local contribution, not a signed source value, complete collar solution, full feedback tail, or five-moment closure.

Let `A=-2 chi K0`, with `Kchi = sum(j>=0) A^j K0` and `D=A K0`. The bounded term is `(Kchi-K0-D)g5 = sum(j>=2) A^j K0 g5`. The initial physical `K0` source integration occurs once; `j` counts potential factors only. This local source group is disjoint from the prior `G0g`, core `Dg`, direct core tail, and Q/P inlet-tail groups. The saved first local word is independently referenced and not included in this remainder.

For `c=sc/8`, `E=exp(sc/16)`, `m0=2 c^2 E^2`, let `C0,C1` bound `|chi|,|chi_Z|` and set `k0=8 c^2 E^2 C0`, `kz=8 c^2 E^2 C1`, `kx=4 c E C0`, `kxz=4 c E C1`, `T=exp(k0/4)`. For whole-collar source radii `S0,S1`, the directed bounds are

```text
R_F   = S0 hb^6 m0 k0^2 T/2
R_FZ  = hb^6 m0 T [S1 k0^2/2 + S0 kz k0]
R_Fx  = S0 hb^5 m0 kx k0 T/2
R_FxZ = hb^5 m0 T [S1 kx k0/2 + S0 (kxz k0+kx kz)/2].
```

`Fx` is the physical radial derivative. These bounds follow from the causal physical `dx` kernel, fixed Z-independent geometry, an ordered-simplex estimate and the certified `0<hb<1/2`; `hb` and the other source powers remain formal. A separate static derivation checked the extra source `K0`, the outer derivative kernel, and both ordinary Z1 insertion locations. The calculation does not convert an error radius or an interval midpoint into a function value.

The returned result has 16 `Z0` and 32 `Z1` sectors for each of `F1` and `d_x_F1`: 96 symmetric error-only sectors. All 32 `chi_Z*g5_Z0` sectors retain both chi Z0 and Z1 source IDs. Each source sector keeps its signed input coefficient and origin as provenance, while its output is a symmetric tail enclosure. The independent synthetic fixture passed 40 checks of source powers, error separation, geometry guards and both derivative routes; it did not perform numerical quadrature. These checks do not make the source enclosures exact functions.

## Reproducible evidence

- Complete compressed receipt: `experiments/root_st073/n1_collar_local_g5_j2_tail_v1_receipt.json.gz`, raw SHA-256 `67894bac58f4fbbcd5437bf970dadafeb918fe2d9638425f8e5165d21af4530b`.
- Source checkpoint: `85604a94d4df5cb5fe107028885634140a91a1bba37451a1bc536ffc69f6f9b1`; entered record: `78c0e2d837776a1e6dd9bdf052b11a56c75e08a0214d90faaeba986da9011a5d`; complete result: `fcc6044d8ae06b8f04233ee90c79f34ea6d9aa42fffadabf1001ce936a1e158c`.
- Leaf SHA-256: `c409a0ed354be4b5149dc3986949a5531db0403f89d83cb61d49a68ece3bde10`; fixture: `740eef403763476d041f1f350ad9f943e4f6a8f870983c028dc3bd4898727265`; runner: `286b22947fa3850ea329b76380eb819d9522fd18571c9c59c68c7636bee20bfd`.

The receipt preserves the complete result, source and entered records, and exact code bytes. Authenticate raw gzip bytes before decoding; do not reenter the consumer or prior mathematical families.

## Remaining dependencies

Endpoint-only Q/P data do not supply the physical phase-cell profiles needed to integrate the reserved m15 density. The E00 exact fixed endpoint, local row-4/6 and pressure source histories, common physical phase-cell data, and complete feedback/error ownership still need source-bound treatment. In particular, full feedback ordinary Z1 contains `B1 W_ZZ`; iterated words need higher Z control or a certified common analytic tube. This local potential bound is not that feedback bound. Next build the missing same-source physical phase-cell interface and no-duplication owner binding, keeping `P1` distinct from the `g4/N1p` pressure contribution; then address the independent own-swirl/pulse comparison without dropping inherited errors.
