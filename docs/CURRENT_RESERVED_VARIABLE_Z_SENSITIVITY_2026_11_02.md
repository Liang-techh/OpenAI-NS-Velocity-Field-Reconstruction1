# Reserved partial response: variable-Z sensitivity — 2026-11-02

The long-term goal remains **active and incomplete**. This checkpoint extends the saved, conditional compact `x=1.27` comparison; it does not evaluate the physical `Rp` response, establish a common live source owner, close the five moments, or validate Navier–Stokes. The original 212 cumulative terms, six signed inner/lower contributions, reserved-once accounting and all historical errors remain untouched. `P1` is not `m15`.

## General ordinary-Z derivative

For the previously defined signed `J=-nonpulse`, signed complete pulse `P`, `own=-A(r0*theta0+r1*theta1)`, and `D=own+r2*J+(1-r2)*P`, retain

```text
A = exp(-17/10-5*mu/2)/(1+Z^2)
A_Z = -2*Z*A/(1+Z^2) - (5/2)*mu_Z*A
D_Z = D_Z_baseline
    + (-A*theta0)*r0_Z + (-A*theta1)*r1_Z
    + (J-P)*r2_Z - (5/2)*own*mu_Z.
```

`D_Z_baseline` includes the geometric derivative of `A`, ordinary `theta0_Z`, `theta1_Z`, `J_Z`, and `P_Z`. The four additional coefficients are kept separately; unknown derivatives are not replaced by zero. The `mu_Z` coefficient is `-(5/2)*own`, not `-(5/2)*A`. For an actual source, the caller must supply compatible derivatives, normalization and shared source semantics; interval inputs and a coefficient formula do not prove those properties.

The new independent leaf is `src/openai_ns_reconstruction/reserved_m15_variable_z.py`. It accepts an injected interval context, interval inputs, and optionally all four derivative intervals. If any of those derivatives is absent, its full `Z1` is `None`; it returns a baseline and sensitivity intervals instead. The input need not be an actual field evaluation. The leaf does not authenticate the owner or provide a separate error certificate.

## Computation and limits

The frozen runner consumed the four saved v3 source cells once each at `x=1.27`, using the saved response, `mu`, theta, signed `J` and signed `P` enclosures. It supplied no values for the four unknown derivatives and returned four sets of directed sensitivity intervals. It did not use old producers or add anything to the historical ledger. The `J-P` sensitivity interval in each cell is contained in `[-0.305,0.305]`, checked by exact rational comparison of the saved binary interval endpoints. This is an enclosure, not a function value or a proof of cancellation; `J` and `P` share an old full-source hash and their live-owner relation is not established. Other coefficients can be much smaller than binary64 can represent; an underflowed floating display must not be read as an exact zero.

The one-entry run returned successfully. Its independent synthetic case supplied all four nonzero derivatives; all 13 recorded checks passed. These checks compare directed synthetic enclosures for intersection and check nonzero coefficients, so they are consistency checks, not a rigorous error certificate or a proof of source identity. The saved v3 four-cell coverage is for the open domain `(-1,1)`; endpoint/seam ownership and closed-domain certification were not established here. `mu_Z` and `r_i_Z` also remain unbound to the actual source.

The metadata-only pre-freeze preparation initially stopped because indented serialization exceeded the brief's 97,000-byte limit. It wrote the source manifest and preimage but did not enter the consumer or alter the brief. The pending brief was then written in compact JSON. This did not run an ancestor or a mathematical consumer.

## Provenance

- Frozen source manifest SHA-256: `ad792f32be6782e18be8cdd96113bec2f90177846cdad40491b59fce9452bb7c`.
- Exclusive entry SHA-256: `862fe69bfa40e4cd35eccb4b88676b273df1567f35a3b576419fda846d8aed1d`.
- Saved result SHA-256: `72134617f9e191abdc51cc9be67047313bcd8ebb94e0c7dcc38ef51f5b3d3638`.
- [Complete receipt](../experiments/root_st073/reserved_m15_variable_z_receipt.json.gz), SHA-256 `e11ce84c7487a9f6dfff2a878949af9b0e7f871bbfa68eab69ed9176d8f5baae`; it embeds the manifest, entry, result, and frozen code bytes. The earlier v3 receipt and separately saved result were raw-hash authenticated before decoding and their contents compared.

Bounded agents used the `luna_max_worker` role with `fork_turns:none`; the role definition specifies GPT-5.6 Luna/max, but spawn responses returned only task names and nicknames, not independent runtime model metadata.

## Next executable work

The saved first-collar `finite_prefix_channels` packet already contains separate axial, swirl and pressure signed interval sectors at ordinary `Z0/Z1`, with a common formal atom ordering. A new unverified local, untracked draft `src/openai_ns_reconstruction/first_collar_m15_sector_group.py` proposes grouping only identical powers and uncertainty/error classes while retaining every source record. Review its actual input schema, freeze a new consumer and fixtures before first entry, and compute the saved finite-prefix sector sum. Do not apply the physical weights a second time. The result would still omit the core m15 inlet and all-word tail and must not be called a complete canonical inner moment or used to delete old errors.

Separately, current physical pulse and frame graph handles are not field values. Binding the frame field-radius and endpoint owners, common units/measure and ordinary-Z transport remains necessary before physical `Rp` or pulse cancellation can be claimed. In the active/compact comparison the `Rc*Delta Q_y(1)` lower-endpoint term remains until a same-source identity proves otherwise.
