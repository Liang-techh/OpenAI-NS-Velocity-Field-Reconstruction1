# Current Reserved Partial Whole-Z Handoff

Date: 2026-11-02

## v3 scope and construction

The new v3 extends the previously saved three-component partial response at fixed compact `x = 1.27` from `Z = 0` to four source `Z` cells:

`[-1, -.5]`, `[-.5, 0]`, `[0, .5]`, and `[.5, 1]`.

Together these cells cover `(-1, 1)` by an exact endpoint sweep. Source-cell data were saved as normalized `theta0`/`theta1` `C0` data and ordinary `Z1` data, with `mu`, the saved whole-Z `J = -nonpulse`, and an independent signed pulse `P` from the same old source family under positive physical normalization. `T = -P`.

Do not call `P` positive by sign: `P` is the original signed complete pulse. `P1` is not `m15`.

The applied factors and decomposition are:

```text
A(Z)    = exp(-17/10 - 5*mu/2) / (1 + Z^2)
A_Z(Z) = -2*Z*exp(-17/10 - 5*mu/2) / (1 + Z^2)^2
own     = -r0*A*theta0 - r1*A*theta1
D       = partial - T = own + r2*J + (1-r2)*P
```

The ordinary `Z1` product rule is used, with the premise `r2_Z = 0`.

Each cell receives separate signed interval products and `Z0`/`Z1` handling. Raw uncertainties remain within the inputs; this is not a separately certified error decomposition. There is no cell hulling, ledger append, midpoint evaluation, or exact function-value evaluation.

## Status and limits

This is a conditional four-cell enclosure at fixed compact `x = 1.27`, not physical `Rp`. There is no same live owner, no proof that the response or `mu` is independent of `Z`, and no terminal match or five-moment closure. Prior `212`/six/errors remain inherited.

The v3 checks included 19 fixture checks, including signed `P`, the derivative, and an `r2 = 1` synthetic case. Two coverage test cases passed; the actual new leaf was entered once and exited cleanly. The v3 work used an interior independent probe and fresh synthetic inputs.

## Provenance

- v3 source SHA: `2a321ccd405b062b63998e02294bcb8f48db07af038b55318c565196f6e1737e`
- v3 entry SHA: `5ae312985c9b4c51c6797481109f9df7f21151861c209d8fd7749d56822dffd2`
- v3 result SHA: `5bd627d0744d36804547c6ebbf02884a68b4163ee4eb795ee29b35ecd561a2cb`
- [v3 receipt](../experiments/root_st073/reserved_partial_wholez_v3_receipt.json.gz)
- v3 receipt SHA: `967fd2f9fe3fb8361b2bf632bac2987213c350056b1e3aa7ad5c83056cbf1308`

The v1 synthetic failed its `A`/`A_Z` check from an independent 90-digit probe at the exact boundary versus a 130-digit interval. The exact failure cause is unconfirmed; boundary precision or rounding is only a possible explanation. No actual receipt was decoded. Failure artifact: `experiments/root_st073/reserved_partial_wholez_v1_failure.json.gz` (SHA `a1dc1f5ea4c1e613ac86a2061aea9ecbe1ca8ade8cf9c1996247d8ac179ef2ca`).

The v2 predecessor failed its pre-entry guard mismatch; no actual code was entered. Startup archive: `experiments/root_st073/reserved_partial_wholez_v2_startup_failure.json.gz` (SHA `47be62777ad2b4d11317891cd1a3ecccb277622e0015ad372c74c297018b92be`).

No old producers or helpers were invoked. Raw hashes were checked before decode.
The v1 and v3 archives were checked against the saved source, entry, result and executable bytes. Staged inputs match all three manifests. Frozen files retain trailing blank lines; the staged whitespace check disables only blank-at-EOF checking.

## Next required work

Obtain a source-bound partial response at the actual `Rp` over the same `Z` domain and establish the graph owner and physical normalization. Then compare forward `212` with an independently defined pulse while retaining all inherited errors. Do not infer cancellation from interval overlap.
