# ST073 signed lower-collar and first-micro endpoint comparison — 2026-11-18

This continues the [whole-axial canonical signed primitive](CURRENT_ST073_SIGNED_COLLAR_WHOLE_Z_2026_11_18.md). It links the singular endpoint and width remainder to the old forward ledger, keeping the old approximations and all old error rows intact.

## Endpoint and finite coefficient

The lower-collar source ends at `S=sc/2`, with `rho=4 exp(hb S)` and `R=epsilon rho`. The shared inlet coefficient is `A0=32 Z Psi_core,rho(4,Z)/L`. Its ordinary first derivative uses the quotient rule. The canonical endpoint has weight `+epsilon² A0 sigma(S)` in zero-based slot 4 and `−epsilon² A0 sigma(S)` in slot 8, where the first micro increment starts. Their common radial endpoint is independent of `Z`. The two recorded finite producers use the same radius-four inlet recipe and the same positive `sigma` source. The native leading collar identity identifies the analytic quotient on this interval; the finite Taylor inlet still carries its own truncation and directed input error.

The frozen comparison uses a common exact scale with `odds=(1−S)^−2−S^−2` and `sigma(S)=exp(odds)/(1+exp(odds))`. It compares coefficients in that scale. The printed logarithm and the old coefficient's printed value are independently serialized approximations; rounding a log at this scale makes exponentiating its printed value unsuitable as a source identity. The complete original rows and their error fields are retained in the receipt. The signed canonical sixth quantity is compared to its serialized zero model and audited against the old unsigned width bound. It is not identified with that bound.

## Checked results

The accepted run followed the full first-micro archive at endpoint `1` into each saved quiet `Rp` packet. It verified the unchanged first six terms and the entire first-micro prefix. In each 212-row ledger, the principal is the same canonical endpoint with two opposite orientations. The two shared principal references cancel at source level, with the orientation and physical unit fixed before any coefficients are compared. The prior raw forward endpoint rows and both width caps at slots 5 and 9 remain in place.

| Z | Canonical principal minus the stored model, value | Canonical principal minus the stored model, ordinary Z1 | Canonical signed sixth quantity, value / ordinary Z1 |
| --- | ---: | ---: | ---: |
| −0.3 | +1.64902736e−19 | −7.93034566e−19 | −21.45063174 / +103.15833936 |
| 0 | exact 0 | −2.52955607e−100 | exact 0 / strictly negative |
| +0.3 | +1.64902736e−19 | +7.93034566e−19 | −21.45063174 / −103.15833936 |

These are rounded directed coefficient diagnostics. The fifth-row unit uses `epsilon² exp(odds) Pstar²` at `±0.3` and `epsilon² exp(odds)` at `0`. The sixth-row unit is `epsilon² hb exp(odds) Pstar²` throughout. All six signed sixth-row comparisons are inside the stored width cap. The archived fifth-row ordinary Z1 intervals at `±0.3` do not contain the full new canonical interval. This diagnostic does not establish whether an exact function value lies outside the old interval. No containment or interval overlap was used to prove a source identity.

[Passing v2 receipt with exact binary endpoints](../experiments/root_st073/original_singular_forward_compare_v2_receipt.json.gz). Raw SHA-256 `66510cffb63e748591ed61c4152f8628c15ee2a14e10cf527dfd5a4130d07144`. The [first frozen attempt](../experiments/root_st073/original_singular_forward_compare_receipt.json.gz) rejected the shorter first-micro archive because it lacked the actual terminal prefix. The corrected entry uses `first_micro_full_m15_numeric`, without replaying its producer.

## Remaining work

Bind the fourth forward regular-pressure term to the actual `N1p` partition with the already shared principal and signed remainder. Keep the raw cutoff response, core histories and finite-`N` tail terms explicit. The whole six-term equality and the independent terminal comparison remain open. The reserved own-swirl correction and independent negative pulse target are each retained exactly once in their existing roles.
