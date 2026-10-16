# Original repaired Rp m15: reserved matrix response and signed endpoints

Tasks14/15 remain **IN PROGRESS**; the full goal is **ACTIVE / INCOMPLETE**.
`OriginalRepairedRpM15` reads the accepted forward checkpoint without running any ancestor producer. It preserves all212 physical sectors, appends the distinct reserved positive-order own-swirl contribution exactly once, and compares the resulting213-sector state with an independently read negative complete-pulse m15 target. Ordinary Z0/Z1 and physical scales remain explicit. No interval midpoint, finite-N cap, P1 datum, or terminal-derived forward value is selected as a function value.

## Source correspondence and physical units

The immutable input is `lei_ren_part1_paper_compliant_current_original_compact_quiet_m15_numeric.json`, SHA256 `7ccb30dbb3e440a3e598eb30deaf03509cf8624e429478aeeb364671baaac071`. Its zero-based source partition is:

| Old component | Preserved forward indices | Physical interval |
| --- | --- | --- |
| inner+lower, joint replacement | 0..5 | axis to rminus; core plus full cutoff/lower collar |
| active | 6..206 | rminus to Rc, all original intervening charts |
| compact | 207..210 | Rc to2Rc |
| quiet | 211 | 2Rc toRp |

The first six sectors are compared verbatim with the accepted signed inner/lower report. Indices0..1 are the core only, ending at rho=4. Indices2..5 are cumulative collar terms on q=0..2, crossing Rin at q=1. Thus assigning0..1 to old inner and2..5 to old lower would be wrong. The actual coordinate is `R(q)=4epsilon*exp(hb*(sc/4)*q)`: old inner is core plus collar_prefix(1), and old lower is collar_prefix(2)-collar_prefix(1). Saved q=1 packets verify the same core and exact phase; their shared endpoint cancels as a source function. These extra q=1 evidence terms are never appended to the forward state.

The old aggregate inner/lower/active/compact/quiet boxes are never added again. `n1_preRp_moment_functions` and `n1_preRp_numeric_bridge` must have identical actual piece/integral handles and checked original-source replays. The old inner/lower physical graph pairs are53161/53162 and36592/36593, their sum53171/53172; active is30707/30708, compact31155/31156 and quiet31215/0. The complete preRp pair is31221/31220.

This is a checked domain, density, source-handle and implementation correspondence. The new direct numeric sectors have **not** been registered as identical archived graph nodes, and correspondence alone does not resolve their signed numerical difference. Each sector already contains its own radial weight and native measure; no additional R, hb, x, or fast-phase Jacobian is applied.

The common positive, Z-independent pressure unit is `reference*A*Rv^2`, where `reference=10.99*sqrt(2*pi)*(exp(5)-exp(3))`. Both independent archives must have the same original Rc/Rp/Rv and A recipes and the m15 combined inverse-mu numerator37/2. This is the p=2 m15 unit, not the p=1 m13/P1 unit. The unit stays as directed logarithmic atoms; enormous radii or exponentials are not materialized.

## Reserved response before independent beta boxing

`OriginalReservedM15Response` contracts the partial weighted pressure row with the original common3x3 angular inverse before multiplying independent input enclosures. With compact s=log(R/Rc), the physical density is `-R*theta0*Utheta1 ds`; its bump exponent is1/2-mu. The five-bump family uses angular beta1/beta2/beta3 here, separately from the background compact beta0/beta2 control.

For the original positive bump normalizer and centers, the partial row has components

```
w_j(s) = exp(-mu*center_j)/normal
         * integral[-1,min(1,(s-center_j)/ell)]
             beta(t)*exp((1/2-mu)*ell*t) dt.
```

Forward and backward directed integrals enclose the same row. Empty support gives exactly(0,0,0). Full support is the **defining third matrix row**, so `w*T^-1=(0,0,1)` by the checked selected inverse identity, rather than a numerical overlap. The original rational A_rc factor cancels its identical denominator in the third theta input before boxing. Consequently, the full own-swirl response is the negative original complete unrepaired m15, with its ordinary Z1 row. Other inputs retain `A_Z*theta + A*theta_Z`.

At Z=0 the terminal C0 enclosure width decreases from about12973.26958 under independent beta boxing to0.6086607263, and C1 from171684.0680 to3.451744527. These are enclosure-width improvements, not function values or seam accuracy. Partial-support response/input components still use independent enclosures; a joint representation of all their correlations is not installed.

## Six exact internal endpoint cancellations

The raw213 terms remain available unchanged. A separate201-term view removes only these source-proved leading endpoint pairs:

| Original indices | Shared source endpoint |
| --- | --- |
| 4,8 | epsilon^2*A0(Z)*sigma(sc/2), A0=32Z*Psi_rho(4)/(1-delta Z^2) |
| 124,143;125,144;126,145 | restoration upper R^2F versus postrestore lower R^2F, epsilon degrees0..2 |
| 133,151;134,152 | restoration upper RQy versus postrestore lower RQy, epsilon degrees0..1 |

The singular pair requires identical exact sc/phase tuples, the same rho=4 core source and the same full q=2 inlet. The restoration pairs require the unchanged full predecessor state, identical atomic interface radius, opposite endpoint orientations and exact local offsets `-8+1=-7+0`. Their coefficients are replayed as exact directed negatives as an additional guard. The functional identities come from the common source frames and fixed physical endpoints, not from interval equality. Both ordinary Z rows cancel because the endpoints are Z-independent.

**Every finite-N and width-error sector remains.** In particular, terms5/9 and all accumulated endpoint/volume/viscous error caps survive. No assertion about a zero terminal relative correction deletes an earlier integrated error.

## Independent target and the remaining signed difference

The target is read from `n1_pulse_numeric_bridge.samples.*.normalized_true_pulse_C1_enclosures.m15`; its physical source handles are35556/35557 and normalized handles53428/53431. It is not computed from the new Rp prefix and is never appended to that prefix. The checked terminal pairs remain:

```
terminal_sources.m15 = 53832/53833
patch_endpoint.m15  = 53905/53906
reduced_patch.m15   = 53921/53922
reduced_terminal.m15= 53931/53932
```

The original complete defect graph explicitly has

```
full_m15 = (inner + lower + preRp) + complete_pulse,
reserved_full_support = -full_m15,
reserved_full_support - (-complete_pulse) = -(inner+lower+preRp).
```

The graph traversal verifies the actual pulse handle occurs exactly once in each outer sum and follows the exact nonpulse inner/lower/preRp handles. This permits cancellation of the **same pulse source** before interval arithmetic. The independently computed forward212 sectors and their errors still remain in the residual ledger. The code does not replace their difference with the functional terminal theorem or declare zero from broad overlap.

The dominant individual normalized error log cap is still `first_micro_finite_N_viscous_error`, about `6.287445994e+408906090034569676` at these queries. This is a bound on an individual sector in the logarithmic comparison, not a measured residual. It is far too broad to certify a useful numerical seam. The current implementation deliberately leaves this obstruction visible.

## Reproduction and next work

The production run exited0 in6.256 seconds, including a2.600-second independent finite-mu fixture. Saved-report inspection also exited0 without recalculation. Three saved axial packets at Z=-.3,0,.3 retain213 raw and201 reduced sectors each; three actual reserved partial-row queries at x=1,1.27,2 and three invalid/absent axial rejections pass. The fixture independently reconstructs the smooth-bump matrix and directly integrates the physical C0/C1 density; it detects transpose, sign, missing-R and missing-product-rule errors. No accepted ancestor constructor is run.

The12,445,568-byte report `lei_ren_part1_paper_compliant_current_original_repaired_Rp_m15_numeric.json` has SHA256 `bd9fb2fc23072c743eeeadd300041ba187aff908a3bff4a3399f7222b72e6871`, with1857 input pins and131 implementation pins. Scoped Git attributes preserve both new hash-pinned scripts byte for byte on other hosts.

Use `python scripts/run_original_repaired_Rp_m15.py --inspect-report` to inspect the saved report and all source/implementation pins. `--write-report` runs only the new saved-state assembly, bounded reserved-row queries and independent finite synthetic fixtures. It does not replay accepted ancestor calculations. Unsupported axial queries are rejected rather than silently reconstructing earlier stages.

Two bounded workers used the fixed `luna_max_worker` role (GPT-5.6 Luna/max, fork_turns:none): source/endpoint audit and independent smooth-bump/direct-physical integration checks. Spawn results exposed task IDs/nicknames, with no separate returned model field; the tool schema fixes the role model and effort. No GPT-6 Astra child was created.

Next, bind **signed actual finite-N RQy traces** on the first-micro/second-micro interface to common original source IDs, with both C0 and ordinary Z1 and the exact physical radius. Prove the two traces are the same source before replacing their two absolute endpoint budgets by a telescoping pair. Extend that evidence across subsequent charts only where the same-source interface proof exists. Keep finite-N volume errors, numerical source errors and all unrelated endpoint errors. Do not delete the large aggregate caps merely because this telescoping is expected.

Independent signed finite-N centers/useful cancellation and the global all-five/common-N functional seam remain open, as do exact heat/finite energy, admissible stress/flat remainder, genuine n-dependent temporal recursion and smooth summation, realizable mean/oscillatory corrections and independent corrected Cartesian NS validation. Three saved axial examples are evidence for this implementation, not a functional five-moment closure or completed physical field.
