# Original repaired preRp weighted history transport - 2026-10-16

`OriginalPreRpM15Numeric.evaluate(Z,chart,coordinate)` now transports the original repaired weighted moment m15 across the compact repair and the complete quiet segment. Use `chart='compact'`, `coordinate=x=R/Rc` in[1,2], or `chart='quiet'`, `coordinate=q` in[0,1], where s=q*log(Rp/(2Rc)). Inputs are exact decimal strings, with -1<Z<1. The output contains source-bound ordinary C0/C1 enclosures as factored exponential terms. It consumes the original selected compact inlet/source family and never uses the terminal-derived numerical Rp inlet to drive this forward history.

This is an implemented cumulative enclosure evaluator. It does **not** resolve the signed numerical centers of the inherited inner/lower/active components or the compact finite-N pressure correction. Their accepted directed source enclosures remain attached to their actual functions; their midpoints are not selected as values. The full objective remains **ACTIVE / INCOMPLETE**, tasks14/15 remain IN PROGRESS, and an independently resolved signed numerical Rp seam/global five-moment numerical admission remains open.

## Incoming history and actual compact repair

The evaluator consumes m15 inner/lower C0/C1 sectors from the accepted full-five numerical bridge and the matching active m15 rows from the preRp bridge. It adds neither the whole preRp component again nor the pulse forcing. The retained incoming/cumulative function handles are `own_field_bindings.incoming.m15` and `own_field_bindings.cumulative.m15` from the original reserved own-field source.

On the compact band, the exact fixed leading profile Q=-K/R gives

```
Delta m15_leading = -K^2*log(R/Rc)/4,
K=Rc*k,    ordinary Z1=0.
```

This follows the original R-weighted S1 identity. Its linear boundary/integral terms cancel; RQ_y is the same at both finite endpoints, and its quadratic integral supplies the nonzero logarithm. The evaluator retains this separate negative term in unit Rc^2*k^2/4. It vanishes only at the empty prefix x=1.

The finite-N correction uses the accepted `error_m15` and the genuinely differentiated `normalized_first_Z_upper.m15`, each normalized by its own p=2 `scale_ratios.m15`. These positive endpoint/integral majorants apply to every partial prefix because x<=2, the boundary factor is at most5 and integral x*dx is at most2. Both finite RQ_y and RQ_yZ endpoints are already included; no endpoint is added a second time. The exact empty prefix has zero correction. Other prefixes retain a nonzero directed error enclosure and the sharper source log bound, without calling the cap an evaluated signed correction.

## Weighted own-swirl history

The original density is

```
dm15_own = -theta0*Utheta1 dR = -R*theta0*Utheta1 ds,
theta0 = Pstar*A_rc*exp(-(1/2+mu)*s).
```

For the same normalized bump shape_j(s), the new positive scalar primitive is

```
integral_0^s exp((1/2-mu)*u)*shape_j(u)du
 = exp(-mu*center_j)/normal
   * integral_-1^clipped_endpoint beta(t)*exp((1/2-mu)*ell*t)dt.
```

The negative sign, extra R and distinct rate are essential. The former unweighted m13 kernel cannot be reused unchanged. Each term retains the original selected beta function enclosures and actual angular amplitude rows. Ordinary Z1 uses the full product rule a_Z*b+a*b_Z; its unit is independent of Z.

The physical scale is combined analytically before interval boxing:

```
Utheta1_unit + logP + logSigma + logRc = fixed_m15_unit,
fixed_m15_unit = log(Phi*A*Rv^2).
```

The production fixed-unit method and the existing production Utheta1-unit method are both executed symbolically in the focused driver. The equality includes the actual shared angular scale logSigma=-logP/2+1/2 exactly once. Adding another logSigma atom after this cancellation would count it twice. No hyper-exponential is materialized.

## Quiet continuation and relationship to the pulse

At2Rc, compact history remains unchanged. The accepted quiet source `actual_quiet_source.m15` has negative sign, radial scale power2 and exact zero ordinary Z1. The new local increment is

```
Delta m15_quiet = -K^2*s/4,    s=log(R/(2Rc)).
```

The full original quiet length and the unchanged K are retained. Vanishing n1 velocity after2Rc does not erase this fixed forcing or the incoming moment. At quiet q=0, every compact term is identical to its2Rc value and the new quiet increment is exactly zero.

The same checked terminal identity proves

```
inner+lower+active+compact+own_swirl+quiet+pulse=0.
```

The reader checks the actual m15 terminal handle against the reserved own-field handle and consumes the accepted translated-integral, transported-partition and selected-inverse proofs. Thus the forward history and negative complete pulse are the same source function atRp. This proof is separate from numerical cancellation accuracy: enclosure overlap and an interval containing zero are not independent signed seam residual evidence. The existing pulse evaluator remains terminal-derived, and its recorded numerical seam limitation is not silently changed.

## Focused evidence

The final focused run exited0 in78.628s; saved-report inspection also exited0. It pins1839 scientific inputs and77 implementation files. The4,520,830-byte report SHA256 is `34b2aca4182bf2e9cbaaee18912efc02c33843f8eb74a4ac3cdf318bae38106e`.

- 6 independent symbolic checks execute the actual p=2 unit, existing Utheta1 unit, combined own scale and original S1 baseline/density/finite endpoints.
- 18 independent C0/C1 comparisons integrate the signed physical own-swirl density against dR at3 finite positive mu values, with27 positive prefix comparisons and3 nonzero negative complete cases. The reference uses the physical R coordinate and a separate quadrature.
- 39 actual original compact/quiet packets,702 source/scope conditions and104 partial finite-N log-bound conditions.
- 138 exact nonzero K-squared-history conditions,18 empty-prefix zeros,9 compact/quiet seam or signed-source partition conditions, and7 rejected invalid requests.

These are scoped construction checks. They do not establish resolved signed component centers or a small independent numerical Rp seam residual.

Use `python scripts/run_original_preRp_m15.py --inspect-report` to inspect the saved evidence without rerunning the construction. Use `--write-report` when implementation changes require fresh focused evidence. Existing read-only worker `/root/n1_atlas_map`: **GPT-5.6 Luna / max**; no new agent. The worker confirmed rate/sign/endpoints and withdrew a shared-scale objection after the combined-unit identity was checked.

## Required next work

1. Replace the inherited symmetric pressure-component bounds with directed evaluations of their original signed radial-free S1 sources, preserving actual inner datum/lower/active histories and finite endpoints. Evaluate the compact finite-N correction with its true signed fields and ordinary Z derivatives. Do not select interval midpoints or differentiate magnitude caps.
2. Retain the selected inverse's source correlations when combining the reserved coefficient response and complete defect. Report independent signed Rp seam residual/error separately from the existing function identity and broad enclosures.
3. Assemble all five original numerical histories and required mixed joins on the same source family. Signed centers/relative cancellation and global five-moment numerical admission remain open until their evidence exists.
4. Continue original flatten/analytic preheat pressure/exact heat and finite-energy tails, common global N, admissible stress/flat remainder, genuine temporal n>=2 recursion/smooth sum, oscillatory correction and independent corrected Cartesian NS validation.
