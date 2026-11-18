# ST073 signed collar m15 over the real axial domain — 2026-11-18

This checkpoint extends the [selected-input radial correction](CURRENT_ST073_SIGNED_COLLAR_AVERAGE_2026_11_18.md) to closed axial cells. It uses the same canonical pressure, core, and prescribed-shear collar as the [signed radial FTC](CURRENT_ST073_SIGNED_COLLAR_FTC_2026_11_18.md).

## Pressure slope and analytic axis value

Write `L=1−δZ²`, `d=1−Z²`, `x=1+Z²`, and `γ=1+δ`. The checked flattening integral is a positive measure of `x^(−β)`, with `β=2(1−σ(t/100))` in `[0,2]`. The integrals of `x^(−β)`, `β x^(−β−1)`, and `β(β+1) x^(−β−2)` are bounded using the saved positive mass interval and factors `1,2,6`. Differentiation under the integral is covered by the saved holomorphic source proof. This preserves the exact factor of `Z` in the pressure derivative, including at `Z=0`.

The ordinary C1 core-inlet pair encloses `Psi_core,rho(4,Z)/Pstar²`. Its pressure part is `Z B/L`, where `B=M2(2d/x³+γ/x²)+γ M0+d∫βx^(−β−1)+γ Fflat`. Its second part is the nonpressure axial slope and the nonlinear `X_h` error, both scaled by `Pstar^−2`. The positive reciprocal factor is bracketed by powers of two derived from the saved `logPstar` interval. The nonlinear corrections at `(radial, axial)=(1,0),(1,1)` use the coefficient embedding at radius `4`; the first derivative is ordinary `∂Z`.

The primitive in units `Q=epsilon² hb Pstar² sigma(S)` is represented as `Z*g(Z)`. The separately factored correction in units `Q*S` is also represented with an explicit `Z` factor. Its remainder uses the saved uniform radial bound on the common analytic tube and the prefactor `2/L`, so the implementation never divides by a sampled `Z` or uses axial reflection symmetry of the hydrodynamic core. The shift `j>0` is retained in the ordinary axial derivative at zero.

## Verified output

The [frozen run](../scripts/run_original_collar_whole_z_once.py) returned 64 contiguous closed cells covering `[-1,1]`. All cells carry a directed ordinary `Z0/Z1` pair. On 62 cells the signed primitive is strictly negative, and the radial correction is strictly positive. The two cells touching the origin retain their full intervals, which include the exact axial zero. Derivative signs are unresolved on some cells near turning points. These are enclosures of one canonical source function evaluated on entire cells.

| Z | Signed primitive | Ordinary Z derivative | Separate radial correction | Ordinary derivative of radial correction |
| --- | ---: | ---: | ---: | ---: |
| −0.3 | −21.45063174315333 | +103.1583393581326 | +42.90126348630667 | −206.3166787162653 |
| 0 | exact 0 | strictly negative | exact 0 | strictly positive |
| +0.3 | −21.45063174315333 | −103.1583393581326 | +42.90126348630667 | +206.3166787162653 |

The displayed nonzero entries are rounded endpoints in the respective units `Q` and `Q*S`. In those units, the axis slopes are enclosed by approximately `[-2.140,-1.069] × 10^(-204453045017284587)` and `[2.139,4.280] × 10^(-204453045017284587)` for the complete primitive and radial correction, respectively. They differ from zero despite the symmetric pressure contribution. The earlier probe at `371/1000` also appears in the receipt. No historical evaluator or quadrature was rerun.

[Exact endpoints and frozen input hashes](../experiments/root_st073/original_collar_whole_z_receipt.json.gz), raw SHA-256 `28268c3f2ef89941c02d792faf859cfbe3cbdf3b4816b7a336298fe9b20c884a`.

## Next source comparison

Use this same-source inlet and the whole-cell signed primitive to compare the fifth and sixth forward channels. The fifth channel uses the principal endpoint `epsilon² A0 sigma(S)`. The sixth was stored as a finite-`N` approximation with an unsigned width error, so numerical enclosure agreement cannot establish common source identity. Bind the exact endpoint orientation and source functions first. The 212 forward entries, all inherited errors, the single reserved own-swirl registration and the independent negative pulse target are retained. Five-moment terminal matching remains open.
