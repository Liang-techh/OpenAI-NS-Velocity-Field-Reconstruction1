# ST073 signed lower-collar m15 radial correction — 2026-11-18

This checkpoint continues the [signed lower-collar FTC source](CURRENT_ST073_SIGNED_COLLAR_FTC_2026_11_18.md) at the saved input `Z=371/1000`, `S=sc/2`. The old packet enclosed the whole signed primitive around its inlet model. The new consumer encloses the signed radial correction separately, in a unit that retains its additional factor of `S`.

Let `P` be the signed primitive `epsilon² hb ∫₀^S sigma′(s) [(A−A0)/hb](s,Z) ds` evaluated at `S`. Keep the same `A0` as in the inlet model. The radial FTC gives

```
P/(epsilon² hb sigma(S)) + A0
   = (1−hb)(2Z/L) ∫₀^S [1−sigma(t)/sigma(S)] K(t,Z) dt .
```

The integration variable on the right runs from `0` to `S`. The identity follows by interchanging the two integrals at fixed `Z` and using `sigma(0)=0`. For the admitted flat logistic cutoff, `(log sigma)′ ≥ 1/S³` on `0<t≤S`. If `M` is the mass of `1−sigma(t)/sigma(S)`, this gives `0≤1−M/S≤S²`. The exact positive atom `sigma(S)` stays in the physical scale.

Write `K=2T0+(K−2T0)` with `T0=16 Psi_core,rho(4,Z)`. The code uses the saved actual core and ratio bounds on their common complex axial tube to enclose the second term. The inlet coefficient comes from the signed directed C1 model already saved for the same source. The ordinary first `Z` derivative uses the saved Cauchy gap on the error only. The code keeps separate enclosures for the inlet model, the radial correction divided by `S`, and the resulting whole primitive.

## Checked result

The [frozen consumer](../scripts/run_original_collar_signed_average_once.py) executed once. The coefficients below are rounded displays of directed intervals. The full binary endpoints are in the [receipt](../experiments/root_st073/original_collar_signed_average_receipt.json).

| Ordinary axial order | Complete signed primitive in `epsilon² hb Pstar² sigma(S)` | Radial correction in `epsilon² hb Pstar² sigma(S) S` | Relative width upper bound for the correction |
| --- | ---: | ---: | ---: |
| 0 | −28.38182369938787 | +56.76364739877575 | 1.644e−29 |
| 1 | −90.11025533330633 | +180.2205106666127 | 2.717e−29 |

Both complete coefficients are strictly negative. Both radial corrections are strictly positive. The factor `S` is about `3.3049e-204453045017284839`, so isolating it is essential to resolving the correction. The positive interval for this correction uses the same inlet C1 data and the actual collar source identity, with each norm bound used as a radius. The run did not evaluate the exponential source atoms or re-execute an ancestor calculation.

## Next source work

The signed collar primitive is certified at the saved `371/1000` input. The legacy forward six-term ledger uses different axial samples, so it cannot be identified numerically with this interval. Extend the same inlet C1 source provider to axial cells using the actual pressure-slope function and its nonlinear core enclosure. Apply the radial operator over those cells, then compare the matching source terms, their endpoint directions and their existing finite-`N` error leaves. All 212 forward entries and the independent negative pulse target are retained. Neither the cross-owner forward identity nor terminal five-moment matching has been established.
