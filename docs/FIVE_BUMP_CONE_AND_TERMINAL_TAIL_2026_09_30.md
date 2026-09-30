# Actual-source five-bump interval and finite terminal tail

The actual common-source reference background has now been joined to the finite
degree-three bump response throughout a 38-point sample of `1 <= x <= e` at
`Z = .3`. The sample contains all three support boundaries, nine interior points
per support, gaps, and reference-interval endpoints. It reuses only the internally
generated actual source snapshot; the cache cannot provide arbitrary Z data.

## Nominal relaxed-cone sample

Receipt: `experiments/root_st073/lei_ren_part1_paper_five_bump_cone_scan.json`.
Reproduce after generating the local reference snapshot with:

```powershell
python experiments/root_st073/lei_ren_part1_paper_five_moment_reference_background_check.py --resume
python experiments/root_st073/lei_ren_part1_paper_five_bump_cone_scan.py
```

The evaluator uses the same P0, its first axial derivative, all five moments,
and corrected fields. Each sample retains 69 source labels, a reference label,
and a correction label. No pressure reset or moment replacement occurs.

Nominal metrics evaluate the retained P9/W2 jets at pressure = width = 1:

- `a = 1 - 2 Utheta_y / Utheta`;
- `b = 2 Uz_y / Utheta`;
- `D = I_theta / F`, `E = I_z / F`, `w = E / D`;
- `kappa = a + b^2 / a`;
- for the relaxed branch `kappa <= 2`, margin `D(a - b w) - 2a`.

All 38 samples satisfy positive a, positive D, kappa <= 2 and positive margin.
The smallest nominal margin divided by D is approximately
`0.7999999994055487`; maximum absolute b is about `1.186e-16`.
The raw margin has an enormous exponent because the construction uses an
enormous reference radius. It is not a physical residual or error norm.

This is finite sampling of a relaxed criterion, not full-cone certification,
continuous interval enclosure, uniform-Z control, or a proof of the source
quadrature, jets, and coefficient remainders. Very small same-atom contributions
can disappear from nominal totals; source labels and existing monomial receipts
remain authoritative for their hierarchy.

The distinction is paper-faithful: the five-bump region claims only the relaxed
cone, explicitly not the admissible cone (cached Lei--Ren text lines
11042--11045). The sampled branch is Eq. (9.18)--(9.19), with the five-bump
specialization Eq. (10.11)--(10.15). Diagnostics also record F > 0 and the
expected D >= 8, and check that the returned shear agrees with the derivative
definitions of a and b. Other regions still require their own stronger tests.

## Actual finite-response terminal tail

Receipt: `experiments/root_st073/lei_ren_part1_paper_five_bump_terminal_tail.json`.
Command: `python experiments/root_st073/lei_ren_part1_paper_five_bump_terminal_tail.py`.

The script first nominally evaluates the actual P9/W2 input defects, then builds
the scalar degree-six polynomial `d + A h_3 + Q(h_3,h_3)`. Its nonlinear scalar
composition does not apply an additional P9/W2 product truncation; equivalence
to the nominal truncated field jet is not certified.
It evaluates each defect monomial before summing, retaining all finite quadratic
products. This avoids subtracting large physical moment totals. A scalar
aggregate replay agrees to absolute error below `1e-400`, using the same declared
quadrature map; this is not an independent quadrature validation.

| Centered moment row | Degree 4--6 tail at Z=.3 |
| --- | ---: |
| 1 | 0 (linear row) |
| 2 | +8.93108408848e-85 |
| 3 | 0 (linear row) |
| 4 | -8.32484942515e-71 |
| 5 | +5.55929683187e-71 |

The first three formal degrees are cancelled up to finite coefficient arithmetic.
The remaining nonzero tail must not be discarded or called exact terminal
closure. In particular, the fifth-row tail is much larger than the extremely
flat fifth-row input defect. Small absolute values do not establish the
source-relative hierarchy or the terminal identity as a function of Z.

## Next construction work

1. Establish an actual uniform C1/C2 defect bound for the common source, rather
   than adopting the Z=.3 sample as a bound. Recover second-Z data where needed.
2. Couple that bound to the five-bump analytic inverse / contraction majorant,
   with finite coefficient and quadrature errors separately enclosed.
3. Control the full infinite response tail and its axial derivatives; retain
   nonlinear source products separately wherever rounded nominal sums would
   obscure flat rows.
4. Enclose the continuous support intervals and extend cone diagnostics to the
   other matching, flattening, and heat regions.
5. Complete the exact heat and radial-energy matching before implementing
   temporal n-dependent coefficient recursion and oscillatory stress correction.

No temporal recursion, heat certificate, global cone certificate, or full
Cartesian Navier--Stokes residual certificate is claimed by these receipts.
