# Remote patch radius and width against the relaxed cone

The [OpenAI paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf), Appendix A, uses a large radial scale in its outer construction. We tested whether moving or widening our simpler two-moment bump is sufficient to create the normalized relaxed-cone margin required before its Appendix-C shear modulation.

The correction field now accepts a remote patch interval and refits its eta-dependent coefficients. All moment and cone quadratures split the long heat interval at successive factors of four; without those splits, the largest-scale two-moment defects were visibly contaminated by quadrature error. With the corrected partition, the sampled maximum two-moment defect stays near `9.88e-4` across the tested scales.

Moving a fixed-ratio patch gives **zero relaxed-cone passes in all 45 samples**: each of `[4,16]`, `[16,64]`, `[64,256]`, `[256,1024]`, and `[1024,4096]` was tested at three eta slices and three relative radii. At the first sampled radius, `v_s - U(P_c,J_c)` falls from about `863` for `[4,16]` to about `24` for `[1024,4096]`, but remains positive. At the outer sampled radius, `P_c` is negative at every scale.

Widening a patch anchored at `X=4` also gives **zero passes in 36 samples** for ends `64,256,4096,16384`. At the midpoint and `eta=0`, the strength and allowed upper bound approach one another:

| Patch | `v_s` | `U(P_c,J_c)` | Gap |
| --- | ---: | ---: | ---: |
| `[4,64]` | 4357 | 497 | 3860 |
| `[4,256]` | 3679 | 1620 | 2059 |
| `[4,4096]` | 3571 | 3320 | 251 |
| `[4,16384]` | 3610 | 3533 | 77 |

The last gap is numerically resolved: doubling the radial Gauss order from 20 to 40 and halving both finite-difference steps changed the `X=8194,eta=0` values from `(3609.984,3532.895)` to `(3610.007,3532.918)`, leaving the gap essentially unchanged. Near the inner side of that patch (`X=2734`), the refined values are `(37.10584,36.75056)`. The outer side still has `P_c<2`. These are finite-window point samples, not a continuum claim or an asymptotic limit.

The widening trend identifies a more promising profile-design region, but scale or width alone does not produce an admissible input for the paper's shear loop. The five outgoing moments and physical pressure remain unmatched; the exterior is not axially localized, and the full momentum and volume-L2 gates are far from acceptance. A next profile search must alter the **shape and pressure datum while enforcing five-moment compatibility**, rather than rely on bump placement alone.

Run `python experiments/root_st073/remote_patch_scale_screen.py` and `python experiments/root_st073/remote_patch_width_screen.py`. Their JSON reports retain `accepted: false`.
