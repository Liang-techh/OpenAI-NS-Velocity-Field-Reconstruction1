# Spatial carrier and viscosity scale audit

Run `python experiments/root_st073/curl_wave_scale_audit.py` from the repository root. The detailed 5 x 5 patch grid and integer multiplier sweep are in `experiments/root_st073/compact_potential/curl_wave_scale_audit.json`.

The current two-mode exact-curl wave uses angular modes `m=1,2` at `tau=0.0078125`, with radial and axial envelope halfwidths `0.0025` and `0.00075`. Its radial carrier phase changes by only `0.046` and `0.004` radians across a radial halfwidth; the axial changes are `0.121` and `0.170` radians. The phase is therefore almost constant across each local radial/axial envelope, rather than rapidly oscillating there.

For each mode, the script separates the velocity obtained by taking the curl of the compact potential into its principal carrier derivative and the remaining envelope/cylindrical-frame derivatives. The measured norm ratios are:

| Angular mode | At envelope center | RMS over 5 x 5 patch | Maximum over patch |
| ---: | ---: | ---: | ---: |
| `1` | `0.47` | `10.04` | `48.46` |
| `2` | `0.22` | `4.88` | `24.87` |

Thus the terms treated as curl remainders in the intended oscillatory picture dominate the principal carrier term on much of this patch. This provides a more fundamental explanation for the large supported-wave residuals and rapidly changing time-correction coefficients. It does not prove that every possible wave construction fails.

An illustrative frequency scaling multiplies both integer angular modes and radial/axial wavevectors by an integer `L`, and divides the potentials by `L` to preserve the leading velocity scale. The curl-remainder ratio falls as `1/L`, while the physical viscous exponent `nu |k|^2 dt` rises as `L^2` for the current `dt=2.5e-6`. Even with deliberately lenient thresholds of remainder no larger than carrier and viscous exponent no larger than one per element, mode `1` needs `L>=49` but damping allows `L<=26`; mode `2` needs `L>=25` but damping allows `L<=15`. No integer satisfies both for this unchanged envelope and time step. At the minimum spatial multipliers, the viscous limit would require steps at most `7.3e-7` and `9.1e-7`, respectively. Strong asymptotic separation would demand more.

This is a **physical scale diagnostic**, not a check of the paper's normalized function classes or a proof of a general obstruction. In the [OpenAI Navier--Stokes paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf), the amplitude-pressure equation (7.5) explicitly retains damping and the exact-curl formula (7.38) controls a remainder by scale separation. The present ST073 patch does not demonstrate that separation. The next construction should redesign the wave location/support, carrier, and pulse time scale together before further tuning high-degree time polynomials; then it must recompute the full momentum and mean-flow corrections.

## Fixed-radius axial relocation screen

Run `python experiments/root_st073/curl_wave_location_screen.py`. It keeps the current radius `0.005689` and evaluates the same **physical local cone analogue** and residual primitive at seven axial locations. The compact baseline's positive axial support boundary at this time is `0.004309`. A symmetric envelope centered at the original `z=0.003433` has only `0.000876` of room before that boundary. Moving toward the midplane increases this room, but the sampled local cone does not survive:

| Axial center | Symmetric halfwidth ceiling | Local cone result |
| ---: | ---: | --- |
| `0.00250` | `0.001809` | Fails: target points the wrong way (`T·N=+39.6`) |
| `0.00300` | `0.001309` | Fails: `T·N=+176.1` |
| `0.00325` | `0.001059` | Fails narrowly: cone ratio `1.031` |
| `0.003433` | `0.000876` | Passes locally: cone ratio `0.067` |

The midplane and `z=0.0015` also fail with negative local `lambda_squared`. This is a coarse fixed-radius screen; another radius or a changed baseline might admit a wider patch. It shows why simply shifting this same wave inward is not yet a supported redesign. The next search should map the cone-feasible region jointly in `(r,z)` and compare its usable support width with the carrier/viscous scale constraints.
