# Two-dimensional local cone and candidate-envelope screen

Run these commands from the repository root:

```powershell
python experiments/root_st073/curl_wave_cone_region.py
python experiments/root_st073/radial_peak_cone.py --radius .003 --z .0025 --output-name radial_peak_cone_r003_z0025.json
python experiments/root_st073/curl_wave_cone_patch_screen.py
```

The coarse screen evaluates the current baseline's **physical local stress-cone analogue** at seven radii and eight axial locations at `tau=0.0078125`. It integrates the full baseline momentum residual radially using shared quadrature segments. At the original point, its stress target agrees with the independently generated reference within relative `1.1e-8`. Eleven of 56 sampled points pass the strict local cone test. The broadest geometric axial room among passing points occurs at `(r,z)=(0.016,0)`, but its stress magnitude is only about `0.525` and local momentum residual about `500`, so it does not target the large collar defect.

One more relevant alternative is `(r,z)=(0.003,0.0025)`: its local cone ratio is `0.396`, and two positive Kelvin covariance weights (`182.87`, `126.22`) reproduce its stress target exactly in the local frozen-ray calculation. Its geometric axial support room is `0.001809`, about twice the original point's `0.000876`. This is a **pointwise** opportunity, not yet a supported wave.

The refined line screens expose the missing uniformity. At `z=0.0025`, sampled radii `0.0015` through `0.00325` pass, while `0.0035` through `0.0045` fail. At `r=0.0015`, sampled axial locations `0.002`, `0.0025`, and `0.003` pass, while `0.0015` and `0.0035` fail. Thus a large envelope centered at the alternative point would cross local cone failures. A smaller support would revive the carrier/cutoff scale problem recorded in `ST073_WAVE_SCALE_AUDIT.md`.

All cone calculations here use the current full residual primitive and are **not** the normalized leading stress cone of Section 7 of the [OpenAI Navier--Stokes paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf). Local positive weights do not verify the phase/amplitude equation, exact-curl remainder, mean correction, pressure, or global momentum. This candidate is not accepted. The next design decision is whether to modify the baseline collar so a strict cone holds on a larger connected patch, or to use a smaller patch with a correspondingly shorter and properly damped pulse.
