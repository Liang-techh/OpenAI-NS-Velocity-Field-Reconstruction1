# Transported phase experiment and the remaining wave obstruction

The [OpenAI Navier–Stokes paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
uses a phase with an evolving normal, a projected amplitude/pressure
equation that retains viscosity (Section 7, Equations 7.3–7.6), and a
full curl of a compact potential (Equation 7.38). The prior ST073
prototype used a fixed radial/axial wavevector and frozen amplitude.

For the current 10% localized-poloidal candidate, a local physical
characteristic phase was implemented as

`Phi=m*theta+kr(tau)*(r-rc(tau))+kz(tau)*(z-zc(tau))+phi0(tau)`.

The angular integer `m` remains fixed. The moving center follows
backward-time characteristics, and `kr,kz,phi0` are advanced to cancel
`(-d_tau+u·grad)Phi` through first spatial order at the center. On a
`5 x 5 x 5` moving sample box, the maximum eikonal defect fell from
`164.2` to `16.5` for `m=1`, and from `963.4` to `87.9` for `m=4`.
The RMS defects fell from `71.2` to `5.70` and `433.5` to `26.6`.
This is a physical-coordinate local transport calculation, **not**
the paper's normalized-chart estimate or its full amplitude equation.

The transported phase and moving support were then inserted into a
spatial exact-curl vector potential. Its ring covariance at the
reference time is inherited from the original local fit. The complete
physical momentum residual on a 16-angle ring did **not** improve:

| Backward time | Base max | Frozen-wave max | Transported-wave max |
| ---: | ---: | ---: | ---: |
| 0.0084 | 54,626 | `9.78e8` | `9.87e8` |
| 0.00846 | 38,232 | `3.71e8` | `2.46e10` |

At the later time, the transported wave speed on that ring reached
`591.8` versus `9.20` for the frozen packet. The narrow moving axial
envelope makes the exact-curl derivative remainder large off its
center; evolving the phase alone does not control it. Finite-difference
divergence estimates for the transported wave at steps `4e-6,2e-6,1e-6`
decrease `390.9,98.1,24.5`, consistent with second-order truncation
of an analytically divergence-free curl. The previous carrier/damping
scale audit also has no overlap for this simple packet.

The next wave attempt must solve a transverse **amplitude and pressure
evolution** compatible with the time-dependent phase, while changing
the spatial/time support or background scale separation enough to
control the curl-envelope remainder and viscosity. A different carrier
frequency alone cannot repair the tested packet. This experiment does
not establish a global wave construction or a Navier–Stokes solution.

Reproduce:

```powershell
python experiments/root_st073/transported_phase_screen.py
python experiments/root_st073/transported_curl_wave.py
python experiments/root_st073/transported_wave_diagnostics.py
```

Reports are under `experiments/root_st073/compact_potential/` with
matching file stems.
