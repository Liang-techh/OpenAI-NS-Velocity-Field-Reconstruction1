# Fixed-pressure spacetime cone screen

The previous radial pressure patch passed a `7 x 3` spatial grid at
`tau=0.0084`. Holding its amplitudes fixed, 58 of 63 nodes passed on
`tau=0.00825, 0.0084, 0.00855`. A **constant** additional combination of
the same 16 compact radial/axial pressure shapes is feasible on all 63
nodes; no time-dependent pressure coefficient was needed. Direct physical
finite-difference residual evaluation at those three times and two
intermediate times gives **105/105** passing local cone nodes on
`r=0.0075..0.013`, `z=0.00355..0.0037`,
`tau=0.00825..0.00855`. The maximum sampled cone ratio is `0.794315`.
This is a sampled box only, not a continuum certificate.

The new pressure combination **worsens the global momentum objective**.
Its physical-volume L2 on Gauss-6 increases from `183.81` to `193.02`
at `tau=0.00825`, `182.21` to `191.54` at `0.0084`, and `180.87` to
`190.32` at `0.00855`. Gauss-8 likewise increases from `126.25` to
`131.75`, `129.66` to `135.02`, and `132.97` to `138.21`. The maximum
residual is near `1e5`–`2e5` in these samples. The quadrature orders
disagree, so neither is a converged global error estimate. Keep this
pressure combination as a **local cone experiment**, not an accepted
Navier–Stokes field.

At box center `(r,z,tau)=(0.01025,0.003625,0.0084)`, the local stress
primitive has cone ratio `0.6534` and a two-ray nonnegative covariance
fit with zero numerical relative error. Yet the existing two-carrier
compact exact-curl prototype still has no lenient envelope/damping scale
over the box's `0.00015` time halfwidth. With radial halfwidth `0.00275`
and axial halfwidth `0.000075`, its required carrier multipliers are at
least `411` and `140`; damping exponents at most one allow at most `3`
and `1`. These are diagnostics for this simple packet, not an
impossibility theorem for the construction in the
[OpenAI Navier–Stokes paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf).
The paper's transported phase, amplitude equation, and exact-curl
potential require a different quantitative construction than this frozen
packet.

The next step is a volume-aware pressure/velocity optimization subject
to a **margin** on the sampled spacetime cone, then a wave design whose
transport and viscous damping fit inside that margin. Pressure-Poisson
compatibility, radial moment closure, full spacetime residual bounds,
and the critical-time scale recursion remain open.

Reproduce the screens and volume comparison with:

```powershell
python experiments/root_st073/radial_pressure_time_screen.py
python experiments/root_st073/radial_pressure_time_fit.py
python experiments/root_st073/radial_pressure_time_validate.py
python experiments/root_st073/radial_pressure_time_volume.py
python experiments/root_st073/radial_peak_cone.py --field radial-pressure-time --radius 0.01025 --z 0.003625 --tau 0.0084 --output-name radial_pressure_time_source.json
python experiments/root_st073/curl_wave_scale_audit.py --source-name radial_pressure_time_source.json --output-name radial_pressure_time_wave_audit.json --radial-halfwidth 0.00275 --axial-halfwidth 0.000075 --time-step 0.00015
```

JSON results reside in `experiments/root_st073/compact_potential/` under
the corresponding command names.
