# Projected local Kelvin amplitude along the transported phase

The [OpenAI Navier–Stokes paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
retains viscous damping, background shear, and the normal-pressure
constraint in its principal pulse equation (7.5)–(7.6). The preceding
ST073 experiment transported the phase only and found that a moving
exact-curl packet still had a very large full momentum residual.

This experiment advances a **centerline physical-coordinate analogue**
of the homogeneous transverse amplitude equation together with the
moving phase. In backward time `tau=T-t`, the Cartesian amplitude has
`a_tau = J a + nu |k|^2 a - 2k(k·J a)/|k|^2`; the rotating cylindrical
frame adds its connection term. The last term is the projected
pressure contribution that preserves `k·a=0`. This is not the paper's
normalized chart equation, supported slow-box amplitude, or pressure
field solution.

For the current two integer angular modes over
`tau=0.00828..0.00852`:

| Mode | Norm at early end / center | Norm at late end / center | Maximum relative `|k·a|` error |
| ---: | ---: | ---: | ---: |
| `m=1` | 0.958 | 0.922 | `3.9e-11` |
| `m=4` | 0.186 | 2.335 | `6.9e-10` |

At `tau=0.00846`, where the transported exact-curl wave had a much
larger momentum residual than the frozen packet, the `m=4` centerline
amplitude is already `1.669` times its reference value; `m=1` is
`0.978`. Because `tau` increases toward **earlier physical time**, the
positive `nu |k|^2` term here represents ordinary forward-time viscous
decay. At the reference time the backward viscous rates are about
`641` and `6092`, respectively. The high-mode amplitude evolution
cannot be assumed to suppress the observed off-center envelope
velocity peak (`591.8` versus `9.20` for the frozen wave).

The simple exact-curl packet's separate scale audit requires carrier
multipliers at least `411` and `165` to make the sampled envelope
remainder no larger than the carrier; at those multipliers its
viscous exponent reaches one after only about `9.2e-9` and `6.0e-9`
units of time, far shorter than the tested `1.5e-4` halfwindow.
These are lenient diagnostics of the current packet, not a theorem
excluding the paper's multiscale construction.

The next wave construction must first obtain wider effective
transverse scale separation or a different supported potential, then
solve a spatially dependent amplitude/pressure equation and test the
full momentum. Merely inserting this centerline amplitude into the
current narrow packet is not a plausible residual repair.

Reproduce:

```powershell
python experiments/root_st073/kelvin_amplitude_transport.py
```

The report is `experiments/root_st073/compact_potential/kelvin_amplitude_transport.json`.
