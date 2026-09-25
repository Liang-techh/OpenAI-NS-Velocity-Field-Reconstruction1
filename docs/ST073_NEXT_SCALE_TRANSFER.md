# One new scale knot does not reduce the momentum defect

`experiments/root_st073/separated_moment_next_scale.py` extends the
exploratory time registration from `k=20` to `k=24`, retains all 36
coefficients at `k=11,15,19`, and fits only 12 new swirl/poloidal radial
window coefficients at `k=21`. This tests an actual one-knot continuation
of the current finite ansatz. The default time registration and earlier
reports are unchanged.

| Scale | Frozen-last-row momentum max | New-knot momentum max |
| --- | ---: | ---: |
| `k=19` | `1.879e9` | `1.879e9` |
| `k=20` | `5.279e9` | `5.274e9` |
| `k=21` | `1.483e10` | `1.480e10` |
| `k=22` | `4.166e10` | `4.158e10` |

The four sampled physical outer moment components at `k=21` have
normalized maximum `0.4266` when the `k=19` coefficients are held fixed.
The new knot brings this to `8.53e-14`, with a successful constrained
optimization. Yet the `k=21` direct Cartesian momentum peak drops by only
`0.198%`, and the `k=22` holdout peak continues to grow. From `k=19` to
`k=22`, the new-knot peak grows by a factor of about `22.1`, corresponding
to roughly `2^1.49` per dyadic step. The correction therefore repairs
sampled moment leakage but does not improve the residual exponent.

This experiment extends a registered exploratory time slab; it does not
validate the order-12 core there. Its constraints use two axial positions,
two physical moment components per position, and a 15-point momentum grid.
The divergence figures in the JSON are finite-difference diagnostics. No
continuous moment, stress cone, volume norm, global support, or uniform
recursion bound follows.

The [OpenAI paper, Proposition 9.6](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
uses a full cycle: supported nonzero-harmonic amplitude inverse, primary
stress update, compact mean correction, and radial moment repair. It
recomputes the complete residual and gains a positive power at each cycle.
Our new-knot repair supplies only the last piece in a sampled analogue.
The next implementation should couple a supported transported wave and
mean/stress correction to moment repair, then require a reduction of the
*complete normalized momentum defect* at a new scale and its holdout.
