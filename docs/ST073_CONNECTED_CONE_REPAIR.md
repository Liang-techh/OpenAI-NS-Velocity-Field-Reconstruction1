# Joint moment and support-cone mean repair

**Moment-closure claim superseded by quadrature audit.** The small moment
values below used 12 Gauss points over the entire bridge and do not survive
radial refinement. Splitting at the correction supports and using 48 points
per panel gives maxima about `213,3372,53240` at `k=11,15,19` for this same
candidate. It is not moment closed. See `ST073_MOMENT_QUADRATURE_AUDIT.md`.
The recorded cone samples remain useful diagnostics; they do not rescue
the failed moment condition.

The widened wave's inner/lower support corner failed the physical stress
cone on both tested scales. A joint constrained mean solve now removes
this sampled obstruction at `k=11,15,19` while retaining small sampled
outer moments. This supplies a new candidate mean, not a residual-
contracting wave or a continuous cone certificate.

The optimization uses the existing 36 separated mean coefficients.
Velocity and its differential jets are affine in the coefficients, so
their cached full momentum and radial stress integrals are quadratic.
SLSQP minimizes scaled coefficient displacement subject to twelve
outer-moment equalities and cone, sign, and positive-growth inequalities
at nine physical support nodes per scale. The radial and axial training
offsets are `(-0.99,0,+0.99)` times the support half-widths. The axial
half-width remains `1.4` times the original wave width. The target cone
ratio is at most `0.95`.

The solve converged in eleven iterations. Direct evaluations of the
resulting field reproduce the predicted cone ratios within `7e-10`.

| Scale | Training nodes passing | Independent nodes passing | Largest independent cone ratio | Local mean momentum max / prior mean |
| ---: | ---: | ---: | ---: | ---: |
| `11` | `9/9` | `16/16` | `0.815162` | `1.09444` |
| `15` | `9/9` | `16/16` | `0.793527` | `1.11213` |
| `19` | `9/9` | `16/16` | `0.781688` | `1.10871` |

The independent grid uses offsets `(-0.95,-0.3,+0.3,+0.95)` in each
direction. Cached normalized moment error is `9.53e-14`. Recomputing
the moments from the actual field, using the same 12-point Gauss panels,
gives maximum absolute error `7.69e-6` and normalized error `1.60e-10`.
The difference reflects finite-difference/cancellation effects and is
reported rather than substituting the cached result for a direct check.
These values are sampled moments, not full spatial momentum norms or
quadrature-independent identities.

An earlier fit only constrained offsets through `0.8`; it passed its
nine training nodes but failed one of sixteen independent edge nodes
at each scale. That result is retained in the non-`edge` JSON files.
Moving the constraints to `0.99` and including the middle scale fixed
the reported holdouts. The geometry improvement costs roughly `9–11%`
in local mean momentum, so it must earn that cost through a supported
wave correction before acceptance.

Reproduce:

```text
python experiments/root_st073/midplane_connected_cone_repair.py --edge
python experiments/root_st073/midplane_connected_cone_holdout.py --edge
```

The candidate coefficients are in `midplane_connected_cone_edge_repair.json`.
Next check intermediate scale support and interior pulse times, rebuild
the stress-matched wave consistently with this mean, and demand lower
direct complete momentum after a coupled amplitude/pressure/mean update.
Continuous support, endpoint control, arbitrary scale recursion, and the
requested maximum/volume-L2 residual tolerances remain unestablished.
