# Stress target after remote moment closure

`ST073_REMOTE_MOMENT_MEAN.md` constructed a solenoidal radial-patch mean field whose two Appendix-A-style moments have small held-out defects. The [OpenAI paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf) next requires an admissible stress and a supported oscillatory realization; moment closure alone does not supply either.

`remote_moment_stress_screen.py` integrates the **full physical momentum residual** radially to three radii inside the patch (`X=6,10,14`) at `eta=0,.2,.35`, `tau=.0084`. It compares the original pure-heat outer field with the moment-patched field and evaluates the same physical-field analog cone used in previous local screens. Piecewise radial quadrature splits at the inner/bridge/patch interfaces. The result is diagnostic: the paper's cone belongs to its normalized leading background, not automatically to this corrected finite-window physical field.

All nine patched sample locations fail the strict local analog cone. At `X=10`, `target·N` is positive at all three axial positions:

| `eta` | Heat outer `target·N` | Moment-patched `target·N` | Patched point residual norm |
| ---: | ---: | ---: | ---: |
| 0 | 0.041 | 3.499 | 2.62e3 |
| 0.2 | 0.040 | 6.240 | 4.61e3 |
| 0.35 | 0.036 | 10.663 | 4.38e3 |

The patched local `lambda²` remains positive at these samples; the main obstruction here is stress orientation, with some other radii also failing the cone-ratio bound. A bounded sample failure is neither a theorem of infeasibility nor a statement about a constructed wave. It does mean the new mean field cannot be handed to the existing pulse screen as if moment closure had solved the stress stage.

A concrete next shape degree of freedom is a **moment-null swirl bump** in `X∈[4,16]`: start with `B(y)(2y-1)` and subtract its weighted projection onto `B` so its integral against `√(2X)` is zero. Its coefficient then changes the local stress without changing the angular moment `I_2`. For each coefficient, re-solve the quadratic meridional amplitude to keep `I_1=0`, construct the streamfunction field, and search for a strict cone over an axial/radial region. Only after that succeeds should supported phase/amplitude waves be built. The full-volume momentum and smooth compact-force gates remain separate requirements.

Reproduce with `python experiments/root_st073/remote_moment_stress_screen.py`. Results are in `experiments/root_st073/remote_moment_stress_screen.json`, marked `accepted: false`.

The first moment-null swirl-shape trial is documented in
`ST073_REMOTE_NULL_SHAPE.md`; its five screened amplitudes did not fix the
stress-cone orientation.
