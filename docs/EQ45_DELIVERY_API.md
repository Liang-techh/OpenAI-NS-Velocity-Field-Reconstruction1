# Named Eq. (4.5) velocity delivery API

The integrated Eq. (4.5) candidate has a named public wrapper that keeps it
separate from the older package default field:

```python
from openai_ns_reconstruction.eq45_delivery import velocity, u, v, w

uvw = velocity(0.1, 0.0, 0.1, 0.5)  # ndarray shape (3,)
ux = u(0.1, 0.0, 0.1, 0.5)
vy = v(0.1, 0.0, 0.1, 0.5)
wz = w(0.1, 0.0, 0.1, 0.5)
```

`velocity(x,y,z,t)` broadcasts NumPy-compatible coordinate/time inputs and
returns one final component axis `[...,3]`. `Eq45DeliveryField.at_points(...)`
uses the same frozen candidate on arbitrary point batches, and
`Eq45DeliveryField.grid(...)` returns `[time,x,y,z,component]` samples suitable
for downstream Python/MATLAB export.

The default wrapper is exactly `Eq45VelocityCandidate.seed()`. A checked JSON
candidate can instead be loaded with
`Eq45DeliveryField(candidate_path=...)`, and the wrapped candidate can be saved
again with `save_candidate(...)`; canonical candidate SHA256 is exposed as
`field.sha256`.

This named API does **not** replace
`openai_ns_reconstruction.velocity_components:velocity`, whose current default
family remains `coupled_velocity_v1`. It also does not change the Eq. (4.5)
field numerically. The Eq. (4.5) candidate remains an autonomous bounded
visualization/research candidate: candidate-local velocity export is ready,
while physical support, visual correspondence, PDE validation, paper-exact
identity, hidden OpenAI-field identification, and blow-up claims remain false.
