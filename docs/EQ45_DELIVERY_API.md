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

## Support-connected candidate and portable files

Use the current support-connected public API:

```python
from openai_ns_reconstruction.eq45_supported_delivery import velocity, u, v, w
uvw = velocity(.1, 0, .1, .5)
```

Create matching candidate JSON, Python NPZ and MATLAB MAT samples:

```sh
python -m openai_ns_reconstruction.eq45_export_bundle --output artifacts/delivery/eq45_supported
```

The committed bundle contains9 times and17 points per spatial axis. Arrays
u/v/w use (time,x,y,z), with a combined velocity shape(9,17,17,17,3).
The manifest binds canonical candidate identity and all file hashes.
The generated evaluate_velocity.m supplies a MATLAB griddedInterpolant example;
it interpolates saved samples and is not exact off-grid candidate evaluation.
MATLAB itself was not available/run; MAT dimensions and sample identity were
checked by Python loadmat against the public evaluator. Scientific readiness
flags remain false. No separate animation work is needed to evaluate the field.

## New central-motion candidate (explicit opt-in)

The live static seed does not reproduce the source's central opposite axial
outflow at the tested symmetric points: w(.1,0,+/-.1,.5)=+1.23219 on both sides;
the radial component reverses sign across z=0. The existing eta-even Phi(1,0)
temporal schedules alone do not remove that parity constraint.

The new candidate uses the existing Phi(0,1)=1 and zeros other poloidal
coefficients, preserving all swirl coefficients, coefficient bounds, geometry,
Eq45 map and physical-support transform. This is an explicitly chosen qualitative
seed, not an inferred hidden coefficient or final candidate selection.

```python
from openai_ns_reconstruction.eq45_supported_delivery import Eq45SupportedDeliveryField
field = Eq45SupportedDeliveryField.load_candidate(
    'artifacts/delivery/eq45_bipolar/candidate.json')
u, v, w = field.components(.1, 0, .1, .5)
# (-0.06730245983355401, 0.13998343506571057, 0.17340970258309657)
```

At z=-.1, radial and swirling velocities are the same and w changes sign.
At96 fresh central similarity-coordinate points, at3 off-keyframe times,
all96 show inward radial flow, positive swirl and axial velocity directed away
from z=0. Fourth-order sampled central divergence max is below4e-10. This is
local structural evidence only, not whole-domain support/PDE/visual acceptance.

Candidate canonical SHA:
8027075948fc4cbbf0c79cf080729b543c7a7d00b011a4b35b582b95deab5970.

The eq45_bipolar bundle includes candidate.json, velocity.npz, velocity.mat,
evaluate_velocity.m, file manifest, the base-versus-child central probes, and
fresh-point central_motion_holdout.json. The unchanged default seed has a
separate bundle at artifacts/delivery/eq45_supported. Do not mix their hashes
or inherit default-candidate validation results for this changed field.

Source interpretation: docs/VISUAL_TARGET.md (official announcement and paper
Section2.1). No animation was generated; velocity functions remain the priority.

## Energy-normalized source-aligned candidate (named direct API)

The separately serialized source-aligned child with the recorded positive
energy-normalization scale can now be evaluated without knowing its artifact
path:

```python
from openai_ns_reconstruction.eq45_normalized_bipolar_delivery import velocity

uvw = velocity(0.1, 0.0, 0.1, 0.5)
```

This wrapper reconstructs and identity-checks the already committed candidate
with SHA-256
`c0e27269adfb6f305c2c0b5a2483f9eddd54f79691f69cacd7bb735ba592a702`.
It binds source bipolar SHA
`8027075948fc4cbbf0c79cf080729b543c7a7d00b011a4b35b582b95deab5970`
and the previously measured common scale `1.8097870818686452`; no coefficient
is fitted or changed by the delivery wrapper. `default_field().at_points(...)`,
`grid(...)`, `save_candidate(...)`, and `load_candidate(...)` use the same
Cartesian `[u,v,w]` contract as the support-connected API.

This is an **experimental source-aligned, energy-normalized candidate**, not a
canonical or visualization-selected field. Existing momentum/energy-balance
obstructions remain in force, so `visualization_ready=false`,
`visual_correspondence_verified=false`, and `pde_validated=false`. The
serialized candidate's existing provenance vocabulary is preserved unchanged;
its separate schema/provenance audit is outside this delivery increment.
