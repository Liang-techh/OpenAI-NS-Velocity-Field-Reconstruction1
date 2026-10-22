# CR-A9-038 — retained ST006 radial circulation profile

## Task and direct delivery value

This increment measures the **actual retained ST006 public `[u,v,w]`** on fixed center-plane rings. It is intended to expose, before expensive rendering, the candidate-side radial profile of circulating speed, circulation, inward radial motion, and their finite-window time trend. This directly lowers the effort needed to judge whether a delivered callable field has the qualitative swirl/radial structure visible in public material.

It is deliberately distinct from the existing cylindrical swirl/poloidal energy and pitch fingerprint: that lane summarizes global component balance, while this lane preserves a radius-resolved center-plane profile and binds it to the published ST006 identity.

## Fixed measurement contract

The contract is autonomous and fixed before inspecting ST006 values:

- physical plane: `z=0` in the registered Cartesian frame;
- times: `0.25, 0.50, 0.75`;
- radii: 31 equally spaced values from `0.05` through `1.50`, kept away from the singular cylindrical axis and the `r=2` support face;
- angle samples: 64 uniformly spaced endpoint-excluded samples on each ring;
- no camera rotation, spatial registration, hidden-time alignment, target fitting, component rescaling, or visual acceptance threshold.

For each ring the public Cartesian velocity is converted to cylindrical `u_r,u_theta,u_z`. The report records ring means, mean absolute circulating speed, mean speed, `Gamma(r)=2*pi*r*mean(u_theta)`, ringwise angular-deviation RMS, peak radii/speeds, inward/outward radial-sign counts, and a radial swirl-weighted radius. Endpoint deltas are descriptive candidate-side observations only.

## External result screened / migration

- source repo: `scipy/scipy`
- screened commit: `7253df3f99e564fdf3978506bfe51bd9de0d1c73`
- relevant public API: `scipy.integrate.simpson`
- source location: `scipy/integrate/_quadrature.py`
- license: SciPy BSD 3-clause terms
- classification: **direct migration / public API only**
- copied source: none
- dependency delta: none (`scipy>=1.10,<2` is already required)

Only composite Simpson integration is delegated to SciPy, for the radial integrals used in the swirl-weighted radius. Ring construction, Cartesian-to-cylindrical conversion, candidate identity binding, report hashing, nontrivial sampled-swirl guard, provenance, and truth-state semantics are repository-local.

The `128*eps` sampled-swirl floor is a floating-point zero guard relative to sampled speed and radial span; it is **not** a visual or PDE acceptance threshold.

## Constraint governance / truth boundary

No canonical constraint changes. The registered problem remains `nu=0.01`, physical `R^3`, evaluation box `[-2,2]^3`, smooth compact support `r<2, |z|<2`, time `[0.25,0.75]`, preregistered restricted two-parameter forcing only, `E(0.25)=1+/-0.001`, separated optimization/validation samples, and unchanged divergence/PDE thresholds.

This diagnostic does not compute PDE derivatives and may not be used as PDE acceptance evidence. It does not modify velocity, pressure, forcing, normalization, seeds, support, or thresholds. It does not infer OpenAI numerical velocity, physical radius scale, hidden frame times, camera pose, streamline seeds, or hidden parameters. Visual similarity cannot promote PDE validity.

Hard false in the report: `visualization_ready`, `visual_correspondence_verified`, `pde_validated`, `paper_exact`, `openai_field_identified`, and `blowup_proved`.

## Public-observable relationship

The profile can later provide candidate-side evidence for the qualitative public statement that circulating speed varies with radius. Its endpoint summaries can also participate in a separately governed review of contraction/speed evolution. They are not by themselves evidence that ST006 corresponds to the OpenAI field, because no quantitative OpenAI radial data or physical frame-time mapping is public here.
