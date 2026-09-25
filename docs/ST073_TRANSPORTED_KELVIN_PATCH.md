# Transported Kelvin amplitudes on a moving patch

`experiments/root_st073/transported_kelvin_patch.py` evolves the two
existing wave modes' projected Kelvin amplitudes on a 3 by 3 moving
radial/axial patch near `tau=.0084`. It interpolates those amplitudes in
a compact vector potential and takes its analytic cylindrical curl.
The resulting velocity is divergence-free by construction. A scalar
pressure proxy comes from the local Kelvin projection; this is not the
coupled pressure PDE of the OpenAI paper's Section 7.

Complete Cartesian finite-difference momentum was evaluated at 16
interlaced patch locations and eight angles at each of two times:

| Remaining time | Mean without wave | Frozen wave | Transported-amplitude wave |
| ---: | ---: | ---: | ---: |
| `.0084` | `1.14e5` | `1.45e10` | `1.45e10` |
| `.00846` | `1.23e5` | `6.04e9` | `1.31e10` |

These are sampled full-momentum maxima, not volume L2 or continuous
suprema. The transported amplitudes and pressure proxy fail to offset
the narrow envelope's curl/viscosity cost; at the second time they
worsen the frozen wave. This packet is rejected (`accepted:false`).

The experiment supplies a spatially varying amplitude field that the
earlier centerline Kelvin ODE lacked, but does not solve the paper's
phase-amplitude-pressure system, the quadratic mean correction, radial
moment constraints, or a scale-recursive family. A viable wave attempt
first needs an admissible stress cone across the defect-bearing
transition region and a support/carrier/time-scale design with actual
cutoff and viscous separation. The local ST073-V core's dyadic transfer
is documented separately in `ST073_CORE_DYADIC_TRANSFER.md`.
