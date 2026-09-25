# Five-moment interface audit and the shear-loop route

The [OpenAI Navier–Stokes paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf), Section 4.2 and Appendix C, separates two operations that our first remote patch had conflated. Equation (4.15) tracks **five** cumulative radial quantities:

| Name | Integral from the axis to `X` |
| --- | --- |
| `M` | `U` |
| `I` | `H=sqrt(2X) E` |
| `J` | `U H` |
| `S` | `U^2-E^2/2` |
| `Cp` | `E^2/(2X)` |

The prior `[4,16]` patch targets only global angular and kinetic identities analogous to `I` and `S`. It does not match all five values at the outgoing interface. A direct 80-node patch audit at `tau=.0084, eta=0` gives the following differences between the two-moment field and the *unchanged* joined field at `X=16`:

| `delta M` | `delta I` | `delta J` | `delta S` | `delta Cp` |
| ---: | ---: | ---: | ---: | ---: |
| `~0` | `+0.74085` | `-0.03395` | `+3.58165` | `+0.001671` |

These are **differences from the unchanged exterior**, not defects relative to a fully specified paper target. They show that the present patch cannot be regarded as a complete interface match to that exterior. In particular, its physical pressure is currently copied unchanged from the joined field. Integrating the changed centrifugal pressure relation with the same axis datum would require an outer pressure increment of approximately `+0.00209` at this slice. The unchanged pressure is therefore not a dynamically consistent continuation of the changed swirl.

A second even meridional streamfunction mode was added, with `B(y)(s^2-beta)` and `s=2y-1`; `beta` makes it orthogonal to the original bump in the `B^2` weight. The four sampled odd/even amplitude pairs `(-5,-12),(-5,12),(-3,-12),(-3,12)` produced **zero strict local physical-cone passes out of 36 points**. At `X=10,eta=0`, pair `(-5,12)` gives a cone ratio near `0.25` but stress-normal projection near `+7.28`, while `(-5,-12)` gives a negative projection near `-0.77` but cone ratio near `34.4`. The large `-5` pairs also retain a fitted kinetic moment defect up to `0.00177`; they are not accepted moment closures. The two-shape field changes `J` at `eta=0` by about `-0.237` relative to the unchanged exterior, while `Cp` remains about `+0.00167`.

The Appendix-C construction does **not** require the unmodulated joined shear to satisfy the strict admissible cone everywhere. It first requires a **relaxed cone** on the transition interval, then builds a rapidly varying periodic shear loop whose mean is the old shear. Integrating at high radial frequency changes the profiles and five moments only slightly; a reserved correction interval restores all five moments. Our physical-field cone is an analog, not this normalized relaxed-cone test. Further coarse amplitude grids of the current bumps will not implement that mechanism.

The next implementation should compute the normalized profile variables and all five target moment differences on the same slices, including a pressure profile consistent with `Cp`; test the relaxed cone `(4.21)` and its boundary collars; then construct a shear loop and five-moment restoration on reserved intervals. Only after that should the oscillatory stress and full physical momentum gates be assessed. The current heat exterior still lacks axial localization and finite total energy.

The first normalized snapshot screen is now in `ST073_NORMALIZED_RELAXED_CONE.md`: the baseline and both remote mean variants fail the relaxed gate at every sampled point. A shear loop is therefore premature for these specific inputs.

Reproduce the five-moment audit with `python experiments/root_st073/remote_five_moment_audit.py` and the two-shape screen with `python experiments/root_st073/remote_two_shape_screen.py`. Both JSON reports are marked `accepted: false`.
