# Centrifugal pressure correction and its interface datum

The remote moment field now has an optional pressure increment satisfying the [OpenAI paper's](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf) leading radial relation `Pi_X=E^2/(2X)` for the *change* in swirl on the heat patch. It integrates the difference from the unmodified heat profile. The two available integration constants preserve either the inner pressure (`pressure_datum='inner'`) or the outer pressure (`'outer'`). Velocity and divergence are unchanged by this option.

At `tau=.0084, eta=.2, X=10`, the finite-difference radial derivative of the added physical pressure is `0.2059888920`, while the change in `u_theta^2/r` is `0.2059888924`; the discrepancy is about `4.3e-10` for both datums. This checks the intended radial identity at a nontrivial point.

Full momentum tells a different story. At `eta=.2`, comparing an inner probe, three patch probes, and one exterior probe:

| Pressure option | Inner residual norm | Patch sampled maximum | Exterior residual norm |
| --- | ---: | ---: | ---: |
| Unchanged | `1.03e-6` | `5603.68` | `2.49e-8` |
| Preserve inner datum | `1.03e-6` | `5603.59` | `0.17049` |
| Preserve outer datum | `0.17049` | `5603.68` | `2.49e-8` |

The same transfer is seen at `eta=.35`, where the unwanted quiet-region residual is about `0.23477`. At `eta=0`, parity makes that transfer tiny, which is why off-midplane probes are essential. The patch residual itself remains thousands and is predominantly axial/poloidal; satisfying the centrifugal pressure derivative does not solve it.

The pressure offset at the outgoing interface is forced by the nonzero `Cp` difference measured in `ST073_FIVE_MOMENT_AND_SHEAR_LOOP.md`. A single radial integration constant cannot preserve both inner and outer pressure data when that difference varies with eta. This is a concrete reason to solve the five-moment interface and axis-pressure datum together, rather than tune the pressure only after fixing the velocity patch. The correction is optional and the existing uncorrected diagnostic remains the default. Neither option is an accepted field.

Reproduce with `python experiments/root_st073/remote_pressure_datum_screen.py`; see `experiments/root_st073/remote_pressure_datum_screen.json` (`accepted: false`).
