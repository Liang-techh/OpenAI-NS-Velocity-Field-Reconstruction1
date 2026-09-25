# Broader axial return-flow collar for the finite-energy field

The first full-space ST073-V candidate used an axial cutoff from
`|eta|=0.30` to `0.45`. This made a finite-energy solenoidal field, but
produced a very steep late collar. `axial_cutoff_width_screen.py` changes
only the cutoff geometry; the wide local core, radial join, pressure
coefficients, viscosity, and unforced residual operator stay fixed.

At `k=5.5`, the complete Cartesian momentum was sampled at two radial
locations and three *relative* axial positions inside each collar:

| Unchanged interior to `|eta|` | Zero by `|eta|` | Largest sampled residual | Vector RMS |
| ---: | ---: | ---: | ---: |
| 0.30 | 0.45 | 212636.58 | 78121.07 |
| 0.30 | 0.49 | 94670.27 | 36520.13 |
| 0.25 | 0.49 | 52336.75 | 19416.50 |
| **0.20** | **0.49** | **33111.75** | **11838.02** |

The last geometry lowers this six-point maximum by a factor of about 6.4
relative to the original cutoff. At the **same physical slice**
`eta=0.375`, its late inner/bridge residuals are `13310/31101`, versus
`43664/94693` for the original. At `k=0.4` those same inner/bridge values
are `68.7/161.7`, versus `225.3/494.3`. The radial heat-exterior sample
at `eta=0.375`, however, rises from `190.8` to `338.4` late, so the new
geometry is not a uniform pointwise improvement. A separate `k=3` six-point
screen of the selected collar has maximum `2504.37`.

The selected cutoff is now the default in `AxiallyCompactField`. It still
keeps a nonzero unchanged central slab, vanishes smoothly by `|eta|=0.49`,
is solenoidal by the streamfunction identity, and has finite energy at each
registered time. Its finite-difference divergence at one late inner collar
point falls `6.58e-7 → 4.11e-8 → 2.59e-9` when the spatial step is
halved twice, consistent with fourth-order differencing.

These samples do not establish a full-domain maximum, physical-volume L2,
moment match, or acceptable Navier–Stokes residual. Tens of thousands are
still far above `1e-3`; the remaining axial pressure and return-flow
dynamics cannot be supplied by cutoff-width selection. The next
construction should use the paper's coupled inner/outer transport and
moment constraints rather than tune this scalar cutoff further.

Run `python experiments/root_st073/axial_cutoff_width_screen.py` and
`python experiments/root_st073/axial_compact_join.py`. The JSON reports
beside those scripts retain all sampled values and the rejected PDE status.
