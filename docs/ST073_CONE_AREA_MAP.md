# Annular candidate cone-area screen

The coupled annular velocity/pressure candidate was screened at backward time
`tau=0.0084` on 9 radii and 9 heights. The candidate and physical radial
primitives are the same as in `ST073_ROBUST_CONE_COUPLING.md`. A node is marked
passing only when the strict local cone test passes and its cone ratio is below
`0.8`.

Only **10 of 81** nodes pass. No sampled rectangle with two distinct radii and
two distinct heights has all four corners passing. Passing nodes cluster partly
near `r=0.0055`, `z=0.0027–0.0032`, but adjacent radial samples fail; other
passing nodes are isolated. Thus this candidate has no *demonstrated* area in
which the compact exact-curl packet can be widened in both spatial directions.
The coarse scan does not prove that no smaller continuous patch exists.

This strengthens the scale obstruction in `ST073_WAVE_REENTRY.md`: the narrow
packet there needs carrier multipliers at least `529` and `175` to control its
envelope remainder, while its viscous damping permits only about `21` and `6`
over the chosen time window. Merely increasing the carrier frequency cannot
close that gap. The next construction step is to modify the underlying
velocity/pressure candidate to make a connected, spacetime-persistent cone
region, then repeat the packet-scale screen. The present candidate and this
map are experimental, not a Navier–Stokes solution certificate.

Reproduce with:

```powershell
python experiments/root_st073/annular_cone_area_map.py
```

Machine-readable point diagnostics and the pass grid are in
`experiments/root_st073/compact_potential/annular_cone_area_map.json`.
