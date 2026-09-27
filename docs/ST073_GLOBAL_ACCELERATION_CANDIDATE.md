# Unified global acceleration candidate

`global_acceleration_candidate.load()` now returns one callable full field:
`u, p = field.fields(points, tau)`, with physical time t=-tau. Points are an
N by 3 Cartesian array; u is N by 3 with components (u,v,w). This is a
research candidate, not an accepted Navier--Stokes trajectory.

The construction combines the globally localized balanced mean/wave,
compact exterior collar tangent, and compact inner acceleration correction.
It verifies the acceleration parent hash, common reference time, and strict
disjointness of the three collar boxes from the fixed wave support. Compact
support and curl construction are retained. The external correction is
based on the enriched mean report, but the balanced candidate differs only
inside the wave support, so that source substitution leaves the collar
field unchanged.

At reference time and k0+1e-6, sampled inner velocity/pressure agree exactly
with the acceleration candidate; sampled outer values agree exactly with
the collar candidate. Axis values are finite, and an off-support axial point
returns zero. These are assembly checks, not a complete momentum replay.
The shape response of the acceleration term is still being checked.

Do not add or average the constituent residual norms: their reported domains
and quadratures differ and leave unmeasured regions. Whole-domain and
multi-time momentum acceptance, terminal regularity, and scale recursion
remain open.

```python
from global_acceleration_candidate import load
field, _, _, snapshot, _, _ = load()
velocity, pressure = field.fields(points, snapshot['inputs']['mean']['tau'])
```

The completed actual shape replay now shows that this full acceleration
candidate fails relative axial elongation: endpoint aspect changes by
-7.44786e-7 relative to reference. Contraction and weighted angular-speed
increase pass. The unified callable is therefore a diagnostic assembly,
not the selected goal-consistent field. Constrained acceleration fitting
must repair the aspect direction before adoption.

## Full candidate support bounds

The unified support helper now accepts additional collar boxes. At reference
and the checked endpoint the complete bounding cylinder has radius 0.0105
and z in [-0.001,0.001]. The earlier mean-plus-wave bounds alone omitted
small parts of the added collar support. Assembly probes beyond the corrected
union bounds return exactly zero velocity and pressure. Full-domain residual
integration must use these enlarged bounds; the older narrower domain is
not sufficient for the combined field.
