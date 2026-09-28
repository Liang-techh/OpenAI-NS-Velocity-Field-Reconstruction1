# Full-support momentum budget

The unified diagnostic field is now replayed at reference and k0+1e-6 over
a single fixed cylinder containing the moving mean, compact wave, exterior
collars, and a two-spatial-step stencil margin. The initial narrower domain
omitted collar support and is retained as failed evidence in the report.
The corrected domain has radius 0.0105015623 and z approximately
[-0.0010015623,0.0010015623]. Outside probes return zero velocity and pressure.

Order-4 panel quadrature with 20 angles has 32,000 points. Its disjoint
wave/collar/rest partitions exactly cover the quadrature. Actual Cartesian
momentum gives:

| State | Volume L2 | Sampled maximum |
|---|---:|---:|
| Reference | 1,952,778.4253 | 7.3489982e10 |
| Endpoint | 1,881,817.1077 | 6.9469369e10 |

At reference, the wave patch contributes 98.50% of squared L2, the collars
1.40%, and the remainder 0.098%. Endpoint shares are similar. The next
construction should prioritize inner residual structure and shape-constrained
acceleration rather than further exterior-only fitting.

This coarse rule has only 640 wave-patch points, so its sampled maximum misses
peaks visible on the denser local grids. Quadrature and continuum maximum
convergence are unverified. Do not compare these maxima with different-grid
numbers as evidence of improvement. The unified diagnostic uses the old
unconstrained acceleration and fails relative axial elongation. Lower
endpoint residual alone is not a successful NS step or scale recursion.

Reproduce with `python experiments/root_st073/global_full_momentum.py`.
