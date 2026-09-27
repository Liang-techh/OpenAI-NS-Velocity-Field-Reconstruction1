# Local value/slope restoration of interscale moments

The smooth interpolated outer repair had moment defects 0.239 and 3.541
at k=13 and 17. A local coefficient correction now separates state
values from derivatives with respect to k=-log2(2 tau).

For a velocity-linear correction basis V_j, physical-time differentiation
adds V_j a'_j(k)/(tau ln 2) to momentum. Integrating this term gives an
explicit 4-by-12 moment response matrix D. Its singular values are
(63.63,38.76,2.34e-15,1.83e-15) at k=13 and
(1032.28,628.88,4.39e-14,3.01e-14) at k=17: only two directions survive.
The axial response entries are at roundoff, consistent with
integral r u_z dr = psi(outer)-psi(axis)=0 for compact poloidal bumps.
Thus slopes can restore angular transport but cannot independently fix
axial moments. A coefficient-value solve is required for the latter.

The experiment first solves the two axial moments with outer poloidal
coefficient values, then solves the two angular moments with outer swirl
slopes. Values and slopes are realized by a single local callable field
with coefficients affine in k, and full field derivatives are replayed.

| k | Value correction norm | Slope correction norm | Max moment, order 96 | Max moment, order 128 | Inner cone |
| --- | ---: | ---: | ---: | ---: | --- |
| 13 | 0.0016209 | 0.0086413 | 3.50e-5 | 7.10e-5 | 9/9 |
| 17 | 0.0013005 | 0.0069609 | 3.66e-4 | 2.14e-3 | 9/9 |

The local state/slope correction strongly reduces the intermediate
moment defects without disturbing the sampled inner cone. At k=17,
however, the higher-order replay still fails the absolute 1e-3 moment
gate. Differences between orders must not be hidden by the near-zero
algebraic prediction. Full sampled momentum remains about 7.79e7 and
4.91e9; no full momentum acceptance is claimed.

These are two local affine trajectories, not an integrated differential-
algebraic system or a common corrected trajectory across the interval.
Next use integration-by-parts moment identities or an independently
checked derivative strategy to separate quadrature and finite-difference
cancellation. Then evolve angular coefficient moments while enforcing
the axial algebraic constraints along one continuous trajectory. Preserve
support separation and verify full interval residuals and contraction.
Neither more independent knot fits nor a local slope correction alone
establishes scale recursion.

Reproduce:
```text
python experiments/root_st073/midplane_outer_moment_dae.py --k 13
python experiments/root_st073/midplane_outer_moment_dae.py --k 17
```
All outputs remain accepted=false and scale_recursion_established=false.
