# Swirl time slopes restore all sampled outer stress directions

Pressure-only diagnosis first rules out a misleading next step. At the
fixed shear-stage velocity, six pressure modes remain infeasible after
removing coefficient bounds. Adding a farther pressure window (.93,.99)
also fails, even with moment constraints removed. Independent single-node
linear tests allow arbitrary axial stress and still find all three
fraction-0.75 nodes infeasible when tangential stress is fixed. This is
a fixed-state axisymmetric-pressure obstruction, not a general PDE
impossibility result: axisymmetric pressure cannot change tangential
momentum stress.

The constructive correction therefore adds compact swirl coefficient
slopes in k alongside pressure modes. The swirl correction is multiplied
by (k-k0), k0=11, so instantaneous velocity, spatial derivatives and
lambda geometry remain unchanged, while the physical time derivative
changes tangential momentum stress. Three outer radial windows times
three axial powers provide nine pressure and nine slope directions.
The farther window compensates total moments without entering the
local cone stress primitives. All four moment constraints and six
outer cone inequalities are linear in these corrections at k0.

Both bounded and unbounded 18-variable linear feasibility checks pass.
A norm initializer and full-momentum cost minimization then converge.
The complete callable pressure/velocity field is independently replayed:

| Check | Result |
| --- | ---: |
| Max moment, 96-point replay | 4.71814e-5 |
| Max moment, 128-point replay | 3.65450e-5 |
| Inner cone nodes | 9/9 |
| Outer cone nodes | 6/6 |
| Largest outer cone ratio | 0.949966 |
| Held-out full momentum peak | 1.42202e6 |

The full sampled momentum peak is unchanged from the preceding
shear-stage field and remains far above 1e-3. The two farther off-midplane
nodes have small negative normal stress projections, about -0.101 and
-0.124, and ratios near the fitted 0.95 margin. Pointwise success does
not establish a useful finite-width pulse support. Time-slope coefficients
reach magnitude 11.125; their effect over a time interval must be evolved
and tested, not inferred from the unchanged instantaneous velocity.

Next screen spatial neighborhoods and nearby times on the same callable
field, then formulate the slope/pressure solve as a state-dependent
continuous evolution with moment and cone margins. Do not transplant
the earlier old-width DAE trajectory. Before pulse realization, require
positive margins on its actual support and use the complete curl and
amplitude/pressure dynamics from the Section 7 route. No finite-energy,
full maximum/volume-L2 or recursive contraction result is claimed.
All outputs remain accepted=false and scale_recursion_established=false.

Reproduce:
```text
python experiments/root_st073/midplane_remote_pressure_feasibility.py
python experiments/root_st073/midplane_outer_slope_pressure_repair.py
```
The `nodewise` array in the combined report records the pressure-only
fixed-tangential-stress diagnostic; `outer_cones` records the final combined
field. They answer different questions and must not be conflated.
