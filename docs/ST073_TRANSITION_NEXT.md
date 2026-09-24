# ST073 continuation: local kernel to matched scale recursion

The current user objective requires nonzero divergence-free finite-energy 3D
velocity, shrinking/elongating/winding cores, dynamically matched inner,
transition and outer regions, and complete physical momentum max and spatial
volume L2 below 1e-3 on an explicit domain/forcing contract. Python/MATLAB and
interactive delivery remain required. This is not the old visualization-only goal.

Imported frozen user delivery under experiments/root_st073; all 173 manifest
payloads verified and original ZIP SHA retained. No bundled sync script executed.
Current main integration at 6a293b3c does not contain this complete bundle; work
is isolated on codex/st073-transition-next without disturbing prior local edits.

Host progress: archived <f16 NPZ arrays cannot load on this Windows NumPy.
host_interface.py regenerates float64 interface traces from unchanged ST073-V
model at k=0,3,6 (17 eta nodes each). Maximum computed complete momentum residual
is 5.315e-7; pressure-gradient difference from supplied JSON <=2.701e-13.
This is same-model replay, NOT a new independent operator or whole-domain test.
Original files stay byte-identical. Windows longdouble epsilon is2.22e-16.

Next construction tasks:
1. Use actual corrected velocity/pressure jets at X=1/64, |eta|<=.5; rebuild
   transition compatibility rather than reusing ST068 leading-only moments.
2. Derive moving-interface normal velocity, mass flux and stress from source
   coordinate map, including time-dependent boundary motion; include axial caps.
3. Construct a solenoidal transition with matched velocity and stress; evaluate
   its full residual, not only continuity or divergence. No zero-extension claim.
4. Match an outer finite-energy solution with explicit forcing; do not select
   arbitrary residual-canceling force. Add nonaxisymmetric corrections only with
   an explicit closure/realizability check.
5. Distinguish prescribed k-dependent sampling from dynamically generated scale
   recursion; test material winding and shrinking/aspect metrics over the full
   registered window. No extension toward tau=0 without new evidence.

Local scope remains X<=1/64, |eta|<=.5, tau in[.5/64,.5], nu=.01 unforced.
No global finite-energy field, exterior match or scale-recursion acceptance yet.

## Moving interface mass budget completed

`experiments/root_st073/moving_interface.py` integrates the curved side and both
flat caps at k=0,3,6 using orders12 and24. Side inflow balances cap outflow;
net fluid flux at order24 is below2.8e-20. The relative outward flux is NOT zero:
it equals minus the shrinking domain volume rate. At k6 it is1.702657768e-5.
Order12/24 results agree; this is numerical integral consistency, not a rigorous
uniform bound or momentum acceptance. Interface labels move at b_r=-r/(2tau),
b_z=-(.5-h)z/tau. Treating this interface as a material no-through-flow surface
would contradict the frozen kernel. The transition must transport the measured
side/cap flow and match momentum stress; a closed impermeable shell is unsuitable.
Report: `experiments/root_st073/moving_interface/mass_flux.json`.

## Stress and moving momentum interface data

`interface_stress.py` constructs the physical Cartesian gradient directly from
frozen radial coefficient jets, sigma=-pI+nu(grad u+grad u^T), and moving flux
u((u-b).n)-sigma n. Both side and caps are included, orders12/24, k0/3/6.
Portable nodewise arrays include geometry, gradients, stress and oriented flux.
At k6/order24, total axial momentum outward flux=2.515468846e-7 and angular
momentum outward flux=1.080789491e-8. These are integrated fluxes, NOT residual
norms or admission thresholds. Gradient trace max5.68e-14, curl/archived-evaluator
vorticity replay difference1.14e-13; same-model algebra consistency only.
Next use these stresses together with mass/velocity traces to build a transition;
check volume momentum-rate plus boundary flux before selecting an exterior.
Artifacts: `experiments/root_st073/moving_interface/stress_flux.json` and
`stress_k0.npz`, `stress_k3.npz`, `stress_k6.npz`.

## Radial continuation constructed

An explicitly separate order10 continuation keeps ST073-V axis data and nu,
extends Xmax from1/64 to3/64 (radius factor sqrt3), and preserves solenoidality
through the same full recurrence. Frozen model unchanged. At the original
interface k6, order8/10 velocity difference4.37e-13 and pressure6.02e-12.
The new annulus X in[1/64,3/64], |eta|<=.5 has k6 sampled boundary max5.748e-5;
12x18 volume quadrature max5.463e-5 and physical L2=4.668e-9, volume1.780e-7.
Orders8x12 and12x18 were evaluated at k0/3/6. This is a finite local extension,
NOT decay to an exterior or proof across all points/times. Tiny volume explains
part of L2. Independent FD and separate directional convergence remain pending.
At X=1/16 (twice original radius), order10 boundary max1.0296e-3 fails the gate;
recorded wider-radius failures must remain. Next use the new outer boundary
for a dynamical transition, or assess higher-order radial continuation with
independent derivatives; do not mask outer errors by cutting the field off.
Reproduce: radial_continuation.py then annulus_check.py in experiments/root_st073.
Artifacts and explicit new model: experiments/root_st073/radial_continuation/.

## Independent annulus derivative check

annulus_fd.py reuses the independent Cartesian/time operator, not recurrence
residual assembly, at six off-calibration points (X=.026,.044; eta=-.22,.37;
k=.4,2.7,5.5; nonzero azimuth). Spatial and temporal steps varied separately.
All FD momentum maxima remain below1.031e-5, with finest late value1.021e-5,
well below1e-3 at these points. This supports using the extended boundary as a
matching datum; it does not prove a uniform bound or an exterior solution.
Artifact: radial_continuation/independent_fd.json under experiments/root_st073.
The next unresolved construction is still dynamical outer decay/axial closure,
not further repetition of these local residual checks.
