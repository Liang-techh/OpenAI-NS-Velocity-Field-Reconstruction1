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
