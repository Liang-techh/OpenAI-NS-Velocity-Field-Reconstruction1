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
