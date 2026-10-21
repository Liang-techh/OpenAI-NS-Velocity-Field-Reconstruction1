# Actual compact pressure density difference - 2026-10-21

Goal ACTIVE / INCOMPLETE. Saved source result passed=True.

The two actual compact fields are now pulled back to the common physical
coordinate x=R/Rc in [1,2], with the accepted original Rc, Pstar and delta.
Native U/Q/Q_y jets come from the C5 physical source. The other field uses
the original C3 band V/m histories, one Pstar, and its exact Q quotient.
Neither incoming history nor compact control branch is identified or reset.

For d=1-Z^2, L=1-delta*Z^2 and a=(1-delta)/2, the physical density is
S=(Q+Q_y+a*Z*Q_Z)/L + Q^2/2 + Q*Q_y
  + U*(-2*Z*Q+d*Q_Z-2*Z*Q_y)/L.
The new source defines its exact difference DeltaS by polarization, retaining
all cross terms and ordinary Z1. Real Q_ZZ, Q_yZ and U_Z ports are used.
New exact polynomial certificates verify both polarization rows, the
continuity Q+Q_y identity and the axial integration-by-parts coefficient.
Only the fresh expression graph is normalized; archived functions stay opaque.

The compact m15 difference is now the explicit C1 function
D_compact = Rc^2/2 * integral_1^2 x*DeltaS dx
            - Rc*(2*DeltaQ_y(2)-DeltaQ_y(1)).
Its physical measure is R^2 dy = R dR = Rc^2*x dx. Both finite endpoint
fluxes retain their source functions and signs. The old reduced S1 formula
represents the full density integral; its factor 1/2 and RQ_y subtraction
are applied once. C3/C5 regularity on the original compact domain justifies
fixed-endpoint ordinary differentiation; no new quantitative norm is claimed.

The saved terminal function is consequently represented as
F-T = (native_Rc-corrected_active) + D_compact.
This is an explicit correlated source-function definition. D_compact and
the active difference are not proved zero, and no numerical defect value,
midpoint, interval overlap or error bound is substituted for either function.
P1 is not m15; physical Q is not loop q; ordinary Z is not temporal n.

Code and inputs were frozen, entry persisted before invocation, and the
complete result saved before review. New consumer calls=1; source checks=
45; execution seconds=1.4188613000442274.
Static preflight and independent review passed. A static unentered draft
had an incorrectly populated decoded-hash map and one AST formatting needle;
the draft and repair receipt are preserved. No mathematical attempt was
entered on that draft. The unused identity-runner draft is also unaccepted.
Runtime model metadata is unavailable; configured child role was Luna/max.

All 212 cumulative / 213 raw / 183 reduced terms, six signed channels, old
errors and one reserved contribution remain frozen. Old operators were not
replayed and the legacy numerical chain was not reconfigured. The earlier
failed source attempt remains pinned through the parent. Full five moments,
original-scale bounds, global physical joins, exact heat, stress/remainder,
true temporal recursion, corrections and independent Cartesian NS remain open.
Numerical baseline remains 8f36f19e3c524b424a9286c46cc8e0b855f75307.

Next: Use the explicit DeltaS and finite RQ_y differences to resolve native versus original compact U/m inlet and control source functions. First test the 2Rc endpoint against the already bound common physical K, then retain the Rc flux correlation with the active difference. Bind the native first-micro initial RQ_y to the corrected original lower endpoint before any active cancellation. Do not reset histories, change the reserved slot, identify source functions from bounds, or remove old212 errors. Global five moments, exact heat, stress, temporal recursion, corrections and independent Cartesian NS remain open.
