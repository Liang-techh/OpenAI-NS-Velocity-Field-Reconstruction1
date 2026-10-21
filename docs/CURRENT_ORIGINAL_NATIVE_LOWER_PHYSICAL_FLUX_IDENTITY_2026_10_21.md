# Lower physical flux C1 source identity - 2026-10-21

Goal ACTIVE / INCOMPLETE. Saved source result passed=True.

The actual native first-micro initial RQ_y ports 129567/129568 now equal
corrected original lower ports 53450/53462 as ordinary C1 source functions
on the source-owned Z cell [0,1/2]. The original source ports 30442/30443
pass through the saved physical patch 31228/31240 and its recorded relocation
once. The old extra Pstar multiplication is not applied again.

Both functions use the same actual physical U, U_y, own m, delta and radius.
At the radial sc/2 endpoint, the prescribed lazy branch and original Duhamel
lower limit set finite-N corrections to zero. The leading field, own mean,
pressure and nonzero incoming histories remain source functions.

For d=1-Z^2 and L=1-delta*Z^2, N_Q=2ZU-(1-delta)Zm-d*m_Z,
N_y=2Z*U_y-(1-delta)Z*(U-m)-d*(U_Z-m_Z), and
N_plus=(1+delta)ZU+2Z*U_y-d*U_Z.
The new exact polynomial witness proves N_plus=N_Q+N_y through ordinary Z1.
Consequently native Q_y=N_y/L equals corrected q_plus_q_y-q; multiplying
by the same physical Rminus gives the actual lower flux identity.
The mass ODE is m_y=U-m. Physical Q and B are not the loop q and B.

The cell-wide source identity comes from the defining source predicates and
symbolic C1 algebra. The old quantitative packet encloses only Z=0.371 in
an extremely narrow interval; it supplies no new whole-cell magnitude bound.
No midpoint, interval overlap, small error cap or sampled agreement is used
as a function value or as proof of this identity. Global axial closure is false.

The terminal function remains
F-T=(native_Rc-corrected_active)+Rc^2/2*integral_1^2 x*DeltaS dx
    +Rc*DeltaQ_y(1).
The accepted native Rc flux join is retained. The corrected active upper Rc
flux has not yet been identified with the original compact lower flux, so
that inner term is not cancelled. Active/compact cumulative m15 and field
identities remain unproved.

Code/input were frozen, entry persisted before the one new consumer call,
and the complete result saved before review. Source checks=
37; execution seconds=0.7225684999721125.
Only the new universal C1 polynomial ran; no ancestor operator was replayed.
All 212 cumulative / 213 raw / 183 reduced terms, six signed channels, old
errors and one reserved contribution remain frozen. P1 is not m15; ordinary
Z is not temporal n. Runtime model metadata is unavailable; children used
fixed Luna/max roles. Full goal remains incomplete.

Next: Recover the corrected active upper Rc flux and its same-source relation to the original compact lower flux, preserving the accepted native Rc join. Use each actual pressure cocycle before combining the density differences and shifted lower inlet. Keep the compact/active field and history differences explicit until their actual source identities are proved. Do not promote old point enclosures to cell-wide bounds, replace functions with error caps, reset histories, duplicate the reserved slot, or delete old212 errors. Global five moments, exact heat, stress, genuine temporal recursion, corrections and independent Cartesian NS remain open.
