# Controlled actual source switches through profile R=110

The fresh axial family [0.49,0.51] now extends the actual exit endpoint at R=100 through the two Section9.4 source switches and the constant-power segment to R=110. It uses the same accepted analytic pressure and amplitude data. This is a profile construction, not a completed time-dependent NS solution.

## Two source switches

With x=log(R/100) and hb=0.005, the first interval [0,hb] keeps a=epsilon_exit D while the axial source is multiplied by 1-sigma(x/hb). The second interval [hb,2hb] has zero axial source and blends a from epsilon_exit D to 4/5. The frozen comparison fields are unchanged, while their moments are recovered as functions of the current radius. C2 comparison moments supply C1 drivers. No off-center finite difference or midpoint projection is used.

Each switch has 16 directed cells. On each complete cell, g_x=-a/2 supplies the log-amplitude range; that range supplies the exponential axial RHS range. The resulting actual phi and U ranges enclose all six cumulative moment integrands. Both values and first axial derivatives are propagated by the same Taylor interval algebra. The cell-range integral argument controls the discretization, rather than treating a nominal RK endpoint as an enclosure.

## Exact constant-power segment

After R=100 exp(0.01), phi scales as exp(-2x/5) relative to the switch endpoint and U remains constant. All six moment increments are integrated with exact exponential primitives. Pressure remains the original P0 plus the physical pressure moment. Ur value is reconstructed from the same axial moment and its first axial derivative. Ur_Z remains unavailable.

Independent arbitrary-precision quadrature checks 16 field/moment value and axial derivative coefficients on a known fixture. Seven cutoff checks confirm the increasing sigma convention against the preexisting source definition. These fixtures do not independently validate the complete switch ODE implementation.

## Stress result and scope

Of 32 switch cells, 31 certify a cone branch. The final first-switch cell, x in [0.0046875,0.005], remains unresolved because its kappa bounds cross 2 and have wide dependence. In the second switch Sz=0 identically, so kappa=a and the relaxed margin is tt-(2-a); cancelling the common angular shear algebraically removes spurious interval division. A dedicated fixture checks this cancellation with a spanning 1e-1000 to 0.8. The terminal R=110 packet and the complete constant-power range certify the relaxed kappa<=2 cone. Neither result certifies the stronger kappa>2 cone; subsequent shear modulation and matching are still required. An unresolved cell is not evidence of actual cone failure.

This extends the actual field/moments through R=110. It does not complete functional terminal five-moment identities, heat-exterior matching, whole-axis stress control, parameter/contraction estimates, temporal n-dependent recursion, flat-remainder estimates or oscillatory corrections.

## Next work

1. Resolve the remaining first-switch cell with controlled subdivision or correlated algebra, preserving the tiny nonzero shear and common amplitude.
2. Use these fresh family R110 moments to construct functional terminal defects and independent moment repairs. Do not import the old scalar-center repair as a family identity.
3. Carry the same pressure datum into flatten, collar and heat exterior matching; quantify relaxed-to-strong stress recovery where required by the paper.
4. Continue actual temporal coefficient recovery and oscillatory corrections only after the preceding dependencies are satisfied.
