# Automatic axial jets of the prescribed joint exit

The local Section 9.25–9.26 exit now differentiates its initial data, auxiliary drivers and prescribed finite RK state together. `AxialDual` wraps the existing pressure/width ring, so each retained pressure-tail and width atom has both a value and its first Z derivative. This does not use finite differences in production.

The component core retains a second axial Taylor jet at every radial row. This is necessary because the first derivative of the auxiliary stress involves the derivative of a first-Z field or moment. The shared nonlinear core evaluator and stress formulas are reused; only their scalar algebra and axial coordinate are changed.

`AxialPressureWidthExitBridge` exposes all five moment Z derivatives and recovers Ur from Mz and Mz_Z with the original similarity-coordinate formula. Angular and axial velocities, pressure, stress and a local physical/Cartesian chart use that same exit data. Analytical shear expressions preserve tiny increments that subtracting nominal velocities would erase.

## Checks and scope

- The independent component/scalar core fixture still passes with maximum scaled difference 2.5729667790e-101.
- The auxiliary dual fixture agrees with a separate resolved centered-difference replay of fields, pressure and drivers to maximum scaled difference 6.6538938460e-16.
- Run `python experiments/root_st073/lei_ren_part1_paper_pressure_width_axial_bridge_fixture.py` for the independent resolved prescribed-exit derivative check.
- That resolved prescribed-exit replay passes at s=0.5,1,1.5,2: maximum scaled first-Z derivative discrepancy for F, Uz, P and all five moments is 1.16837500572e-18. Independently reconstructing Ur using the prescribed Mz finite difference gives a maximum scaled value discrepancy of 1.60095798564e-13, including the finite ring evaluation effects. The nonzero Cartesian chart smoke passes. This does not validate a second Z derivative or Ur_Z.
- Run `python experiments/root_st073/lei_ren_part1_paper_pressure_width_axial_bridge_check.py` for the actual degree-18, pressure-order-9/width-order-2 source receipt, with delta=1e-200 and h_b=exp(-100-1e154). The receipt records boundary agreement, retained tail derivatives, first-width moment derivative increments, radial velocity and stresses.

The actual-source run completed: Uz_Z and P_Z retain nonzero pressure-tail atoms, and every one of the five moment Z derivatives has a nonzero first-width endpoint increment. At the initial boundary F and Uz match the core exactly at working precision; P and moment-derived Ur differ only by tiny arithmetic terms, with scaled errors below 1e-200. The nominal (pressure=0,width=0) algebraic divergence numerator is about -3.56e-261; this number has only the identity scope described below.

The automatic derivative is that of the finite RK algorithm and retained finite jets. Pressure/width truncation, core/axial truncation, analytic datum aggregation and ODE errors are not enclosed. The reported divergence numerator uses the moment derivatives prescribed by the ODE; it is an algebraic identity, not an independent radial derivative of a finite RK interpolant or a Cartesian divergence certificate.

## Next construction dependencies

1. Use the same derivative-aware state in continuation, switching and long reshape; recover radial velocity consistently there. The legacy ExitContinuation freezes nominal values once their positive-epsilon changes fall below working precision. That policy cannot be transferred directly to this atom-preserving route: retain the epsilon-driven field/moment changes and their Z derivatives as separate width atoms, rather than reporting their positive shears on frozen values alone. Use a finite log-radius coordinate after the collar, not astronomical s=y/h_b.
2. Close the five terminal identities as functions of Z and restore compatible angular/axial/pressure data together.
3. Resolve the nonzero radial tail and establish physical finite energy and exact/controlled heat exterior.
4. Establish admissible cone margins and flat remainder before true n-dependent temporal recursion and oscillatory cancellation.

No complete joined field, functional moment closure, finite energy, admissible cone or temporal recursion is claimed by this local adapter.
