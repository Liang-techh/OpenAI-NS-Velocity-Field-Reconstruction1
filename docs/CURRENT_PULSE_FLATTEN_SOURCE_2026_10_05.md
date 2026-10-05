# Current native pulse terminal to O5 flatten - 2026-10-05

CurrentPulseFlattenSourceAssembly adds the original 100-unit flatten to the accepted current source graph: twenty-one profile source owners, while physical Cartesian/time ownership remains at twenty. CurrentFlattenMixedC4 inherits the original flatten ODE implementation and CurrentPowerInlet inherits the original canonical incoming functions.

## Current defining inputs

The new inlet uses the current native pulse object's U constants, live buffer-derived Hp/Pin/Xp and analytic pressure datum. Xv is evaluated from the native terminal source formula, including its signed nonzero memory. For every Z, future_energy calls the exact same native selection.future.future callable and takes complete_future_energy_Taylor/2. No saved terminal Xv, energy sample or fitted polynomial supplies the new owner.

The original exact positive amplitude has log((Ev0/Pstar)^2) = 2 log(U) - 13/mu - 26. The source radius is log(Rv) = log(Rp) + 13/mu. Absolute pressure retains P0 plus the original FTC Mp primitive. Small positive caps only enclose these exact source functions; they do not replace them.

At pulse s=0, the two future end-bump supports are empty, so the two linear moments and the axial/radial velocities vanish by their actual defining integrals. The angular moment remains Xv and energy remains the completed positive future integral divided by two. Incoming nonzero Rp reference histories have not been reset to enforce matching. The inherited incoming() object remains Rp reference data used only for u, Mp and P0 in the original integration formula. terminal_histories(Z) separately publishes all five actual Rv primitives plus absolute pressure from the same native end(Z,0) call. Every new flatten packet includes that adapter; no unused Rp linear moment is presented as an Rv terminal input.

## Admission and output

The canonical pulse_end_flatten_join_check receipt retains 155 functional, 34 inlet/datum and 24 endpoint source identities. The new owner consumes that theorem after binding the unchanged native methods, current object identities, actual endpoint/flatten assignments, class-scoped Xv/amplitude construction and the exact current future_energy return. Numeric interval overlap is not the source join proof.

The new report evaluates the whole real Z domain at flatten inlet, across t in [0,100], and at the exit. These three packets contain 180 ordinary mixed logR/Z derivative source rows through total order four for velocity and pressure. The focused checker also requests a fresh unsaved Z and checks empty supports, positive future energy, original domains and current dependency hashes.

Run `python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage currentpulseflatten`. The checked default CurrentPulseFlattenSourceAssembly exposes provider(chart) and evaluate(chart,Z,coordinate); existing twenty routes remain with the prior dispatcher. The new flatten route uses velocity/Ev0 and pressure/Pstar^2 units. It does not install a twenty-first physical Cartesian owner.

## Next work

Extend the new source through post-flatten power/angular/steep/waiting/heat owners with the same five histories, datum and source radii. Add the current flatten Cartesian/time adapter separately. Uniform quantitative native pulse mixed-four interfaces, complete nonlinear production point fields, global tensor/cone/flatness/required-domain energy, true n-dependent temporal recursion and corrected dynamics remain open.
