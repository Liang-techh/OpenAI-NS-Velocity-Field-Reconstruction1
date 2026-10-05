# Current Rp native pulse source join — 2026-10-05

The current terminal O3 power histories now join the original native CompliantPulseMixedC4 entrance. CurrentRpPulseSourceDispatcher retains all fourteen accepted downstream owners and adds pulse_entrance as owner fifteen. It reuses the original pulse implementation and the same current pre-pulse history object; no saved power-inlet sample is promoted to a source definition.

## Implemented

The join binds the actual slope, axial turnoff, slope-mu and power source assignments and their scalar kernels. For q=1+Z² it derives the terminal functions u=U/q, m=MZ, h=H/q, k=KZ/q, e=E_Z Z²+E_Q/q² and p=P_in/q² from those assignments. Twenty-eight stage coefficient-function identities and five canonical unit identities are exact symbolic source relations.

The native inlet obtains H and P from its live SharedOuterBuffer.power(0,1) call. The actual axial-high/fifth-order selection, publication and native data path are replayed, establishing 54 inlet Taylor coefficient identities through order five. The current and native sides retain the same analytic P0(Z), its defining datum, the Cstar record and original parameter definitions.

Actual pre-pulse and native algorithms establish 135 velocity, pressure and five-primitive mixed derivative identities through total order four at Rp. Full radial velocity factors, Pstar normalization and pressure derivatives are retained. The original pulse entrance function and its ordinary jets through order four vanish there. The original logRp and positive mu/Tw source definitions are bound separately.

The generated report evaluates 135 current left-boundary rows, 60 native right-boundary rows and 60 native entrance rows over the full original Z∈[-1,1], t∈[0,.02/mu] domain. The upper endpoint enclosure covers the exact endpoint. These are interval source enclosures; they are not selected production point values.

## Reproduce and use

Run `python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage actualrpjoin` on branch codex/st073-transition-next. The stage runs the producer and focused checker. The default CurrentRpPulseSourceDispatcher consumes the checked receipt before reporting current_Rp_external_pulse_join_certified=true.

Use `evaluate("pulse_entrance", Z, t)` for ordinary log-radius/axial mixed source grids. The same native pulse object is returned by `provider("pulse_entrance")`. Its coordinate is t=log(R/Rp), with xi=mu*t. The existing fourteen-chart physical adapter remains in its accepted historical scope until an explicit extension consumes this new dispatcher.

## Scope and next work

This completes F42-Rp-source-recipes, F42-Rp-canonical-units, F42-Rp-neighbor-and-radius and F40-current-Rp-pulse-join in the native entrance source scope. Earlier pending F41/F42 source-transfer items are superseded here. The saved CompliantPowerInletC4 diagnostic sample path is not used for this admission and is not retroactively certified.

Next, expose the other five original pulse charts through this same native object and consume their functional internal interface certificate. Then extend the current Cartesian/time adapter to the resulting twenty source owners, keeping original radius and amplitude units. Shared full leading-source/nonlinear-remainder admission, production point histories/fields, global completed stress/cone/flatness/required-domain energy, true n-dependent temporal recursion and oscillatory corrected dynamics remain open. Radial Taylor coefficients and coordinate maps are not temporal recursion.
