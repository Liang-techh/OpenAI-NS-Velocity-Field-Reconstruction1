## 2026-09-30: actual five-moment bundle installed; pressure mismatch isolated

Shared profile pressure and partial quadratic APIs are installed. Nominal bounded partial-energy and pressure receipts are available. Integrated actual five-moment/pressure/velocity stress replay at Z=.3 (post-Rv flattening and heat) completes but yields enormous stress/shear ratios; coupled matching is not closed.

The post-support pressure combination H=2(1+delta) Z P-(1-Z^2) PZ is approximately -0.01936304114, with P=-0.01208803828 and PZ=0.01330793206 at both points. N_z equals the pressure contribution to I_z when axial velocity/jets vanish; the saved I_z and N_z agree to their serialized precision. This identifies pressure compatibility as an immediate actionable obstruction, while smaller moment terms remain unresolved. Do not remove it with a Z-dependent pressure gauge reset.

Artifacts: experiments/root_st073/lei_ren_part1_paper_continuous_{partial_axial_energy,pressure_moments,moment_bundle,bundle_stress_check,pressure_tail_defect}.{py,md,json} as applicable. Finite energy, admissible stress, full residual and scale recursion remain uncertified. The previously materialized nonzero radial-energy tail remains open.

## 2026-09-30: quantitative terminal radial-energy obstruction audit

Added continuous_radial_energy_tail.py/md/json for the current regenerated candidate. It retains Mz/Mz_Z and C=[(1-delta)Z Mz+(1-Z^2)Mz_Z]/(1-delta Z^2), giving physical ur=-nu*C/r. At fixed tau, radial kinetic-energy density per dZ dlogR is pi*nu^2/2*(dz/dZ)*C^2 with the explicit physical axial Jacobian.

At Z=.3, nu=1 and tau=exp(-4), materialized C remains nonzero and the pointwise logarithmic energy density is positive. Post-Rv and heat mass/transport agree nominally to about 1.35e-292; density to about 2.70e-292. These are arithmetic consistency results, not accuracy certificates for the residual. No Z interval is integrated, no terminal value is forced zero, and the exact functional coefficient/integral closure remains unproved.

This quantitative audit keeps the full finite-energy requirement visible while partial-energy and pressure worker checks continue. Next accept their separate-component results, run unified moment/stress diagnostics, and enclose actual integral inputs/terminal transport. The current numerical radial tail must be resolved before any global finite-energy or recursion claim.

## 2026-09-30: analytic radial velocity jets prepared for stress

Added velocity_radial_jets(profile,logR,Z) using the same installed angular slope, incoming cutoff derivative, live pulse product jet and owned end-bump derivative. It supplies Utheta_y and Uz_y in logR coordinates alongside the existing Z jets, without changing coefficients or terminal moments.

Independent shared-field point derivative checks passed in incoming, pulse, end-bump and heat regions; maximum local relative difference is 1.61e-16. These are nominal functional checks, not input uncertainty or stress-cone certification. The helper is available for the next same-candidate stress adapter after actual five-moment/P providers are accepted.

Partial axial energy and preheat pressure workers remain active. Their Python processes were confirmed live during this turn. No incomplete worker dependency or prepared local quadratic dispatch is published as completed work.

## 2026-09-30: shared pressure heat cumulative and infinite-tail atoms

Added pressure_increments(t,Z) and complete_pressure_heat_integral(Z) to the common heat-moment owner. They use the pressure weight Utheta^2/(2R), exponent -(1+delta), and the same retained heat polynomial/collar atoms as the installed angular field. Reference, tiny correction and correction Z remain separate. Infinite-tail analytic polynomial truncation is bounded; quadrature/arithmetic are not.

Shared-field finite pressure-integrand checks at t=.6 and 4 passed with maximum relative difference 9.50e-14. Infinite-tail Z difference passes nominally at working precision, and Rtail increments vanish. Actual preheat/inner pressure provider remains in progress; full pressure matching, five moments, finite energy and recursion are not claimed.

The partial-axial worker's Python PID44536 was confirmed live during this turn. Parent has prepared lazy quadratic dispatch locally but will publish it only with the accepted provider, keeping uploaded code free of unfinished dependencies.

## 2026-09-30: incoming swirl reference regenerated coherently

Replaced inherited incoming I_swirl by the same normalized ContinuousAngularMoments propagation used by actual angular cumulative moments. Reference energy is generated before actual inner offsets are applied once; energy target and live axial coefficients are then recomputed together. Reference swirl/Ep^2 remains Z independent, preserving the analytic tangent cancellation for this input.

Shared-field test passed: actual Rp swirl and Z primitive match regenerated reference plus measured inner offsets to about 1.17e-84 nominal relative difference. The change from inherited I_swirl is about 7.92e-19; this is not a global accuracy certificate. Updated linear-equation replay errors are about 6.87e-174 and 1.09e-173, and the live coefficient tangent evaluates. Quadrature/input uncertainty, terminal radial transport and complete five moments remain open.

Pressure cumulative implementation and partial axial-square integration are assigned to bounded Luna/max workers. Their unfinished files are not accepted or uploaded as completed work. Next integrate their verified results, complete full pressure/heat targets and stress diagnostics.

## 2026-09-30: cumulative heat angular moments and infinite swirl component

Extended actual seeded angular and swirl-square moments/Z jets beyond Rtail using the installed common heat polynomial. Collar Gauss atoms and analytic exterior exponential atoms retain reference, heat corrections and Z terms separately. Actual quadratic moments now dispatch through finite heat radii. Six forward-integrand checks at t=.6,4,30 pass with maximum relative difference 7.69e-14; Rtail increments are zero and corrections remain nonzero.

Added complete_swirl_heat_integral for Rtail..infinity Utheta^2 dR, retaining its reference/correction/Z atoms and an analytic integrated polynomial-truncation bound. Positive delta makes this swirl radial integral converge. This does not include physical Z weights or radial kinetic energy; total finite energy and recursive closure are still unproven. Quadrature/arithmetic and full heat-target uncertainty remain open.

Next regenerate exact heat-defect/angular/pressure targets from these shared atoms, regenerate inherited incoming swirl coherently, complete pressure moments/Z, finish partial axial energy and resolve terminal radial transport.

## 2026-09-30: continuous heat point jets installed

Installed the new common heat functional in ContinuousAngularSchedule and propagated separate heat deficit/log-amplitude/slope correction atoms to actual angular field adapters. Preserved nonzero heat Z derivatives beyond Decimal exponent storage using arbitrary-exponent MP values. Shared-field collar and exterior checks passed: retained-log radial derivative maximum relative difference 1.83e-12; Z derivative and actual field jet nominal agreement at working precision. The live nominal future-energy target is regenerated on construction; exact heat-defect integral targets still retain their explicitly declared Taylor approximations. Complete heat cumulative moments, input bounds, finite energy and recursion remain open.

Evidence: continuous_heat_install_check.py/md/json. Next integrate the retained heat atoms into full cumulative angular/quadratic/pressure moments, regenerate correction targets coherently, and complete partial axial energy.

## 2026-09-30: separate heat deficit/jet provider prepared

Added continuous_heat_kernel.py/md/json: one Gamma integral for H/H-prime/H-second, normalized positive deficit quadrature, separate logH and tiny derivative correction, plus a common small-x quadratic Taylor approximation with explicit analytic truncation bounds. A first run exposed absolute-tolerance loss in the tiny deficit integral; normalization by h fixed it. Functional derivative checks now pass with maximum relative difference 4.28e-19; exp(-1e6) argument retains nonzero deficit. Arithmetic/quadrature are not enclosed. Provider is not yet installed in the shared field or targets.

Next install using separate heat atoms, regenerate angular/pressure/energy targets consistently, and continue cumulative heat moments. Never interpret this standalone provider as full heat or finite-energy closure.

## 2026-09-30: mixed moments and complete exterior axial energy installed

Added actual inner-seeded mixed moments/Z jets and complete axial-square/Z jets, using the live continuous incoming, pulse and end atoms. Installed linear_axial_moments_jet and quadratic_moments_jet on ContinuousIncomingProfile. The quadratic provider currently supports Rv through Rtail. Actual inner offsets are applied once; materialized terminal residuals remain nonzero. Pulse/end supports are disjoint by xi = 13 + mu*t_v > 11 on the end band.

Fixed the tiny positive pulse startup primitive using a positive scaled endpoint integral. The startup checks retain a nonzero value at xi=1e-5; full source energy integration runs again. Mixed Rp matching is about 2.28e-81; local mixed integrands about 2.08e-23 and 1.19e-17. Exterior quadratic integrand checks are about 7.34e-15 and 1.88e-14; axial energy is constant after axial support. These are nominal consistency checks, not certified quadrature accuracy or NS momentum residual bounds.

Next: partial axial-square/Z primitives inside incoming/pulse/end; consistent heat H and derivative primitives; complete pressure/Z; regenerate inherited incoming swirl from the installed angular primitive and update targets coherently. Enclose inputs and coefficient tangents, resolve nonzero terminal radial transport, then finite energy, stress/remainder and scale recursion. No full closure claim.

## 2026-09-30: actual inner-seeded angular moments extended to Rtail

Added continuous_angular_moments.py/md/json and installed ContinuousIncomingProfile.angular_moments_jet. It transports actual Mtheta, full swirl-square integral and both Z derivatives from Rh through every preheat stage to Rtail. Actual inner raw swirl is half the square integral; the API uses twice that atom. Reference baseline, signed compact-bump contributions and actual inner offsets remain separate; offsets are reapplied once. Constant stages use exact exponential transfers and finite stages declared Gauss nodes on the installed continuous schedule. Absolute checkpoint guards avoid rejecting recorded endpoints after subtraction of enormous log radii.

Actual field checks passed: Rh seed nominal relative agreement about 3.21e-291 (not an input accuracy certificate); four forward moment integrands agree with actual velocity/velocity-Z in flattening with maximum relative difference about 6.83e-14. Separately retained tiny bump primitives agree with point correction integrands to about 1.11e-17. Samples include Rh, flattening, angular bump and Rtail. Quadrature and inherited inner uncertainty are still unenclosed. Complete five moments, heat continuation, finite energy, stress and recursion remain open.

Next: combine these angular rows with the same live axial row-2 and quadratic primitives plus actual inner seeds, producing mixed and z_theta moments/Z jets. Continue angular/energy primitives through heat using consistent H values and derivatives, then enclose inputs and propagate coefficient tangents. Preserve nonzero terminal radial residuals; do not infer finite energy from these local checks.

## 2026-09-30: continuous angular correction and pressure Z jet installed

Completed continuous_angular_correction.py with one continuous beta definition for point values, coefficient equations, partial pressure and full bump energy atoms. The actual field now includes amplitude*h_Z after flattening; tiny h/h_Z/h_t are retained separately. Coefficient Z derivatives use the implicit two-equation Jacobian, keyed by full MP Z strings. Parameter parsing preserves declared precision; atoms use 100 working digits, not a claim of 100-digit accuracy. Pressure exposes the backward bump contribution and its Z derivative separately from baseline pressure. Existing correction/schedule identities are retained and stale caches cleared. The energy target is regenerated from the current live tail nominal contribution instead of preserving the previous future term.

Actual shared-field run passed: five angular installation samples agree with the provider and the first bump retains nonzero Utheta_Z. Bump pressure Z derivative fourth-order difference relative error is 2.329932752163308432839506125565243891769e-20; Rp mass and Z matching remain about 6.92e-83. Pressure and energy providers confirm shared bump atoms. Standalone bump derivative, parity, support and endpoint checks passed. These are numerical consistency results, not tiny pressure-target cancellation, full pressure/moment closure, finite energy or recursive matching.

Next: continue full angular/pressure moments and their Z jets using these atoms; bound the inherited preheat baseline and heat H contributions; enclose actual inner/swirl/future energy inputs; propagate complete input intervals through coefficient tangents; resolve the actual nonzero terminal radial tail before any finite-energy or scale-recursion claim. Keep unrelated scale_reference edits untouched.

## 2026-09-30: continuous angular schedule installed in the shared field

Installed continuous_angular_schedule.py in the existing schedule object. Full pre-heat log A, switch slopes, Z flattening and public MP Z entry points now use continuous MP primitives. Interior primitive working precision is 100 digits, arithmetic 443; accuracy is not inferred from either. Angular field adapters consume the returned jets. Heat normalization and interpolation retain epsilon in MP while H remains inherited. Incoming log Ep is regenerated before solving. Five transition Z differences at 1e-4 and 1e-30 pass; heat interface log jump is about 4.97e-292; Rp mass and Z matching still pass around 6.92e-83. Shared runtime identities and nonzero terminal residual retention pass.

Next tasks: implement continuous relative angular correction and its coefficient Z jets; update angular cumulative and pressure moment providers together; enclose H and its energy contributions; bound actual inner offsets and inherited swirl/future energy; propagate full coefficient tangent intervals and resolve radial tail. No full five-moment, finite-energy, stress or scale-recursion certification.

## 2026-09-30: continuous incoming angular primitive and second-row bounds

Completed continuous_incoming_angular.py and continuous_incoming_angular_enclosure.py with exact J(1)=1/2 and full incoming mixed-row bounds. At 4096 panels mixed relative width is about 5.858e-4. regenerate_incoming now computes the nominal mixed factor with MP J instead of float-backed schedule._log_A; the input changes about -9.57e-17 relatively, at 100 working digits/order96 (not certified accuracy). Both reference matching rows now pass their quadrature uncertainty into the coefficient report, preserving inner offsets once. Updated field Rp mass and Z matching pass at about 6.92e-83; shared runtime identity checks pass and both terminal residuals remain nonzero.

Next tasks: install continuous J and consistent angular jets throughout the complete exterior schedule; update angular cumulative/pressure/energy atoms together; enclose measured inner offsets, inherited swirl and future heat contribution; propagate all Z jets through coefficient tangents. Incoming angular bounds are for the ideal source, not a certificate of the legacy complete angular velocity. Finite energy, stress and scale recursion remain open.

## 2026-09-30: reference incoming axial intervals installed in coefficient report

Completed continuous_incoming_enclosure.py with full-support bounds and analytic reference Z derivatives. At Md=.5 and Z=.3, 4096-panel relative widths are about 2.221e-4 for Iz and 2.156e-4 for Iuz2; the installed values are contained. The incoming provider exposes full_incoming_enclosure. The coefficient report propagates mass uncertainty into base1 and squared-velocity uncertainty with the correct negative sign into the target, preserving actual stored inner offsets once. Nominal amplitude and coefficients remain contained. Energy-target relative interval width from this component alone is about 9.49e-43; coefficient width remains pulse-dominated at about 2.49 percent.

Conditional gaps remain: mixed angular row2, inherited swirl energy, inner offsets, future heat energy, normalization-source uncertainty and their Z jets. Reference axial derivative intervals are available but not yet propagated through coefficient tangents. Next tasks: enclose the mixed angular primitive and its factor; enclose inner moment offsets and future energy; transport each source separately into base/target intervals and tangent equations; tighten pulse range bounds. No global mean, finite-energy, stress or scale-recursion certification.

## 2026-09-30: pulse energy uncertainty now propagated

Completed continuous_pulse_energy_enclosure.py, JSON and proof notes. Startup uses monotonic switch bounds to enclose its nested primitive; plateau is analytic; cutoff uses outward positive rectangles. Entire support is covered. At 4096 panels Kp relative width is about 4.10e-9 and the installed nominal value is contained. The coefficient interval solver now consumes this Kp range; amplitude relative width is about 2.05e-9, coefficients about 2.49 percent. Incoming rows and target remain fixed conditional inputs. No coefficients or terminal residuals were overwritten.

Next tasks: bound actual continuous incoming rows and squared-velocity energy; separate mixed angular uncertainty and target uncertainty; pass those intervals into the coefficient solve. Tighten correlated coefficient numerators where needed, then add incoming/target Z-derivative enclosures. Basis, pulse rows and Kp have bounds; full mean, finite energy, stress and scale recursion remain unproved.

## 2026-09-30: conditional coefficient interval solve

Completed continuous_solve_enclosure.py with exact endpoint JSON and notes. It propagates basis and full pulse bounds through the factored determinant and the weighted energy quadratic. The live runtime exposes coefficient_enclosure. Conditional unique positive root and nominal containment pass at 1024 and 4096 panels; coefficient relative widths narrow from about 10.36 to 2.49 percent. Incoming rows, target and Kp are still fixed stored dyadic parameters, with originating uncertainty explicitly unbounded. Fixed-coefficient residual intervals containing zero establish compatibility only; nonzero materialized terminal residuals remain.

Next tasks: (1) Enclose continuous incoming cutoff mass and squared-velocity integrals using a scaled endpoint transformation and explicit positive tails. (2) Enclose pulse energy Kp rather than using refinement differences as bounds. (3) Feed actual input intervals into solve_atoms and preserve each uncertainty source. (4) Tighten correlated inverse numerators if row dependency loss dominates. (5) Carry base/target Z derivative intervals into the tangent solve before certifying physical radial tail closure. Finite energy, all five moments, stress and scale recursion remain open.

## 2026-09-30: complete continuous pulse interval bounds

Added experiments/root_st073/lei_ren_part1_paper_continuous_pulse_enclosure.py and its JSON/notes. Both weighted rows have outward bounds for the full functional pulse support, including three positive omitted pieces. At 4096 panels the centered relative interval width is about 1.89 percent. The live shared runtime exposes pulse_enclosure and its provider passes independent containment. Basis and pulse integration are now bounded for declared decimal parameters; incoming/angular uncertainty, pulse energy and complete coefficient closure remain open. Terminal mass and Z residuals are still retained and nonzero; finite energy and scale recursion are not certified.

Next agent tasks: tighten the pulse interval width with analytic curvature or adaptive range panels; preserve exact dyadic endpoints and positive omitted support; propagate correlated basis and pulse bounds into the coefficient solve without subtracting nearly equal matrix products; enclose incoming rows and pulse energy before claiming full mean closure. Keep user changes in scale_reference files.

## Interval enclosures for continuous basis atoms - 2026-09-30

New continuous_basis_enclosure.py/md/json uses outward interval arithmetic
and an analytic midpoint remainder for the actual continuous bump definition.
It encloses I0, both weighted B rows/matrix atoms, and squared-bump energy Gram.
Exact dyadic endpoints are serialized; displayed decimal digits are not the
certificate. Real declared mu/ell are covered, inherited parameter errors not.

For 4096 panels: I0 width 1.216e-6, B1 width 5.706e-6; determinant width 1.477e-33.
The correlated analytic determinant identity (with positive sinh bounds)
proves this candidate's real basis matrix invertible despite tiny mu.
Nominal I0/matrix/determinant/Gram all lie inside independent intervals.
1024-to4096 panel refinement decreases the analytic remainder by 16.

SharedContinuousAxialRuntime now exposes optional basis_enclosure; its fixture
checks the actual owned basis against these bounds. Complete pulse quadrature,
incoming moments, coefficient errors and source target remain unenclosed.
The broad full-axial/global-energy certificates stay FALSE, and materialized
terminal residuals remain visible. This bounds basis integrals, not the whole
NS field or exact terminal cancellation.

Next: enclose pulse centered quadrature and incoming primitives; propagate
correlated atom/input errors through functional coefficients; establish both
terminal constraints without replacing materialized residuals. Full outer
moments/pressure, global energy, recursive matching and stress/remainder remain.

## Shared live continuous axial atoms installed - 2026-09-30

New continuous_axial_runtime.py/md/json owns ONE live basis/pulse, complete
matrix and pulse rows, numeric p, bump energy atoms and solved a/c. The point
component, algebra and complete primitive share objects by identity. Partial
integrals retain the same basis; complete end atoms use the stored matrix.

The actual ContinuousIncomingProfile now uses this owner for point values,
means and coefficient Z tangents. No coefficients or full integral atoms are
reconstructed from JSON inside this installed path. Kp retains its inherited
100-digit atom precision and open accuracy bounds; algebra uses 200 digits.

Full-field reruns: Rp mass/Z-mass matching 6.92e-83. The actual primitive and
independent full-atom replay (including Z derivative) now differ by zero at
reported precision. The SHARED materialized residual remains nonzero:
normalized mass balance 8.50e-202; Z balance 1.12e-201. Global finite energy
and recursive closure are NOT certified. No mass or derivative is reset.

Next tasks:
- [x] Share live complete atoms between solve, values and cumulative means.
- [x] Differentiate coefficients with the same live algebra and actual inputs.
- [ ] Enclose complete and partial integral errors, coefficient sensitivity and
      signed pulse omitted pieces without conflating them with exact identities.
- [ ] Represent and establish BOTH terminal constraints for the exact continuous
      functional coefficients, preserving materialized residual diagnostics.
- [ ] Finish angular jets, all outer moments/pressure and physical energy;
      recursive matching, stress/remainder and oscillatory corrections remain.

## Direct remaining end integrals installed - 2026-09-30

ContinuousAxialBump.tail uses reflection of the same even bump to integrate
remaining mass directly. Flat endpoint integrals factor out the endpoint
maximum before quadrature, avoiding absolute-error stopping on tiny values.
At precision100, s=.1499 retains log(tail)=-763.2991417 while full-minus-partial
rounds to zero. This is a numerical representation improvement, not closure.

ContinuousAxialCorrection now caches and RETAINS its evaluated terminal
balance rho. End-region cumulative rows use rho minus the direct remaining
end integral. The installed mean/Z-mean jet path uses the same representation
and retains rho_Z; no terminal mass or derivative is replaced with zero.

Actual full-field reruns: Rp mass and Z-mass matching remain 6.92e-83.
Installed normalized terminal mass balance is 1.11e-200; Z balance6.29e-201;
physical radial tail remains nonzero. Incoming scalar serialization now keeps
its declared100-digit quadrature precision, with separate source/algebra and
mixed-integral precision provenance; extra printed digits are not accuracy.

Next tasks:
- [x] Direct end-tail primitives and retained-residual cumulative representation.
- [x] Install representation in actual outer mean AND Z-mean jets.
- [x] Preserve incoming atom precision in regenerated receipts.
- [ ] Share complete continuous atoms with functional coefficient definitions,
      keeping evaluated coefficient/atom residuals and quadrature bounds apart.
- [ ] Establish both terminal constraints under that complete integral model.
- [ ] Finish angular jets, outer five moments and pressure, global energy,
      recursive matching, stress/remainder and oscillatory corrections.

## Actual terminal balance and energy-tail obstruction - 2026-09-30

New continuous_terminal_balance.py/md/json evaluates the INSTALLED continuous
incoming/exterior and actual input/coefficient Z jets, not a formal quadrature
graph. Installed terminal mass/sumabs balance is 1.73e-200; Z mass balance is
6.62e-201. Both remain nonzero. Independent full-atom replay differs slightly
from the cumulative primitive; that summation discrepancy is retained.

Recovered physical radial tail: ur=-nu*C(Z)/r, with
C=((1-delta)*Z*Mz+(1-Z^2)*Mz_Z)/(1-delta*Z^2).
The radial kinetic energy density per dZ dlogr is pi*nu^2*C^2*dz/dZ.
At Z=.3 the nominal C/Rh is nonzero, logabs about 1.1504127226124e28.
Thus increasing ordinary coefficient precision does not establish finite
GLOBAL energy: any surviving continuous 1/r tail makes that integral diverge.
Pulse omitted-piece bounds are transported separately and remain conditional
on nominal coefficients; other quadrature/input/coefficient errors are open.

Next: tie complete continuous integral atoms, coefficient constraints and
partial primitives through a cancellation-preserving definition; close BOTH
terminal mass and its Z derivative while preserving actual incoming seeds.
Do not assign mass zero, add an arbitrary radial mask, or use a favorable
alternate term replay as physical-field closure evidence. Angular primitive,
all five outer moments, stress/remainder and recursive matching remain open.

## Shared continuous incoming axial field installed - 2026-09-30

New continuous_incoming.py/md/json provides the same MP smooth cutoff for
Uz, Uz_y, Uz_Z and weighted cumulative mass/axial-energy integrals. It reuses
the pulse switch and evaluates the reflected cutoff to retain positive flat
tails. Incoming mean_Z is analytic; quadrature accuracy is not enclosed.

New continuous_incoming_outer.py/md/json installs that provider before Rp,
regenerates I_z, I_theta_z and I_uz2, reapplies measured inner offsets ONCE,
updates the prior-energy target and re-solves with retained continuous atoms.
Both incoming point values and mass/Z-mass now feed the joined field.
Rp mass and mass_Z relative matching: 6.92e-83 each. New linear row replay:
1.27e-173 / 1.24e-173; energy replay 1.63e-201 (numerical algebra only).

A float-rounded amplitude derivative produced a 2.26e-17 derivative jump.
Using the exact MP derivative of (1+Z^2)^-1 removes that discrepancy.
Zero-Z tangent factors are taken from the regenerated incoming receipt.

Remaining limits: the mixed integral retains the float-backed angular
primitive, swirl energy is inherited, future angular target derivative is
still finite-differenced, and integral/coefficient uncertainty is unenclosed.
No terminal mean is assigned zero. Global finite energy, full five moments,
NS stress/remainder and scale recursion remain incomplete.

Next actionable tasks:
- [x] Install continuous incoming axial values and matching mass/Z-mass.
- [x] Regenerate axial incoming rows/target and install updated coefficients.
- [ ] Replace retained angular primitive and cache/jet precision coherently.
- [ ] Resolve terminal mean AND its Z derivative with explicit integral bounds.
- [ ] Continue all five moments and pressure jets through outer/heat layers.
- [ ] Establish finite global energy and cross-scale matching before claiming
      recursive background closure; then stress/remainder and pulse corrections.

## Actual axial coefficient and cumulative mean Z derivatives - 2026-09-30

New continuous_axial_tangents.py/md/json differentiates the same continuous
linear rows and quadratic energy constraint. Independent declared-direction
fourth-order replay refines 2.19e-17 -> 1.37e-18; actual nominal input jets
replay derivative rows around 1e-201. Arithmetic replay is not a certificate.

New seeded_input_tangents.py/md/json transports actual corrected Rh momentsZ
and reference/source-factor derivatives into normalized incoming rows/target.
It does not recover tiny reference inputs by subtracting a large seed. Swirl
reference energy derivative cancels analytically after Ep normalization.
Unresolved inner entries [3,5], subtraction and quadrature uncertainty remain.
Future angular-energy derivative is explicitly float-backed fourth-order FD;
its relative two-step change is about 3e-12, not an error enclosure.

New continuous_seeded_outer_jets.py/md/json installs post-Rp Uz_Z and mean_Z
from those input/coefficient jets. Local radial divergence refines 6.24e-24
-> 3.90e-25 without neighboring full coefficient solves. Pulse integrals are
separate from cumulative means; each seed/amplitude is counted once.
Schedule/angular/tail/pressure identities and nonzero terminal rows remain.
Incoming before Rp retains a declared finite-difference fallback; angular
velocity Z jets remain incomplete. Signed pulse omission bounds are separate
from incoming/coefficient/quadrature uncertainty.

Next: analytic angular coefficient/heat target derivatives; shared continuous
incoming primitive and Z cache precision; actual mean and mean_Z closure;
outer five moments/jets and global energy; recursive matching, stress/remainder
and oscillatory correction. Global finite energy/scale recursion are not done.

## Continuous axial correction installed in joined exterior - 2026-09-30

New continuous_seeded_outer.py/md/json adapter installs continuous Uz AND
its cumulative seeded mean together. Original schedule/angular/tail/pressure
objects are preserved. Seven fixed-Z positions cover Rp, startup, plateau, cutoff,
first end bump, Rv and after-Rv. Post-Rv constant mass replay is 2.10e-443.
Actual Ur is recovered from the same mean. Independent radial divergence
relative cancellation refines 6.24e-24 -> 3.90e-25 as step halves.
These are local diagnostics; float-Z derivatives remain uncertified.

Startup weighted primitives now use endpoint scaling and integration by
parts, retaining the shared smooth sigma definition. Truncating a SUBTRACTED
sigma integral produces a negative correction: absolute omission bounds and
sign are propagated by the component and installed mean. Actual-mu startup
xi=.015 bound transport passes. Numerical loss of a positive primitive raises;
no nonzero terminal mean is replaced by zero. Quadrature and inherited input
uncertainty are not covered by omitted-piece bounds.

Next derive shared continuous incoming and Z derivatives (avoid float cache
loss), certify actual terminal mean/Z-mean closure, and finish exterior jets,
five moments and energy. The installed finite numerical tail remains nonzero;
full finite energy, recursive matching, stress/remainder and oscillatory
closure are still open. Do not equate tiny local divergence with NS closure.

## Continuous cutoff primitive and shared cumulative component - 2026-09-30

Completed cutoff partial integral branches for xi in (10,11): endpoint,
saddle and after-window evaluation, with positive omitted-piece bounds.
Four actual-candidate positions run; independent endpoint derivative replay
refines 3.19e-12 -> 7.98e-13 when step is halved. Startup (0,.02] remains
unsupported explicitly. Quadrature error is not enclosed by omitted bounds.

ContinuousAxialCorrection now uses the solved pulse/end coefficients for
point values and normalized cumulative rows, with the SAME continuous atoms.
Actual terminal nominal values remain nonzero; no mean/tail reset occurs.
This component is not yet installed in the global physical field.

Injection audit: preserve original schedule/angular/tail/pressure identities;
carry the entire seeded incoming receipt, not just linear base rows. Its
row_normalization is seeded but dimensionless_integrals.I_z remains unseeded.
Pre-Rp mass seed must be added once. Convert stored full-row pulse weight to
current-point weight by exp[lambda*(13-xi)/mu]. Replace velocity and primitive
providers together. Details in continuous_axial_solve.md.

Next complete startup weighted primitive and shared incoming/Z provenance;
install consistent axial values/means, then establish terminal mean AND its
Z derivative closure. Global finite energy and scale recursion remain open.

## Continuous pulse energy and re-solved axial correction - 2026-09-30

Completed: continuous_pulse_energy.py/md/json and continuous_axial_solve.py/md/json
under experiments/root_st073. Kp uses the same continuous pulse definition,
MP startup/cutoff nodes, and analytic plateau. Same-order 70/100 precision
refinement is 7.54e-72; order80/112 refinement is 2.76e-21 (not an enclosure).
Fixed an import-time decimal-boundary precision leak. Old float Kp differs
by 9.84e-15 relative. Working digits do not equal certified integral accuracy.

Re-solved actual seeded base rows/energy target with continuous bump, pulse
and energy atoms. Old coefficients have 6.68e-15 relative continuous-row
defects; new row arithmetic replay is about 1e-173. Energy algebra replay is
about 1e-201. Consistent solved end-bump value/derivative/weighted primitive
provider is available. These coefficients are not installed globally yet.
Continuous pulse partial integrals now cover the analytic plateau in addition
to the saddle window, retaining an explicit positive startup bound.

Next: complete startup/off-window cutoff and incoming primitives; derive
shared Z derivatives; install continuous coefficients in both point velocity
and cumulative means together. Establish actual mean/Z-mean closure before
any tail removal. Then exterior energy, outer jets/five moments, recursive
transitions, stress/remainder and oscillatory correction. Full goal stays open.

## Continuous integral providers and bounded core energy - 2026-09-30

New artifacts: lei_ren_part1_paper_continuous_axial_basis.py/md/json,
lei_ren_part1_paper_continuous_axial_pulse.py/md/json,
lei_ren_part1_paper_core_energy.py/md/json, under experiments/root_st073.
The continuous bump shares values, derivatives and full/partial primitives.
Its matrix differs from the old float-node canonical matrix by 2.21e-19.
The MP pulse shares pointwise values and full weighted rows; saddle-window
partial weighted primitives now include positive omitted-piece bounds and
an independent derivative diagnostic. Off-window partial evaluation raises
an explicit error and remains to be implemented. Pulse full-row log inputs
change by 6.68e-15 from the float-centered evaluator. No global installation.

Actual bounded core energy integration confirms declining kinetic energy
under the selected shrinking time scales, despite growing velocity amplitude.
This is core-only; the exterior nonzero radial tail remains the global energy
obstruction. Do not infer recursive scale closure or full NS residual closure.

Next tasks in order:
1. Complete off-window continuous pulse partial primitives with bounds.
2. Add continuous pulse energy and shared incoming primitive/Z derivatives.
3. Re-solve with the same continuous bump/pulse/incoming atoms used by velocity.
4. Establish actual terminal mean AND its Z derivative cancellation before
   modifying any 1/r tail; keep nonzero arithmetic tails visible.
5. Complete outer jets/five moments, annular/exterior energy and measured
   core-width profiles; then uniform stress/remainder and oscillatory layers.
The full goal remains active; local core energy is not global finite energy.

## Exact quadrature dependency graph and measured multi-time core geometry - 2026-09-30

Run lei_ren_part1_paper_axial_closure_graph.py for actual shared candidate
inputs. Canonical matrix and row atoms are now exported by the axial solver.
Exact rational Cramer expressions eliminate both quadrature row terms;
their formal Z derivatives also cancel for a fixed matrix and differentiable
shared atoms. This is a quadrature-model identity, not continuous closure.
Independent basis perturbations produce nonzero residual terms, not zero.
Rounded graph coefficient/energy replays remain separate. The physical
nonzero tail is unchanged; continuous integral/Z-derivative provenance and
the exact energy root still need work. Read closure_graph.md for the next
shared continuous/conservative basis requirements.

Run lei_ren_part1_paper_core_scale_geometry.py. It uses actual regular-core
physical callable roundtrips at logq=-alpha/delta, alpha=0,2,4,6. Measured
aspect gains are 1,e,e^2,e^3; swirl/axial amplitudes, local winding density
and axial vorticity component grow with the mapped scale laws. Log-aspect
errors are about 1e-244 and relative velocity roundtrip errors about 1e-247.
These extreme MP chart times are not ordinary simulation times. One core
point diagnoses geometry/coordinate scaling; it is not a measured vortex
core width, an integrated streamline, full vorticity, global energy, or a
completed recursive transition. Both Python/JSON/Markdown artifacts exist.

Next connect a shared continuous basis for values, partial/full primitives,
and Z jets; prove the actual mean and its derivative close before removing
the 1/r tail. Independently integrate energy over shrinking core/annuli and
exterior, measure core-width profiles across time, and finish outer jets and
five moments. Uniform stress/remainder and oscillatory correction remain
open. Do not treat formal row cancellation as a finite-energy certificate.

## Installed seeded axial field; mean conditioning exposed - 2026-09-30

Read lei_ren_part1_paper_seeded_outer_field.py/md/json and the regenerated
seeded_shared_candidate.json. New coefficients are now installed in both
axial velocity and the actual cumulative mean. Schedule/angular/pressure
identities are retained. Direct seeded transport avoids subtract/add loss.
Fixed Decimal stage-offset loss in axial_incoming._row, MP bump normalization,
and stored full pulse row reuse. Prior joined_outer.json is labelled historical.
At Z=.3 Rp relative mean jump=3.541e-260; Rv reported mean jump=0 at precision;
pulse independent divergence relative cancellation=3.899e-25. These are local
numerical diagnostics, not full momentum residuals or uniform certificates.

IMPORTANT: tail coefficient ratio is -2.590e-183, but log(|tail/Rh|) remains
1.1504e28. The numerical tail is still nonzero and the implemented radial
1/r energy obstruction remains. Increasing ordinary MP precision is not the
next solution: derive an exact closure identity for the same actual field,
with continuous integral/conservative-model provenance and uncertainty.
Do not set a nonzero terminal mean to zero or hide an arithmetic residual.
Then restore corrected outer jets, full five moments, energy/support tests,
actual multi-time scale diagnostics and stress/remainder. Oscillatory layers
and full residual <=1e-3 remain later requirements.

## Shared inner-through-heat callable and actual seeded exterior solve - 2026-09-30

Run lei_ren_part1_paper_joined_outer.py. JoinedOuterField keeps the SAME source
profile and the actual five-bump corrected Rh mass/pressure offsets. It now
provides a physical velocity callable from the regular core through heat.
The Z=.3 heat-chart roundtrip passed at working precision. This is a
provisional continuation: sampled nonzero radial transport still prevents
claiming finite energy. Float-Z derivative errors, missing corrected outer
velocity jets and full five-moment exterior stress remain explicit.

Run lei_ren_part1_paper_seeded_axial_inputs.py --shared to reproduce the actual
Rh-seeded Section 7.31/7.34 coefficient solve (seeded_shared_candidate.json).
Canonical linear row replay is 5.924e-416 and 1.945e-416; energy replay is
2.870e-444. These are finite-quadrature coefficient equations, NOT full NS
residuals or a continuous mean closure proof. New coefficients are not yet
installed in the joined velocity/primitive. The default command is a clearly
labelled historical fixture regression, not the actual shared candidate.

Next implement a seeded axial pulse/primitive provider together: before Rp,
retain old outer cumulative mass plus actual Rh offset; from Rp onward use
updated incoming m1 and the newly solved pulse/end-bump coefficients. Avoid
double-counting the inner offset in the joined radial transport. Recompute
Z jets from the same actual primitive, update incoming integral provenance,
measure terminal mean by replay rather than assign zero, and quantify
quadrature/precision refinements. Preserve shared swirl and axis pressure.
Then resolve finite energy and axial support before claiming temporal scale
recursion. Global norms, admissible stress/remainder and oscillatory layer
remain open.

## Actual inner seed transport and serialization precision - 2026-09-30

New lei_ren_part1_paper_seeded_axial_inputs.py maps explicit raw inner offsets
(z, theta_z, z_theta) into the two normalized Section 7.31 linear rows and
the Section 7.34 prior-energy input. Use the same schedule log_Rp and log_Ep.
No absent offset defaults to zero. This is preparation for a new outer
coefficient solve, not a completed repair or a replacement axial primitive.
The copied dimensionless_integrals remain historical reference inputs;
update the actual primitive consistently before using a seeded receipt in
CorrectedSourceProfile. A dimensional replay at logRp=5e152+17 passed.

from_signed_log now preserves arbitrary_exponent_value when available,
validates its sign, and retains the legacy log-only fallback. The stored
precision receipt shows relative error 2.24e-161 versus 2.23e-8 for a
log-only roundtrip. The existing actual-tail solve/replay at Z=.5 passes
(linear differences below 5.80e-133, energy difference 6.85e-160).

Next: complete the same-object inner/outer join, measure terminal transport,
then feed actual offsets into the pulse/end-bump coefficient solve and replay
the same primitive. Do not zero the residual mean. Finite energy, uniform
source constants, temporal scale recursion, and oscillatory residual remain
unestablished.

## Five-bump numerical correction and physical inner callable - 2026-09-30

Read experiments/root_st073/lei_ren_part1_paper_inner_corrected_field.py/md/json
and inner_moment_map.py/json. Equation (10.8) now solves all five representative
nonlinear moment equations with analytic coefficient Z derivatives. Actual
partial moments, pressure with unchanged P0, and continuity-derived Ur are
connected to the shared core-to-Rh candidate. Unknown defect entries3/5 use
explicit midpoint representatives and conditional value/Z intervals; they
are not claimed zero. Exact normalized bump mass identities are retained.
At Z=.3 eleven center/flank probes pass the relaxed cone. Higher quadrature
replay of the terminal representative equations gives normalized residuals
around 1e-33 or lower, distinct from the much smaller algebraic solve residual.
The local mapped-divergence flank replay gives q div(u) about 1.57e-24.
CorrectedInnerField.physical_field().velocity(x,y,z,t) now covers the actual
regular core through Rh and raises outside the constructed domain.
Next audit terminal fields/moments over Z, resolve or enclose input uncertainty,
then append the supplied corrected outer/heat profile with actual mean
identities. Tiny nonzero mass tails cannot imply finite energy. Uniform
source constants, inner admissible collar, global moment closure/energy,
temporal scale recursion and oscillatory-corrected residual remain open.
Detailed ordered tasks are in inner_corrected_field.md.

## Axial reference restored; actual defect inputs recorded - 2026-09-30

Read experiments/root_st073/lei_ren_part1_paper_axial_restore.py/md/json.
Section 9.38 transports actual moments and Z jets from Rz through Rh.
At Z=.3, V reaches 4Z=1.2, V_Z reaches 4, and V_y vanishes at the end.
Five relaxed-cone probes pass. Largest 32/64 moment quadrature difference
is 1.51e-32. Independent radial mapped-divergence replay at restoration
midpoint gives q div(u)=2.14e-24, relative cancellation 3.26e-25.
Centered defects d1,d2,d4 are about 2.234e-15,5.690e-16,5.915e-41.
Angular/pressure entries d3,d5 are unresolved by subtraction: nonzero
roundoff artifacts exceed their conditional physical bounds. Do not set
them to zero or use them as measured defects. The receipt preserves
conditional bounds and flags the missing uniform C1/e_star certificate.
Next implement the fixed Section10.8 five-bump matrix and nonlinear
moment correction, with stable signed-log/interval treatment of tiny
defects and actual partial bump moments. Detailed checklist is in the
axial_restore.md. Source requires relaxed cones in these later intervals;
strict admissibility is needed in the inner collar. Finite energy,
heat exterior matching and temporal scale recursion remain unfinished.

## Long reshape now reaches the axial restoration entry - 2026-09-30

Read experiments/root_st073/lei_ren_part1_paper_long_reshape.py/md/json.
Section 9.30 is connected to actual shared-candidate R=110 moments and
Z jets. Endpoint-normalized MP quadrature reaches Rsh and analytic
reference-power primitives continue to Rz=exp(-8)Rref. Five Z=.3 probes
pass the relaxed cone; full admissible stress remains unestablished.
Largest 16/32 normalized quadrature difference is 6.57e-25. The angular
reference log-value difference at Rsh/Rz is about 5.59e-263. Inherited
inner moments are preserved, and conditional value/Z tail bounds are
recorded separately from uncertified quadrature error.
Next implement axial restoration on [Rz,e Rz] with actual five moments
and Z jets, continue to Rh, compute actual Section 10 repair defects,
and implement moment repair. Complex A_Omega, mixed core A/K, pressure
Z-tail bounds, heat exterior, finite energy, admissible stress and temporal
scale recursion remain unfinished. Provisional source budgets are unchanged.

## New shared candidate reaches R=110 - 2026-09-30

Read experiments/root_st073/lei_ren_part1_paper_shared_candidate_1.md/json.
The reproducible Python script recomputes pressure and degree 18 core at
j=1e-14, Lambda=1e36, logCstar=5e151, logPstar=14, delta=1e-200.
All five necessary input gates and six sampled core signs pass. The four
new-parameter connection probes pass the relaxed cone; stricter admissible
cones are not established. Explicit phases resolve distinct switch shears
although the physical radius offsets round away. R=110 has a=.8,b=0.
Real-axis norm bounds are available; mixed core A, full K and complex
A_Omega bounds remain open. A=1e150 and logK_upper=1e152 are provisional.
Next implement long angular reshaping with actual moments/Z jets while
bounding these source constants, then axial restoration, moment repair
and heat exterior matching. Finite energy, temporal scale recursion and
full oscillatory-corrected residual remain unfinished. The detailed ordered
checklist is in shared_candidate_1.md; mark tasks complete only with receipts.

## Short source switches now reach R=110 - 2026-09-30

ExitSwitches implements both exact short switch intervals, actual moment
and Z-jet transport, and analytic a=.8,b=0 continuation to R=110. At Z=.3,
R=100 field/moment/jet/slope matching is 1.60e-160 relative. All six sampled
relaxed cones pass. The 16/32 switch endpoint difference is 7.03e-12; the
largest interior difference is 2.44e-9. The recorded b=0 cone margin gives
D(R=110) about 5.68e37, above 3. Read exit_switches.md/json under
experiments/root_st073. These remain development-fixture diagnostics.

Next priority is the shared parameter rebuild, not blindly extending
this fixture to Rsh. Its necessary parameter gate explicitly fails.
Follow connection_scale_gate.md for selecting j, controlling A/K,
choosing Cstar/Rref/delta and a common source collar, and recomputing
pressure/core. Then replay the reusable exit/switch modules, reshape
angular velocity, restore axial velocity and repair actual five moments.
Full matching, finite energy and temporal scale recursion remain open.

## Shared source parameters must be rebuilt - 2026-09-30

The current development fixture fails necessary Section9 input conditions:
Cstar/A, long angular-shape radius, axial-radius separation, j and short
collar width. Read connection_scale_gate.md/json under experiments/root_st073.
This is an explicit necessary-condition rejection, not just missing proof.
Real endpoint max(G) bounds A; it is not the complex A_Omega bound.

build_source_core now accepts shared j/logC/logPstar/delta and recomputes
outer schedule, pressure anchor, axis and core together. Its nondefault
wiring receipt passes local checks, but does not certify a replacement.
Do not reuse old fixture receipts after changing parameters. Next select
one compatible shared candidate with controlled A/K, then rerun pressure,
core, switches and long reshape. Existing local source modules remain
useful algorithm implementations; full matching and recursion remain open.

## Source exit now reaches R=100 - 2026-09-30

The conditional actual exit propagates five moments toR100; six finite
samples pass relaxed cone tests. Comparison extends to R=110. Next implement
the source short shear switches and satisfy the shared parameter gates
before reference angular/axial restoration and Section10 repair. Read
exit_continuation.md and connection_next.md under experiments/root_st073.
Full matching, finite energy and temporal recursion remain open.

## Actual three-component initial exit - 2026-09-30

Actual Z-jet transport, continuity-based Ur and Section 3 stress are
implemented with a physical core/exit callable. The one independent
divergence sample has relative error 3.95e-20. Midpoint admissible and
endpoint relaxed cone samples pass; near-join interior cone remains unverified.
See exit_field.md/json under experiments/root_st073. These local results
leave full matching, finite energy and temporal recursion open.

## Source exit update - 2026-09-30

The same high-degree core now feeds Section 9.23 comparison and actual
Section 9.25 initial exit. Five moment increments are jointly accumulated.
See experiments/root_st073/lei_ren_part1_paper_exit_bridge.md for y<=.01
scope. Next transport actual Z derivatives, recover Ur, evaluate actual
stress/cone, extend to the reference region and repair moment defects.
The baseline future-pressure Z envelope leaves angular-correction
coefficient derivatives open. Temporal scale recursion remains open.

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

## Compact solenoidal localization control rejected

compact_control.py constructs velocity by tapering the actual poloidal
streamfunction (degree9 C4 cutoff) and swirl, not multiplying velocity alone.
It is bounded and compactly supported at every registered time, so finite
energy follows; analytic solenoidality follows from the streamfunction.
Only X<=1/64, |eta|<=.3 is preserved (NOT the full frozen axial domain).
Pressure is explicitly tapered. No residual-canceling force is added.
Independent FD at radial/axial/corner points gives late k5.5 momentum norms
1.104e5,4.045e4,2.111e5, stable after halving steps. Divergence FD defects
shrink about16x, consistent with fourth-order truncation; they are not blanket
numerical divergence acceptance. This localization is rejected as an NS field.
The failure quantifies the missing transition dynamics. Do not repeat simple
cutoff sweeps or present compactness as matching. Next solve a stress/flux-driven
transition correction or a decaying exterior with evolving interface data.
Artifact: experiments/root_st073/compact_control/report.json.

## Updated source priority

Read docs/ST073_PAPER_ROUTE.md before new construction. Official PDF was reopened;
actual heat-exterior direct-join screen rejects value-only matching.
Evidence: experiments/root_st073/heat_join/screen.json.

## C2 radial swirl component constructed

swirl_bridge.py uses actual corrected coefficient jets at X=3/64 and the frozen
heat-exterior amplitude to build an autonomous quintic on r in[ri,2ri]. Values,
first and second physical radial derivatives match at both ends. At nine k/eta
probes endpoint absolute errors are <=4.22e-15,3.09e-12,3.05e-9 respectively;
101 radial samples per bridge retain positive swirl. This removes the tangential
viscous traction jump at those interfaces without fitting force. The polynomial
is OUR interpolation, not a formula claimed from the paper. It is only a swirl
component, not the complete transition: no poloidal velocity, pressure, dynamic
residual, axial caps or finite-global-energy acceptance is supplied by this step.
Next construct compatible poloidal streamfunction transport and derive the full
transition residual stress; do not interpret endpoint matching as NS acceptance.
Artifact: experiments/root_st073/swirl_bridge/coefficients.json.

## Poloidal streamfunction bridge constructed

poloidal_bridge.py builds a septic psi(r,z,t) matching corrected inner radial
jets through order3, and zero jets at outer radius2ri. Combined u_r=-psi_z/r,
u_z=psi_r/r is solenoidal when moving endpoints are differentiated consistently.
This is an autonomous interpolation, not a paper-derived dynamic solution.
Its zero outer streamfunction forces return transport: at k6,eta=.4 core axial
flux5.22367e-5 is balanced by annular flux-5.22367e-5. Sampled annular axial
velocity ranges -2.3740 to2.1197; return flow cannot be ignored in stress fitting.
Endpoint residuals are recorded; full z/time derivatives, pressure joining,
complete momentum and axial closure remain pending. No arbitrary force added.
Next evaluate the combined swirl/streamfunction bridge and integrate its actual
stress defect, then test the paper-inspired moment/covariance correction.
Artifact: experiments/root_st073/poloidal_bridge/coefficients.json.

## Callable full radial joining field

joined_field.py now evaluates all velocity components and pressure across the
inner field, Hermite transition and heat swirl exterior inside |eta|<=.5.
Poloidal psi_z analytically differentiates coefficient jets AND moving endpoints;
pressure matches value and first two radial derivatives. No forcing is fitted.
Independent Cartesian FD at four transition points, two times and two steps:
early momentum norms38.88--62.31, late7851.6--12562.9, stable under refinement.
Divergence errors decrease about16x on step halving (late finest<=1.054e-6).
This is a callable divergence-structured background with a large dynamic defect,
NOT a completed transition, axial closure, finite-energy global field or theorem.
Next compute actual integrated tangential residual stresses of THIS background
and test realizability/correction following paper Sections4,7-9; do not transfer
leading-only cone values or call smooth endpoint joins a momentum solution.
Artifact: experiments/root_st073/joined_field/report.json.

## Actual transition residual moment targets

transition_stress.py integrates the full FD tangential residual at fixed physical
z,t. With zero inner stress, the radial inverses give nonzero terminal stresses:
at late k5.5 rtheta is about-1.52 to-1.66 and rz ranges-2.24 to2.19 across three
axial locations. Thus this restricted inverse cannot vanish at both interfaces
without moment repair. Orders8/12 quadrature are recorded, not a continuum bound.
The sampled maximum reaches6.05e4 on the denser radial nodes; previous four-point
checks were not maxima over the whole transition. Both records are retained.
The pointwise PSD trace lower bound2*sqrt(T_rtheta²+T_rz²) is3.32--5.42, only an
algebraic bound. PSD completion also contributes diagonal/radial momentum and
is NOT a wave realization. No arbitrary forcing or post-hoc stress cancellation
has been counted as a solution. Next solve these measured moment defects with
boundary-jet-preserving background corrections before wave/cone admission.
Artifact: experiments/root_st073/transition_stress/moments.json.


## Boundary-preserving swirl moment repair: rejected as a PDE improvement

Implemented `swirl_moment_repair.py`: an axisymmetric pure-swirl bubble with
cubic zeros at both interfaces, preserving velocity and first/second jets.
Pressure, poloidal velocity and frozen core stay fixed. This is an autonomous
basis choice inspired by the moment-repair route, not a formula from the paper.
Two bounded coefficients fit at k=3, eta=-.3,0,.3 are approximately
[2.94846774, 1.25961520]. Independent holdouts use k=.7,5.5 and eta=+/-.2;
training quadrature uses 8 points, reporting uses 12. No force is introduced.
At late holdouts terminal angular stress magnitude falls from about1.60 to
.0105--.0123 (over99% reduction), but full sampled momentum grows from
4.36e4--4.55e4 to7.02e4--7.08e4. Axial moment remains unchanged.
Thus integral compatibility alone is not a residual reduction. Do not promote
this optional candidate to the default JoinedField or count it as PDE progress.
The final held-out point is also evaluated with half the spatial FD step.
Artifact: experiments/root_st073/swirl_moment_repair/report.json.

Next implementation tasks:
- Replace moment-only fitting with a constrained angular PDE collocation solve:
  use multiple radial cubic-endpoint bubble modes and axial modes, preserve
  interface jets, and minimize the full theta residual while bounding moment
  defects. Freeze core data and retain all rejected baselines.
- Include multiple training scales and disjoint scale/axial holdouts; one frozen
  two-parameter fit does not remove scale-dependent moment defects exactly.
- Couple poloidal streamfunction modes to axial moment and axial PDE repair;
  swirl alone cannot change the current axial residual at fixed poloidal field.
- Recompute radial pressure compatibility after velocity corrections. Angular
  residual reduction alone does not control centrifugal/radial momentum.
- Only then evaluate actual correction-stress compatibility and realizability;
  retain the separate axial closure, finite-energy and scale-recursion work.
Full momentum max and physical-volume L2 gates remain1e-3; none has passed.


## Angular PDE collocation: component gain, no full-field gain

`angular_collocation.py` adds 12 endpoint-preserving swirl modes (four radial
Legendre modes, three axial polynomial modes). It trains at k=1,4 and
eta=-.3,0,.3, using normalized pointwise theta residual plus a soft moment
penalty, bounded coefficients and mild regularization. Base field evaluations
are cached without changing numerical differentiation or the fitted fields.
Disjoint holdouts k=2.5,5.5, eta=+/-.2 use18 radial nodes. Theta maxima fall
about13%, but full momentum maxima rise about0.35--0.38%; terminal angular
stress magnitudes also rise about0.3%. Candidate is NOT promoted to baseline.
At late eta=.2, full max45530.7485 agrees with spatial-step-halved45530.7480.
The coarse divergence estimate8.78e-4 drops to1.79e-6 on refinement; stencils
near the C2 join affect coarse finite-difference errors. No supremum or volume
L2 certificate is claimed. Axial/poloidal correction is now the priority.
Artifact: experiments/root_st073/angular_collocation/report.json.


## Coupled poloidal / pressure collocation

`poloidal_collocation.py` fits six streamfunction and six pressure bubbles to
FULL nonlinear momentum at k=1,4 and eta=-.3,0,.3. Quartic streamfunction
endpoint zeros preserve velocity and all spatial jets through order2; physical
psi_z includes moving-interface and source-to-physical derivative factors.
The compact correction adds zero net axial flux through each annular slice.
It is analytically divergence-free; finite-difference divergence is separately
reported. Pressure bubbles have cubic endpoint zeros. Frozen core, original
swirl and heat exterior are retained, and no force is fitted.

`affine_momentum.py` precomputes affine velocity/derivative jets; evaluating
advection retains the exact quadratic coefficient interactions. The resulting
surrogate differs from a fresh full-field Cartesian FD evaluation by at most
2.50e-8 at the checked training slice. This is a discretization consistency
check, not independent proof of the PDE. The bounded nonlinear fit converges
in10 evaluations. These are autonomous basis/optimizer choices, not a direct
implementation of the paper's oscillatory stress realization.

At disjoint k=2.5,5.5 and eta=+/-.2 holdouts, full sampled maxima decrease
about8.8--9.1%. Late eta=.2 falls45373.40 ->41259.86; halving spatial FD step
gives41259.8573, divergence2.96e-6. Training decreases12--29%.
Angular moment compatibility slightly worsens; this is an improvement candidate,
NOT an accepted NS solution. Finite global energy, axial closure, scale
recursion, non-axisymmetric stress realization and global gates remain open.
Artifacts: experiments/root_st073/poloidal_collocation/report.json;
physical-volume subdomain audit: poloidal_collocation/volume_audit.json.

Derivative convention: q_z and eta_z from source_coordinates must be divided
by sqrt(nu). For y=r/ri-1, y_z=-(1+y)q_z/(2q). For increasing physical time,
q_t=-q_tau and eta_t=-eta_tau; the full FD evaluator already implements this
minus sign. Do not confuse remaining time tau with physical time t.


Physical-volume audit over eta in[-.4,.4], ri<r<2ri, all azimuths, uses
12 radial by8 axial Gauss nodes and the exact coordinate volume Jacobian.
At k2.5, sampled maximum2989.51 ->2606.71 and volume L2 5.48547 ->5.12251.
At k5.5, sampled maximum68136.50 ->59387.68 and volume L2 26.3981 ->24.6533.
These are approximately12.8% maximum and6.6% L2 improvements on this subdomain;
BOTH remain far above1e-3. Quadrature convergence and continuum bounds are
not established. These denser axial checks exceed the earlier slice maxima,
which must not be reported as whole-transition maxima.
Next combine poloidal, pressure and angular modes in one full-momentum solve,
with actual moment penalties and scale-dependent modes; preserve disjoint
holdouts and volume metrics before considering any candidate promotion.


## Joint full-momentum and moment fit: tradeoff, not a new accepted baseline

`joint_collocation.py` composes the12 poloidal/pressure and12 swirl modes,
retaining nonlinear cross interactions in the affine-jet momentum evaluator.
It uses an analytic optimizer Jacobian, bounded coefficients, six training
slices at k=1,4 / eta=-.3,0,.3, and soft weight3 penalties on normalized
r-squared angular and r-weighted axial residual moments. Core and interface
jets stay fixed; no external force is fitted. Source files include direct
full-field FD holdouts and the physical-volume audit for reproducibility.

At disjoint k=2.5,5.5 / eta=+/-.2, angular terminal stress magnitude improves
about16--17% against the original JoinedField; full maxima improve about9%.
Against the previous poloidal/pressure candidate, the additional slice maximum
improvement is only about0.07--0.10%. Axial moment changes are mixed by parity.
On the volume subdomain at k5.5, joint max59244.28 is below prior59387.68,
but L2 24.7024 is ABOVE prior24.6533. Thus the joint candidate does not dominate
the previous candidate on the user's two gates. Preserve both, promote neither.
Original baseline on the same quadrature: max68136.50, L2 26.3981.
All remain far above1e-3 and no finite-energy whole-space field exists yet.
Artifact: experiments/root_st073/joint_collocation/report.json.

Next avoid endless fixed-coefficient tuning: introduce explicit scale-dependent
streamfunction/swirl/pressure amplitudes, differentiating those amplitudes in
space and physical time, then test a cross-scale solve. Inspect actual moment
transport and outer matching freedoms before imposing a radial stress inverse.
If gains stay small, advance the paper-inspired stress realization route rather
than equating increasingly complex Hermite bubbles with dynamical matching.


## Explicit log-time amplitude comparison

`scale_collocation.py` extends24 joint coefficients to a0+s*a1, where
s=(k-3)/3 and k=-log2(2*tau). The amplitude depends only on time, so it
preserves the spatial divergence-free streamfunction/swirl construction and
all endpoint jets. The optimizer explicitly includes s_t*delta_u in momentum,
s_t=1/(3*ln(2)*tau), for increasing physical time. Product-rule jets are
compared against full-field FD with varying amplitudes at every time stencil.

A24-coefficient constant model and48-coefficient varying model are fitted on
IDENTICAL k=1,3,5 and eta=-.3,0,.3 data, with the same pointwise/moment objective,
bounds and regularizer. Initializing the varying model from the fitted constant
model avoids attributing extra training data to amplitude variation. Intercepts
and slopes are individually bounded[-4,4]; effective coefficients can exceed
that interval, which is reported explicitly. Holdouts remain k=2.5,5.5 and
eta=+/-.2, with a separate physical-volume subdomain audit. This is a finite
log-time ansatz, NOT a closed scale-recursion law or singularity construction.
Artifact: experiments/root_st073/scale_collocation/report.json.


Matched-data outcome: at k5.5, varying coefficients slightly worsen subdomain
max59244.12 ->59266.70 while slightly reducing L2 24.7044 ->24.6954.
At k2.5, max2600.444 ->2600.241 while L2 5.133039 ->5.133393 worsens.
No consistent improvement on both gates; do not promote this variant or
continue increasing temporal polynomial degree without a new mechanism.
The next route is the actual shear/phase/amplitude dynamics of paper Section7,
with the approximation limits recorded in ST073_PAPER_ROUTE.md. Global energy,
axial closure, recursive scale control and the1e-3 gates remain unfulfilled.


## Frozen phase/amplitude implementation

`frozen_pulse.py` extracts physical angular velocity F=u_theta/r and radial
tangential shear g=(d_r u_theta-F,d_r u_z) from the actual JoinedField and
poloidal/pressure candidate. It implements a frozen cylindrical local amplitude
ODE with the pressure projection needed to preserve n dot a=0, and viscosity.
The frozen phase has constant tangential wavevector and n_r_dot=-g dot n_tan.
This is an autonomous physical-coordinate approximation inspired by paper
Section7, not substitution of full ST073 into the leading normalized theorem.
Radial strain, axial gradients and background variation along trajectories are
omitted; angular mode1 is not a high-frequency justification.

Viscous damping is factored exactly as exp(-nu*(|n0|^2*t+(n0 dot n_dot)*t^2+
|n_dot|^2*t^3/3)) before integrating the undamped matrix ODE, avoiding relative
transversality errors once damped solutions fall below absolute solver tolerance.
Each sampled time uses the matrix singular value for optimal amplification;
the covariance uses the initial seed that maximizes the sampled peak gain.
Nonnegative covariance fitting is recorded separately and is NOT PDE admission.
Artifact: experiments/root_st073/frozen_pulse/report.json.


Screen result:156/192 sampled points have positive inviscid frozen growth
parameter, but only8/192 pass the diagnostic stress-direction inequalities.
The late candidate point with largest growth rate has the WRONG target-dot-N
sign. Its45 phase choices (angular modes1,2,4; five axial ratios; three radial
ratios) achieve peak optimal gain1.49859 over0.1*tau with viscosity retained.
The initial restricted radial-direction family missed this transient growth;
no absence-of-growth claim is retained. Nonnegative time-averaged covariance
fits the two-component target algebraically, but this cannot certify compact
pulse construction, mean cancellation, or full residual reduction.

Next integrate rays/amplitudes at points that actually satisfy the directional
screen, compare with full local velocity-gradient evolution, and construct a
supported vector-potential field only after checking its own full residual.
Keep compact moment tails and global axial/energy closure as separate failures.


## Full-gradient rays and actual collar residence

`affine_pulse.py` uses the full physical Cartesian velocity gradient J at the
two late poloidal/pressure points passing the previous directional screen.
It solves n_dot=-J^T*n and the incompressible projected amplitude equation,
retaining viscosity through an integrated scalar damping factor. The45 initial
wavevectors per point include angular-like modes1,4,16,64,256 and axial/radial
ratios-1,0,1. None of these sampled frozen-affine evolutions amplifies above its
initial norm over0.1*tau. Transversality relative errors stay below1.84e-9.
This is not a no-growth theorem for the changing background or other phases.
The average covariance fits algebraically with nonnegative weights, again
showing why an algebraic fit alone is not a dynamical construction.

The points sit only1.788e-5 physical length from the moving inner interface.
The dimensional cutoff rate nu/distance^2 is3.128e7, about3.57--3.59e5 times
the reduced inviscid growth rate; this is a scaling warning, not an operator
bound. A0.1*tau frozen horizon also predicts radial travel about30 times that
distance, making that local approximation unsuitable for a supported pulse.

`ray_residence.py` therefore integrates particles in the ACTUAL time-dependent
candidate. Both cross the moving inner interface in about3.84e-5 time units,
only0.00348*tau, far shorter than the assumed0.1*tau pulse interval. The event
uses q(z,tau), not a fixed cylindrical boundary. The y=.02 outer probe was
explicitly defined and is not claimed to be a cone boundary; neither path
reached it. Artifacts: affine_pulse/report.json and affine_pulse/residence.json.

Next redesign the transition to admit an interior region with appropriate
stress direction and transport residence BEFORE adding localized oscillations.
Do not attach the earlier short-time1.50-gain wave to this narrow collar: it
came from a different point that fails the direction screen. Preserve the
core, moment accounting, full residual gates and global-closure requirements.


## Variable-width transition construction

`width_field.py` generalizes the radial join to ro=ratio*ri while preserving
the same ST073 inner velocity/pressure jets and fixed heat amplitude. It
rebuilds all Hermite polynomials using width=(ratio-1)*ri. For
Y=(r-ri)/width, Y_z=-(1+(ratio-1)*Y)*ri_z/width; width_z/width=ri_z/ri.
These terms are retained in analytic psi_z, so changing width does not replace
the solenoidal streamfunction construction with a velocity cutoff.
The ratio2 evaluator is compared directly with the original JoinedField.

`width_screen.py` compares ratios2,3,4,6 at k=2.5,5.5 using16 radial and6 axial
Gauss nodes over eta in[-.4,.4]. It records physical-volume L2 on each ACTUAL
annulus, full sampled maxima, divergence, velocity magnitude and the same
approximate stress-direction screen. Larger annuli have larger volumes;
no volume-normalized score or smaller evaluation region hides that cost.
No fitting is performed, so this is an exploratory geometry comparison.
Artifact: experiments/root_st073/width_screen/report.json.


Width comparison result: at late k5.5, ratio2 -> ratio6 lowers sampled max
66741.68 ->4514.92 (93.2%) and physical-volume L2 26.3981 ->5.7862 (78.1%),
INCLUDING the larger volume. At k2.5 max2928.42 ->200.61 and L2 5.48547 ->1.21796.
Ratio6 is an exploratory candidate, not accepted by either1e-3 gate.
The directional screen passes14/96 nodes instead of4/96; some new passes are
near the OUTER edge, where compact-stress moment tails are still unresolved.
No continuous interval of admissibility or pulse residence is certified.

Use `WidthField(6)` / width_screen/candidate.json for the next bounded trial.
Do NOT blindly reuse old AngularModes or PoloidalModes: their normalized
radius and derivative formulas assume width=ri. Generalize them to width=5ri,
then recompute actual moments and transport residence with the new background.
The fixed heat amplitude, original core and original ratio2 implementation
remain available for matched comparisons. Larger width alone does not solve
axial closure, finite global energy, pulse realization or recursive control.


## Width-aware coupled correction basis

`wide_modes.py` generalizes the existing24 swirl, streamfunction and pressure
modes using y=(r-ri)/((ratio-1)*ri). Streamfunction radial derivatives divide
by the ACTUAL width, and y_z=-(1+(ratio-1)*y)*q_z/(2*q*(ratio-1)) uses physical
q_z. Endpoint zeros and analytic divergence are preserved. Original narrow
modules remain unchanged. Nonzero-coefficient ratio2 compatibility differs
by at most3.1e-16 in velocity and zero in pressure on the recorded four points.

`wide_collocation.py` fits the generalized basis on WidthField(6), starting
from zero coefficients rather than transferring the narrow-layer optimum.
Full nonlinear momentum and actual tangential moments enter the objective;
training k=1,4 / eta=-.3,0,.3 and holdouts k=2.5,5.5 / eta=+/-.2 stay distinct.
All residual, moment and volume quadrature radii now cover ri<r<6ri; using
old ri<r<2ri sample helpers would incorrectly omit most of the new transition.
Artifacts: wide_collocation/report.json and wide_collocation/compatibility.json.


Wide coupled outcome: at late k5.5 on MATCHED12x8 volume nodes, sampled max
4629.00 ->4261.34 (7.94%) and physical-volume L2 5.78633 ->4.89297 (15.44%).
At k2.5, max205.326 ->189.200 and L2 1.21799 ->1.02850. Earlier16x6 width
screen maxima differ due to node placement; compare paired values, not maxima
from different grids. Disjoint axial slice checks show about2% full maximum
reduction and42% angular terminal-stress reduction, with modest axial-moment
improvement. Neither moment is zero, and neither1e-3 global gate is passed.
This is a useful wide-background candidate, not an accepted global NS field.
Next recompute the directional screen and actual transport residence for THIS
corrected wider field; earlier narrow-field cone or trajectory results do not
transfer. Keep explicit stress tails, axial closure and finite-energy tasks.


## Corrected wide-background pulse geometry

`wide_pulse_screen.py` recomputes full residual prefix stresses, local shear and
the diagnostic direction conditions on the fitted Width6 field. It uses20
radial nodes at k=2.5,5.5 and eta=-.3,0,.3; the selected late passing point
maximizes physical distance to either radial interface. No old narrow-field
stress direction or trajectory is transferred. Actual particle transport uses
the changing background and moving boundaries, stopping at a radial boundary,
|eta|=.48 probe, or0.1*tau. The initial directional condition is NOT asserted
along the trajectory. A separate frozen full-gradient amplitude calculation
uses the observed residence duration, with27 explicitly recorded initial
wavevectors. This remains preparatory geometry/dynamics, not a supported
non-axisymmetric field or a full momentum cancellation.
Artifact: experiments/root_st073/wide_pulse_screen/report.json.


Wide pulse screen result:14/120 nodes pass the approximate direction test.
The selected point lies at eta=0,y=.956117 (near the outer boundary), with
physical distance7.0616e-4 to that boundary. Its actual particle reaches the
moving outer interface after0.0718045*tau, about20 times the previously tested
narrow-collar residence, although these are different selected positions.
All27 sampled frozen-affine perturbations have peak gain<=1 up to rounding.
This does not certify absence of growth along the changing trajectory.
The target radial stress remains nonzero near the outer boundary, so moment
tails and compact support are still unresolved even where direction passes.
Next evolve phase/amplitude against J(t) along the actual trajectory and
re-evaluate stress-direction compatibility there before constructing a
supported vector-potential wave. Do not treat longer residence as PDE progress
or as a substitute for reducing the remaining full residual and global energy.


## Evolving phase along the corrected wide trajectory

`evolving_pulse.py` evaluates the actual Cartesian velocity gradient at25 saved
particle positions, using fourth-order spatial differences, then integrates
the Kelvin wavevector and pressure-projected viscous amplitude with a smooth
interpolated J(t). It compares27 initial wavevectors with a fixed-start J on the
SAME observed residence time. In both models the sampled peak gain is1.0;
the best evolving final gain is0.97759. The evolving gradient trace is at most
1.22e-7 and wavevector/amplitude orthogonality errors stay below5.23e-11.
This rules out gain for this explicitly sampled local family and horizon only;
it is not an impossibility statement for all waves or a global pulse.

At the seven late points passing the approximate direction test, a mode-one
azimuthal wave has minimum local damping nu/r^2 larger than the reduced
inviscid growth rate: the ratios range about1.30--12.1. This dimensional
comparison explains the observed damping but is not a spectral bound for the
changing, inhomogeneous field. The positive-direction region still occurs
only in narrow inner/outer collars. The next constructive step should change
mean swirl/shear/pressure in the wide transition to create an interior cone
margin with growth that overcomes viscosity, while retaining the full residual
and volume scores. Then re-solve transported phase and stress tails.
Artifact: experiments/root_st073/evolving_pulse/report.json.


## Interior cone/viscosity fit on width6: rejected

`cone_fit.py` uses the exact nonlinear momentum residual of the24-mode wide
field and radial-prefix stress from it. At interior radial Gauss nodes, soft
penalties favor the approximate stress direction and inviscid growth exceeding
the mode-one local viscosity rate nu/r^2. It keeps the existing residual and
moment terms, training k=1,4,5.5 and eta=-.3,0,.3. The diagnostic geometry is
in PHYSICAL cylindrical coordinates and must not be equated with the paper's
normalized leading-cone theorem.

This fit creates three interior passing nodes at late eta=0, whereas the prior
had zero, but creates none at eta=+/-.3. On held-out late physical-volume nodes,
maximum momentum grows4261.34 ->6218.60 and L2 grows4.89297 ->5.39810;
angular terminal mismatch grows from about.307 to1.60--1.63 at eta=+/-.2.
It is rejected as a full-field candidate. The simple cone penalty cannot
replace a radial/axial profile satisfying BOTH stress direction and moment
identities. The current poloidal/pressure basis has axial modes1 and eta only;
add an even eta^2 mode before assessing whether the side failures are structural.
Artifact: experiments/root_st073/cone_fit/report.json.


## Even axial mode extension for the wide background

`even_modes.py` adds eta^2 terms to both streamfunction and pressure radial
bubbles, giving30 coefficients in all (18 poloidal/pressure plus12 swirl).
The physical z derivative of eta^m includes m*(eta/.3)^(m-1)*eta_z/.3;
this is required for exact streamfunction divergence cancellation. With the
new six coefficients zero, the extended evaluator agrees exactly with the
old24-mode evaluator at three nonzero-coefficient sample points.

`even_cone_fit.py` repeats the bounded full-momentum/moment and approximate
interior direction/growth fit with that expanded parity basis. It starts from
the previously accepted width6 fit, with all new eta^2 coefficients zero.
The aim is to check whether a symmetric axial correction can supply interior
pulse conditions at BOTH eta=+.3 and eta=-.3 without sacrificing full momentum
and physical-volume L2. Fitting penalties remain diagnostics, not acceptance.
Artifacts: even_cone_fit/report.json and even_cone_fit/compatibility.json.


Even-mode result: the new eta^2 terms leave the former field EXACTLY unchanged
when set to zero. The cone-penalized fit increases late midplane interior
passing nodes from0 to4/10, but STILL gives0/10 at eta=+/-.3. Its held-out
late slice maxima improve about5%, yet the late matched-volume maximum grows
4261.34 ->4965.24 and physical-volume L2 grows4.89297 ->5.70224.
Angular terminal mismatch at eta=+/-.2 grows about.307 ->1.72--1.85.
Therefore it is rejected. A low-degree axial polynomial did not solve the
side stress orientation and made the user's actual gates worse.

Next use the paper's radial shear modulation idea as a separate controlled
candidate: alter local shear while keeping mean velocities and radial moments
nearly fixed, then compute the FULL viscous residual and pulse compatibility.
High radial frequency can greatly increase viscosity; test this cost explicitly
before claiming any stress cone progress. Do not replace the width6 best
full-momentum candidate with either cone-penalized field.


## Radial shear loop and coupled mean repair: rejected as full-field candidates

`radial_shear_loop.py` adds a compact two/four-cycle axisymmetric swirl loop
with quartic endpoint zeros, keeping exact divergence freedom and C2 interface
jets. At k5.5, a two-cycle amplitude-.3 creates one approximate direction+mode1
growth passing node at y=.729 on BOTH eta=+/-.3 slices; the existing width6
mean field had none there. The late side sampled maximum grows only about0.2%,
and angular terminal stress falls about22%. A four-cycle -.3 loop creates
passing nodes too, but its momentum maximum grows markedly from viscosity.
These are finite-node diagnostics, not a wave or a strict paper cone certificate.

The two-cycle candidate fails the actual physical-volume gate comparison:
late max4261.34 ->4294.04 and L2 4.89297 ->5.90882; early L2
1.02850 ->1.23678. `loop_repair_fit.py` refits all24 mean coefficients in the
presence of the same loop. It recovers part of the L2 cost (late5.45039) but
remains worse than the width6 baseline, and the two side passing nodes disappear
on direct re-evaluation. Thus shear orientation, moment compatibility and full
PDE residual must be solved together. Neither raw nor refitted loop replaces
the accepted best sampled baseline; no supported non-axisymmetric pulse exists.
Artifacts: radial_shear_loop/report.json, candidate_audit.json, and
loop_repair_fit/report.json.

Next work should follow the paper's linked moment-preserving profile continuation
and oscillatory stress program. A sinusoidal shear loop by itself is too weak a
substitute. Preserve failed coefficient sets and full-volume metrics for
comparison; do not count approximate cone passes as an NS residual reduction.

## Spatially compact solenoidal closure: finite-slab construction only

`compact_potential.py` puts the current width6 field into an axisymmetric
vector potential: A_theta=psi/r and A_z=-integral_0^r u_theta(s,z,t)ds.
The new velocity is curl(chi A), where chi is a smooth physical-r/physical-z
cutoff. On the cutoff plateau it replays the entire width6 velocity and pressure;
outside the finite support it vanishes. Thus it is nonzero, spatially compact
and divergence-free by construction for every registered time, and its spatial
kinetic energy is finite at each such time. This is only defined for
tau in [.5/64,.5]; neither uniform energy control nor continuation to tau=0
is established. The radial joins inherited from the width6 field are C2.

`compact_potential_audit.py` checks two times. Four plateau samples match the
width6 field exactly. A coarse 6x6 physical cylindrical quadrature gives
kinetic-energy estimates 6.64e-5 at k=.4 and 1.18e-5 at k=5.5; these are
neither convergence checks nor uniform bounds. Finite differences at four
cutoff-collar points give momentum-residual norms 343--1674 at k=.4 and
5.73e4--3.22e5 at k=5.5. Sampled divergence is <=2.85e-6, consistent with
the curl identity and numerical differentiation. These collar defects exceed
the 1e-3 target by many orders and are not global max/L2 measurements.
The pressure is explicitly tapered with chi; no cancelling force was added.
The outer closure therefore solves the energy/support property on the finite
slab but worsens the full NS momentum gate. Next solve the cutoff-induced
momentum and pressure balance together with the mean/pulse stress construction,
then audit physical-volume norms and critical-time scaling. Report:
`experiments/root_st073/compact_potential/report.json`.

## Weighted residual moments at the compact closure

`compact_moment_audit.py` integrates the *full physical* unforced momentum
residual at k=5.5 on two fixed-z slices. The radial pieces are split at the
inner interface, width6 outer interface, radial cutoff plateau, and support
edge. These are diagnostics analogous to the weighted radial primitives in
paper Section 8, not that section's normalized auxiliary-mean defects.

At z=0, the 12-node-per-piece moments are integral r^2 R_theta dr =
+0.00143895 and integral r R_z dr = +0.000913628. The outer radial cutoff
alone contributes +0.00158829 to the angular moment; the width6 transition
contributes -0.000149333. The inner core is below 2e-6 pointwise in the
sampled residual, so the outer cutoff dominates this slice's angular defect.

At z=0.00394501, inside the axial cutoff collar, the total axial moment is
-0.525289. The width6 transition contributes -0.553372 and the core adds
+0.0264002. This is a distinct axial closure problem, not just a radial
edge-layer defect. Six-versus-twelve-node moment changes are <=2.4e-5 for
these four totals; this is limited quadrature evidence, not a norm bound.

Section 8 requires the relevant weighted source moments to be repaired before
a compact radial stress primitive can cancel the mean residual. A direct
primitive of the current defects would leak beyond its annulus. Next construct
moment-preserving velocity and pressure corrections, and an actual realizable
non-axisymmetric stress if needed, while keeping full momentum max/L2 as the
acceptance gates. Artifact: `experiments/root_st073/compact_potential/moment_audit.json`.

## Conservative balance of the two weighted defects

`compact_balance_audit.py` evaluates the exact axisymmetric conservative
identities behind the preceding residual moments. For the compact field,
Mtheta=int r^2 utheta dr and Fz=int r uz dr (the latter is analytically zero
because uz=r^-1 d_r(chi psi)). The angular balance is
int r^2 Rtheta dr = d_t Mtheta + d_z int r^2 uz utheta dr
                   - nu d_zz Mtheta.
The axial balance is
int r Rz dr = d_t Fz + d_z int r(uz^2+p)dr - nu d_zz Fz.
All derivatives use physical time t=T-tau and physical z.

At k=5.5, z=0, angular time, axial transport and axial viscosity contribute
+0.00107881, -0.000052400 and +0.000411138. Their sum 0.00143755 differs
from direct integrated angular residual by 1.40e-6. The axial kinetic and
pressure flux derivatives contribute +0.000728821 and +0.000184808; their sum
agrees with the direct axial moment to 2.23e-10.

In the axial cutoff collar z=0.00394501, the axial kinetic-flux derivative is
-0.0867909 and the axial pressure-flux derivative is -0.438498. Their sum is
-0.525289109, within 3.28e-9 of the independently integrated axial residual.
The pressure term supplies ~83.5 percent of this signed moment. Fz is ~1e-18
numerically, so its time/viscous contributions vanish at this resolution.
This identifies the main scalar closure defect as a pressure/kinetic axial flux
imbalance, not merely a large pointwise derivative of the radial cutoff.

A moment-only pressure adjustment could cancel the integrated axial defect,
but may create a large radial pressure-gradient error. The next candidate must
reconstruct pressure jointly with radial momentum and velocity/stress fluxes;
it is admissible only if full pointwise and physical-volume L2 momentum improve.
Artifact: `experiments/root_st073/compact_potential/balance_audit.json`.

## Pressure-only axial-moment repair: constructed and rejected

`pressure_moment_repair.py` adds a compact pressure bubble outside the preserved
inner core. At each fixed (z,tau), it subtracts the full integrated axial flux
H=int r(uz^2+p)dr times a radial beta bubble normalized to int r b dr=1.
At strength1 this makes H exactly zero up to quadrature; velocity and analytic
divergence are unchanged. This is an actual correction candidate, not a fitted
force, and it directly tests whether the large axial moment could be the only
remaining obstruction.

At k=5.5, axial collar z=0.00394501, the repaired H is 6.78e-20 and direct
int r Rz dr is -3.76e-9, down from -0.525289. Yet on the 5x5 physical-volume
sample grid over the compact support, full momentum maximum moves
189000.81 -> 189001.45 and sampled L2 moves 209.05091 -> 209.10133.
A least-squares strength -3.05459 lowers that grid L2 only to 208.98472 and
moves opposite to the moment correction. Neither is an accepted field. These
coarse grid numbers are candidate screens, not full-domain maxima or certified
volume norms. The baseline worst residual and the angular component remain
orders above1e-3.

Thus cancelling one weighted pressure/kinetic flux identity is insufficient.
The next correction must jointly handle radial and angular momentum, the
axial pressure profile and velocity/stress covariance, with the full physical
residual retained after each step as in the paper's Sections8-9. Artifact:
`experiments/root_st073/compact_potential/pressure_repair.json`.

## Joint axial-collar velocity/pressure candidate and scale screen

`joint_collar_fit.py` adds six compact axisymmetric modes in the axial collar:
two streamfunction modes (hence exactly solenoidal velocity), two regular
swirl modes, and two pressure modes, each with even/odd axial parity. All vanish
on the preserved z plateau and outside the compact support. A linearized full
Cartesian momentum fit at k=5.5/Gauss5 supplies coefficients; the subsequent
screen retains the exact quadratic velocity term. Strength1 overfits: training
L2 falls209.05 ->141.29, while independent Gauss6 L2 rises314.77 ->1432.30.
At strength.1, training L2 is197.95 and independent Gauss6 L2 is215.36,
with independent max181939 ->91679. This is a bounded candidate, not a new
accepted NS field. Another independent Gauss8 late grid gives max321680 ->
210291 and L2 383.97 ->288.38. The strong grid-order dependence prevents any
claim of converged physical-volume L2.

The constant coefficient set fails at earlier k=.4 (Gauss6 max946 ->2436,
L2 22.56 ->58.82). `joint_collar_scale.py` therefore tests the SAME fitted
coefficients with amplitude multiplier (tau_ref/tau)^1.5 and strength.1. Full
time derivatives of that multiplier are included; k6 uses fourth-order forward
time differences at the registered lower-tau endpoint. On independent Gauss6
samples, baseline -> candidate (max; L2) is:
  k=.4: 946.16 ->932.34; 22.56 ->22.29
  k=4: 38740 ->31923; 144.97 ->125.12
  k=5.5: 181939 ->91715; 314.77 ->215.33
  k=6: 304685 ->127480; 407.60 ->262.11.
This is genuine sampled improvement across the registered finite window, but
still misses both1e-3 momentum gates by many orders and establishes nothing
as tau tends to zero. The angular maximum remains unchanged at these grid
nodes; axial maxima sometimes increase. The sampled quadrature is not a norm
certificate. `ScaledJointCollarField` is a callable full velocity/pressure
candidate. An independent fourth-order Cartesian check at three late collar
points found momentum norms124396 ->122894, 95717 ->95622, and55512 ->55501;
divergence remained at the ~3e-6 finite-difference level.

Next fit across several times and radial/axial quadrature orders, target the
remaining angular channel and interface/cap stresses, then evaluate full
physical-volume gates. Paper Sections7-9 still require a realizable
non-axisymmetric pulse/covariance and iterative mean repair; these six slow
axisymmetric modes do not implement that program. Artifacts:
`experiments/root_st073/compact_potential/joint_collar_fit.json`,
`joint_collar_scale.json`, and `joint_collar_direct.json`.

## Radial-collar angular correction added to the joint field

The remaining angular peak at k6/Gauss8 lies at r=0.02154 and |z|=0.002264,
inside the radial cutoff collar but below the axial cutoff onset. Thus the
previous axial-collar modes are exactly zero there. `radial_swirl_fit.py` adds
two C4 compact pure-swirl radial bubbles, with the same (tau_ref/tau)^1.5 time
weight and the existing smooth axial taper. Their velocity is analytically
solenoidal. Fitting their linearized FULL momentum response at k6/Gauss8 and
retaining the quadratic self-transport gives coefficients7.12829 and6.07674.
`RadialSwirlRepairedField` exposes the complete sum as a callable field.

At fitted k6/Gauss8, angular max drops179634 ->29815 and sampled physical-
volume L2 drops342.08 ->194.37. Total maximum remains276850 because a radial
peak outside these bubbles dominates. On an independent k6/Gauss6 grid, total
max falls127480 ->94315 and L2 262.11 ->191.48. The same coefficients also
lower sampled L2 at k5.5/Gauss6 from215.33 to161.13, at k4 from125.12 to
115.55, and at k=.4 from22.29 to22.17; total maxima at these earlier grids
are unchanged because other components control them. These are finite samples
with strong grid-order dependence, not converged maximum or volume norms.

`radial_swirl_direct.py` independently evaluates the full fourth-order
Cartesian momentum at six physical collar points at k5.5. At the midplane
collar center, the norm falls111121 ->8486; the largest post-correction norm
among these six points is11773, versus111927 before. Sampled divergence stays
around1e-10. This establishes a real local angular improvement, while leaving
large radial/axial errors elsewhere and the 1e-3 target far unmet.

The next work is a coupled correction for the *radial* maximum and interface
stress, together with non-axisymmetric realizability/moment constraints from
the paper. Do not promote this finite-slab screen to an accepted NS solution.
Artifacts: `experiments/root_st073/compact_potential/radial_swirl_fit.json`
and `radial_swirl_direct.json`.

## Late radial-peak incremental joint repair

After the radial-swirl correction, the largest k6/Gauss8 point is at
(r,z)=(0.005689,0.003433), inside the axial cutoff collar. Its residual is
mostly radial (+272322), with angular -29815 and axial -39970. This explains
why the radial-collar swirl modes cannot change the total maximum there.
`radial_repair_refit.py` reuses the six compact streamfunction/swirl/pressure
modes as an *increment* with a steeper (tau_ref/tau)^2 weight. The fit uses
full Cartesian momentum at k6/Gauss8 and retains nonlinear self-transport;
the fitted strength1 overfits the independent k6/Gauss6 grid (max94315 ->
107841; sampled L2 191.48 ->218.15). It is rejected.

Strength.1 is the conservative screen: at k6/Gauss8, total max276850 ->266160
and sampled L2 194.37 ->188.34; independent k6/Gauss6 max94315 ->91544,
L2 191.48 ->190.13. At k5.5/Gauss6 max91715 ->88740 and L2 161.13 ->
158.02. At k4 and k=.4 both sampled measures also decrease slightly.
`IncrementalRadialRepairField` exposes the combined compact field with this
increment; analytic divergence freedom follows from its streamfunction and
swirl construction. `radial_repair_direct.py` independently checks six late
peak-neighborhood points using fourth-order Cartesian derivatives and a
one-sided time stencil at the registered endpoint: the worst of these norms
falls276850 ->266160, with finite-difference divergence around1e-6.

The candidate remains massively above1e-3, the quadrature is not converged,
and no critical-time extension or non-axisymmetric stress is supplied. More
small slow-mode fits will not substitute for the paper's realizable oscillatory
stress and exact moment correction cycle in Sections7-9. Preserve this
candidate as a numerical benchmark, not an accepted Navier-Stokes solution.
Artifacts: `experiments/root_st073/compact_potential/radial_repair_refit.json`
and `radial_repair_direct.json`.

## Exact-curl two-harmonic wave: local covariance passes, full momentum fails

At the current k6 radial peak (r,z)=(0.00568908,0.00343263),
`radial_peak_cone.py` computes the actual corrected velocity gradient and the
physical local stress primitive needed for the tangential residual:
(sigma_rtheta,sigma_rz)=(45.5950,138.9583). This primitive is integrated only
from the axis to this radius; nonzero global weighted moments prevent its
unmodified compact extension. The physical analogue of paper Eq(7.1) gives
lambda^2=3919.85, target dot N=-145.02, and cone ratio0.0673<1. Frozen
Cartesian Kelvin rays give two positive time-averaged covariance weights.
These facts are LOCAL diagnostics, not the paper's normalized leading-cone
hypotheses or supported amplitude solution.

`curl_wave_prototype.py` periodicizes the two selected modes to distinct
integer angular harmonics m=1,2, projects their polarizations transverse to
the corrected normals, solves positive covariance weights23.1621 and75.8107,
and defines an actual compact C4 vector potential. Its analytic curl is
exactly solenoidal. At the selected point, 16-angle sampled covariance of the
exact velocity matches the target within2.3e-15 relative error; different m
values remove cross terms under angular averaging. This completes a concrete
local stress realization, but not the wave amplitude/pressure PDE of paper
Section7 or its mean correction cycle in Sections8-9.

The decisive full momentum test rejects this frozen wave. At target amplitude,
the local 16-angle maximum rises266160 ->42781680 and the angular-mean radial
residual rises262085 ->11826714. `curl_wave_amplitude_screen.py` retains exact
linear and quadratic residual terms: among 0 and101 positive amplitudes from
1e-5 to1, zero has the lowest sampled momentum RMS. At amplitude0.001 the
covariance is only1e-6 of target, yet maximum residual already rises to281524.
The finite-difference divergence of the target-amplitude wave is ~2e-4 while
analytic divergence is zero; this reflects differentiation of a steep compact
wave, not a certified numerical divergence bound.

The axial support halfwidth is0.00075, so its inverse scale is1333, versus
carrier magnitudes239 and419. The dimensional cutoff diffusion scale
nu/h_z^2 is17778, versus the local physical analogue growth rate62.6. These
are only scale diagnostics; the paper's normalized estimates cannot be
replaced by this comparison. They explain why freezing a locally favorable
covariance into a narrow compact pulse has a large omitted curl/viscous cost.
Next derive evolving periodic phases and solve the coupled amplitude-pressure
equation before any wave is added to the benchmark; preserve exact radial
moments and re-evaluate the complete residual. Artifacts:
`experiments/root_st073/compact_potential/radial_peak_cone.json`,
`curl_wave_prototype.json`, and `curl_wave_amplitude.json`.

## Separate the explicit wave's nonzero and mean errors

`curl_wave_amplitude_screen.py` now keeps the exact affine expansion of full
momentum for base + a*wave. At the late 16-angle sample, angular means of the
LINEAR wave residual are below1.2e-5 in all cylindrical components (as
expected from distinct nonzero harmonics). Its pointwise linear maximum is
1.6386e7. The QUADRATIC wave residual has angular mean
(1.15646e7,3.14050e6,-1.92154e4) in (r,theta,z). Thus the enormous radial
mean in the target-amplitude wave is quadratic self-interaction, not merely
a poor mean of the linear phase equation. The nonzero harmonic linear error is
also enormous, so both the paper's amplitude-pressure equation and its mean
correction cycle are independently needed.

`curl_wave_mean_pressure.py` adds a compact axisymmetric pressure bubble with
radial derivative -a^2*1.15646e7 at the selected point. This cancels the
wave-induced mean radial residual LOCALLY: at a=1, its angular mean falls
1.18267e7 ->2.62085e5, back near the original background value. The full
16-angle maximum is still3.14580e7 and RMS1.94002e7 (versus baseline
2.66160e5); at a=.1 the maximum is1.94510e6. The azimuthal mean and nonzero
harmonics remain. This pressure is a local proxy, not the exact compact radial
primitive, and no wave amplitude or harmonic pressure equation was solved.
Do not add either wave to the candidate field. Artifacts:
`experiments/root_st073/compact_potential/curl_wave_amplitude.json` and
`curl_wave_mean_pressure.json`.

## Harmonic pressure and local amplitude time slope: local gain, spatial failure

`curl_wave_harmonic_pressure.py` Fourier-decomposes the exact-curl wave's FULL
linear momentum residual at the selected k6 point. For integer angular modes
m=1,2, only14.9% and11.1% of the respective harmonic norms lie along the
local phase normal. Compact pressure harmonics remove those projected parts.
At target wave amplitude, full local maximum changes42781680 ->42719654;
with the earlier local mean-pressure proxy it changes31457954 ->31342811.
Pressure alone cannot cancel the mostly transverse linear error.

`curl_wave_taylor_amplitude.py` adds an exact-curl time-slope vector potential
whose local Fourier curl amplitudes replay the transverse targets to relative
errors below3.2e-16. At the registered endpoint, its added velocity vanishes
but its physical-time derivative cancels the targeted local linear harmonic.
With wave amplitude.1 and both pressure corrections, the 16-angle maximum
falls1963009 ->571391 and RMS1066886 ->299568. Yet the unperturbed field's
local maximum and RMS are both266160. A dense amplitude scan over[-.3,.3]
finds its best maximum at ZERO amplitude. Best RMS occurs at amplitude-.0205,
266160 ->266094, but the maximum worsens to278588 and covariance is only
0.00042025 of the stress target. The polynomial residual reconstruction agrees
with direct full-field evaluation at that amplitude within3.9e-6.

The local Taylor correction is especially fragile in space.
`curl_wave_taylor_neighbor.py` evaluates complete momentum at eight angles on
five meridional positions. At amplitude.1 the center maximum is571391; moving
radially by +/-0.0005 gives about3.58e6/3.61e6, and moving axially by
+/-0.0002 gives about1.29e7/1.35e7. Each exceeds its no-wave baseline by
large factors. These are sample maxima, not global norms, but they reject the
local Taylor wave decisively. Exact curl and a single-point amplitude derivative
are insufficient. A supported solution of the coupled phase/amplitude/pressure
PDE across the whole patch, followed by the paper's mean moment correction,
is required before revisiting the full-field acceptance gates. Artifacts:
`experiments/root_st073/compact_potential/curl_wave_harmonic_pressure.json`,
`curl_wave_taylor_amplitude.json`, and `curl_wave_taylor_neighbor.json`.
