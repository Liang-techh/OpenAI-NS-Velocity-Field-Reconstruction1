# Current angular terminal closure (2026-10-05)

The current exact-repair heat view now closes the angular terminal constant as a function on the whole axial domain. Its runtime directly consumes the new repair coefficients, prescribed waiting root and radius source. The absolute-pressure constant Cp is recomputed from the same forward history and full future, and remains retained. Full physical graph installation and temporal scale recursion are separate work.

Implementation: `experiments/root_st073/lei_ren_part1_paper_compliant_current_angular_terminal_closure.py`.

Receipt: `experiments/root_st073/lei_ren_part1_paper_compliant_current_angular_terminal_closure_check.json`.

Focused execution:

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage currentangularterminal
```

The controller shares one producer/checker object and consumes accepted upstream evidence. It does not regenerate the earlier full physical-chart reports.

## Exact angular identity

Use r=1-mu, k=1-delta/2, L=-30 log(mu), W the prescribed waiting time, and logone=log(1-epsilon). The waiting root is computed from the **uncorrected** reference XT0. The actual angular history includes both supported bump corrections:

```text
A_bump = Aweight*(exp(-3r)*d1 + exp(-r)*d2)
A_corr = exp(-rL)*(Xf(0)-Xf(Z))
       + exp((r+k)/2+kW-logone)*S*Theta(Z)
Xtail_actual-Xtail_raw0
       = exp(-(r+k)/2-kW)*(exp(-rL)*(Xf(Z)-Xf(0))+A_corr)
(1-epsilon)*Xtail_raw0 = 1/k + epsilon*J
(1-epsilon)*Xtail_actual = 1/k + epsilon*J + S*Theta(Z)
Dtheta = 0.
```

Production AST bindings identify the actual bump integral, forward XR/XS/XQ/XT/Xtail formulas, original inverse/quadratic first equation and physical coefficient scale. Common-ball uniqueness identifies the exact coefficient branch. The proof also identifies the same original sigma/phi integrand and measure in J and the current full collar JW integral. Different directed quadrature boxes do not prove function equality.

The correction at Z=0 is retained. It is generally nonzero because the full heat Gamma repair is nonzero there. Using the corrected XT in the prescribed raw waiting-root equation would remove that distinction and invalidate the proof.

## Runtime source installation in this view

`current_forward_terminal()` consumes `exact.angular5.angular(Z)`, `exact.repair.weights` and native `flatten.flatten(Z,100)` angular/absolute-pressure data. It replays both forward terminal histories. All thetaR/S/Q/T factors are recomputed from the same parameters and power length; theta_base uses the newly replayed W/logone. Pressure keeps the original bump-square contribution and the required one-half normalizations.

A restricted heat view exposes only the original shape, local_Gamma and collar_tails methods. Its exact S recipe is the new repair's logS terms. The native stable factor callable is actually invoked on their sum, then intersected with the already proved strong cap. This is an enclosure of the positive exact inverse radius, never a replacement radius or cap-defined field. Gamma derivative bounds use Scap only as a bound. The view exposes no old outer, steep, data, collar or exterior terminal methods.

The returned stress and pressure use this view throughout. They do not call the old companion's evaluate or old outer.angular. This closes the local source-wiring gap while leaving the admitted upstream graph intact.

## Pressure and stress scope

The retained pressure constant is

```text
Cp = new_Ptail + new_full_collar_pressure_future * new_pressure_scale
new_pressure_scale = Ev2 * new_theta_base^2.
```

Both terms are in Pstar-squared units. The original analytic axis datum is preserved. No pressure correction is appended after propagation. The angular stress constant contribution vanishes; the canonical full Gamma theorem then makes the angular exterior stress zero. The axial pressure-constant stress, whose physical radial factor is R^(1/2), remains.

The three scoped gates are:

- `current_native_angular_repair_function_identification_certified`
- `current_heat_angular_terminal_constant_eliminated`
- `current_heat_angular_stress_mixed4_after_terminal_closure_recovered`

The receipt covers 63 production AST bindings, 27 exact identities, four collar/exterior views, 60 angular stress mixed4 rows and 24 stable angular axial5 coefficients. A separate moderate-parameter fixture integrates the original normalized beta and sigma/phi J, constructs both bumps, and verifies the raw/corrected distinction including the nonzero axis repair. Source identities establish zero; numerical interval overlap only checks consistency afterward. Checked loading and an additional unsaved Z=.479, t=.73 evaluation pass. The read-only reviewer remained GPT-5.6 Luna at max reasoning effort and found no material remaining gap in this scoped adapter.

`current_heat_pressure_terminal_constant_eliminated`, both-constant elimination, full stress-free exterior, all physical-chart installation, global admissibility, physical energy, flat remainder, complete Cartesian NS validation and temporal recursion remain false.

## Next tasks

1. Bind the common analytic pressure datum and native Mp/P0 history through pulse/flatten to the corrected pressure branch. Keep every Pstar, Ev and Rrel factor explicit.
2. Replay the quadratic pressure moment equation together with full collar/Gamma pressure future and epsilon atoms; prove whole-Z Cp=0 from those exact equations.
3. Add a dedicated current pressure-terminal receipt and then a separate full current stress-free-exterior receipt. Do not reinterpret this angular-only receipt as either one.
4. Install the new branch in full future C4/C5 and selected ap/c1/c2 owners, with matching prefixes, scalar factors and fresh caches; then complete quantitative native physical interfaces.
5. Continue physical stress/residual composition, point evaluation, global cone, independent flat remainder, prescribed-domain energy and genuine n-dependent temporal recursion before oscillatory correction.
