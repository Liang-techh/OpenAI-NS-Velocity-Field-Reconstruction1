# Current angular and heat-tail cone

Implementation commit: pending pin after the source/receipt commit.

## Constructed result

**F57C-cone1b-angular/tail is complete.** Six additional nonzero regions now satisfy the original strict two-vector cone on the same checked current graph: angular s[-4,0], steep entry/power/exit/waiting phase[0,1], and heat collar t[0,3). The previous seven regions are inherited, giving **13 nonzero regions**, plus the separate **exact zero heat exterior**. The other 19 registry regions have no current whole-domain cone admission.

The full current A/E/P/K, variable axial K, selected signed C5 controls, quadratic terms, original absolute pressure, tensor/divergence/Cartesian remainder records, seven adjacent joins and four angular support traces remain attached. The heat zero edge is excluded from strict inequalities. These are bounds on the original source functions, not resolved numerical point fields.

This advances admissible-stress stage 4. Genuine n-dependent coefficient recursion, actual two-family oscillatory velocity corrections, global completed-tensor/lift, temporal flatness, corrected NS/energy and measured material winding remain open. The nonflat leading origin remainder has not been cancelled.

## Current source reductions

Set a=delta/2, r=1-mu, k=1-a, b=(1-2a)/2 and g=(mu-a)/r. In angular coordinates the CURRENT normalization is G=exp((a-mu)s), K=G*F, A=K*X, F=1+h; G>=1. There is no second KR factor.

The actual angular numerator is W=k*N-b*Z*N_Z-F, N=F*X. The previous current flatten/power proof supplies the same Hf(Z)>0, so the native base term is g+Hf*exp(-r*(Lrel+s))>=g. The two angular repairs are bounded by

`exp(r*(4+c))*A_full*(k*abs(d_c)+abs(b)*abs(d_c,Z)), c=-3,-1`.

The actual partial A integral lies between 0 and the same full A. This follows from its defining nonnegative density `exp(r*u)*beta(u)` and the nonnegative complementary integral; it is not inferred from a zero interval lower endpoint. The proof binds the original normalized beta, width 3/20, same normalization, A rate/power, current selected controls and both original support branches. Full E/P are retained in the axial bound.

The whole angular W lower exceeds both repair bounds and abs(h). Bounds on h_s/F retain kappa-2=2mu-2h_s/F>0 and a negative actual shear. An exact source logR comparison bounds the negative inverse-radius shear relative to W. The positive Hf is kept before enclosing; no unnecessary flatten error is subtracted.

At steep entry the same actual angular terminal gives M0>0, M=kX-bZX_Z-1. The native ODE is

`M_t+r*(1-sigma)*M=mu-a+r*sigma`.

Thus M>=M0_lower*exp(-r) throughout [0,1]. Current steep power and exit already have positive whole-domain theta source enclosures. Their exact log slopes give kappa-2=2 and kappa-2 in[delta,2], respectively. Current source B maxima use the actual inlet positions and decrease along the original log-radius coordinate.

For waiting, cancel stationary full moment baselines BEFORE interval evaluation. With q=wait*(phase-1) in[-wait,0], Kloc=1-eps,

`theta_loc=exp(-k*q)*(I0-2*S_exact*(1+a)*Kloc*exp(-a*q))`,

`axial_loc=Je0*exp(delta*q)+Jp0*exp(p*q)`.

The relative shear cap uses exp(a*wait); the outside exp(-k*q)>=1 is retained. Both axial modes contract backwards. This is not a cap on the un-factored absolute shear.

Waiting and heat use C=exp(KT-logone)=KR^-1: current theta=C*local_theta and current axial=C^2*local_axial. The cone therefore uses current logB+logC in local units. The proof binds the actual current waiting/heat endpoint K/A/E/P inputs, the complete Gamma/epsilon-atom energy identity, current Dtheta=0 function, original pressure suffix, and actual O7 `entry_stress_rows`. It replays the actual `collar_defect_rows`, `collar_tails`, flat inlet shape and `collar_stress_rows`, proving their common full input and mixed3 operator identity with `waiting_stress_rows`. The source S is a single exact positive inverse radius; caps enter only the final bounds.

## Heat collar and zero edge

The two heat-collar proofs have different lower bounds:

- On t[0,1], the actual sigmoid/Gamma difference gives Gmin=eps/e-Dcap>0, with Dcap=2a(1+a)Scap. Full future and shear contributions have positive margins; an absolute full axial bound yields the strict cone.
- On t[1,3), phi=exp(-4/(3-t)^2)>0. Extract eps*phi before enclosing. With y=3-t, theta_loc/(eps*phi)>=Hmin>0 and abs(axial_loc)/(eps*phi*y^3)<=Cz. The exact flat future integral satisfies `integral_0^y phi(t+u)/phi(t)du<=y^3/8`, yielding the directional cone uniformly.
- As t->3-, the stronger same-source ratio is bounded by `exp(finite_log_prefactor)*y^6`, uniformly Z[-1,1], so the direction tends to e_theta. The finite constant may be enormous; it is not replaced by a small numerical cap. This is a regional spatial edge result, not global temporal flatness.
- At t=3 and throughout the unbounded Gamma exterior, the checked current full terminal functions give exactly zero tensor/divergence/remainder/momentum. A finite anchor supplies layout only. The strict cone is not applied to a zero tensor.

There is **no uniform positive unscaled theta lower on the entire collar**; the record stores `theta_uniform_lower_on_whole_collar=None` and separate `sigmoid_part`/`heat_part` proofs. A query box containing t=3 is not admitted as strictly nonzero.

## Entry points and evidence

Prefix: `experiments/root_st073/lei_ren_part1_paper_compliant_`.

- `current_angular_tail_cone_operator.py`: generic cone algebra, actual current source bindings, current inlet/full moment/stress AST replays and continuous source bounds.
- `current_angular_tail_cone.py`: `CurrentAngularTailCone(flattencone=checked_current_flatten_power_cone)`. The constructor obtains all seven whole views; earlier admissions are nested and scoped separately.
- `current_angular_tail_cone.json` and `_check.json`: 76 generic algebra identities, 79 current AST identities, 46 strictly positive directed bounds, source hashes and scoped gates.
- `current_angular_tail_cone_views.json.gz`: deterministic gzip with all seven complete unpruned signed tensor/history views. Decode with `json.loads(gzip.decompress(path.read_bytes()))`.
- Focused controller: `currentangulartailcone`. Reuse the warm source graph; do not rerun inherited cold producers when their defining inputs are unchanged.

Producer/checker, focused controller and checked scoped/zero-edge API passed. All 799 dependency hashes matched the working tree and Git index before publication. Read-only reviewer: GPT-5.6 Luna/max. The review led to explicit current waiting/heat full-input operator bindings and separate heat subdomain lower bounds; the initially questioned waiting relative shear bound was confirmed correct.

Inventory remains 33 regional tensor routes, 32 adjacent and 14 internal tensor traces; primitive atlas remains 14 adjacent / 8 internal. No global completion gate is promoted by a regional cone receipt.

## Next executable tasks

- [x] **F57C-cone1b-angular-1:** continuous original angular W estimate, current Hf and partial-A positivity/complement, same signed controls/normalization, current K/A/shear/log units.
- [x] **F57C-cone1b-angular-2:** all original E/P/K_Z retained; four support traces, power/angular and angular/entry functional tensor joins preserved.
- [x] **F57C-cone1b-tail-1:** whole current steep entry/power/exit/waiting cone; actual current heat endpoint full moments and O7/collar/waiting operator identity.
- [x] **F57C-cone1b-tail-2:** separate sigmoid/flat heat bounds, exact zero endpoint/unbounded exterior, uniform edge direction and original absolute pressure.
- [ ] **F57C-cone1b-entrance-1:** start from checked `registry.owners['entrance']`, original pulse entrance domain and Z[-1,1]. Recover actual incoming and local velocity/shear contributions in original cone units. Retain all signed incoming stress and moment sectors, source amplitude logs, full absolute pressure and entrance/main join. Output complete source views and a continuous sign/correlation bound; interval crossing zero is an inconclusive bound, not a proven sign failure.
- [ ] **F57C-cone1b-entrance-2:** port original pulse entrance comparison using the SAME current terminal/incoming data. Keep exact exponential amplitude ratios before applying caps. Bind full source shear and numerator to native AST; prove T dot S<0 and the directional inequality for every point of the original domain. Do not extrapolate from main/exit samples.
- [ ] **F57C-cone1b-O3:** original O3 slope-to-mu, power and pulse entrance joins. Use actual variable log slopes and source primitives, ordinary log-radius derivatives and complete forward moments. Close each full domain separately; identify any actual wrong sign with chart, source sector and interval.
- [ ] **F57C-cone1b-O2:** Rh reference, O2 slope/axial/buffer and O2/O3 join. Keep nonzero axial source and actual full pressure/energy. Derive shear stress correlation before widening intervals, then publish same-current full-domain cone or a precise analytic obstruction.
- [ ] **F57C-cone1b-inner-1:** actual patch, restore and reshape. Preserve variable K_Z and all mixed histories, source-coordinate joins and axis restrictions. Derive original required shear loop, reserve independent moment-repair intervals and determine a finite uniform N. Do not substitute an unrelated field to force cone signs.
- [ ] **F57C-cone1b-inner-2:** switch power, both microswitches, first/second/macro bridges and positive-radius core. Use original periodic shear-loop construction where necessary; independently repair five moments on reserved intervals. A defining-family change requires rebuilding affected closure, pressure, tensors, joins and cone receipts. Keep the analytic core-axis extension separate from positive-radius cone coordinates.
- [ ] **F57C-cone1b-global:** only after all 33 source regions and relevant zero-edge limits have been handled, assemble the completed admissible stress with a smooth global cone lift. Verify derivative bounds and sign throughout joins. Current two-vector cone certificates alone do not prove the completed diagonal tensor admissible.
- [ ] **F57C-cone1c-pulses:** implement the TWO actual homogeneous oscillatory pulse families from the paper. Specify their velocity/potential, support, divergence identity, stress covariance integrals, frequency and finite errors. A positive reference covariance matrix is not an actual pulse family.
- [ ] **F57C-cone1c-amplitudes:** bind covariance to current background stress; prove positive squared amplitudes uniformly. Construct flat edge weights and smooth square-root extension with derivative bounds, using source zero edges rather than clipping negative boxes. Keep mean and oscillatory correction contributions explicit.
- [ ] **F57C-cone1c-cancellation:** implement the signed linear lift/mean correction and actual averaged quadratic momentum-flux cancellation. Compare uncorrected and corrected stress on common current physical coordinates, reporting finite-frequency and derivative errors.
- [ ] **F57C-recursion-n1:** recover the paper's actual n=1 equations and coefficient sources on common inner domains, with independent moment repair. Address the certified nonflat leading origin Ez before declaring a flat remainder.
- [ ] **F57C-recursion-nge2:** implement distinct n-dependent recovery for n>=2, finite-order estimates and smooth summation. Apply truncation to streamfunction/potential before curl. Coordinate rescaling of n=0 is not coefficient recursion.
- [ ] **F57C-final-fields:** resolve u/v/w/p with tolerances at multiple physical times; independently evaluate corrected Cartesian divergence and NS residual, prescribed-domain energy, radial tail, measured scale recursion and true material-line winding. Report interval bounds separately from resolved point values.

The long-term goal remains active. Completed flatness, actual oscillatory correction, corrected NS/energy and genuine coefficient-recursion gates remain false.
