# Current original-pressure terminal closure (2026-10-05)

The prescribed original analytic preheat datum is identified with the negative complete native raw pressure integral. Together with the admitted second-quadratic/full-Gamma pressure balance, this proves current Cp=0 as a function of Z in [-1,1], through axial order five. This closes the pressure terminal in the restricted current heat view; all physical charts and the full exterior stress still require separate source installation/admission.

Implementation commit: [146fe1b1](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/146fe1b1).

## Exact functions and directed enclosures

`ExactOriginalPreheatPressureOperator` exposes the fourteen original Section 6.1 velocity densities, exact stage integral expressions, M2, M0 and Fflat(Z), and a callable `selected_datum(Z)`. The latter substitutes the prescribed raw waiting root. The symbolic parameters retain their checked common original constructor definitions; huge positive exponential factors are not replaced by numerical caps.

Native densities are independently translated from production pressure primitive AST. The reference value and radial primitive ODE identify the complete incoming Mp prefix. Stage density identities precede integral additivity and the P0 identity. There are no unbound replacement symbols M2, M0 or Fflat in the admitted cancellation, and the native raw integral does not use P0 to define itself.

The datum's immutable source definition remains unchanged. `CompliantPressureDatum.m2/m0`, flatten bounds, `normalized_jets()` and native numerical integral outputs remain directed enclosures of the exact functions. No midpoint, endpoint or arbitrary interval value defines P0. Old forward Cp and P0+raw-total enclosures are retained only as consistency diagnostics.

## Raw waiting, amplitude and convergence

The exact axis chain starts at X(0)=5/8 and passes through the live buffer/pulse sources and the complete flatten. At Z=0 its flatten exit is `2*(Xv + integral_0^100 2^(-sigma(v/100))*exp(r*v)dv)*exp(-100*r)`. This is axis-only. The following XR, XS and XT are uncorrected and determine the unique original raw waiting root. Corrected XT never enters this equation. The older coarse waiting box is bound to the same equation, rather than identified by overlapping endpoints.

The actual amplitude identity is `exp(2*logBase)=exp(2*logEv)*theta_base^2=forward['pressure_scale']`, with the original Pstar units and both half factors retained. Infinite-tail domination uses monotonic J and the actual negative remainder `-mu*(Jd-Jrel)-(Jrel-Jq)-delta*Jq/2`, the guarded raw collar bracket in [0,1], and the exact tail amplitude. A complete holomorphic-strip majorant justifies integral differentiation through axial order five. Integrating an unrelated envelope alone is not used as a proof.

## Evidence and use

- 123 exact source identities, including 14 production density translations and the complete raw angular chain.
- 84 independent original global-to-segmented pressure density rows over all fourteen stages, including six nonzero infinite-tail rows. Moderate fixture lengths test algebra, not final parameter admissibility.
- Seven unilateral incorrect-source mutations rejected: native density factor, missing flatten exponent, flipped/perturbed P0, missing axis flatten factor, perturbed Xv and perturbed raw waiting root.
- Four current heat views retain 60 pressure mixed4 rows and 72 source-zero axial coefficients. Checked loading and fresh Z=.709, collar offset .83 acquire the admitted view from the live common graph.
- No full upstream graph/report regeneration was needed.

Producer/checker: `experiments/root_st073/lei_ren_part1_paper_compliant_current_pressure_terminal_closure.py` and `_check.py`; scoped receipt: `current_pressure_terminal_closure_check.json`. Focused stage: `python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage currentpressureterminal`. Producer and checker share one object; a live session can pass the already admitted balance/raw graph directly.

## Remaining work

Issue a separate full current stress-free-exterior receipt only after the canonical full-Gamma stress theorem and all required angular/pressure/energy/meridional histories are tied to one current source. Identify/rebuild selected and complete future C4/C5 owners on the new branch where needed; consistently update ap/c1/c2 and caches. Install source views in physical chart owners and finish quantitative native interfaces, points, global cone, independent flat remainder and energy on the prescribed physical domain. Genuine n-dependent temporal recursion and oscillatory stress correction remain unimplemented. Spatial Taylor coefficients do not establish temporal recursion, and this scoped pressure result is not a blow-up theorem.
