# Finite-width Ra-to-R100 signed source integrals and own moments

The new `actual_bridge_integrals` provider preserves the original known comparison direction from Lei-Ren v2 (4.34)-(4.36), (9.23)-(9.26), while integrating the actual field's own histories from the current six core atoms. It does not replace the comparison direction with an actual-stress feedback equation.

Run the bounded stage:

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage actualbridgeintegrals
```

The producer/checker and their receipts are `experiments/root_st073/lei_ren_part1_paper_compliant_actual_bridge_integrals{,_check}.{py,json}`.

## Source field and finite-width terms

Here `ell=log(F/f)`, where `f` is the actual core exit value. The axial quotient is `(f/barF)*exp(ell)`. This is equivalent to the paper's `log(F/barF)` representation; the two logarithms must not be confused.

The first microscopic control is `(1-sigma)+hb*sigma`. Its full endpoint mass is `hb*(1+hb)/2`, so both width powers are retained. The second microscopic interval contributes at `hb^2`, and the macro interval uses `chi=hb`. All source coefficients remain signed before any final interval bounds. Intermediate cutoff masses use directed quadrature remainders; only the complete first interval has the exact reflection-symmetry mass `1/2` for each component.

The finite-width comparison field keeps its admitted `hb^2` polynomial and separately bounded `hb^3` remainder. The macro comparison starts at **R0=Ra*exp(2hb)** with **L=log(100/Ra)-2hb**. Its exact moment modes have exponents `0,-1,-2`. Radius kernels for the angular/hydro/pressure terms and the swirl term are integrated with powers `R` and `R^2`, simplifying the fixed endpoint to `R100` before arithmetic.

The actual axial source is decomposed into signed modal integrals and a controlled nonlinear angular-prefix error. The weighted axial jet norm bounds every prefix, using integrals of absolute modal norms. A signed endpoint value is never used to bound earlier prefixes after cancellation. Pressure and swirl factors remain logarithmic and formal.

## Actual own histories

H/M/K/A/B/C use the fresh actual core atom initial values, with rates `2,1,2,1,2,1`. Their sources are respectively `2phi,V,2phi*V,V^2,phi^2,phi^2`. Thus B has equilibrium target `phi^2/2`. Raw axial `V=Uz` remains distinct from the pressure primitive `V_pressure=4C`.

The microscopic own moments use rigorous positive Volterra-kernel enclosures of the actual prescribed F/V prefixes, including every mixed and quadratic change. They are not selected point histories.

The macro propagates that actual microscopic endpoint and integrates signed linear feedback with

```text
int_0^L exp(-rate*(L-t))*int_0^t R0^p*exp((p-j)*s) ds dt,
p=1,2; j=0,1,2; rate=1,2.
```

Both resonant branches are retained. The full `exp(ell)-1-ell` remainder is controlled at width power2. Axial `exp(ell)-1` errors and all `delta_phi*delta_V`, `delta_V^2` and `delta_phi^2` terms remain. This restores source-defined signed feedback and conservative finite-width enclosures, rather than replacing actual moments with the comparison moments.

## Whole axial domain and amplitude bounds

The report supplies the three original source charts (Z=0, Z=.5, exact shared root), both microscopic joins, and a complete Z[-1,1] R100 packet. The API accepts real source Z intervals.

Whole-axis interval division loses the shared `H^2` in `chi=H^2/(H^2+sigma^2)`. The provider locally replays the unchanged original atom/profile methods with a final model-tail callback that intersects only `chi0` with the exact real range `[0,1]`. Every derivative row and original finite-field expression is unchanged; no module or existing source file is patched. The exact positive source width admits the stronger enclosure-only cap `1e-50000`, needed for the large whole-axis derivative bounds. This is neither a new width nor a selected field parameter.

**Gbar is a bound, not G(Z).** The new provider encloses the full exact factor `F0(Z)^2=Cstar^-2*exp(-2Lambda*G(Z))`, from the lower log bound `-2logC-2Lambda*Gbar` through the upper bound `-2logC`. A lower bound alone cannot define or enclose the swirl amplitude. The older macro/first-switch/switch signed receipts require this amplitude-range correction before being composed as actual finite-width source bounds.

## Focused evidence and remaining work

The stage passes 298 current hashes, 5 original source/control AST bindings, 813 logarithmic source-product bounds, 108 fresh atom coefficient checks, 72 microscopic source joins, 1349 finite actual-field/own-moment coefficients, and 240 signed source terms. Thirty-three independent moderate-parameter direct integrals check the radius kernels, both resonant double Volterra kernels, nonlinear axial prefix bound and own-moment kernels. Tolerance is `1e-55`; maximum positive enclosure miss is `7.982e-85`. Fixtures are algebra evidence, not production parameter choices. The read-only GPT-5.6 Luna/max review accepts the double-kernel and nonlinear-feedback algebra.

`coordinate_R100_endpoint` asserts only the exact radius. The actual R100 functional join to the existing switch provider remains false. Production point histories, installation into existing mixed4 providers, corrected R100-to-R110 composition, implicit leading input/Jacobian recomputation, completed global tensor admissibility, global temporal-flat/volume/energy bounds, true n-dependent recursion and oscillatory correction remain open.

Next actions:

1. Correct the existing macro/first-switch/switch F0-squared logarithmic scale range and refresh only the affected signed-integral receipts.
2. Bind and transfer this arbitrary-Z actual R100 field, all six own moments, the canonical P0 and axial jets to the two short switches. Retain their distinct angular/axial controls and all incoming histories.
3. Integrate/refine the microscopic own-history feedback where the subsequent implicit input requires sharper coefficients; preserve the same source field and positive kernel.
4. Propagate the R110 actual history into reshape/restoration and the implicit five-defect functions, including Jacobian and omitted-term bounds. Do not mark this complete from endpoint overlap or a narrow-width leading coefficient.
5. Continue the original global gates and actual n-dependent recovery/correction sequence.
