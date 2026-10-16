# Signed m15 from the core through cutoff and lower collar - 2026-10-16

`OriginalN1InnerM15Numeric.evaluate(Z,q)` now carries the signed original core m15 through the selected velocity cutoff and lower collar to rminus. Inputs are decimal strings, -1<Z<1 and0<=q<=2. The same cumulative function has q=0 at rho=4, q=1 atRin and q=2 at rminus. It returns six factored ordinary C0/C1 terms with numerical approximations and complete directed source/error enclosures. The signed core inlet from the preceding milestone is retained unchanged once.

The full goal is **ACTIVE / INCOMPLETE**; tasks14/15 remain IN PROGRESS. This implements the signed inner/lower boundaries, but they are not yet substituted into the external preRp enclosure API. Signed active/compact finite-N pressure, independent numerical Rp cancellation, global five-moment admission/mixed joins/common N, heat/energy/stress/flat remainder, genuine temporal n>=2 and corrected NS remain open. A useful relative error for the combined multiscale moment is not claimed.

## Exact original phase and two different cutoffs

Write a=sc/4, h=hb, t=a*q and rho=4*exp(h*t). The physical source geometry is

```
R=epsilon*rho,  d-rho=h*a*rho*dq.
```

The leading prescribed-shear field uses sigma(t). The own n1 velocity cutoff is instead chi(q)=1-sigma(4q-2): it equals1 through q=1/2 and0 from q=3/4 onward. On the lower q in[1,2], t<1/8; leading sigma(t) is not1. The original a is about1.65e-204453045017284839. Materializing sigma(a*q) is infeasible, so its actual nonzero amplitude remains a log atom.

The original cut velocities f_cut and u_cut are used only inside their admitted raw-inner domain. They vanish outside the own support. Leading N1p remains active after that support, and the accumulated moment is not reset at q=3/4 orRin.

## Pressure derivative source and exact inverse-width cancellation

For u=Uz1/epsilon and p=P1/(epsilon*Pstar^2), the exact local moment increment is

```
Delta m15 = epsilon^2*h*a*integral
    (2*rho*U_lead*u_cut-Pstar^2*rho^2*p_rho) dq,
p_rho = Np + 2*S*Phi*f_cut/(epsilon^3*Pstar^2).
```

Here Np is the normalized leading pressure source from the unchanged first-collar forcing body. The actual original cut-F coupling expression is executed independently in the focused driver to check the plus sign, factor2 and physical powers. U_lead includes the axis j and epsilon*Psi.

Np has one h^-1 Laurent term. After multiplying by the physical h measure, its complete contribution is

```
epsilon^2*integral_0^t A(v,h) d-sigma(v),
A=2*rho^2*Z*(1-h)*(Phi_actual/Phi_comparison)*Psi_comparison_rho/L,
A0=32*Z*Psi_core_rho(4)/L.
```

The Pstar^2 factor cancels exactly in this singular term. The entire1-h factor stays in A; its h part is included in the width remainder and is not added again to the regular source. A private regular leading-frame calculation removes this complete sigma-prime source after the split. Independent symbolic evaluation of the actual production forcing body checks both its regular row and singular C0/C1 factors.

Original `WidthChange` bounds retain rho variation, the actual/comparison Phi ratio, comparison Psi_rho and1-h. They prove |A-A0|_norm<=h*E uniformly on the complete first/lower phase. The comparison-change methods return coefficients before h; the original caller applies h through its source-width bound. No a<=h assumption is made. Positivity of d-sigma then gives

```
singular increment = epsilon^2*sigma(t)*A0
                  + epsilon^2*h*sigma(t)*[-E,E].
```

The second term has the correct ordinary Z derivative norm factors. This computation never subtracts a rounded core pressure endpoint from a nearly equal collar endpoint. Constant pressure cancels algebraically; all original core m15 terms stay in the cumulative result.

## Nonzero flat amplitude and regular integral errors

For t>0,

```
odds=(1-t)^(-2)-t^(-2),
sigma(t)=exp(odds)*kappa,  kappa=1/(1+exp(odds)).
```

The original phase certifies odds<=log(t)-1000, so exp(odds)<=t*exp(-1000). Kappa is enclosed directly near1, while odds remains a source log atom; neither exp(odds) nor h^-1 is materialized. The q=0 branch is the exact empty prefix. Source log-atom approximations have their own directed enclosures; they do not select an implicit original scalar as an exact numeric value.

Regular axial and pressure integrals retain their separate epsilon^2*h*a and epsilon^2*Pstar^2*h*a units. Cells split at q=1/2,3/4,1,2. Local leading sigma values use the original flat derivative provider; actual/comparison field errors remain in the cell ranges. Directed cell ranges include quadrature and field uncertainty; midpoint values are numerical proxies with their full distance to the integrated source enclosure charged. The original full core fields and first-collar sources are reused without changing accepted implementations.

## Focused evidence

The final run exited0 in114.667s; saved-report inspection exited0. The1,751,167-byte report SHA256 is `0615f83f6a0c67ee506ca143c1d4ba698d06a92a63e018992b608c15dc930d95`; it pins1719 scientific inputs and18 implementations.

- 6 symbolic C0/C1 comparisons execute the actual original regular/singular forcing and cut-F coupling, including the complete1-h factor.
- 6 independent flat-log/kappa conditions and6 direct positive d-sigma C0/C1 comparisons.
- 24 original signed core/first/lower packets: Z=-.3,0,.3; q=0,.5,.625,.75,1,1.5,2 with16 cells per fixed piece, plus the3 terminal8-cell refinements.
- 399 source/serialized-error conditions,21 unchanged signed-core inlet conditions,24 exact empty-prefix zeros and36 nonzero flat-log atom conditions.
- 12 conditions distinguish stopped velocity increments from continuing nonzero lower pressure history,24 directed cell-refinement C0/C1 comparisons and4 rejected invalid requests.

Use `python scripts/run_original_n1_inner_m15_numeric.py --inspect-report` to inspect saved evidence. `--write-report` reruns the focused construction when implementation changes require it. Existing read-only worker `/root/n1_atlas_map`: **GPT-5.6 Luna / max**; no new agent. It confirmed the source/units and corrected the phase-map and width-bound objections after the exact two-cutoff and pre-h error conventions were checked.

## Required next construction

1. Consume this q=1 inner boundary and q=2 inner+lower cumulative result in the external preRp assembly. Replace only the old inner/lower broad components; retain active/compact/own/quiet once, in their original physical units and function bindings.
2. Directly evaluate signed active and compact finite-N S1 pressure sources with actual controls, ordinary Z1 and finite RQ_y/RQ_yZ endpoints. The prior finite-N error majorants remain bounds, not selected signed corrections.
3. Preserve selected inverse/defect correlations to quantify the independently transported Rp value/derivative cancellation, then complete all five original numerical histories and required mixed joins.
4. Continue original flatten/analytic preheat pressure/exact heat/energy/common N, admissible stress/flat remainder, true temporal n>=2 recursion/smooth sum, oscillatory correction and independent corrected Cartesian NS validation.
