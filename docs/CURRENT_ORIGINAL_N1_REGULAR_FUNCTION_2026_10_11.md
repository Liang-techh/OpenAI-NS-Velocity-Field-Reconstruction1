# Original n1 common-inner regular series function and first radial jets

The original n1 coefficient now has a shared function definition and a callable enclosure interface over the accepted common inner domain. All six components use the same original leading core, actual first collar, exact coefficient equations and remaining feedback bound. The interface supplies function Z jets through order 3, first radial mixed Z jets through order 2, and its own mean/Q/V recovery. Six nonaxis queries each pass all **12 six-component value/Z1 width targets**. The full reconstruction goal remains **ACTIVE / INCOMPLETE**.

Read this checkpoint before [the natural-component checkpoint](CURRENT_ORIGINAL_N1_NATURAL_COMPONENTS_2026_10_11.md). The numerical centers still compute the zeroth source integrals; the rest of the convergent series is enclosed by an explicit tail. Higher Z entries that lack signed centers receive analytic bounds. This is an interval function interface, with no claim of an explicitly summed feedback prefix or a completed repaired n1 coefficient.

## Function definition and shared source path

The natural unknown is ordered as

```text
W = (F1, Uz1, K1, P1, d_x_F1, d_x_Uz1)
W0 = Gprincipal*g
V(f) = Gprincipal*(C0*f + B1*d_Z(f))
W = sum_{k>=0} V^k(W0).
```

`exact_operator()` restores the genuine natural-normalized B0/B1 and source rows from the accepted leading own-mean packet. It includes every nonlinear/pressure sector and the remaining F0 gauge. Positive majorants are separate data used to bound the function; they do not replace these matrix callbacks.

The principal inverse has source columns 4, 5 and 6. Its pressure kernel is 1; the axial kernels are `s*log(x/s)`, `-s*(1-s^2/x^2)/2` and `s/x`. The swirl kernel is the entire potential series

```text
Kchi = sum_{j>=0} (-2*chi)^j * T^j(K0)
K0 = s*(1-s^2/x^2)/2
T = (d_xx + 3*d_x/x)^(-1), with regular zero-axis data.
```

The finite kernel adapter returns value/Z1 kernels and a factorial tail. It rejects unsupported kernel Z2 requests. The exact function definition uses the infinite kernel. Every word uses the same path from 0 to x and inherits the core/actual-collar source join; no coefficient or cumulative history is reset at the join.

The regular space uses `W_i/x^p_i`, with `p=(2,2,2,2,1,1)`. First four components are even at the axis; last two are odd. A conservative Green factor 64 controls the supported residual columns on the original holomorphic chi tube, where |chi|<=1. For source row 5 and input power 1, the combined zeroth kernel column is `t^2*(1-t^2)/2+t^4<=1`; the remaining columns sum to at most `q*exp(q)`, since `(j+1/2)/(j!*(j+1)!)<=1/(j-1)!`. Input power 2 adds at most a. The uniform cover `a*(1+q*exp(q))` is below 35 and therefore below 64; pressure and axial columns are smaller. The near-axis initial bound explicitly recovers uniform core source rows from the admitted production majorants:

```text
core_L1 = 2 * core_uniform_source
g4_uniform = 2*a*N1p/Pstar^2
bound(g4/x) = core_L1[g4]/(2*a)
bound(g5)   = core_L1[g5]/2
bound(g6)   = core_L1[g6]/2.
```

Their directed source norm is checked against M0. For x<=1 the zero-axis kernels give the weighted initial bound; for x>=1 the existing unweighted M0 is sufficient. An arbitrary L1 bound is never inferred to be a pointwise source bound. Source M0 uses eta/8, the value series uses eta/16, and final Z derivatives use eta/32 with a separate eta/32 Cauchy gap. The original physical eta, hb, F0, pressure and Lambda parameters remain unchanged.

| Regular natural quantity | Conservative upper bound |
| --- | ---: |
| Weighted initial M0 | 6.884152e74 |
| All-word ratio | 8.242090e-282 |
| Remaining value tail | 3.255715e-205 |
| Ordinary Z3 remaining tail | 5.835879e-97 |

These are coefficient-function bounds in the stated normalization. They are not physical Navier–Stokes residuals or a percentage of the full project.

## Radial derivatives and own recovery

At x>0, the first derivatives use the genuine six-row equation:

```text
d_x W = -D*W/x + B0*W + B1*d_Z W + g
D = (0,0,2,0,3,1)
d_rho W = d_x W/(2*x), rho=x^2.
```

At the axis the code uses finite limits instead of dividing by x:

```text
d_x W = (0,0,0,0,g5/4,g6/2)
d_rho(first four W components) = (g5/8,g6/4,-g6/8,(g4/x)/2).
```

The odd last two components have no finite rho derivative at the axis and those entries are omitted there. The n1 mean and radial fields are recovered from its own coefficients:

```text
mean1 = Uz1 + K1
Qnorm = [(1-delta)*Z*Uz1 - (1+delta)*Z*K1 - (1-Z^2)*d_Z(mean1)]/L
Vnorm = rho*Qnorm
d_rho Vnorm = Qnorm + rho*d_rho Qnorm.
```

Physical Q1 uses epsilon and V1 uses epsilon^2. Swirl physical Z derivatives add the original F0 gauge once. Output grades are W:Z3, d_x/d_rho W:Z2, Q1:Z2 and d_rho Q1:Z1. An available enclosure grade does not imply a narrow relative error at that grade.

## Callable entry and checked samples

Production module:

`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_regular_series_function.py`

```python
owner = CurrentOriginalN1RegularSeriesFunction(accepted_natural_component_owner)
packet = owner.evaluate('.371', '.625', chart='first')
report = owner.report(packet)
contract = owner.contract(N=0)
```

The constructor requires the accepted live natural-component owner and its original source graph. `core` queries accept exact real Z in [-1,1] and rho in [0,4]. `first` queries accept the selected actual-collar fraction in [0,1], mapped to the corrected original source phase. It does not use the old canonical half-collar exit.

Checked samples: core axis, core rho=1.5, core exit, actual first inlet, actual first interior, corrected selected inner exit, and another core query at Z=-0.4. All six nonaxis samples pass 12/12 component value/Z1 rows. A checked-owner query at an additional first-collar fraction 0.625 also passes 12/12. The zero axis is not counted as a relative-width success. Core exit and actual inlet agree in function, first radial jets, own recovery and exact equation packets.

The checker supplies 52 exact polynomial-kernel ODE/boundary identities, arbitrary-forcing pressure/axial/K FTC identities, an independent physical theta-source Z2 fixture, source-join bindings, near-axis source-bound bindings and packet/unsupported-scope guards. A read-only GPT-5.6 Luna / max review examined the kernels, regular bounds and recovery formulas. It prompted the explicit near-axis pointwise ledger above.

## Completed scope and next work

- [x] Define one genuine six-component regular Volterra series over the common core/selected actual collar.
- [x] Expose the true coefficient/source callbacks separately from their norm covers.
- [x] Supply an exact principal inverse definition, finite value/Z1 adapter and infinite factorial kernel tail.
- [x] Bind near-axis pointwise source factors before applying the weighted initial bound.
- [x] Enclose all remaining feedback words and ordinary Z derivatives on explicit nested tubes.
- [x] Publish six-component function queries and first radial mixed Z2 jets, with finite axis limits.
- [x] Recover n1's own mean/Q/V and their first radial derivatives from the same coefficient function.
- [ ] **Next: second radial jets.** Differentiate the genuine equation, including actual-collar source and matrix radial derivatives. Preserve inverse hb factors until their appropriate integral/Jacobian cancellation; do not set unsupported derivatives to zero.
- [ ] Bind the leading higher-radial mixed derivatives needed by N1p_r and N1z_r. List the required grades before extending the source owner; keep accepted ancestors immutable and use a new adapter.
- [ ] Establish finite axis Taylor/remainder compatibility and a common function extension at the derivative grades used downstream. The current finite limits alone do not close this task.
- [ ] Establish full derivative joins at the core/actual-collar seam. Existing value/first-radial equality is a scoped result.
- [ ] Compute cumulative n1 moments from the genuine coefficient function, including core history and actual collar measure. Carry the all-word tail into each integrated error; do not reuse leading-order defects as n1 defects.
- [ ] Build n1's own moment target and independent repair map. Specify target units, pressure compatibility, correction support and terminal identities as functions of Z. Keep the n1 repair gate false until this is checked.
- [ ] Extend the repaired coefficient to all required completed-leading charts with compatible support, pressure and heat-exterior data.
- [ ] Only after n1's own completion: implement genuine n-dependent n2 sources, its own regular solve and own moment repair. Copying n1 or changing coordinates is not temporal recursion.
- [ ] Then finite-order recursion errors, smooth temporal sum, admissible stress, two-family oscillatory correction and full corrected NS residual/physical scale diagnostics.

Still false: `n1_full_inner_interval_solution_certified`, `n1_radial_Taylor_remainder_certified`, `coefficient_F1_Uz1_P1_solved`, `positive_order_five_moment_repair_closed`, `n1_Volterra_iteration_prefix_computed`, `temporal_recursion`, global admissibility, flat remainder and full corrected NS solution. The new scoped acceptance gate is `current_original_n1_regular_common_inner_series_and_first_radial_function_enclosures_installed`.

## Concrete handoff for second radial jets

Use a new adapter; do not change a checked ancestor. For x>0 the equation to implement is

```text
W_xx = D*W/x^2 - D*W_x/x
       + (B0)_x*W + B0*W_x
       + (B1)_x*W_Z + B1*(W_x)_Z + g_x.
```

Start with mixed Z1. The existing W:Z3 and W_x:Z2 grades are sufficient on the coefficient side. Differentiate source and matrix formulas analytically with rho_x=2x; do not finite-difference independently boxed point queries. In particular,

```text
(g4)_x = (2*N1p + 4*rho*N1p_rho)/Pstar^2
(g5)_x = 2*x*(g5)_rho
(g6)_x = 4*x*N1z_rho
          - epsilon*(8*Z*x*N1p + 8*Z*rho*x*N1p_rho)/L.
```

Extend the factored leading rows through rho3, retaining at least Psi/Phi:Z2 and meanPsi:Z2 there, so Q_rhorhorho:Z1 is genuine. Derive the actual collar rows from its prescribed-shear ODE with chiSpatial_rho and chiSpatial_rhorho, preserving both hb^-1 and hb^-2 sectors. The own-mean identity supplies `meanPsi_rhorhorho=(Psi_rhorho-3*meanPsi_rhorho)/rho`. Reuse the original core's admitted total-order-five inputs where available, and explicitly reject a missing source grade rather than pad it.

At x=0 derive separate finite limits/Taylor coefficients and an error remainder before using W_xx there. At the core/actual inlet require equality of the new source and coefficient derivative packets. Only this new bounded scope should be admitted; own five-moment repair and full n1 completion still need their independent work.
