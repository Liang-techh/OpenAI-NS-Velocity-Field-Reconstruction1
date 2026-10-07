# Original native signed averaging and mixed slow density adapter

Checked implementation: [8c481507](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/8c4815071ce603e6b2f7b401e670a873b915cffc). `NativeSignedAveraging.route(Z)` rebuilds the original inlet-to-Rc source integrals on **24 true cells /17 charts**, for Z=[-1,1] and Z=[.49,.51]. It uses the original zero inlet, uniform global fractional phase and true radial widths. The five correction histories and first Z derivatives now have **uniform N^-2 covers for every integer N>=160**, with endpoint terms and inherited memory retained.

This is an additional asymptotic certificate, not measured small error or an installed velocity/repair field. The old chart-wide slow y/yZ caps are coarse. The new averaged certificate is **not sharper at N=160 for any of the five normalized rows in either Z domain**; both direct and averaged certificates are retained. The first-bridge source/mixed-derivative conditioning is the next production task. Five actual controls, terminal identities and genuine scale recursion remain incomplete.

## Exact signed split

E=Utheta/Pstar, V=Uz/Pstar and B=original loop primitive/Pstar. A below is the phase primitive, distinct from terminal repair amplitude Ac/S. Original reflection gives A(1-phi)=-A(phi), B(1-phi)=-B(phi), using t0=-b/a. With uniform fractional phi, the leading vector and its fixed-phase C0/y/Z/yZ derivatives have zero mean:

`f1=(B,E*A,E*(V*A+B),2*V*B-E^2*A,E^2*A)`.

Use `deltaE=E*A/N+R_E/N^2`, `deltaV=B/N`, `R_E=E*A^2*R2(A/N)`, `R2(x)=integral_0^1(1-s)*exp(s*x)ds`. Its Z rule is `R_E_Z=E_Z*A^2*R2+2*E*A*A_Z*R2+E*A^2*A_Z*R3/N`, with `R3(x)=integral_0^1 s*(1-s)*exp(s*x)ds`. The signed density remainder vector is `(0,R_E,V*R_E+F_N*B,-E*R_E+B^2-F_N^2/2,E*R_E+F_N^2/2)`. All coefficient functions still depend on N; directed covers are uniform. Original |A|<=159 bounds the exponential remainder without clipping field values.

## Mixed source definitions and units

The separate compressed source graph exposes inverse_yZ, A_yZ, normalized loop B_yZ and all20 leading-density slow rows on each of17 charts. It differentiates original free-angle function nodes, then applies the exact implicit inverse, T1 and A/B chains. Generic partial derivative instructions have formal function semantics; a numeric point inverse/jet evaluator is not installed here. Accepted original fixed-phi y/yZ caps are used only as absolute bounds. Original E/V are already normalized by Pstar with ordinary log-radius derivatives; no extra width, Pstar or radial half-shift is applied. Input source B=2*V_y is distinct from the loop B used in f1.

## Endpoint-retaining integration

Set G_j=integral_0^phi f1_j ds. Zero mean makes G periodic, G(0)=G(1)=0, and |G_alpha|<=sup_phi|D_slow^alpha f1|/2 for alpha=0,y,Z,yZ. With K=exp(-lambda*(b-y)),

`integral_a^b K*f1(y,Z,frac(N*y+phi0))/N dy = ([K*G]_a^b-integral_a^b K*(G_y_slow+lambda*G)dy)/N^2`.

Use G_Z and G_yZ for the Z row. The cell coefficient bound is `(1+decay)*Gcap+mass*(Gycap+lambda*Gcap+remainder_cap)`, with true mass applied exactly once. Adjacent endpoints are retained; no seam cancellation or chosen phase is assumed. Each new increment is transported with the original decay from the actual zero inlet. Quiet power has exactly zero local increments, and pressure rate0 retains upstream C0/Z memory. Original P0/P0_Z remains unchanged.

Normalized targets use the same original Rc A/A_Z, mu and signed divided-row numerator as the preceding stage. Both direct and averaged bounds apply to the same functions; the normalized N*r C1 certificate is `min(prior_uniform_cap, averaging_coefficient/N)`. Sufficient crossover log-N thresholds compare conservative certificates and never materialize an astronomical integer. They do not measure actual error or admit global N/cone/controls.

## Focused evidence

- 340 mixed leading source graph rows, 960 C0/Z contribution/transport rows and 480 inherited rows.
- 40 normalized target rows, 40 exact quiet rows and 8 retained quiet pressure rows.
- 150 independent original-density remainder C0/Z references at N=160,1024,10^20, plus 6 independent finite-interval IBP references with nonzero endpoints. General mixed inverse/T1/A/B chains checked symbolically.
- Producer 72.250s; focused checker 77.500s on the same warm original source graph. 1090 working/index dependency hashes match. Read-only review: **GPT-5.6 Luna / max**.

## Agent tasks and acceptance criteria

- [x] **AVG1 original zero-mean split:** bind the current original primitive reflection theorem and t0=-b/a, use the uniform fractional-phi measure, expose all five signed leading density functions and their fixed-phase Z derivatives. Keep the nonlinear F_N out of zero-mean claims.
- [x] **AVG2 exact second remainder:** expose R2/R3, the complete R_E_Z rule and all five signed cross/quadratic density terms. Bound them for every integer N>=160; do not rescale selected-frequency histories.
- [x] **AVG3 mixed source interface:** expose exact formal inverse_yZ, A_yZ and normalized loop B_yZ definitions and all C0/y/Z/yZ leading-density product rules on17 charts. Bind accepted original slow-jet caps. These are function definitions and bounds, not a point inverse evaluator.
- [x] **AVG4 source integral C0/Z:** integrate each original cell with its true log-radius Duhamel kernel, retain both endpoints and G_y/G_yZ, transport all inherited memory, and keep quiet rate-zero pressure. Export same-source N^-2 Rc correction and target covers.
- [x] **AVG5 preserve competing certificates:** report min(prior direct bound, averaged coefficient/N^2) for the same correction function; publish normalized N*r crossover log thresholds and conditional exact repair log conditions. The new certificate is not sharper at N=160 in either current Z domain.
- [ ] **COND1 first-bridge budget decomposition:** use active_first_bridge on the original [3sc/4,1] coordinate interval. Split each weighted terminal C0/Z bound into right endpoint, decayed left endpoint, G_y or G_yZ, kernel lambda*G, quadratic remainder, and inherited route. Normalize with the same original Rc A/A_Z and mu; retain the joint divided-row numerator. Report which term and source derivative dominate each of five rows.
- [ ] **COND2 native slow-source conditioning:** extract same-source E/V/a/b/p2 ordinary y,yZ rows from packet.velocity and original signed roots. Bind ordinary axial derivatives with factorials, the original half-shift convention and exactly one Pstar normalization. Retain width correlations, signed p2 and all source/pressure/energy terms. Bounds must not become new field coefficients.
- [ ] **COND3 parameter-uniform mixed primitive bounds:** extend the existing signed-r Mobius/Fourier estimates to genuine fixed-phi y and yZ on the first bridge. Avoid large artificial Poisson-denominator powers, retain implicit inverse cross terms, and handle small-|r| and outer-|r| sectors plus original flat/body/transition cutoff branches. Add a separate adapter; keep accepted older files immutable.
- [ ] **COND4 actual source partitions:** if conditioning still dominates, partition the original first bridge and Z domain using source-defined branch/positive certificates. Keep true microscopic width/coordinate Jacobian and one global phase. Prove every whole cell and include no-gap endpoints; isolated favorable samples do not certify a partition.
- [ ] **COND5 seam cancellation:** compare G on the same physical interface with the same original E/V, a/b/p2, phase and inverse function. Cancel adjacent endpoint terms only after a source-function identity; otherwise keep them. Include the initial flat collar and all16 cross-chart joins. A shared radius alone is not a primitive identity.
- [ ] **COND6 quantitative smaller targets:** recompute only changed source descendants, weight each contribution through its actual remaining route, compare certificates at identical N and over all required Z. Obtain a useful five-row C1 budget; record upper-bound improvements separately from measured actual errors. Keep both direct and averaged certificates whenever their dominance crosses.
- [ ] **CONTROL1 five original controls:** use exact five-bump matrix and Gram functions, same original mu and Ac/S, and equation B_exact(mu)*h+N*r+Q_exact(mu,h)/N=0. Construct function-domain controls; do not select interval midpoints or set a downstream inlet to zero. Keep controls and target repair increment -r distinct.
- [ ] **CONTROL2 terminal identities and admission:** prove five terminal Z identities and first derivatives on the required whole Z interval. Carry controls through original Rc->2Rc geometry, matching and right collar. Actual post-repair2Rc->Rb source/cone admission remains separate.
- [ ] **HIGH/OUTER/ENERGY:** genuine higher slow/fast/physical derivatives and C4 seams, one global finite N and full stress-cone margins, compatible analytic pressure, exact heat exterior, physical tail energy and flat remainder. Current fixed-phi1x1 bounds do not certify these.
- [ ] **REC/WAVE/PHYS:** correct n=1/n>=2 coefficient recovery and independent moment repair, finite-order remainder/smooth sum, two actual oscillatory stress-cancelling pulse families and mean lift, then corrected Cartesian NS and measured scale recursion/material winding.


Scoped gate: `current_original_native_zero_mean_C0_Z_Nminus2_source_integral_covers_executed`. Uniform source-integral C0/Z N^-2 covers and formal mixed-source interface are true. Actual controls/terminal closure, point inverse, higher jets, global finite N/cone, heat/energy, coefficient recursion and full NS remain open.
