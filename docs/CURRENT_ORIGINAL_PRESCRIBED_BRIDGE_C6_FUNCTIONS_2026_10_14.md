# True prescribed background bridge C6 functions, 2026-10-14

The genuine rho4 inlet is now propagated through the two microscopic smoothing windows and full macro bridge to R100 by original integral functions. The known comparison data has C7 functions, the prescribed background field and all six own moment histories have C6 functions, and full Q/inertial/p2 has C5 functions. This is the prescribed background before the finite-N loop modification. Modified fields, R100-to-Rm propagation and all17-chart sources remain unfinished. Full goal **ACTIVE / INCOMPLETE**.

Inspect `python scripts/run_original_prescribed_bridge_C6_functions.py`. It reads source-pinned archives and does not rebuild ancestors. On the original host reuse session55770/cache26, `public_inlet6`, `bridge6_owner`, `bridge6_receipt`, `hier` and existing lower/pulse/terminal owners. The [preceding inlet milestone](CURRENT_ORIGINAL_CORE_INLET_AXIAL6_GENERIC_SOURCE_2026_10_13.md) remains valid.

The extra order is necessary: a C6 prescribed bridge field uses a C6 known comparison direction, whose first axial derivatives consume C7 mean/energy/P0 input functions. The same unique twenty-term analytic core, root-anchored G/F0 and six true core moment kernels provide that seventh row. Earlier C6 handles remain. P0_Z7 is an explicit derivative of the original independent fourteen-atom pressure function; its executable projection uses factorial5040 once. These function definitions do not select numerical enclosure endpoints or certify new C7 numerical mixed norms.

The comparison smoothing is exactly the original one:

```text
alpha(y)=1-sigma((y-hb)/hb)
log(phi_bar(y))=log(Phi_core(4,Z))
                  +integral_0^y alpha(t)*d_t log(Phi_core(4exp(t),Z)) dt
V_bar(y)=V_core(4,Z)+integral_0^y alpha(t)*d_t V_core(4exp(t),Z) dt
```

Only the comparison profiles freeze after y=2hb. Their own six histories keep evolving by the original history ODEs. The known direction is built from those comparison histories; actual histories never replace it.

The actual field keeps the original coupled controls and full physical source sectors:

```text
chi(y)=1-(1-hb)*sigma(y/hb)
phi=Phi_core(4,Z)*exp(-1/2*integral chi*R*(Dbar/R) dy)
V=V_in-integral chi*(phi/phi_bar)
                    *(R*hydro+R*Pstar^2*pressure+R^2*F0^2*swirl) dy
```

The first microscopic chart has chi=(1-sigma)+hb*sigma; the second and macro charts have chi=hb by the original flat sigma predicate. In micro coordinates y=hb*s, the measure is hb*ds and R=Ra*exp(hb*s). In the macro chart:

```text
Y=log(100/Ra), length=Y-2hb
y=2hb+length*fraction, dy=length*d_fraction
R=Ra*exp(y), fraction in[0,1]
```

The true nonzero micro endpoint phi/V and all six actual moment functions are passed directly into the macro field and history functions. Endpoint R100 functions preserve the finite physical domain.

Each actual history has its complete Duhamel body:

```text
X(y)=exp(-rate*(y-y0))*X(y0)
       +integral_y0^y exp(-rate*(y-t))*density(t,Z) dt
(rate,density): H=(2,2phi), M=(1,V), K=(2,2phi V)
                A=(1,V^2), B=(2,phi^2), C=(1,phi^2)
```

The B forcing is phi^2; its constant-profile equilibrium is phi^2/2. No incoming history is reset to zero. The modified-loop correction's separate zero inlet is not substituted for these background memories.

All physical F0 derivatives are dressed before differentiating the comparison direction. The graph's normalized swirl drive divides by full F0(Z)^2; physical reconstruction multiplies that same function once. The independent check compares the full physical sector to the original Taylor producer, whose F0(base)^2 is a fixed external scale. F0 itself depends only on Z. A lexical source-DAG check proves that F0 and all bridge phase bodies contain no free core rho after the Phi/Psi substitutions.

Generic E/V/m/h/k/e/p/P0 normalization is the same as the inlet adapter, now at the current physical radius. Full signed Q and p2 through Z5 use genuine C6 source rows, all pressure terms, one Pstar quadratic weight and one radius factor. Radial E_y and V_y through C6 follow the original actual field ODEs.

Focused validation covers actual source/coordinate bindings, original comparison-direction Taylor equations with nonzero seventh input rows, independent exp/log calculus, signed physical axial transport with nonunit scales, and all six C6 history ODE/incoming identities. The result proves function construction. Numerical C6 bridge containment, derivative norms, the finite-N modified branch and global recursion remain separate requirements.

Completed:

- [x] Extend actual analytic core/moment/P0 source functions to C7 where the comparison direction consumes an order.
- [x] Construct original true comparison smoothing and six complete histories through C7.
- [x] Preserve frozen comparison macro profiles and continuing comparison histories.
- [x] Construct actual coupled prescribed micro/macro phi/V through C6 with signed full physical drives.
- [x] Carry six actual C6 histories and nonzero micro-to-macro memory to R100.
- [x] Recover current normalized generic fields, radial E_y/V_y C6 and full Q/inertial/p2 C5.
- [x] Save explicit R100 endpoint functions, source acceptance and an archive inspection entry point.

Next tasks:

1. [ ] Transfer same-source quantitative C6/C7 enclosures through the actual comparison and prescribed field integrals. Use complete analytic tails and original physical log bases; existing comparison axial6 boxes do not define actual C6 fields.
2. [ ] Propagate this prescribed background through the actual R100 switch to R110, reference extension/restoration and long reshape to Rm, retaining pressure datum, all six incoming histories and original coordinate/radius maps.
3. [ ] Implement versioned Rm axial6 source exports and six-row Q/inertial/p2 axial5 outputs. Preserve the old axial4 outputs and never promote Qpad.
4. [ ] Extend the actual loop/inverse/A/B and finite-N five-density functions through C5/C6 as required. Preserve original flat support predicates, phase origins, selected N and positive Jacobian. Compose the modified velocity and signed correction-history transport with this true background.
5. [ ] Complete all17-chart actual sources, normalized C4/C5 targets and source-derived D4/D5 bounds. Extend the same compact-control branch with the original Jacobian/N and factors4,3 and5,10.
6. [ ] Build full pre-Rp pressure C3, complete five C3 defects and actual reserved coefficient derivatives; recover repaired mean_ZZ and Q/V Z1 with weight1+delta.
7. [ ] Complete required mixed/radial norms, global joins, compatible preheat/finite-energy/exact heat matching, admissible stress and flat remainder.
8. [ ] Implement true n>=2 coefficient equations with independent moment repair, finite remainders/smooth sum, oscillatory averaged-stress cancellation and independent corrected Cartesian NS Linf/volume-L2 validation below1e-3.

Artifacts: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py`, same stem `.json.gz`, `_check.py`, `_check.json`, and the inspection script above. This is construction progress toward the original goal; it does not yet establish genuine global scale recursion.
