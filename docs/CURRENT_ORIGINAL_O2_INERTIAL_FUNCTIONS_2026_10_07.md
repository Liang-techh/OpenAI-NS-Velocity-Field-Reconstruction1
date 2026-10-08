# Original O2 full inertial functions and exact pressure Z jets

Checked source: [8c619506](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/8c6195065a5515119d976e848be8f00d7206690a). Predecessor: [original O2 slope point profiles](CURRENT_ORIGINAL_O2_SLOPE_POINT_PROFILES_2026_10_07.md). Focused independent checker PASS; 607 direct dependency hashes audited against the Git index. Producer 14.672s; checker 28.422s. Existing read-only reviewer: **GPT-5.6 Luna / max**; no new workers.

The O2 point-profile service now feeds the full original inertial inputs p1,p2 and their ordinary Z derivatives. The original absolute pressure is bound to the exact fourteen-stage preheat integral, its prescribed raw waiting root and its first two Z derivatives. This supplies source-function expressions needed by the Section 11 loop and the five-history correction route.

Radial quadrature coefficients remain approximate. Original radius and scale parameters remain formal; exact pressure integrals have not been numerically evaluated at the native parameters. The exact conditional C1 control family remains available, but native phase/integral evaluation, certified numerical errors, installed five controls, global common N, genuine coefficient recursion and corrected NS velocity are still open.

## Source equations and units

At fixed original O2 coordinate y in[0,1], reuse the original f,H,D,P defining-integral functions from the predecessor. Put C=1/(1+Z²) and E=C*f=Utheta/Pstar. In this replay the physical axial velocity is Vraw=4Z. The physical stress program receives the original moments with their correct radius/Pstar factors:

```text
h=C*H; k=4Z*C*H
e=16Z²/Pstar²-C²*D; cumulative_p=C²*P
m_theta=sqrt(2)*R^(3/2)*Pstar*h
m_z=R*(4Z)
m_theta_z=sqrt(2)*R^(3/2)*Pstar*k
m_z_theta=R*Pstar²*e
absolute_pressure=Pstar²*(normalized_P0+cumulative_p)
F=Pstar*E/sqrt(2R)
p1=Itheta/F; p2=Iz/F.
```

This replay's m and k use Vraw. The predecessor common-basis point service uses V=Vraw/Pstar and k/Pstar; its basis has not changed. The new template lifts the original axial modes before division, so p2 retains both Pstar inverse and Pstar positive contributions. The positive axial energy term, negative angular energy term, transport, mixed histories, nonlinear meridional terms and absolute pressure all remain. Only the two named radial shear sectors are excluded from I/F; using total stress T/F here would be incorrect.

The full original raw stress AST is replayed unchanged. The generated formulas are independently compared against the original physical I_theta/I_z program, including both ordinary Z rows. p2_Z includes the second derivative of the original pressure datum; pressure cannot be reset to zero or replaced by a range cap.

## Exact pressure function service

The pressure returned here is normalized_P0=P0/Pstar², with H=1 in the original preheat definition:

```text
normalized_P0(Z) = -sum_stage integral (Utheta_pre/Pstar)^2/2 dy
normalized_P0_Z  = -sum_stage integral partial_Z(density) dy
normalized_P0_ZZ = -sum_stage integral partial_ZZ(density) dy.
```

All fourteen original integration bounds and the original uncorrected waiting root are Z independent. Six stages carry the (1+Z²)^-2 factor, seven have Z-independent densities, and the flatten stage retains its original variable exponent. The source-defined waiting root includes the actual incoming angular history and flatten integral; it is not a selected endpoint of the saved waiting box. Existing same-family pressure identification/differentiation receipts bind the exact generator. No m2/m0 mass cap, datum enclosure midpoint or zero-pressure replacement is consumed as a scalar function.

Pressure and inertial formulas share one `original_delta` symbol, including the substituted raw waiting root. Other source scales and stage lengths remain explicit formal parameters pending the original parameter-contract task. In particular the radius relation R=Rref*exp(y), yd=log(Pstar), and all original scale/length dependencies must be installed together before claiming native scalar values.

## API and focused evidence

Implementation/result: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_O2_inertial_point_functions.py/.json`. Independent checker/receipt: same stem plus `_check.py/.json`. Gate `original_O2_full_inertial_point_function_templates_and_exact_P0_Z_jets_connected`.

`OriginalO2InertialPointFunctions(dps=50).functions(y)` returns p1,p2 as value/Z expression pairs. Default `resolve_pressure=True` eliminates the unresolved P0 function placeholders and binds the exact original integral jets. It also returns E,V,a,b and formal R/Pstar/delta. The actual O2 shear inputs returned are a=4/5+(6/5)*sigma(y), b=0; use these returned values rather than the generic template's free positive a symbol. The optional compact form leaves P0 symbols only for algebraic serialization/reference probes. `.pressure_at(Z,order)` returns exact integral expressions at order0,1,2; it is not a numeric oracle.

The focused checker establishes four exact identities against the separate original physical stress program: p1,p2 and both Z derivatives. It independently checks 28 density derivative identities, all fixed limits, raw-root preservation, the one-delta binding and ordinary differentiation of the full integral operator. Different placements of Z-only factors inside/outside an integral are compared by linearity. Actual returned p2/p2_Z contain original integrals and no unresolved P0 function placeholder.

16 finite-difference comparisons use a clearly manufactured pressure probe and two finite R/Pstar/delta unit choices. These test algebra and Z differentiation only. Maximum reference derivative discrepancy is 5.49188165695201670312432764168e-14; it does not certify original quadrature/roundoff or native pressure accuracy. Seven invalid pressure domains/jet orders/radial coordinates are rejected. No source ancestor constructor is rebuilt.

## Executable next tasks

- [x] **O2-PROFILES:** original E/V and five radial-history point coefficients, same-family pressure reference, original four background graph roles and ordinary Z rows.
- [x] **O2-INERTIAL-EXPRESSION:** full original p1,p2 value/Z expression pairs, correctly lifted axial modes and independent physical I/F identities.
- [x] **O2-PRESSURE-EXPRESSION:** original fourteen-stage normalized P0/P0_Z/P0_ZZ integral expressions, prescribed raw waiting root and same delta binding.
- [ ] **NEXT O2-SOURCE-PARAMETER-CONTRACT:** create one source-bound exact/factored parameter frame shared by the profile, pressure, stress and phase services. Bind original Rref and R=Rref*exp(y); logPstar=exp(Md)+11, yd=logPstar, original mu/delta/epsilon and all stage lengths. Bind existing accepted source definitions, not interval endpoints. Preserve correlations and tiny positive quantities without astronomical exponentials. Publish one frame/receipt with each dependency and selected source definition explicit.
- [ ] **O2-PRESSURE-NUMERICAL-JETS:** evaluate actual defining pressure integrals through Z order2 in the parameter frame. Split original fixed stages; keep flatten dependence and tail terms. Return separate quadrature/roundoff/tail errors. Never pass the existing pressure enclosure as a point value. Report unresolved branches/conditioning instead of clipping or replacing small positive contributions by zero.
- [ ] **O2-FULL-POINT-INPUT:** combine the frame, original radial coefficients and pressure jets into a conditioned value/Z object for E,V,a,b,p1,p2. Preserve common units and radius correlations; make casting to a scalar fail while any required source factor remains unresolved. Bind to the accepted all-N original source roles before allowing a full-oracle mode.
- [ ] **O2-LOOP-SCALES:** bind eta and d_star from the original whole-input cone, p1/p2 and boundary margin contract. Distinguish certified conservative bounds used to choose loop scales from actual field point values. Preserve microscopic kappa excess and phase denominators using factored arithmetic.
- [ ] **O2-CONDITIONED-PHASE:** adapt the existing `GenericLoopPointZ` backend to the native conditioned source inputs. Use the original global fractional phase and implicit inverse; retain a_Z,b_Z,p2_Z,E_Z. Track input, angle-inversion and primitive-quadrature errors. The existing finite-input backend is available, but native scale/phase installation is still open.
- [ ] **O2-ONE-SIGNED-INTEGRAL:** produce one genuine original correction density and its own-rate cell integral, including its Z row and all linear/cross/quadratic terms. Use the exact five-density order formulas already implemented; their coefficients remain N dependent. Keep nonzero incoming histories and publish actual value/error outputs, not rescaled fixed-N caps.
- [ ] **ORACLE-ALL-CHARTS:** extend the stateless source recipes to the remaining fixed original C1 cells. Reuse source-bound histories and cached constant coefficients; honor signed origins, chart interfaces and the common fast phase. Point evaluation must not rebuild the complete ancestor hierarchy.
- [ ] **CONTROL-CERTIFIED-INSTALLATION:** combine coefficient, phase, integration and roundoff error bounds with the checked exact C1 control-map stability/tail. Complete whole-frequency constraints and choose one compatible common N. Install all five original controls and terminal functional error bounds; finite approximate Picard iterates alone do not give exact zero terminal defects.
- [ ] **FIELD-JOINS-OUTER:** install corrected field and pressure functions, admit required spatial joins/mixed derivative rows, complete outer stress-cone/heat/energy compatibility and independently bound the full remainder.
- [ ] **REC-COEFFICIENTS:** implement genuine n-dependent coefficients, independent per-level repair and flat summation. Existing conditional control identities are not a recursively produced sequence.
- [ ] **WAVE-CANCELLATION:** construct the two original oscillatory families and their averaged quadratic stress cancellation; prove the required remainder/forcing regularity.
- [ ] **PHYS-UVW:** expose corrected Cartesian u(x,y,z,t),v(x,y,z,t),w(x,y,z,t), and measure scale-to-scale recurrence, core contraction, relative axial elongation, material winding and complete momentum residual. Animation has lower priority.

Long-term goal remains active. This checkpoint completes the O2 source-expression layer; it does not certify an actual recursive solution or final physical velocity field.
