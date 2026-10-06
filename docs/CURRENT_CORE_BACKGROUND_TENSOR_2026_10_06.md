# Positive-radius core full T/E and completed core/first trace (2026-10-06)

The same admitted nonlinear fixed point now supplies the full physical core background stress, its completed tensor/divergence and all three NS remainder components on **0<rho<=4**. The original integrated equations cancel the **total** leading core stress exactly. Every original signed sector is kept in a separate ledger; the nonzero remainder is retained.

The completed tensor inventory is **33 positive-radius regions, 32 adjacent traces and 10 internal traces**. The separate primitive velocity/pressure atlas stays14/8. This milestone excludes the physical axis, global admissibility, independently bounded time-flat remainder, required-domain energy, resolved point fields and actual n-dependent coefficient recursion. Older status snapshots below the current handoff are historical.

## Implementation and use

- Operator: `experiments/root_st073/lei_ren_part1_paper_compliant_current_core_stress_operator.py`.
- Producer/data: `experiments/root_st073/lei_ren_part1_paper_compliant_current_core_background_tensor.py` and complete deterministic `.json.gz`.
- Independent checker/receipt: matching `_check.py` and `_check.json`.
- Focused controller: `python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage currentcoretensor`.
- Warm reuse: `CurrentCoreBackgroundTensor(interior=checked_current_core_interior_moments)`.
- `field.chart('core', Z, rho, log_tau, theta, viscosity)` accepts finite Z[-1,1], any requested positive compact rho sector up to4, finite log tau and nu>0.
- `field.interface('core_bridge', ...)` gives71 source-identified common completed tensor contributions at rho4 / first phase0.

Do not rebuild historical tensor producers for a focused continuation. Reuse the checked bridge/interior/common/first source graph. Do not select interval representatives as physical point values, replace the nonlinear field by its finite polynomial, or treat the published compact rho[.001,4] as an imposed physical cutoff. A fresh rho.000137 sector below that interval is included; arbitrary positive compact sectors remain callable.

## Source units and exact cancellation

R=epsilon_core*rho, F=F0*Phi, V=Uz, D_y=rho*d_rho. Define unit=sqrt(2R)*F0base/Pstar. The direct core adapter avoids the microscopic bridge's inverse phase-width conversion:

- u_j=unit*(D_y+1/2)^j dress(Phi).
- h_j=unit*(D_y+1/2)^j dress(H)/2 and k_j=unit*(D_y+1/2)^j dress(K)/2.
- m_j=D_y^j M.
- e_j=Pstar^-2*D_y^j A - unit^2*(D_y+1)^j dress(B)/2.
- p_j=unit^2*(D_y+1)^j dress(C)/2; absolute pressure includes the separate normalized original P0.
- Radial velocity uses (D_y+1/2)^j Q with the original external sqrt(R/2) factor.

Full source Phi/Uz rows come from the unchanged rebuild.profile; six moments, all product tails, F0 ratios, P0 and Q come from the same checked interior source. P0 Taylor coefficients are not ordinary derivatives. The four fixed factor logs enclose the same original factors and never define their values. The original F0 enclosure expression is AST-bound to core_physical_field.profiles.

The unchanged raw stress operator is replayed against the integrated original core equations, including pressure and swirl units, FTC and the analytic axis integration constants. It yields T_rtheta=Rtheta/(L*R) and T_rz=Rz/(L*sqrt(2R)); the accepted nonlinear core makes both integrated residuals identically zero. This cancels only the total. The original eleven signed source sectors and their mixed3 enclosures are exported under `original_uncancelled_signed_stress_sector_enclosures`; a separate `integrated_core_total_stress_zero` aggregate is sent through the unchanged linear physical completion.

All six remainder sectors remain: radial time, radial viscosity, nonlinear transport and axial viscosity, plus theta and axial axial-viscosity terms. The full remainder is not zero and is not yet an independent time-flatness estimate.

## Core-to-first completed attachment

The same CurrentCoreFirstInterface common owner, true rho4 atom normalizations, original scaled primitives, original pressure identity, exact positive scalar hb pullback and phase0 ODE uniqueness identify all mixed source functions through order4. The core and first bridge then use identical full stress, remainder and physical completion operators. The resulting tensor/remainder traces agree as functions before triangle bounds are formed. Numerical interval overlap is not an identity proof.

## Targeted evidence

- 135 arbitrary smooth-source raw unit/mixed derivative identities and20 full stress-to-integrated-equation identities.
- Five complete compact/inlet/fresh views,555 physical contribution rows and71 common completed tensor contributions at the core/first seam.
- All eleven signed stress sectors preserved, only the total reduced; all six NS remainder sectors retained.
- Producer source hashes and complete output match independent replay after clearing only the interior adapter's local caches.
- Fresh low-radius sector below the published compact sector; fresh seam time/viscosity sector; invalid/axis/foreign source requests rejected.
- Focused producer/checker/controller and four changed Python compilation checks pass. Working/index source-hash and staged whitespace checks precede publication.
- Read-only source/math review used GPT-5.6 Luna/max; no material blocker for this positive-radius scope.
- Three regional gates become true. `current_core_full_background_tensor_available` remains false because physical axis tensor/remainder admission is unfinished.

## Completed tasks

- [x] F57C5a-core5/core6: full original core raw units, source-bound integrated total stress and original signed cancellation ledger.
- [x] F57C5a-core7: original full physical remainder with all six source sectors.
- [x] F57C5b-core9: completed same-source core/first tensor and remainder trace.
- [x] F57C5a-core11/core12-positive: positive-radius API, compact/inlet/fresh data, independent receipt and focused controller. Axis portion remains under the tasks below.

## Immediate axis tasks: F57C4e-core-axis10

- [ ] axis1: use the same checked CurrentCoreBackgroundTensor and its interior/common owner. Pin the family, selected source, datum, epsilon, delta and all original callables. Create a separate nonsingular axis adapter; retain rejection of rho0 in the current cylindrical chart.
- [ ] axis2: derive the original radial time/viscosity/nonlinear and axial-viscosity expressions after Ur=sqrt(R/2)*Q. Cancel apparent inverse-R terms symbolically before taking bounds; do not plug rho0 into the current inverse-radius lift.
- [ ] axis3: recover Phi, Uz, Q and all needed radial/axial derivatives at rho0 from the unchanged infinite fixed-point rebuild. Include density/model/correction/product tails. Use the normalized moment-axis identities only for their actual primitive source contributions.
- [ ] axis4: construct Cartesian E_x,E_y,E_z directly from smooth transverse-coordinate source expressions, retaining parity and nonzero transverse derivatives. A zero cylindrical radial/swirl value on the axis does not imply all Cartesian derivatives vanish.
- [ ] axis5: compute full Cartesian stress/completion/divergence and remainder through mixed order2 on the axis. Justify the continuous extension of zero total core stress from the exact source equations; preserve the separate signed-sector ledger.
- [ ] axis6: prove equality of the nonsingular axis formulas and the current positive-radius formulas on rho>0 with arbitrary source jets. Check a near-axis sector and fresh axial/time/nu requests without treating small-radius samples as a limit proof.
- [ ] axis7: publish an independent producer/checker and a focused controller stage. Admit axis/full-core gates only after these source proofs and full derivative/tail checks pass. Do not add another radial region to33 merely for a coordinate extension.

## Remaining ordered work

- [ ] F57C-angular-full1: complete the four angular internal support joins using actual whole tensor/remainder functions. Existing local difference and primitive source gates do not admit completed tensor traces. Preserve original signed memory and complete pressure/energy histories.
- [ ] F57C-global-cover1: assemble all regional tensor and axis APIs on one family/source/datum, with a chart/interface coverage table and no unhandled coordinate gaps. State compact positive-time and radial scope separately from global terminal-time bounds.
- [ ] F57C-cone1: evaluate the original admissible stress cone and margins throughout the inner exit, bridges, pulses, flatten, angular supports and heat collar. Use full stress, radial/shear factors and exact source amplitudes; do not infer a cone from zero core/exterior stress.
- [ ] F57C-flat-time1: derive independent time-dependent bounds for every E component and required derivative. Prove the claimed flat decay as t approaches T, rather than copying fixed-time spatial rows or asserting flatness from a scale formula.
- [ ] F57C-points1: recover convergent numerical representations of the actual shared nonlinear core and cumulative histories at specified points, with truncation/error bounds. Return u(x,y,z,t), v(x,y,z,t), w(x,y,z,t), p and their source coordinate mapping without interval midpoint substitution.
- [ ] F57C-energy1: state and integrate the prescribed physical energy domain, retaining core and full heat tails and any required localization. Existing local compact-time bounds do not certify global/terminal finite energy; the unlocalized whole-space source must not be relabelled finite-energy.
- [ ] F57D-recursion1: implement the original n=1 recovery equations and per-order moment repair on the same core interval, including error controls and pressure/history transport.
- [ ] F57D-recursion2: implement the distinct n>=2 coefficient equations with their actual n-dependent operators; verify at least two nontrivial orders. Repeated coordinate scaling of a leading field is not coefficient recursion.
- [ ] F57D-recursion3: construct finite-order remainder estimates and the required smooth sum, with localization/cutoff at streamfunction or vector-potential level before curl.
- [ ] F57E-oscillation1: construct both oscillatory pulse families and mean corrections; independently compute their averaged quadratic flux and compare/cancel the admitted background stress.
- [ ] F57E-corrected1: restore complete corrected NS momentum residual in independent Cartesian coordinates. Apply the full corrected L-infinity/L2 target only at this stage.
- [ ] F57F-dynamics1: measure radial contraction, relative axial slenderness, swirl/vorticity amplification and actual material-line winding across times. Fit genuine recursive coefficients separately from coordinate-based geometry.

Each task should produce its implementation, complete source evidence, an independent scoped receipt, a task-status update and a commit. Mark only its achieved scope done. The long-term goal remains active.
