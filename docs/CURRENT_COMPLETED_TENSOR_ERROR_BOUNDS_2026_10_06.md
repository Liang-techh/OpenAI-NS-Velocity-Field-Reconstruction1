# Completed tensor and leading remainder error bounds for variable N

Latest successor: [CURRENT_MODIFIED_O3_TRANSITION_CONE_2026_10_06.md](CURRENT_MODIFIED_O3_TRANSITION_CONE_2026_10_06.md), implementation [bdbe19b1](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/bdbe19b13b3024c1a42bdd38514cfa258d06763f). Closed modified O3 signed cone is complete for all N>=22; left O2, quiet, common-N and global/recursion gates remain open.

Checked implementation: [ce9268ce](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/ce9268ce120ebabdcf46e1ca64ed7ed82d38d95e). **BOUND2a is implemented as error bounds for the actual changed O2/O3 modulation and quiet repair sources.** It supplies the full recovery/completion inputs needed for signed directional cone estimates. Whole modified cones, a sufficient common N, finite energy and actual coefficient recursion remain open.

## API and current source domains

CurrentCompletedTensorErrorMajorants in experiments/root_st073/lei_ren_part1_paper_compliant_current_O3_completed_tensor_error_majorants.py:

- query(chart,coordinate,N,Z=(-1,1),log_lambda=(-3,-1),theta=None,viscosity='1'). chart='O2' uses the shared offset[-11,0]; 'O3' uses[0,1]; 'quiet' uses local logx[0,1]. Split boxes at the actual O2/O3 seam. Modulation requires integer N>=1; independent repair requires N>=27,303,666.
- at_physical_time(chart,coordinate,N,Z,log_tau,theta=None,viscosity='1') requires strict |Z|<1, finite log_tau and nu>0. It computes the actual log(lambda) before applying the unchanged physical differential operators.

The native closure Z[-1,1] is a formal uniform source bound. Finite positive physical time uses its strict interior. Returned caps/enclosures bound modified source minus the same original source; they are not resolved signed field coefficients. The underlying original source, current mu/delta, axis pressure datum and radius recipes are unchanged. Only checked packets and pure source programs are read; no ancestor constructor is run.

## Full retained recovery and physical completion

The provider replays modified_pre_stress_rows and modified_pre_remainder_sectors, subtracts the unchanged original source as functions, and separates each additional Ad or Ad^2 factor. It supplies 7 theta stress error pieces, 7 axial pieces, and radial/theta/axial remainder pieces 5/1/1. The original Pstar0/1/2 sectors, radius powers, half normalizations, axial beta modes and local/history derivative programs are preserved. Exactly unchanged sectors cancel by source identity.

Original nonzero Utheta, meridional moment Mz and radial velocity remain in cross terms. They are read only from actual_source.current_original_pre_source inside the checked current_modified_pre_stress_views packet. Saved modified N=10^12 controls/values are never used as another N's input. The absolute pressure correction is its cumulative Cp value; post-cutoff pressure and transported moments remain after local profile changes stop. The unique repair equation yields exact zero errors after the final bump.

The formal_Ad symbol records an additional external amplitude factor; it is not a perturbation parameter. Raw original fields already contain their own original amplitude. Complete normalized dk/de/dp histories carry an external Ad^2, including their retained baseline-swirl cross terms. Cross terms with a raw baseline keep that raw factor and apply only the additional Ad factor. All huge positive scales remain logarithms.

Saved physical moment derivatives Mz/current_R contain (Dy+1)^j m; saved radial derivatives contain (Dy+1/2)^j Q. The loader undoes these shifts before treating m/Q as normalized source functions. Recovered radial error caps already include the +1/2 shift; a positive inverse-shift bound restores normalized Q jets. Original stress and remainder programs then apply their own radius rates exactly once. The source theorem verifies the same derivative rows, including the -1/2 radial-viscosity rate.

Every differentiated expression is collected into a differential polynomial of the current mixed4 inputs. Absolute rational geometry coefficient enclosures and positive monomial products bound each native stress mixed3 and remainder mixed2 grid. The original linear lift then supplies:

- cylindrical off-diagonal stress mixed3 and divergence mixed2;
- completed Ttheta_theta=r*partial_z(Trz) mixed2;
- symmetric Cartesian xx/xy/xz/yy/yz/zz tensor error sectors;
- three-component remainder and momentum error decomposition -div(delta T)+delta E;
- exact zero completed radial tensor divergence.

Only the remainder input and its additional Ad log factor are adapted in the lift AST. All derivative, completed-radius, divergence and Cartesian formulas are unchanged. Scalar absolute values run at precision above the interval context so positive endpoint magnitudes do not round downward at the unrelated default mpmath precision.

## Actual lambda/time relation

The coordinate theorem gives lambda^2*(1-Z^2)=tau. Thus loglambda=(logtau-log(1-Z^2))/2. The legacy tensor lift internally computes gamma*mapper_log_tau/2; the new provider passes 2*loglambda and restores the requested physical-time metadata afterward. No artificial 1-Z^2 denominator or duplicated lambda factor is introduced. At Z=0 this agrees with logtau/2. Boxes touching Z=+/-1 cannot use the finite physical-time wrapper. The separate global velocity log_row convention expects loglambda itself and must not receive 2*loglambda.

## Focused evidence and files

The checker passes 79 exact source/derivative/lambda identities, 10 independent signed original-program fixtures, 1720 native mixed comparisons, and 1720 original physical-derivative comparisons. It checks 90 baseline shift round trips, 1862 physical log rows, 7 current queries and 8 invalid guards. Signed numerical fixtures allow only a working-precision roundoff floor for subtracting two baseline enclosures; production source caps are not enlarged. Exact support/repair zeros and surviving post-cutoff cumulative pressure are retained. Read-only source review: GPT-5.6 Luna / max.

Files use prefix lei_ren_part1_paper_compliant_current_O3_completed_tensor_error_majorants: .py, _check.py, .json, _check.json and _views.json.gz. The compressed views hold complete current queries and are hash-bound by the source report/receipt. Run producer and checker directly. The full goal remains active; controller's last full runtime stage remains currentmodifiedphysicalvelocity.

## Detailed next tasks for GitHub agents

- [x] **BOUND2a1:** replay/subtract the actual signed recovery operators, preserving original nonzero radial and meridional histories and all baseline/correction and correction-square terms.
- [x] **BOUND2a2:** correct native derivative units and publish mixed3 stress/mixed2 remainder bounds, same cumulative pressure datum and variable-N implicit repair family.
- [x] **BOUND2a3:** map all error sectors through the completed physical tensor, diagonal/divergence/Cartesian/residual formulas and exact lambda/time relation.
- [ ] **BOUND2b1:** derive current signed a, b_s, vs-2, D and J from the complete modified source. Use b_s=2*Uz_y/Utheta with the actual physical Pstar sectors. Preserve the accepted correlated primitive identity and both signs of b_s. Error magnitudes alone do not certify orientation.
- [ ] **BOUND2b2:** bound D=Ttheta-(b_s/a)*Tz and J=Tz+(b_s/a)*Ttheta, then Q=2*D^2-(vs-2)*J^2, using all retained pressure/energy/incoming moment terms. Supply separate a>0, vs-2>0, D>0 and Q>0 inequalities, with source-normalized units and domains.
- [ ] **BOUND2b3:** handle both flat modulation tapers with actual chi, chi_y and N*offset phase dependence. The primitive floor vanishes at support edges; preserve cancellations before enclosure and never divide by a zero cutoff. Prove edge equalities/limits as source identities, not a phase grid.
- [ ] **BOUND2c1:** specialize the completed error sectors on each disjoint quiet bump and on all gaps. Retain partial/complement history strips, cumulative absolute pressure, original radial cross terms and C versus C^2 axial derivatives. Use the same unique finite-N implicit control vector.
- [ ] **BOUND2c2:** compare the signed completed tensor directions/alignment with the quiet primitive reserve a-2>=3mu/2. The primitive threshold68,533,403 is only an input; derive all additional tensor constraints and keep the original pressure/energy/radius/amplitude factors.
- [ ] **BOUND2c3:** prove entry, internal bump and terminal/source joins for the signed cone estimates. Flat terminal source equality does not establish strict positivity in its preceding strip.
- [ ] **BOUND3 / COMMONN:** combine repair, scalar/history, pressure/radial, diagonal, orientation and alignment constraints into one sufficient finite N. Publish an inequality ledger with constants, domains and unchanged source hashes; neither existing repair nor primitive threshold alone completes this task.
- [ ] **CONT4f2b:** implement reference_restore/restore_buffer physical interface providers from accepted source germs and original radii, then complete the remaining adjacent/internal providers and compose one global API. Never relabel obsolete full33 numerical views. Use each mapper's documented lambda convention.
- [ ] **ENERGY:** integrate physical kinetic cross terms with the actual volume Jacobian and radial/time tails of the revised source. Terminal energy-moment closure does not erase positive kinetic density.
- [ ] **REC:** implement actual n=1 and n>=2 recovery equations, a common inner domain, independent per-order moment repairs, divergence-preserving truncation, finite-order remainder control and smooth summation. This work has not implemented that recursion.
- [ ] **WAVE / PHYS / DYNAMICS:** two oscillatory/mean correction families, averaged quadratic stress cancellation, resolved u/v/w/p, corrected Cartesian NS residual, measured contraction/aspect ratios and material winding.

On completion of each task, mark only its own checkbox, publish source/receipt hashes and remaining gates, and update the current top-level handoff. Preserve unrelated dirty files and historical sections. Previous report: [CURRENT_RADIAL_PRESSURE_ERROR_BOUNDS_2026_10_06.md](CURRENT_RADIAL_PRESSURE_ERROR_BOUNDS_2026_10_06.md).
