# Named patch interfaces and actual finite-N source bounds

Latest successor: [CURRENT_RSH_INTERFACE_AND_CORRELATED_SHEAR_2026_10_06.md](CURRENT_RSH_INTERFACE_AND_CORRELATED_SHEAR_2026_10_06.md); Rsh physical interface [d962130e](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/d962130e08ff98ea27035ac921ac203f5378ad10), correlated primitive shear [c2296cde](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/c2296cde076b6c01559d3e3fc02d59785ac96f1e). CONT4f2a and primitive BOUND1a5a–b are complete. Independent repair, full signed cone/direction conditions, global composition and recursion remain open.

Patch atlas implementation and receipt: [5a4f43af](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/5a4f43afd02b2c0329f79d31afffc3535f79e046). Actual finite-N parameter/density implementation and receipt: [a6b4c2c2](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/a6b4c2c26fd9263eb9ac00e8e000a580ba0df9e5).

**CONT4f6a4a and BOUND1a4b are implemented in their stated local scope.** Six current patch supports have a named physical velocity/absolute-pressure interface provider, with fresh time/constant-viscosity sectors and actual-lambda queries. The finite-N source formulas now consume the actual mu/log-mu interval and original Rd log amplitude, and export all five physical defect-density majorants. These are new callable implementations; the controller's last full runtime stage remains currentmodifiedphysicalvelocity.

## Six named physical interface queries

Use CurrentPatchVelocityInterfaceAtlas from experiments/root_st073/lei_ren_part1_paper_compliant_current_patch_velocity_interface_atlas.py. The interface names are patch_support_49, patch_support_51, patch_support_59, patch_support_61, patch_support_69, and patch_support_71. They retain the original exact x edges49/40,51/40,59/40,61/40,69/40,71/40 and the already proved independent left/right source germs.

interface(name,Z=(-1,1),log_tau=(-3,-1),theta=None,viscosity=1) regenerates216 common physical spatial4/time1 contribution groups, left/right bound views, the common analytic pressure datum, the source family and the fresh bound-definition hash. It accepts any finite log_tau sector and one finite positive constant viscosity. Source Z endpoints are permitted in this enclosure API. Whole-Z primitive uncertainty remains attached even when a narrower operator Z interval is requested; coefficients are not interpolated or fitted.

physical_interface(name,Z,log_tau,theta=None,viscosity=1) requires |Z|<1 and uses q=log(lambda)=(log_tau-log(1-Z^2))/2. Its physical coordinates retain log(r)=log(nu)/2+q+(log(R)+log(2))/2, and signed axial z=Z*exp(log(nu)/2+(1-delta)*q). The unchanged actual physical operator evaluates the requested lambda and viscosity factors. Lambda is kept distinct from sqrt(tau) at nonzero Z. Returned values are source enclosures and derivative bounds, with unresolved coefficient uncertainty.

Evidence:23 exact registry/coordinate/viscosity identities;6x216 common contribution groups; a fresh Z[-.7,.8], log_tau[-12,-11], nu=.8 sector; a physical request at Z=.237, log_tau=-2.337, theta=.337, nu=1.3; an independent signed physical-coordinate implicit equation;35 spatial multiindex viscosity powers;9 invalid request guards;912 working/index dependency hashes. Existing source equality, physical transport and fresh primitive numeric receipts are consumed. The old full33 cached numerical view remains rejected.

Files use prefix lei_ren_part1_paper_compliant_current_patch_velocity_interface_atlas: .py, _check.py, .json, _check.json, _views.json.gz. Run its producer and checker directly. This named provider closes CONT4f6a4a. **CONT4f6a4b, composition into one complete global velocity-interface API, remains open.** The32 adjacent,8 remaining original internal and12 affected current modified germs must retain their own source evidence and be composed before a global gate is admitted.

## Actual finite-N parameter and density adapter

Use CurrentO3ActualParameterMajorants from experiments/root_st073/lei_ren_part1_paper_compliant_current_O3_actual_parameter_majorants.py. bounds(N) accepts a positive Python integer and returns ordinary logR orders0..4 and axial orders0..5 over logR offset[-2,1/2], Z[-1,1]. Its actual scalar/family inputs are checked receipts, without rebuilding the ancestor field.

The current compliant pressure source has Md=40, log(Pstar)=exp(40)+11 and log(mu)=log(.001)-4log(Pstar). The adapter preserves the exact MP endpoint tuples of the current source_mu and log_mu enclosures, checks their defining equations, and takes the exact original J(1)=1/2 from the same pre-pulse endpoint packet. The original common Rd amplitude is retained as

log(Ad)=-exp(40)/2-26/5.

Ad itself is not exponentiated or replaced by a midpoint. The normalized source grid is Fhat[j][k]=T[j]*k!, and the original grid is Ad*Fhat. The finite-N theta, increment and own-axial majorants are linear in that grid, so each keeps Ad as an external positive factor. The finite-N exponential G, its increment bound and the two normalized periodic-loop shear error caps retain their separate mu/N dependence. The microscopic positive increment is computed by H0*exp(H0), avoiding exp(H0)-1 cancellation.

The original and corrected amplitudes remain distinct: E0/Pstar=Ad/(1+Z^2)*exp(-t/2-mu*J(t)), EN/Pstar=(E0/Pstar)*exp(A/N), and the new axial profile is proportional to (E0/Pstar)*sqrt(mu)*chi*cos(2*pi*N*t)/N. Source norms are independent of N. High derivative bounds can grow as N^(j-1); increasing N is not assumed to make every derivative small.

The adapter also supplies absolute two-variable Leibniz majorants for the five actual defect densities before independent repair. If u is the own-axial profile, e the swirl increment and E the original swirl source, the defining densities are M=u, I=e, J=u*(E+e), S=u^2-E*e-e^2/2, and Cp=E*e+e^2/2. They retain the exact physical units:

| Density | Ad power | Pstar power | R power | sqrt(2) power |
|---|---:|---:|---:|---:|
|M|1|1|0|0|
|I|1|1|1/2|1|
|J|2|2|1/2|1|
|S|2|2|0|0|
|Cp|2|2|-1|0|

The selected logCstar is the actual accepted analytic-family enclosure. Absolute log radii remain logRm=logRref-6, logRd=logRref+logPstar, logRw=logRd+1 and logRp=logRw-60log(mu), with logRref=log(110)+10(logCstar+logPstar). The actual original radius source program is bound by hash. For a density factor Ad^a*Pstar^p*R^r*sqrt(2)^h, correlated giant logs are combined algebraically before interval enclosure; derivatives of R^r are included by the ordinary product rule. The output includes finite physical log upper bounds for150 density derivative rows per queried N. These are absolute triangle bounds, not signed stress errors or repaired moments.

Evidence:169 exact original scalar/factor/finite-N identities;30 independent two-variable polynomial product identities;107 independent moderate-mu symbolic cap enclosures;150 fresh physical density rows at N=37; exact J endpoint and nonzero microscopic increment;7 invalid N guards;915 working/index dependencies. Published examples are N=1 and N=10^12. **Neither selects a sufficient common N. N=10^12 remains repair-only for the previous admitted repair.**

Files use prefix lei_ren_part1_paper_compliant_current_O3_actual_parameter_majorants: .py, _check.py, .json, _check.json. Run producer and checker directly. Read-only mathematical review worker: **GPT-5.6 Luna / max**; no worker files changed or ancestor constructors/tests run.

## Next executable tasks and acceptance criteria

- [x] **CONT4f6a4a:** publish the six named patch providers; re-evaluate requested time/nu sectors; use actual lambda and signed physical z at |Z|<1; preserve source uncertainty and current hashes.
- [ ] **CONT4f6a4b:** compose this provider with the existing velocity atlas. Use explicit registry names and aliases, include adjacent/internal/current modified seams, and return an inventory of missing providers. Do not treat copied common bound boxes as a missing source theorem.
- [x] **CONT4f2a:** adapt reshape_reference first. Consume actual_Rsh_source_join.boundary_identities and typed UT/UZ/UR/P source rows; bind the common original radius and analytic P0; transport through the canonical spatial4/time1 operator and expose the two-sided result.
- [ ] **CONT4f2b–4:** add the remaining original inner/reference/restore/pre-pulse adapters from existing receipts; replace axial_buffer, buffer_transition and transition_power with the current modified germs. Preserve exact-zero rows, the typed core axis and all radial units.
- [ ] **CONT4f6b/7:** compose pulse/outer physical providers, map angular_entry to angular_steep explicitly, include4 pulse support and4 angular support interfaces, and cover all12 current affected seams. Admit global smooth velocity only when every required source join and physical trace is covered.
- [x] **BOUND1a4b:** substitute actual mu/logmu, exact J1, original logAd, selected logCstar and original radius/Pstar factors into the finite-N source formulas; keep one query N and external amplitude factors; publish physical five-defect-density derivative bounds.
- [x] **BOUND1a5a (primitive):** retain pointwise chi, chi derivatives, sigma and phase in the finite-N shear errors at both flat support edges. Derive a local edge bound that vanishes with the same source taper rather than using the whole-support absolute cap.
- [x] **BOUND1a5b (primitive):** the successor retains the exact correlated square and proves a quarter taper reserve for all finite N>=1; independent absolute-error/margin ratios are superseded by this joint inequality.
- [ ] **BOUND1a5c:** carry the correlation through independent repair and completed signed stress/pressure/alignment conditions; the primitive floor alone does not admit the whole modified cone.
- [ ] **BOUND1b1:** use the new five density majorants to bound the signed cumulative defects and the independent moment-repair controls as functions of the same N. Bind the actual inverse/Jacobian and avoid choosing controls from midpoint defects.
- [ ] **BOUND1b2:** apply actual bump derivative/support bounds and partial/full cumulative histories, preserving shrinking strips and terminal implicit closure. Carry pressure and radial recovery errors with the same analytic datum.
- [ ] **BOUND2/3:** assemble the completed signed pressure, radial, diagonal and two-vector stress errors, including all cross terms. Compare their whole-domain margins with the actual local loop/taper margins before choosing a common N.
- [ ] **COMMONN/CONE:** solve the combined inequalities for one sufficient finite integer N and certify the whole modified regions. A source bound, a candidate frequency, or the earlier repair-only N does not suffice.
- [ ] **ENERGY:** bound the true kinetic cross terms, the physical spatial Jacobian and time dependence over all regions, including exact heat exterior and infinite radial tail.
- [ ] **REC1–3:** implement the actual n-dependent recovery equations for n=1 and n>=2 on the common core domain; perform an independent five-moment repair at every order; preserve curl/streamfunction divergence; derive finite-order remainder and smooth sum. The modulation frequency N is distinct from recursion order n.
- [ ] **WAVE/PHYS/DYNAMICS:** construct the oscillatory and mean corrections, averaged quadratic stress cancellation and flat remainder; resolve physical u/v/w/p and independently check corrected Cartesian NS residuals; measure radial contraction, relative axial elongation and accumulated material winding.

The persistent paper-faithful reconstruction goal remains active. Global interfaces, completed repair/tensor bounds, sufficient common N/cones, energy, actual n-dependent recursion and corrected NS are not admitted by these local implementations. Unrelated working files are preserved. Predecessors: [CURRENT_PATCH_PHYSICAL_NUMERIC_BOUNDS_2026_10_06.md](CURRENT_PATCH_PHYSICAL_NUMERIC_BOUNDS_2026_10_06.md), [CURRENT_PATCH_SUPPORT_AND_MODULATION_BOUNDS_2026_10_06.md](CURRENT_PATCH_SUPPORT_AND_MODULATION_BOUNDS_2026_10_06.md).
