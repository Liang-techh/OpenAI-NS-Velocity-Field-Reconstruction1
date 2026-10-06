# Original cone map and whole current pulse-end signed cone

Implementation and scoped original-cone/current pulse-end receipt: commit [fde7e10e](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/fde7e10ef112f4dbe915e896186c656dab2e922d).

This milestone continues the source-bound global physical T/E cover; it does not complete the Navier--Stokes reconstruction.

## Implemented result

The exact original cone now has a source-bound operator and a current-graph signed margin adapter for pulse main, exit and end. Every source sector is retained: 15 order-zero signed contributions for main/exit, 10 for end. The full diagonal, divergence, Cartesian tensor and remainder records stay attached to each view.

Four fresh current main/exit/end requests and the whole pulse-end box s[-4,0], Z[-1,1] pass the strict two-vector cone. The whole pulse-end proof is an interval bound on the complete source box, not a collection of phase samples. The positive common physical factor cancels, so this two-vector sign result is independent of positive lambda and constant nu. This does not prove the whole main/exit domains, the completed tensor condition or a global wave lift.

Five complete views, their signed source rows, shear jets and reference coefficient bounds are preserved in `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_cone_views.json.gz`. Read it with `json.loads(gzip.decompress(path.read_bytes()))`. Producer: `_current_original_cone.json`; receipt: `_current_original_cone_check.json`. Both use the common `lei_ren_part1_paper_compliant` prefix. The receipt verifies the gzip's current SHA256.

## Source map and signs

Original OpenAI: X=R, eta=Z, q=lambda^2, h=delta/2, E=Utheta, U=Uz, F=Utheta/sqrt(2R). The integrated inviscid pair is p_s=I/F, not physical pressure. The shear variables are a=-S_theta/F, b_s=S_z/F, t_s=-b_s/a, v_s=a+b_s^2/a=Lei kappa. T=I+S=F*(p_s-s), s=(a,-b_s). The residual sign stays -div(T)+E.

With D=Ttheta-(b_s/a)*Tz and J=Tz+(b_s/a)*Ttheta, the strict original test is a>0, v_s-2>0, D>0 and 2D^2-(v_s-2)J^2>0. Zero stress is excluded from this interior test; core/exterior zero stress and limiting edge directions need separate treatment. A box straddling a required sign is inconclusive, not a demonstrated failure.

Current physical units: u=sqrt(nu)*lambda^beta*U; p=nu*lambda^(-2-2delta)*P; T=nu*lambda^(-2-delta)*Tprofile. The common positive factor leaves the cone unchanged. The adapter binds the actual source shear and physical-lift ASTs and uses the same velocity packet that generated the signed tensor.

For pulse fields, a-2=2mu is retained directly rather than subtracting 2 from a rounded number. End b_s uses the original positive D factor and full Bhat derivatives; main/exit uses its full local axial derivative, with incoming meridional stress sectors still present.

## Numerical construction changes

The common positive sqrt(R/2)*B source is canceled by subtracting its exact powers from each sector mode before interval evaluation. Huge log intervals are never subtracted from one another to claim cancellation.

A fixed exp(-1024) cap can swamp the astronomically small mu after multiplication by a large coefficient. The new cap uses log(mu)-log(max_abs_coefficient)-1000 when the exact relative source log is more negative. It remains an enclosure of the original exponential. Coefficient signs are retained before summation.

The source equilibrium has the exact positive rewrite

`(C*(mu-delta/2)/(1-mu)+(1-delta)*Z^2*C^2/(1-mu))/(1-delta*Z^2)`, `C=1/(1+Z^2)`.

Its direct Z^2 enclosure intersects the original Taylor coefficient enclosure. This removes the lost correlation from multiplying a Z interval by itself while preserving the original source function and all original records. It turns the former whole-end inconclusive box into a strict cone bound.

## Actual lift versus algebra

The original lift is Proposition 7.5: H's columns are the two actual homogeneous pulse covariances, y=H^-1*T0_star must be positive, and a_sigma=sqrt(y_sigma). The implementation now supplies the matrix inverse and the explicit reference-column solve. Reference columns omit the finite covariance errors and use a local diagnostic phase choice; they are not actual pulses or a uniform phase construction.

For later signed increments, Proposition 7.6 uses delta_a_sigma=[H^-1*(Sigma/epsilon)]_sigma/(2*a_sigma), about the fixed positive amplitudes. The algebra permits negative increments and preserves the factor 2 in the cross-covariance identity. It does not require every later increment to satisfy the leading cone. No physical pulse field or smooth edge amplitude has been constructed by this milestone.

## Paper versions

- Original OpenAI PDF (166 pages): [official PDF](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf), SHA256 `0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`. Profile/shear pp.26-27 (4.8)-(4.11); cone pp.30-32 (4.20)-(4.23); uniform edge hypotheses p.33, Theorem 4.6(iii)-(iv); frame p.74 (7.1); actual covariance pp.82-84 (7.24)-(7.30); signed linear lift p.84 (7.31)-(7.34). References are PDF page numbers.
- Supplied Lei--Ren `2609.35406v2.pdf` (245 pages), SHA256 `8396b998dcf737cd6a7c12ff1019c6b2fd308a7e6f0907c0c6647a9b2ec64e8d`: scales pp.11/18 (2.2)-(2.5), (3.1)-(3.3); shears p.21 (3.18); strict and relaxed cones p.22 (3.21)-(3.23).
- Supplied `2609.17642v1.pdf` (31 pages), SHA256 `a38af23a5873e7671ecf0bcd7f3d4c45851c67689a7551ae0e51203ba166fe95`: Ramani Duraiswami's GD1998 comparison paper, not the original OpenAI manuscript. Use its numerical methods and the shear-freedom comparison; do not use its smooth surrogate as proof of the original cone.

## Next executable tasks

- [ ] **F57C-cone1b-main:** apply the same source-correlated reduction to the entire main xi[.02,10] and exit xi[10,11], all Z[-1,1]. Preserve the selected implicit amplitude, incoming histories, exact partial-kernel integration by parts and complete future/pressure energy cancellation. Use real domain subdivision or an analytic positive rewrite where dependency widens the signed margin. Report each failed or inconclusive bound; fresh point checks alone cannot mark whole-domain admission.
- [ ] **F57C-cone1b-regions:** add same-source F/shear adapters for the remaining 30 registry routes. Keep each original radius/amplitude recipe correlated. Distinguish stress-free core/exterior, annular strict margins and edge direction limits. Bind all 32 adjacent and 14 support traces rather than inventing stress/shear from absolute envelopes.
- [ ] **F57C-cone1b-inner:** wherever the real inner construction only meets the relaxed cone, construct the original periodic shear loop and profiles E_N=E_0*exp(A(X,eta,NlogX)/N), U_N=U_0+B(X,eta,NlogX)/N. Recover the exact radial shear chain rule, choose one sufficient finite N uniformly in eta, and perform the independent five-moment repair. Preserve the analytic collar, reserved positive/mean intervals and absolute axis pressure. A changed source family requires new dependent receipts; do not carry old admission to the changed profiles.
- [ ] **F57C-cone1c-pulses:** solve the actual homogeneous pulse equations in the chosen slow boxes with the original two auxiliary support rectangles and rounded phase frequencies. Compute the covariance integrals (7.27), including the cosine factor 1/2 and torus area. Bound finite column errors and determinant, then use the actual columns to prove positive y on the whole active shell.
- [ ] **F57C-cone1c-edges:** construct the original positive flat weight zeta, uniform limiting direction and derivative estimates (4.26)-(4.27); prove smooth zero extension of sqrt(y) and a common uniform phase parameter. T=0 or a finite number of edge samples cannot establish these requirements.
- [ ] **F57C-cone1c-linear:** instantiate the signed linear lift with those fixed actual amplitudes. Preserve both wave families and the full meridional/diagonal tensor contributions through mean correction and averaged quadratic cancellation.
- [ ] **F57C-recursion:** implement actual n=1 and distinct n>=2 recovery/moment repair on the same core interval; cancel or absorb the certified nonflat leading origin E term, then establish finite-order estimates and smooth summation. The remaining global cone, corrected NS, finite-energy and measured scale-recursion goals stay active.

## Validation and scope

18 exact signed source/map/covariance identities; independent interior/exterior/weak-shear/zero/inconclusive cone fixtures; prescribed matrix coefficients (2,1), boundary/singular rejection, a negative signed cross-covariance increment, and coefficient-scaled exponential enclosures. The five real current-source requests all pass the strict two-vector test. Reference coefficient positivity uses the shared frame correlation rather than assuming an interval matrix inverse is sharp.

New gates: original map available; signed current pulse margin adapter available; whole current pulse-end signed two-vector cone certified; original covariance and signed linear algebra available. Global completed cone/lift, converged u/v/w/p, full NS, energy, flat completed remainder and genuine coefficient recursion remain false. Keep the long-term goal active.
