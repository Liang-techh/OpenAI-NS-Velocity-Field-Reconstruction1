> Successor: [CURRENT_COMPLETE_ORIGINAL_AXIAL_TRANSPORT_2026_10_08.md](CURRENT_COMPLETE_ORIGINAL_AXIAL_TRANSPORT_2026_10_08.md), checked source [4033e3dc](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/4033e3dc6d92bf7bb17ed5c12b86d8cdf8e38bae). Endpoint history/stress/phase/inverse, signed five C0/Z integrals and the complete axial operator with genuine slope inlet are now executed on both tiles. Buffer/new Rc update and functional closure remain open.

# Whole original axial signed b/q source

Checked source [2cb96973](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/2cb96973df6e02f63b0751e2572129a0f2656a51); predecessor [CURRENT_COMPLETE_ORIGINAL_TRANSITION_2026_10_08.md](CURRENT_COMPLETE_ORIGINAL_TRANSITION_2026_10_08.md). The root owns mathematical decisions, code, computation and publication. The existing read-only GPT-5.6 Luna / max worker reviewed the formulas and source bindings; no new child was spawned.

## Implemented result

The original O2 axial selector p in [0,1] now has signed b, correlated Delta/eta and q ordinary y2/Z1 source covers over both microscopic endpoint layers and a proof of strict flatness over the entire intervening interval. This is executed on the actual Z[.36,.38] and [-.38,-.36] tiles, using the accepted original eta, Pstar and Md. Each of four endpoint/tile combinations has seven typed source intervals: closed endpoint, bulk-to-seam, outer seam, active seam, approach to cutoff, smooth cutoff seam and flat connector. All 28 queries are enclosed; 168 ordinary q rows are saved and replayed.

Both exact axial endpoints have b=0, Delta=0 and nonzero q=sqrt(eta/2). The middle has q and all required jets exactly zero. The source retains b's original sign -sign(Z), including ordinary Z derivatives of q; an O3-specific q_Z=0 assumption is not used here. The original a is exactly 2, so the q body has denominator 4, not 2*(2+Delta).

The source partition covers the complete original axial selector on these strict-sign tiles. It is a b/q and geometry interface, not a completed axial velocity/history/stress/density integral. Axis Z=0, full Z, five functional targets/controls, global frequency, exact heat/stress, temporal coefficient recursion and corrected NS gates remain open. The already accepted complete O3 transition remains available.

## Original defining formulas and typed coordinates

With M=Md, y=exp(M*p), the original inlet amplitude and turnoff give

    Utheta/Pstar = exp(-1/5)/(1+Z^2) * exp(-(y-1)/2),
    Uz = 4Z*sigma(1-p), a=2,
    b = -8*exp(1/5)*Z*(1+Z^2)/(M*Pstar)
        * sigma_prime(p)*exp((y-1)/2-M*p),
    rho = Delta/eta = b^2/(2eta).

Here sigma_prime is the derivative in p. The source jets are ordinary derivatives in y=log(R/Rref), after D_y=(M*y)^(-1)D_p exactly once. Utheta and Uz use their original physical/common-Pstar units. Pstar^-1 is the -2 power of the same sqrt(Pstar) arithmetic unit. The original source assignments, J(1)=1/2, parameter definitions and reflected sigma-prime identity are bound.

Set H=-log(eta)-2log(Pstar), A_left=0, A_right=exp(M)-1-2M and Q_e=H+3log(H)+A_e. Let r=p at the left or r=1-p at the right; use xi=r*sqrt(Q_e/2) in the closed bulk and k=2/r^2-Q_e on the seam. The source retains the exact H-Q_e cancellation before directed evaluation. The original log-rho identity is

    log(rho) = -k + 3log((Q_e+k)/H) + 2/(1-r)^2
      + log(16*exp(2/5)*Z^2*(1+Z^2)^2/M^2)
      + (y-y_endpoint)-2M*(p-p_endpoint)
      + 2log(1+(r/(1-r))^3)
      - 4log(1+exp((1-r)^(-2)-r^(-2))).

The original nonzero y increment is kept as exp(x)-1=x*integral_0^1 exp(t*x)dt. Closed-endpoint phase-sigma derivative bounds avoid dividing an interval containing zero; original b products give all normalized-excess jets there. On the seams, collecting rho before square-rooting gives the same signed b=-sign(Z)*sqrt(2eta*rho). Original ordinary log-rho derivatives and their Z factor supply the remaining genuine jets.

The common flat-middle distance is rstar=sqrt(2/(H-log(H)-200)). The endpoint source k_star=-4log(H)-A_e-200 makes Q_e+k_star=H-log(H)-200 exactly on both ends. For rstar <= min(p,1-p) <= 1/2,

    sigma_prime(r) >= 8*exp(1-1/rstar^2),
    log(rho) >= log(H)+202-2M
        + log(2048*exp(2/5)*Z^2*(1+Z^2)^2/M^2) > 0.

This proves the entire middle is strictly flat with the original cutoff. It does not reset incoming moment or pressure memory. Endpoint pieces and middle telescope exactly as rstar+(1-2rstar)+rstar=1; physical y endpoints telescope as exp(M)-1. The right endpoint pieces are reversed before increasing physical-radius transport.

Each recorded physical width is a conservative directed enclosure of the original positive width, using the true phase width and M*y Jacobian once. It is not an exact scalar measure or selected parameter value. Downstream transport must consume the bound together with the original defining endpoint expressions, keeping parameter dependency and the full source phase. The k width is rationalized before subtraction of close microscopic endpoints.

## Evidence and limitations

Six files under experiments/root_st073 use stem lei_ren_part1_paper_compliant_current_axial_endpoint_q_source: .py/.json, two signed lossless .json.gz archives and _check.py/_check.json. Archives total 449,386 compressed bytes and retain the accepted native source/coordinate provenance, original signed b/rho/q jets and full source partition.

Focused checks PASS: 2240 independent original scalar b/rho/q derivative and middle comparisons, 80 phase/physical width comparisons, and 8 finite original-source crossing queries. Scalar b is differentiated independently in ordinary y; q uses analytic scalar derivatives with a proved subprecision cutoff guard. Those finite fixture parameters are not substituted into the actual original source. All 28 actual endpoint queries replay exactly in serialized form, with genuine b signs, active endpoints, smooth seam, flat connector and strict middle checked. 1284 source hashes match the Git index. Producer 2.782s; checker 11.671s. No accepted ancestor producer/check or full legacy suite was rerun.

## Executable next tasks

Only this leading handoff is active. Mark DONE only after committed implementation/evidence, required focused checks and explicit scope limits. DONE, accepted and merged remain distinct.

- [x] **WHOLE-AXIAL-SIGNED-b/q:** original signed b, a=2, correlated Delta/eta, genuine ordinary y2/Z1 q, both nonzero endpoint layers, smooth seams, flat connectors and entire flat middle on two strict-sign tiles.
- [x] **ORIGINAL-AXIAL-TYPED-GEOMETRY:** parameter-derived endpoints, positive rationalized phase widths, retained nonzero y increments, physical Jacobian once and complete reflected source partition.
- [x] **AXIAL-ENDPOINT-AMPLITUDE/HISTORY (next):** bind the original five pre-background histories and separate P0 at each axial endpoint. Left inlet E=exp(-1/5)/(1+Z^2); right inlet E=Pstar^(-1/2)*exp(29/5)/(1+Z^2), before the 11-unit buffer. Recover right endpoint histories from the saved genuine Rd pre-background using exact backward 11-unit zero-Uz ODE transport; these raw histories are separate from pulse corrections. Preserve each actual nonzero microscopic expm1 correction. Bound original turnoff weighted integrals with true widths and original sigma; do not divide through zero or replace kernel integrals by selected caps.
- [x] **ORIGINAL-AXIAL-SIGNED-STRESS-ROOTS:** reuse unchanged original physical_mixed and raw_pre_stress_rows programs on these endpoint sources. Preserve linear/cross/quadratic terms, all history and energy/pressure terms, full P0+p once, and original R=Rref*exp(y). Collect source radius/amplitude correlations before evaluating p1/p2 and preserve their actual signs. Keep a=2, t0=-b/2 and same eta/dstar; an O3 b=0 shortcut is invalid here.
- [x] **ACTUAL-ENDPOINT-PHASE/INVERSE:** bind actual N*log(R/r_minus), with N1024 from the accepted route, the genuine common phase origin and nonzero endpoint y increments. Integrate every actual phase box. Reuse signed-u conditional inverse support if needed; evaluate nonlinear kernels before hulling overlapping branches. Avoid midpoint phase selection or treating a source q cover as an inverse solution.
- [x] **FIVE-AXIAL-ENDPOINT-INTEGRALS:** evaluate all five original signed C0/Z density functions and actual positive-width Duhamel contributions on every nonflat endpoint piece. Keep Pstar conversion and physical Jacobian once. The exact flat connector/middle supply zero own increments while retaining incoming attenuation and pressure memory. Compose the full original axial operator in increasing physical y, reversing the right-endpoint order; bind genuine slope-exit correction. Do not treat a width enclosure as an exact scalar or reset signs.
- [ ] **FULL-BUFFER-PERIOD-TRANSPORT:** derive original period means, slow variation and actual phase-boundary errors for the 11-unit buffer with the accepted period data. Keep leftovers and source correlations; integrate signed densities as functions rather than enumerate every period or choose a primitive cap as a value. Propagate the completed axial correction to Rd and reuse the accepted complete transition to Rc.
- [ ] **CORRELATED-RC-TARGETS/FIVE-CONTROLS:** recover all five continuous-Z terminal defects against analytic pressure/heat targets, with original units, full accumulated errors, retained incoming memory and separate P0. Solve independent original controls, recording conditioning and effects on high derivatives/stress; isolated sampled defects are not five functional identities.
- [ ] **AXIS/WHOLE-Z/HIGH-JETS/GLOBAL-N:** derive a separate smooth source treatment near Z=0, where log|Z| is invalid and the middle is not uniformly flat. Cover required edge domains, mixed derivative orders and one compatible finite N over the full route. Two strict-sign tiles do not certify whole-Z closure.
- [ ] **MATCHING/EXACT-HEAT/FINITE-ENERGY:** enforce original five-moment and pressure compatibility, high-order joins, axis regularity and actual finite-energy radial tail before exact heat exterior admission.
- [ ] **ADMISSIBLE-STRESS/FLAT-REMAINDER:** recover -div(T_B)+E_B, all-chart cone margins and separate high-order flat-decay estimates, preserving a nonzero background.
- [ ] **TRUE-n-DEPENDENT-RECURSION:** implement the original n=1/n>=2 recovery equations, common core interval, independent moment repair per order and controlled streamfunction/vector-potential truncation before curl. These ordinary endpoint source jets are not temporal recursion.
- [ ] **TWO-OSCILLATORY-FAMILIES/CORRECTED-UVW:** implement both pulse families and mean corrections, quantify averaged quadratic stress cancellation, then expose corrected Cartesian u/v/w and independent NS/divergence/dynamics diagnostics under fixed or restricted forcing.

Whole axial b/q is now available on the two actual tiles. Actual endpoint amplitude/history/stress/density transport is the next implementation, followed by buffer transport and functional terminal matching. Full-construction gates remain false.

## Next adapter: exact background recovery formulas

For the original raw pre-background at Y=exp(M), write E_Y=Pstar^(-1/2)*exp(29/5)/(1+Z^2). The saved original Rd histories follow 11 zero-Uz units. Their exact inverse transport is

    mY = exp(11)*mRd,
    hY = exp(33/2)*hRd - E_Y*(exp(11)-1),
    kY = exp(33/2)*kRd,
    eY = exp(11)*eRd + 11*E_Y^2/2,
    pY = pRd - E_Y^2*(1-exp(-11))/2.

These restore the original pre-background, not its pulse correction. For backward microscopic distance d=Y-y, recover

    m(y) = exp(d)*(mY-JV),
    h(y) = exp(3d/2)*(hY-E_Y*(1-exp(-d))),
    k(y) = exp(3d/2)*(kY-JUV),
    e(y) = exp(d)*(eY-JV2/Pstar^2+E_Y^2*d/2),
    p(y) = pY-E_Y^2*expm1(d)/2.

Here JV=integral_y^Y exp(-(Y-s))*V(s)ds, JUV=integral_y^Y exp(-3(Y-s)/2)*(Utheta(s)/Pstar)*V(s)ds, and JV2=integral_y^Y exp(-(Y-s))*V(s)^2 ds, with V(s)=4Z*sigma(1-log(s)/M). The Pstar^-2 factor belongs only to JV2. Bound these original integrals over the true positive microscopic width, retain their signs and nonzero upper tails, and differentiate their defining functions in Z for inherited source jets. Add the same analytic P0 to absolute pressure once. Apply unchanged physical_mixed to generate higher ordinary rows from the original five ODEs. The existing read-only worker checked these inverse identities against the original source; they still require root implementation and focused checks before admission.
