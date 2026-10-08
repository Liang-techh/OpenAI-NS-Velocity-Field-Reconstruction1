> Successor: [CURRENT_ORIGINAL_FACTORED_TRANSITION_SOURCE_2026_10_08.md](CURRENT_ORIGINAL_FACTORED_TRANSITION_SOURCE_2026_10_08.md), source [d3e2ee6f](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/d3e2ee6f94f2f544fc0e90cb24d27cac909c0364). Original velocity/history/root adapter and local five C0/Z integrals now execute on four bulk/seam queries; full prefix remains open.

# Original microscopic transition coordinates and smooth q source

Accepted source [ca0fb994](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/ca0fb99479d436ef5b91d90efc16e38f028f5036); predecessor [CURRENT_ORIGINAL_POSITIVE_LOG_SIGMA_2026_10_08.md](CURRENT_ORIGINAL_POSITIVE_LOG_SIGMA_2026_10_08.md). Root owns math, code, computation and publication. The existing read-only GPT-5.6 Luna / max reviewed typed endpoint/source/phase contracts and the next original-source adapter. No new child, worker edits, tests or compute.

## Concrete result

The microscopic O3 cutoff now has executable same-source coordinates and original q value / ordinary y2,Z1 derivative enclosures. Five queries cover an open bulk piece, two active seam collars, a complete local active-to-flat cutoff seam and a flat collar. Active pieces retain nonzero q. The exact original s=0 endpoint retains q=sqrt(eta/2); it is not reset to zero.

This advances a source dependency for five-moment integration. It does not yet evaluate the microscopic velocity amplitude, conditioned inverse, density or complete prefix integral. The bulk before xi=1/4 and the bulk-to-seam interval are explicitly not integrated. Previously accepted quiet interval [2^-128,1] and all original incoming/pressure memory remain unchanged.

s=log(R/Rd) is the original radial transition coordinate, not physical Navier-Stokes time. Candidate N1024 and the original source family, mu, eta, dstar and defining field are unchanged.

## Typed original geometry

    D = log(mu)-log(eta)>0; W=D^(-1/2).
    Bulk: xi=s*sqrt(D), s=W*xi, ds/dxi=W.
    Seam: k=1/s^2-D, s=(D+k)^(-1/2).
    ds/dk=-(D+k)^(-3/2)/2.

Every endpoint retains its defining expression and the entire directed original parameter cover. Compact mpf exponent tuples provide scalar covers for source/phase inequalities; they are not selected endpoint/midpoint values. No astronomical Fraction denominator or exp(-D) is allocated.

For k_left>k_right>=0, increasing-s width is collected exactly before endpoint subtraction:

    width = (k_left-k_right) /
      [sqrt(D+k_left)*sqrt(D+k_right)*
       (sqrt(D+k_left)+sqrt(D+k_right))].

The implementation factors D^(-3/2) and retains its positive coefficient. Width stays positive even when directed scalar endpoint covers overlap at working precision. Original Duhamel masses/attenuations use this width once. Ordinary y derivatives below do not receive W or ds/dk again.

The exact original affine O3 radius identity and genuine native source phase at s=0 supply the phase origin. Nonzero formal s and N*s components are retained separately before periodic projection. The phase uses s, never xi/k. Original analytic phase uncertainty and any full-period cover remain.

## Correlated original excess and q

Let F=1/(1-s)^2 and L=F-1/s^2. The original sigma and Delta are unchanged:

    sigma(s)=exp(L)/(1+exp(L)), Delta=2*mu*sigma, rho=Delta/eta.
    Bulk log(rho)=log(2)+D*(1-xi^(-2))+F-log(1+exp(L)).
    Seam log(rho)=log(2)+F-k-log(1+exp(F-D-k)).

Common D terms cancel symbolically before giant-log subtraction. The positive logistic denominator retains its directed nonzero exponential tail. Ordinary y derivatives use L_y=2/(1-s)^3+2/s^3 and L_yy=6/(1-s)^4-6/s^4 with the original logistic chain rule. All Z rows of rho/q are exactly zero because the original transition coordinate and mu/eta are Z independent.

Original q remains lazy:

    q=sigma(1-rho)*sqrt(eta*(2-rho)/(2*a)), a=2+eta*rho.
    rho>=1 everywhere: q and every slow jet are exactly zero.
    rho<1: eta<2eta-Delta<=2eta, retained before the square root.

Across the queried seam, rho<2 throughout. The root is real and positive even on the flat part, where the original smooth cutoff and all derivatives give zero. Thus the same formula encloses original q and ordinary y2,Z1 over the whole seam. This is not a newly fitted q or artificial branch selection. Source cells unable to prove a positive body still request subdivision.

## Executed source windows

| Query | Typed endpoints in increasing s order | Result |
|---|---|---|
| Open bulk | xi: 1/4 -> 3/4 | active, nonzero q |
| Near seam | k: 4 -> 2 | active, nonzero q |
| Active seam collar | k: 2 -> 7/4 | active, nonzero q |
| Cutoff seam | k: 7/4 -> 5/3 | smooth active-to-flat q and y2,Z1 |
| Flat collar | k: 5/3 -> 0 | exact q / slow-jet zero |
| Exact original endpoint | s=0 | q=sqrt(eta/2), slow jets zero |

These are source queries, not a full prefix partition or five defect integrals. The seam width is of order D^(-3/2), much smaller than W; neither width measures vortex-core contraction or true temporal coefficient recursion.

## Evidence

Five files use stem lei_ren_part1_paper_compliant_current_transition_normalized_q_source under experiments/root_st073/:

- .py: original typed geometry, collected Delta/eta, original q/ordinary derivatives and genuine affine radius/phase.
- .json: original parameter binding, source branches, retained limits and false full-construction gates.
- _evidence.json.gz: lossless geometry/q/phase/mass evidence; 88,811 compressed bytes.
- _check.py / _check.json: 183 independent direct physical-y q/excess comparisons, 15 independent exact widths, original parameter replay, 30 saved q rows, 25 true mass/decay rows, actual N*s and nonzero offsets.

Focused check PASS; 1266 Git-index hashes match. Producer 58.641s; checker 11.766s. Independent fixtures do not choose native field values. Accepted ancestor numerical integrations were not rerun.

## Executable next tasks

- [x] **TYPED-ORIGINAL-MICROSCOPIC-COORDINATES:** exact xi/k source expressions, positive collected width, true Jacobian once, original N*s phase and whole parameter cover.
- [x] **CORRELATED-ORIGINAL-DELTA/ETA-Q-SOURCE:** original D cancellation, conditional gamma, nonzero active bulk/seam q, ordinary y2,Z1, smooth original seam, lazy flat branch and exact original endpoint.
- [x] **SAFE-ORIGINAL-TRANSITION-KERNELS (next) (queried microscopic source component; see successor):** outer_buffer.transition_kernels calls outer_initial.stable_sigma, directly evaluating exp(-1/s^2). At this native scale the exponent is unallocatable. Add a source-owned evaluator/AST callback with unchanged logistic and integral definitions, directed positive-log tails and exact endpoints. Collect each positive original cell width before subtraction. Replace exp(b)-exp(a) and exp(-a)-exp(-b) with their exact expm1 identities. Retain J, theta/energy/pressure kernels, original mu and all accumulated histories. Do not monkeypatch accepted provider methods or rerun prior integrations.
- [x] **ORIGINAL-MICROSCOPIC-EXPONENTIAL-CORRECTIONS (queried microscopic source component; see successor):** retain the formal nonzero corrections in exp(-s/2-mu*J), exp(-s), exp(-3s/2). Use exact expm1/one-minus-exponential integral identities; preserve source expressions before final directed enclosure. A rounded scalar one cannot replace the original profile factor or erase its history increments.
- [x] **EXACT-ORIGINAL-O3-AMPLITUDE-ROOTS (queried microscopic source component; see successor):** retain log Utheta/Pstar=log u_d-s/2-mu*J(s), Z shape and ordinary source rows as one frozen formal unit. Preserve original u1, signed axial/p2 roots, history signs, Pstar basis and separate P0. Carry E/E_Z/p2 in one source packet/basis/ledger before inverse evaluation. No amplitude cap selection, fitted field or axial sector reset.
- [x] **ORIGINAL-SOURCE-ADAPTER-BINDING (queried microscopic source component; see successor):** compile a copied original pre_pulse_mixed_C4.slope_mu with only safe kernel/sigma/arithmetic callbacks. Assert original source assignments and history formulas remain unchanged. Do not substitute SharedOuterBuffer.slope_mu for the defining native field. Retain original J(1)=1/2 symmetry if the adapter reaches that endpoint.
- [ ] **ORIGINAL-OPEN-PREFIX-AND-BULK-TO-SEAM:** cover s=0/open prefix and the gap between xi=3/4 and k=4 with correlated source bounds. Keep the exact active endpoint and all cutoff seams. The five current queries do not cover this gap or establish integral closure.
- [ ] **ORIGINAL-PREFIX-C1-DENSITY/INTEGRAL:** connect safe source amplitude, refined E/p2, original q jets, actual phase, conditioned inverse and all five signed density C0/Z functions. Integrate with true W / rationalized seam width once and retained source error. Refused cells remain explicit; original incoming cannot be reset. Compose accepted quiet operator after actual prefix input exists.
- [ ] **ORIGINAL-AXIAL-b/FULL-BUFFER:** retain joint b/amplitude/radius signs across axial crossings and source-bound period means/slow variation/actual boundary phase over enormous buffer widths. Accepted one-period buffer data remain a local component.
- [ ] **CORRELATED-RC-TARGETS/FIVE-FUNCTION-CONTROLS:** complete the original route and compatible analytic pressure/heat target functions with units/error/P0; solve independent moment controls and continuous-Z terminal identities.
- [ ] **WHOLE-AXIAL-HIGH-JETS/GLOBAL-N:** cover required Z including zero and signed edges, genuine mixed derivatives and a single compatible finite N. Source-only q_Z identities do not satisfy this gate.
- [ ] **MATCHED-HEAT/STRESS/FLAT-REMAINDER/TRUE-RECURSION/CORRECTED-UVW:** five moments/pressure/joins/exact heat first, admissible cone/flat remainder, original n-dependent recovery/independent repairs/solenoidal summation, both oscillatory cancellation families, actual Cartesian field and independent residual/dynamics.

All full-construction gates remain false. DONE applies to the checked source component, not complete transition, terminal closure or scale recursion.
