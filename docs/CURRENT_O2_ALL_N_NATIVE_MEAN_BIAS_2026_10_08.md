# Actual original O2 all-N native mean-bias coefficients

Checked source [e1e5895b](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/e1e5895bceded1f2a092c44eefda7aea64c87c2f). Predecessor: [CURRENT_O2_TRUE_PHASE_MEANS_FINITE_N_BIAS_2026_10_08.md](CURRENT_O2_TRUE_PHASE_MEANS_FINITE_N_BIAS_2026_10_08.md), source0a710c20/docsb8ff8595. The existing GPT-5.6 Luna/max worker reviews read-only; root implements, computes and accepts. No new child. Full paper-faithful reconstruction goal remains active.

## Concrete advancement and scope

The actual original five nonlinear phase means now have explicit **uniform all integer N>=160 C0/ordinary-Z N^-2 coefficients**, covering original y[0,1], Z[-1,1], all64 source cells and192 regular/positive/negative conditional predicates. This does not rescale the previously computed N160 density or mean outputs: actual N-independent original E,V,A,B source jets are used afresh in native products. The original pressure datum, q dependence, paired flat endpoints and inherited histories are retained.

For N160 the h mean C0 upper bound improves from1.39435588e-4 to1.13868435e-6 (about122 times), and p from8.03224887e-4 to6.32151856e-6 (about127 times). Ordinary-Z mean-component bounds also tighten. These are mathematical enclosure improvements, **not** physical residual reduction, solved terminal identities or achieved scale recursion.

Four files under experiments/root_st073 use stem **lei_ren_part1_paper_compliant_current_original_O2_uniform_mean_bias**: producer, deterministic lossless compressed report, focused checker and receipt. Main API: **OriginalO2UniformMeanBias.integrate / at_candidate**; actual same-owner coefficients and integer N are required. Candidate N is not admission of a common global N.

## Source amplitude, exact pair functions and native coefficients

Original b=t0=0, strict Phi_psi>0 and Phi(pi)=1/2 imply psi/(2*pi) in[0,1/2] on true phi in[0,1/2]. Since original a=4/5+6*sigma/5<=2 and A=a*(phi-psi/(2*pi))/2, **|A|<=a/4<=1/2**. This is a bound for the same actual A function. It is never differentiated as a selected source function.

For x=A/N, z=B/N, the exact reflected pair densities are m=0, h=E*(cosh(x)-1), k=V*E*(cosh(x)-1)+E*z*sinh(x), e=z^2-E^2*(cosh(2*x)-1)/2, p=E^2*(cosh(2*x)-1)/2. Nonzero V and full finite-N quadratic biases remain.

For C=1/2, rho=exp(C/160), rho2=exp(2*C/160), cosh(x)-1<=x^2*exp(|x|)/2, |sinh(x)|<=|x|*exp(|x|), cosh(x)<=exp(|x|). Therefore the following **native positive coefficients** enclose N^2 times pair densities and ordinary Z derivatives:

~~~
H = |E|*Acap^2*rho/2
P = |E|^2*Acap^2*rho2
K = |V|*H+|E|*|B|*Acap*rho
H_Z = |E_Z|*Acap^2*rho/2+|E|*Acap*|A_Z|*rho
K_Z = (|V_Z|*|E|+|V|*|E_Z|)*Acap^2*rho/2
    + |V|*|E|*Acap*|A_Z|*rho
    + (|E_Z|*|B|+|E|*|B_Z|)*Acap*rho+|E|*|B|*|A_Z|*rho
P_Z = 2*rho2*(|E|*|E_Z|*Acap^2+|E|^2*Acap*|A_Z|)
C0: m=0, h in[0,H], k in[-K,K], e in[-P,|B|^2], p in[0,P]
Z: m_Z=0, |h_Z|<=H_Z, |k_Z|<=K_Z,
   |e_Z|<=2*|B|*|B_Z|+P_Z, |p_Z|<=P_Z.
~~~

Two independently valid envelope families use Acap=the original native |A| cover and Acap=1/2, respectively. Their output covers intersect for the same actual function. All products happen on the same native source atlas before conversion, so original positive q factors are retained. The constant source bound is not divided by q. Actual E_Z,V_Z,A_Z,B_Z enter product rules; no C0 enclosure, amplitude envelope or norm is differentiated.

Ordinary Z is exported only after exact positive-unit division by **U(Z)=Lambda0/L(Z)^2**, Lambda0=Pstar^11*Cstar^10. U is Z-only and restored as a formal native scale after y integration. The normalized cap is not a new function of Z to differentiate.

Each branch uses four accepted original true-half-phase cells, exact full-phase weight1/4 each. The branch mean coefficients are then hulled locally, not added over overlaps. Original final-endpoint own-rate masses are applied once over y[0,1]: m/e rate1, h/k rate3/2, p rate0. No extra R Jacobian, angle Jacobian or period-count multiplier.

## Uniform full-window mean-component coefficients

For every integer N>=160, C0 mean contribution is enclosed by column2/N^2; ordinary-Z mean contribution is U(Z)*column3/N^2. Thus increasing N by a factor10 reduces these **envelope widths** by a factor100. This says nothing yet about the unbounded spatial remainder or total terminal error.

| Moment | Native-normalized C0 coefficient | Ordinary-Z coefficient divided by U(Z) |
|---|---:|---:|
| m | [0.0, 0.0] | [0.0, 0.0] |
| h | [0.0, 0.0291503193] | [-3.39080946e+10, 3.39080946e+10] |
| k | [-0.129285774, 0.129285774] | [-1.09636224e+11, 1.09636224e+11] |
| e | [-0.0829754276, 0.456482285] | [-2.31897796e+11, 2.31897796e+11] |
| p | [0.0, 0.161830875] | [-2.05134122e+11, 2.05134122e+11] |

At N160 only, the new all-N mean covers also intersect with accepted direct nonlinear true-phase means:

| Moment | Refined N160 mean C0 contribution | Refined ordinary-Z contribution divided by U(Z) |
|---|---:|---:|
| m | [0.0, 0.0] | [0.0, 0.0] |
| h | [0.0, 1.13868435e-6] | [-1324534.94, 1324534.94] |
| k | [-5.05022554e-6, 5.05022554e-6] | [-4282665.0, 4282665.0] |
| e | [-3.24122764e-6, 1.78313393e-5] | [-9058507.66, 9058507.66] |
| p | [0.0, 6.32151856e-6] | [-8013051.64, 8013051.64] |

Mean m is exactly zero, including ordinary-Z mean derivative; this does not set the actual fast-phase radial integral to zero. Derivative enclosures remain broad. No small physical gradient, terminal norm, global N, total averaging theorem or full NS residual is claimed.

## Focused acceptance

Three exact ordinary-Z product-rule identities.32 compatible finite functions, actual exponential pair and independent Z differentiation at N160,257,1024,10^9;640 raw native N^2 pair C0/Z comparisons, including amplitude boundaries, zero/tiny A and both V signs. These finite diagnostics do not replace the actual large-scale source.

One live same-source coefficient reconstruction revalidates all **768 original phase parts**, **6144 source input rows**, **16896 native positive coefficient rows**, **26250 dual range/phase-weight/branch-union/total checks**,320 independently evaluated positive Duhamel masses. Four candidate frequencies retain native C0/Z units; nine copied-owner, partition and N guards. Exact staged dependency hash audit PASS for **1280 files**. Accepted ancestor suites are reused, not rerun.

Producer 51.984s; focused checker 56.219s. All full-goal/oscillatory flags remain false.

## Ordered tasks and acceptance

- [x] **O2-SOURCE-A-AMPLITUDE/DUAL-MEAN-NORMS:** original half-phase monotonic source theorem |A|<=1/2, native q-preserving and physical amplitude envelopes, actual C0/Z derivatives, tighter mean bounds; source e1e5895b.
- [x] **O2-ALL-N-MEAN-NMINUS2-COEFFICIENTS:** N-independent actual source inputs, fixed uniform exponential floor160, exact source units and five own rates, complete original continuous domain; source e1e5895b. Mean component only.
- [ ] **NEXT O2-COLLECTED-SLOW-Y/YZ-NORMS:** collect actual carrier p2 and correlated q_y before inverse-q division. Regular branch must retain raw original p2/p2_y covers in derivative calculus instead of substituting p2=d_star*u/q; keep separate |u|<=1/4 geometry predicate. Signed branches may use the necessary actual implication q>=(3/16)*d_star/(Lambda0*max|g/Lambda0|), without q floor or declaring the entire enclosing rectangle satisfies the predicate. Retain d_star as its exact formal offset and paired flat-endpoint suppression. Acceptance: finite native normalized actual A_y,B_y,A_yZ,B_yZ and first-order density derivative coefficients on every original branch; same source/ordinary derivatives, meaningful comparison to predecessor, no derivative of covers.
- [ ] **O2-ACTUAL-ZERO-MEAN-FIRST-ORDER-PRIMITIVES:** bind f1=(B, E*A, E*(V*A+B), 2*V*B-E^2*A, E^2*A) to actual same-source functions. Define G(phi)=integral_0^phi f1 ds with G(0)=G(1)=0, and actual slow G_y,G_Z,G_yZ. Full density=f1/N+r_N/N^2; nonlinear r_N is not zero mean. Acceptance: explicit native source coefficients and correct units, both signs/axis, no cap-only generic function substitution, mean bias retained.
- [ ] **O2-SPATIAL-OSCILLATORY-REMAINDER/SHARP-AVERAGING:** literal phase frac(N*y+phi0), unchanged actual phi0_Z=0, own-rate kernel K. Integrate first order by parts: ([K*G]_a^b-integral K*(G_y+lambda*G)dy)/N^2; Z retains G_yZ. Keep every cell/global endpoint, wrap and join, unless cancellation is proved for the same source function. Add the full nonlinear remainder and source/oracle/integration errors. Acceptance: complete original five C0/Z spatial contributions with useful common N dependence, not just phase means.
- [ ] **ACTUAL-INHERITED-HISTORIES/ALL17-24:** source-owned five compatible incoming Z functions, analytic P0 and rate0 pressure memory; original axial-buffer parameters, offsets and derivative joins; continuous actual inverse/primitive/density/integral source coverage. Zero/affine history fixtures and constant-q reference routes do not meet acceptance.
- [ ] **FIVE-TERMINAL-FUNCTIONAL-CONTROLS/GLOBAL-N:** join all routes and all signed endpoint/mean/error terms, solve five whole-Z terminal identities and satisfy original inequalities with one admitted global N. Do not label candidate queries as admission.
- [ ] **JOINS/PREHEAT/HEAT/ENERGY/STRESS/FLAT/n-RECURSION/PULSES/CORRECTED-UVW:** preserve original full scope: same-data background joins/pressure, exact heat and finite physical energy, admissible stress cone and flat remainder, actual n-dependent recovery equations/moment repairs/smooth sum, both pulse families and stress cancellation, Cartesian velocity/forced-NS plus contraction/elongation/winding diagnostics. Actual scale recursion remains open.

Use this successor as the active queue. Mark a scoped task DONE only with accepted evidence and a commit. Full reconstruction remains active; proceed with actual source derivative collection and spatial contribution.
