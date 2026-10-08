Current successor: [CURRENT_REFERENCE_MIDPLANE_2026_10_08.md](CURRENT_REFERENCE_MIDPLANE_2026_10_08.md) adds exact-midplane whole-window C0/Z contributions and identifies the retained large axial carrier. The fixed nonzero-Z C1 result below remains valid; whole-Z closure remains open.

# Original reference mixed source and C1 phase integration

Checked source [bd32c71f](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/bd32c71fcf487750c562e6fc666694086d19b668). Predecessor: [CURRENT_REFERENCE_PHASE_AVERAGING_2026_10_08.md](CURRENT_REFERENCE_PHASE_AVERAGING_2026_10_08.md). The existing GPT-5.6 Luna/max worker reviewed the new mixed-chain and weighted-curvature formulas read-only; root implemented, computed and published. No child was spawned. Unrelated work and ancestor producers remain untouched.

## Implemented change

**OriginalReferenceMixedC1** issues source frames carrying actual original C0,y,Z,yZ roots over entire radial cells at fixed nonzero Z. It differentiates the original ordinary-Z templates in native log radius, including every radius power and profile derivative. The original L_Z terms already in the Z rows are retained once. Shared P0,P0_Z,P0_ZZ functions are fixed under y differentiation. Pressure-affine late errors retain explicit Pstar^-1 factors.

A directed mixed-jet algebra applies all product, reciprocal and trigonometric chain rules in the canonical original source-factor atlas. The public primitive method accepts only its own issued frame identity; saved caps, midpoints or replacement frames cannot enter that route.

At fixed original fast phase, the inverse chain has

~~~
D = 1+t^2
psi_i = -T2_i/D
psi_yZ = -T2_yZ/D
           + 2*t*(t_y*T2_Z+t_Z*T2_y)/D^2
           - 2*t*t_psi*T2_y*T2_Z/D^3

T1_total_yZ = T1_yZ - t*T2_yZ/D
           + (t^2-1)*(t_y*T2_Z+t_Z*T2_y)/D^2
           + (1-t^2)*t_psi*T2_y*T2_Z/D^3
~~~

Both first cross terms, fixed-angle mixed derivatives and phase curvature remain. A_yZ and B_yZ use these actual function enclosures and the full original E product rule. These slow derivatives hold phi fixed; they are not the total radial derivative of phi=N*y+phi0.

The signed reference branch bounds source-correlated curvature before taking absolute values:

~~~
h = hinv = 1/sqrt(1+u^2), s=h^2=1-r^2, K_i=u_i*h
n=r+coschi, Dchi=s+2*r*n, t=2*q*n/h
T2_i=K_i*J
J=2*q^2/r^2*(sinchi*(Dchi-3*s)+s*(chi-psi)/r)
t_i=K_i*(2*coschi-r)*t-2*q*K_i*h
t_psi=-2*q*sinchi*Dchi/h^3
~~~

Under actual source conditions q>=1/2 and strict |r|>0, the complete functions 2*t*t_psi*J^2/D^3 and (1-t^2)*t_psi*J^2/D^3 have finite directed enclosures. The proof retains h and n across the turning region rather than multiplying independently enormous direction and inverse covers. These conservative bounds are error envelopes, never field values. The theorem does not cover r=0.

The finite-N nonlinear derivative is recomputed from actual N:

~~~
Q_N=E*A^2*R2(A/N)
Q_N_Z=E_Z*A^2*R2(A/N)+E*A*A_Z*exprel(A/N)
~~~

All five C0/Z nonlinear remainder product rules remain. Original leading reflection supplies a periodic primitive and its Z derivative. C1 integration by parts uses genuine G_Z and G_yZ bounds. Physical endpoints remain; only internal traces of this one continuous reference source telescope. Own rates m/e=1,h/k=3/2,p=0 and radius Jacobian1 are retained. Pressure memory has mass5; separate P0 and unknown incoming histories are not reset.

Direct and averaged contributions are computed independently and the tighter valid bound retained componentwise. N160 and N320 Z averaging is still weaker than direct integration; direct bounds remain active there. N16384 is independently evaluated as a local candidate, not selected as an admissible global frequency.

Four files with stem **lei_ren_part1_paper_compliant_current_original_reference_mixed_C1_integrals** are published under experiments/root_st073: producer, deterministic compressed manifest, checker and receipt.

## Actual bounds

Original Rh_reference[-5,0], Z=37/100,16 radial cells. Rounded outward absolute contribution bounds:

| Row | Direct Z,N16384 | Averaged Z,N16384 | Averaged C0,N16384 |
| --- | ---: | ---: | ---: |
| h |0.000031256|0.000001104|0.000000001411|
| e |0.000044118|0.000001317|0.000000001502|
| p |0.000161953|0.000005077|0.000000003579|

At the same N and partition, Z envelope reductions are approximately28,34 and32 times. These are local source-window integral enclosures, not signed terminal moment errors or Cartesian NS residuals. Whole axial coverage and the full matching route remain open.

At N160/16 cells, direct Z h/e/p bounds remain approximately0.003206,0.004531,0.016632; this conservative averaging bound does not improve them. N320 similarly retains direct bounds. The accepted C0 method is recomputed on the three common grids. Physical bounds agree within a separately declared relative2^-48 comparison budget; formal anchors and export precision can differ. The comparison budget never becomes an installed integral bound. Directed source and integral exports retain the original high precision. Every N-specific nonlinear remainder is reevaluated from the original source rather than rescaled from a saved result.

Producer32.735s, focused checker95.969s. Evidence:11 exact mixed source identities;10 independent original C0-function mixed derivative comparisons;3 r/s/hinv mixed identities;2 original implicit-chain identities;7 nonnegative rational certificates;4 source-geometry/Cauchy certificates; exact Q_Z identity;56 independent defining-integral/chain comparisons;384 finite curvature diagnostics supplementing the algebraic proof; original frame/domain, endpoint, source, fresh-N, precision and zero-rate-memory contracts. **1212 Git-index dependency hashes PASS**. No ancestor producer or unrelated check suite was rerun.

## Ordered executable queue

Read this queue before historical handoffs. Take the first unfinished dependency, implement a bounded source change, compute necessary evidence, and mark DONE with code, scoped receipt and commit. Preserve the persistent full objective. Keep numerical oracle errors separate from formal constants and Picard tails.

- [x] **REFERENCE-ACTUAL-MIXED-yZ-ROOTS-FIXED-Z:** original Z templates differentiated radially, radius powers and pressure factors retained; issued C0,y,Z,yZ frames over whole radial cells.
- [x] **CORRELATED-SIGNED-REFERENCE-PHASE-yZ:** original inverse mixed derivative, both cross terms and correlated curvature theorem; actual A_yZ/B_yZ enclosures.
- [x] **REFERENCE-C1-ENDPOINT-IBP-FIXED-Z:** genuine leading Z/yZ and nonlinear remainder Z derivatives; C0/Z physical endpoints retained; same-source internal telescoping.
- [x] **FRESH-N/DUAL-C0-Z-INTEGRALS:** actual levels(4,160),(16,160),(16,320),(16,16384); direct Z retained where averaging is weaker.
- [ ] **EXACT-MIDPLANE-WHOLE-CELLS-FIRST:** add a separate exact-Z0 reference adapter, preserving the signed service's Z0 rejection. Reuse the original templates and the accepted small-u slow_Z branch; use p2_Z*q/dstar at p2=0 rather than a p2_Z/p2 ratio. Require exact kernel.u.zero or the accepted regular small-r guard (max|r|<=1/4 and positive denominator). If a non-singleton source box meets that guard, use the Fourier series with its tail; otherwise first refine the original pressure parity/late-error enclosure. Never force p2=0 or declare a large-u interval small. Acceptance: whole[-5,0] C0/Z own-rate contributions at exact Z0, genuine nonzero p2_Z retained and no claim of interval-Z coverage.
- [ ] **FUNCTIONAL-Z-SOURCE-ATLAS:** extend exact-Z atlas to rational closed Z intervals with shared original powers and positive L(Z). Preserve source family, selected logC/logP, pressure datum and phase origin. Do not materialize astronomical R or subtract huge common logarithms independently. Acceptance: issued whole-cell C0,y,Z,yZ enclosures in both variables, exact source-factor identities and isolated finite-unit derivative references.
- [ ] **REFERENCE-MIDPLANE/SMALL-u-C1:** preserve p2=0 with nonzero p2_Z at Z=0. Use stable original primitives or defining-integral series with directed tails; no r^-1 at zero and no zeroing true Z derivatives. Bound mixed derivatives before joining the signed branch. Acceptance: an entire Z neighborhood containing zero, branch overlap and exact phase/derivative compatibility.
- [ ] **SIGNED-OVERLAP/WHOLE-Z-REFERENCE:** cover[-1,1] by exact rational intervals. Keep rho,s,hinv positive as formal source functions. Use signed theorem only with strict |r|>0 and q>=1/2. Compare the same actual source values and derivatives at overlaps, not selected caps. Acceptance: complete reference-window C1 envelopes as Z functions, including midplane; point samples are insufficient.
- [ ] **CURVATURE/DRIFT-SHARPENING-IF-NEEDED:** reduce conservative covers analytically or by certified subcells while retaining the complete-product identity. Never differentiate saved bounds. Acceptance: source-linked reduction at common N without degrading direct bounds. Functional coverage takes priority over this optimization.
- [ ] **SIGNED-NONLINEAR-PERIOD-MEAN:** integrate genuine Q_N and quadratic terms with the original phase measure/Jacobian and controlled y/Z drift when signed terminal error needs it. Leading reflection does not zero nonlinear means. Acceptance: signed mean, endpoint, drift and oracle errors separated, actual N remainder and original phase origin.
- [ ] **O2-MIXED-ROOTS/GENERAL-PARAMETER-CHAIN:** build original O2 C0,y,Z,yZ roots from accepted kernels and inherited pressure histories. Include a,b,t0,q derivatives; reference zero identities require proof before reuse. Acceptance: native domain[0,1], correct offsets/Jacobian, issued source frames and independent mixed references.
- [ ] **O2-C1/AXIAL/11-UNIT-BUFFER:** integrate genuine C0/Z densities with inherited histories and all11 buffer units, native offsets and own rates. Acceptance: typed integral outputs with explicit seam/incoming-history contracts rather than zero initial data.
- [ ] **O3/RESTORE/REMAINING-UNITS:** follow the predecessor detailed source queue for every region. Implement unsupported units from original definitions, keep pressure/phase histories, fail closed on missing source data. Acceptance: all17 charts/24 native cells have actual C1 function/integral services at a common N.
- [ ] **ALL-ROUTE-C1/ACTUAL-FIVE-CONTROLS:** attach actual contributions to the accepted centered all-N/unit-C1-ball bridge. Evaluate terminal five-moment maps as Z functions; separate computed uncertainty from theorem constants. Acceptance: common full-route oracle and functional terminal closure; local N16384 is insufficient.
- [ ] **GLOBAL-N/PICARD:** choose N using all-route bounds and contraction conditions, preserving logarithmic/formal frequency when needed. Replay genuine numerical corrections and report oracle error separately from Picard tail. Acceptance: source-backed controls and functional closure at an admissible global N.
- [ ] **BACKGROUND-JOINS/EXACT-HEAT:** use the same accepted pressure and moments across core, annulus, collar and exterior; preserve divergence through streamfunction/vector potential. Verify appropriate smooth traces and finite energy.
- [ ] **ADMISSIBLE-STRESS/FLAT-REMAINDER:** compute the original stress/remainder split, cone margins in every region and actual flat asymptotics. Local low momentum residual does not replace them.
- [ ] **ACTUAL-n-DEPENDENT-RECURSION:** execute n=1 and n>=2 recovery equations, independent moment repair at each order, controlled remainders and smooth summation. Rescaling a fixed field is insufficient.
- [ ] **TWO-PULSE/QUADRATIC-CANCELLATION:** implement both original oscillatory families, mean corrections and actual averaged quadratic stress cancellation.
- [ ] **CORRECTED-CARTESIAN-UVW/PHYSICAL-DIAGNOSTICS:** assemble u(x,y,z,t),v(x,y,z,t),w(x,y,z,t), independently validate corrected forced NS residual and measure contraction, relative elongation, cumulative winding, scales and finite energy.

Only the fixed nonzero-Z reference mixed source/C1 integration milestone is installed. Whole-Z terminal closure, all17/24 oracle, actual five controls, global frequency, recursion, two pulses and corrected Cartesian field remain open.
