# Original restricted q-flat cone margins and composed local C0 frequency budget

Checked source: [1b301374](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/1b30137409804f786f83893c7e0d6cf0e0ace019). Predecessor: [repair-band functions and local power budgets](CURRENT_NATIVE_REPAIR_BAND_CONE_BUDGET_2026_10_07.md). Producer 0.094s; focused checker 0.171s; 1,148 matching working/index dependency hashes. Read-only mathematical reviewer: **gpt-5.6-luna / max**. No new workers or source ancestor constructors were used.

The original modulation route now has quantitative weighted G1..G4 input margins on every applicable q=0 subset, complementary to the accepted active-loop bounds. The original 24 continuous cells and 35 conditional cutoff branches are retained. Thirteen branches on six source charts receive 52 normalized state-frequency requirements. Eleven other source charts have source-proved empty q-flat subsets. The combined lower bound retains all prior source, repair, active, quiet, correction-band and band-positivity conditions.

This closes **CONE2a/c as a local conditional C0 budget layer**. It is not a selected global N, evaluated original field, actual fixed-point control or terminal functional closure. The correction-band budget remains conditional on the original five controls lying in their accepted common C1 ball. Global cone/exterior admission and genuine coefficient recursion remain open. No percentage of the overall mathematical construction is inferred from this local result.

## Restricted source identity and stability

The original loop cutoff gives q=0 exactly on kappa-2>=eta>0; it does not make every source point a strict-cone point. Write

```text
Delta = a+b^2/a-2
D = p1-(b/a)*p2-a-b^2/a
J = p2+(b/a)*p1
Q = 2D^2-Delta*J^2
G1=a, G2=a*Delta, G3=a*D, G4=a^3*Q.
```

Independent symbolic checks recover these four identities and polynomial gradient coefficient sums (1,6,8,208), of degree at most five. Each q-flat margin uses the original whole-chart caps for **all four a,b,p1,p2 components**, rather than imposing the active-loop special bound |a|,|b|<=3 on original input. A directed upper cover for H=2+max(3,|a|,|b|,|p1|,|p2|), with g=min(G_i lower), gives the sufficient state tolerance min(1,g/(512H^5)). The perturbation path then preserves each weighted margin by at least half. Source caps are majorants, never defining function values.

## Whole bridge margins and the actual first inlet

For bridge_first/second/macro, retain the actual correlated source K and the signed comparison vector (Dbar,Ebar), with Hbar=(Dbar^2+Ebar^2)/Dbar and omega=1-chi:

```text
D = omega*Hbar + e_theta + (Ebar/Dbar)*e_z
J = e_z - (Ebar/Dbar)*e_theta
Eomega(K)<=2*K1*cstar*K^-80
S(K)=40*K^6*Eomega(K)<=80*K1*cstar*K^-74<1 for actual K>=10^6
rho=|e|/(omega*|(Dbar,Ebar)|)<=1/(20K^5)
(kappa-2)*rho^2<=K^10/(400K^10)=1/400.
```

The K cancellation is performed **before** independent ranges. The saved ratio at Kmin is an upper bound because S decreases with K; it is not an error value to multiply by an arbitrary live K. Consequently D/(omega*Hbar)>=1-rho_max and Q/(omega*Hbar)^2>=2(1-rho_max)^2-1/400, where rho_max=1/(20*10^30). These are whole strong-subset reserves, greater than .95 and 1.8; a selected inner-exit collar is not promoted to a whole theorem.

The first original bridge has zero stress at Ra. The actual modification begins later at r_minus=Ra*exp(hb*sc/2), so the native first phase is bounded below by sc/2>0. The accepted original sc and hb<1/2 give

```text
sigma(s)=exp(-1/s^2)/(exp(-1/s^2)+exp(-1/(1-s)^2))
omega=(1-hb)*sigma(s)
log omega >= -4/sc^2-log4 for s>=sc/2.
```

The checker verifies the original cutoff's monotonicity on (0,1). This factored logarithm remains valid through the first chart without materializing the microscopic stress scale, changing the original width, or assigning a positive floor at Ra. Second/macro omega=1-hb>1/2. Hbar>=2 completes absolute D,Q floors.

## O2 axial and O3 signed margins

On the O2_axial subset b^2/2>=eta, use a=2 and the original full signed D/Theta and Q/Theta^2 reserves. The zero-shear edges and midplane are outside this strict subset. To recover the normalized canonical Theta floor, the checker hash-binds and parses the upstream producer's exact construction:

```text
axial_log_D_lower = lower(logRref+1+log(canonical_theta_reserve)+log(full_D_over_theta_lower)).
```

Removing the **same factored coefficient** conservatively recovers this original Theta reserve. This is not an inference from two unrelated lower inequalities. The checker requires the original baseline, signed proof, exact product expression and directed attachment. The read-only review identified this qualification and verified its resolution.

For O3_slope_mu, a=2+2mu*sigma and b=0; restrict to 2mu*sigma>=eta. For O3_power, a=2+2mu and b=0; use the original 2mu excess. The accepted normalized full signed D floor and Q/D^2>=2-direction>0 give G3>=2D_floor and G4>=8D_floor^2*(2-direction). Full energy and absolute pressure remain in direction; V=0 never implies p2=0.

## Empty subsets and continuous coverage

The source-proved empty charts are switch_first/second/power, reshape, inner_reference, axial_restore, restore_buffer, actual_patch, Rh_reference, O2_slope and O2_buffer. Their original kappa upper bounds are <=2, with eta strictly positive. Switch recipes and the decreasing-K small-shear cap are taken from the original exit certificate/ledger. Restoration and patch retain their actual numeric source caps; reference/slope preserve exact zero b and original kappa range; buffer uses original a=2,b=0. These are empty **q-flat strong subsets**, not removed regions. Every original weak/active conditional certificate remains required.

The implementation reuses unchanged signed continuous state-error polynomials for a,b,p1,p2, geometry, cutoff identities and source provenance. Each term has a negative N power. For N>=1 the log-sum cap gives a directed sufficient logN for the tolerance above. All 24 cells/35 original branches remain; the 13 applicable q-flat branches carry four requirements each. The initial collar retains exact zero correction. Incoming five histories and pressure are preserved even in the quiet power interval. The final maximum includes the previous full source/repair/active/quiet/band lower bound.

## Implementation and evidence

Producer/result: `experiments/root_st073/lei_ren_part1_paper_compliant_current_native_q_flat_cone_budget.py/.json`; independent checker/receipt: same stem plus `_check.py/.json`. Gate: `current_original_restricted_q_flat_margins_and_modulation_C0_frequency_budget_bound`. Checked evidence: four exact weighted-G identities, two bridge projections, correlated K power cancellation, flat-cutoff monotonicity, five original axial AST recipe bindings, all seventeen chart cases, 24 positive weighted floors, all four whole source state caps, 24 cells/35 retained branches and 52 negative-power state requirements. Dependency hashes are audited against the Git index before commit. The checker consumes saved native graphs/receipts and does not rebuild expensive original owners or rerun accepted ancestors.

The sufficient logN remains an extremely large formal bound. No representable machine integer, original point-field phase, converged control or practical numerical N has been assigned from it. Keep this distinction explicit in the next work. The former warm Python session 77153 is dead; do not target it. The long-term goal remains active and incomplete.

## Detailed executable continuation

Take one bounded production task, publish the artifact and focused evidence, then check its box with the commit. Preserve the same original source family, analytic pressure datum, actual mu, five histories and full signed sectors. Reuse accepted receipts unless source/dependencies change; avoid repeating ancestor validations. Root owns math/production; the existing GPT-5.6 Luna/max reviewer stays read-only. Prioritize actual velocity functions and recursion over animation.

- [x] **CONE2a-domain:** bind q=0 to kappa-2>=eta and classify all seventeen original source charts; retain eleven empty-case active certificates.
- [x] **CONE2a-inner:** whole strong bridge D,Q reserves with pointwise K correlation; actual sc/2 first-phase logarithmic stress floor; three source-proved empty switches.
- [x] **CONE2a-O2:** original axial signed D/Theta,Q/Theta^2 with source-identity-bound canonical Theta reserve; exclude zero-shear edges/midplane.
- [x] **CONE2a-O3:** original transition/power signed direction, full pressure/energy, quantitative weighted margins on the restricted strong subsets.
- [x] **CONE2c-local-cover:** merge q-flat and active conditional budgets over 24 continuous cells/35 branches, including inherited quiet histories and prior conditional repair-band bounds.
- [ ] **NEXT CONTROL1b-condition-inventory:** enumerate every remaining frequency/positivity/cone-region/seam condition in one artifact; include accepted source, active, q-flat, quiet, band, repair contract and exterior requirements. Separate proven local lower bounds from genuinely missing conditions. Report whether the existing maximum covers each one; do not rename the local maximum a global N.
- [ ] **CONTROL1b-outer-compatibility:** identify exactly which modified joins or outer intervals require new estimates. Bind the original 2Rc..Rb source, pressure/heat data and cone domain. Produce original-source higher/outer bounds where needed; if a condition contains an upper N bound or growing derivative, establish actual compatibility rather than selecting a larger N automatically.
- [ ] **CONTROL1b-frequency-representation:** choose a source-faithful representation for the enormous frequency without decimal/float overflow. Keep exact integer/phase requirements separate from log majorants. Do not replace the actual positive microscopic width or amplitude by enclosure endpoints, enlarge them for convenience, or transfer a frequency from an obsolete modified family.
- [ ] **CONTROL1b-whole-N:** only after the complete condition inventory is covered, prove one finite common N exists and record the chosen representation and every admitted condition. Leave `current_whole_N_selected` false until that evidence exists.
- [ ] **CONTROL1c-source-oracle:** implement exact original parameter, chart coordinate/phase, quotient-source and cumulative-history point services from the accepted native graph. Values must be original functions, not whole-chart caps or source samples. Retain explicit source-family identity, arbitrary Z and directed or certified integral contracts.
- [ ] **CONTROL1c-phase-inverse:** resolve the original phase primitive/inverse at that common N, including seams and q=0 transitions. Bind phase=frac(N*original log-radius displacement) where appropriate and preserve the original tiny geometric offsets. Report unresolved inversions rather than substituting sampled phases.
- [ ] **CONTROL1c-integral-evaluation:** cache by genuine free variables, respecting bound integration variables and endpoint dependence. Evaluate original cumulative functions and incoming memory once when constant across an outer quadrature. Avoid rebuilding expensive ancestor constructors for each point or derivative.
- [ ] **CONTROL1d-control-solve:** instantiate all five controls (a0,a2,e0,e1,e2) with the same original mu, B_mu, targets and Q. Use the original divided D=(J-M)/mu convention; retain joint k-row M+muD. Show finite control residuals before claiming convergence.
- [ ] **CONTROL1d-C1-implicit-jets:** compute the controls' ordinary-Z derivatives from the actual equation, including target/amplitude/parameter product rules. Preserve the accepted common C1 ball and distinguish frozen-phase derivatives from total spatial derivatives.
- [ ] **CONTROL1d-tail-certificate:** certify convergent Picard/fixed-point tails for value and Z rows on Z[-1,1], with the original inverse and contraction bound. Report finite residual and certified tail separately. A finite Picard graph is not a fixed point; mark controls installed only after same-source solved-control evidence.
- [ ] **CONTROL2a-terminal-functions:** apply the exact repair-band endpoint identities to the solved controls and tail. Independently integrate corrected m,h,k,e,p and ordinary-Z rows from incoming Rc functions across the three original bumps; prove functional closure throughout Z[-1,1], retaining rate-zero pressure memory.
- [ ] **CONTROL2b-field-installation:** evaluate the solved compact band E_N,V_N and recover radial velocity from the original divergence-free structure. Recover absolute pressure with unchanged P0; expose original coordinate/physical interfaces. Do not minimize residual by adding a fitted pressure tail.
- [ ] **CONTROL2c-join-functions:** establish all needed value/mixed-derivative joins at r_minus,r_plus,Rc,2Rc and later outer seams with original cutoffs. Terminal function closure must precede resetting defect histories outside support.
- [ ] **HIGH-finite-jets:** derive the finite required total-y/Z derivative rows with correct N powers and source width correlation. Higher y derivatives may grow with N; do not impose artificial decay to reuse the C0 proof. Preserve regularity at q-flat transitions, endpoints and the axis.
- [ ] **OUTER-heat-and-energy:** compose corrected outer annuli, flatten/heat collars and exact heat exterior with the original analytic preheat datum. Bound finite physical energy and radial tail in the intended time/space domain; preserve same-source data at every interface.
- [ ] **OUTER-global-cone:** establish modified full signed cone margins in all required regions and joins; compose a global stress lift only after complete coverage. Local C0 budget booleans do not certify a global admissible tensor.
- [ ] **REC-n1:** implement the genuine first-order source/recovery and its independent five-moment repair on the common inner interval, with correct pressure/exterior compatibility. Produce a nonzero coefficient function and directed finite-order remainder.
- [ ] **REC-nge2:** implement the n-dependent higher-order recovery equations and independent repair for each order; retain exact divergence-free structure. Compare adjacent orders on their common interval rather than rescaling one leading profile.
- [ ] **REC-flat-sum:** choose controlled cutoffs on the streamfunction/vector potential before curl, bound coefficient growth/truncation and prove the flat summed remainder. Keep finite-order evidence separate from all-order summation.
- [ ] **WAVE-pulses:** construct the two original oscillatory pulse families and mean corrections from the admitted stress; measure their averaged quadratic momentum flux and cancellation with signed directional data.
- [ ] **WAVE-corrected-NS:** compose corrected velocity/pressure/forcing, preserve smoothness/flat remainder and quantify the full NS residual after wave correction. Background residual magnitude alone is not the stress-cancellation target.
- [ ] **PHYS-uvw:** expose corrected [u(x,y,z,t),v(x,y,z,t),w(x,y,z,t)] in Cartesian space and similarity coordinates; independently check divergence and NS residual, with explicit units/domain/time scaling.
- [ ] **PHYS-recursion-diagnostics:** measure radial contraction, axial/radial aspect ratio, swirl/vorticity exponents and material winding over approaching-critical times and multiple genuine recursion orders. Label geometry-only scaling separately from coefficient recursion.
