Current successor: [CURRENT_REFERENCE_NEAR_MIDPLANE_2026_10_08.md](CURRENT_REFERENCE_NEAR_MIDPLANE_2026_10_08.md), checked source 28498a53. The three near-midplane tasks below are now DONE on signed zeta windows through[-1e-6,1e-6]. Read the successor for pressure/L interval-correlation qualifications and the next mixed/overlap dependencies. This exact-midplane evidence remains valid.

# Original exact-midplane reference whole-window contributions

Checked source [254c9d2c](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/254c9d2cd2d44feabcbd953236b99c7542caf735). Predecessor: [CURRENT_REFERENCE_MIXED_C1_2026_10_08.md](CURRENT_REFERENCE_MIXED_C1_2026_10_08.md). The existing GPT-5.6 Luna/max worker reviewed source-wide pressure parity and the near-midplane carrier read-only. Root implemented, computed and published; no child was spawned. The persistent full reconstruction objective remains active.

## Implemented change

**OriginalReferenceMidplaneWholeCells** adds a separate exact-Z0 original reference service. The accepted strict signed service still rejects Z0. Two reversible AST projections reuse the original whole-cell programs: remove that signed-only constructor guard inside the separate adapter, and use the actual full-pressure odd-derivative identity to omit only the order1 pressure late-error budget at exact zero. Reversing these changes reproduces the original methods exactly; all inertial formulas, source powers, coefficients, nonlinear terms, phase origins, rates and integration assembly remain.

The original fourteen pressure densities are even in Z and have first Z derivative exactly zero at zero. Their stage origins are independent of Z. The accepted normalized-pressure receipt separately certifies the finite axial end and all late-stage first-derivative budget zero at the midplane. This is a property of the complete original pressure function, not just the baseline approximation -alpha/(1+Z^2)^2. The even order0 and order2 pressure errors remain.

With this original identity applied before term errors, the full-cell source has

~~~
p2(y,0)=0, u(y,0)=0
p2_Z(y,0) != 0
E_Z(y,0)=0
V(y,0)=0, V_Z(y,0)=4/Pstar
L(0)=1
~~~

The accepted regular primitive branch executes its exact-midplane formulas, with psi held fixed for the partials and true original phi held fixed for the inverse derivative:

~~~
T1 = 2*q*sin(psi)
T1_Z = q*u_Z*sin(2*psi)
T2_Z = 4*q^2*u_Z*(sin(psi)+sin(3*psi)/3)
psi_Z = -T2_Z/(1+t^2), t=2*q*cos(psi)
A_Z = -a*psi_Z/(4*pi)
B_Z = -a/(4*pi)*(E_Z*T1 + E*(T1_Z+t*psi_Z))
~~~

The genuine u_Z=p2_Z*q/dstar source remains formal and nonzero. Neither an r^-1 formula nor p2_Z/p2 is evaluated at zero. Period/halfperiod symmetry traces are exact zeros where the original formula requires them; the full derivative function is not flattened.

The unchanged five-density coefficient program computes actual N-dependent exprel(A/N), exp(A/N), both coefficient orders and all C0/Z product rules. Whole radial cells cover[-5,0] and the actual phase unions. All own-rate masses and the radius Jacobian remain; pressure zero-rate memory totals5. Unknown incoming histories and the separate P0 datum are not reset. High-precision records are exported inside the working precision.

Four files with stem **lei_ren_part1_paper_compliant_current_original_reference_midplane_integrals** are published under experiments/root_st073: producer, deterministic compressed manifest, focused checker and receipt.

## Actual result and revised next dependency

Four genuine source levels were computed: (4 cells,N160), (16,N160), (16,N320), (16,N16384). At N16384/16 cells, rounded outward absolute source-window contribution bounds are:

| Row | C0 bound | Z bound divided by Lambda |
| --- | ---: | ---: |
| m |0.000006184|5.648|
| h |0.000010024|1.600|
| k |0.000004074|3.980|
| e |0.000013540|2.067|
| p |0.000051562|5.962|

**Lambda=Pstar^11*Cstar^10*L^-2**, with L(0)=1 exactly. The last column must be multiplied by this actual original positive source factor. It is not a small unscaled Z error. No astronomical exponential is materialized and no finite cap replaces the source factor.

The new result fills exact-midplane C0/Z integration coverage and changes the next action: a nonzero-Z integral bound cannot establish uniform C1 control across the axis. Genuine p2_Z is strictly negative on every computed midplane radial cell and retains the large carrier. Near-midplane coverage needs a factored coordinate that collects this carrier with Z before conversion; uniform ordinary-Z samples cannot replace the missing functional source service. Global N admission must include this derivative factor. N16384 remains only a local candidate.

These are direct interval contribution bounds, not signed terminal moment errors, a selected global frequency or Cartesian NS residuals. Exact Z0 does not establish any open Z neighborhood. The earlier fixed Z=37/100 C1 averaging result remains separately valid.

Producer24.734s, focused checker46.141s. Evidence:14 full-source pressure density parity identities and fixed stage origins;12 exact p2/even-pressure-sensitivity identities;2 reversible source-method projections;10 independent derivatives of original C0 functions;28 independent defining-integral/primitive comparisons with both derivative signs;20 live original point coefficient containment comparisons;4 domain/guard rejections;20 own-rate mass comparisons and pressure mass5; fresh N-dependent coefficients and retained nonzero formal Z carrier. **1204 Git-index dependency hashes PASS**. No ancestor producer or unrelated test suite was rebuilt.

## Ordered executable queue

Read this before historical handoffs. Implement the first unfinished dependency, retain genuine source functions/errors, and mark DONE with code, scoped receipt and commit. Keep the full original goal and the predecessor detailed queue.

- [x] **EXACT-MIDPLANE-WHOLE-REFERENCE:** original C0/Z roots, actual inverse/primitive derivatives and full[-5,0] own-rate contributions at four N/partition levels.
- [x] **FULL-PRESSURE-ODD-JET/EVEN-ERROR-RETENTION:** all14 original densities and fixed origins; exact P0_Z zero including remainder; order0/order2 errors and nonzero p2_Z retained.
- [x] **ORIGINAL-SOURCE-PROJECTION/LIVE-POINT-CONTAINMENT:** original mathematics unchanged outside the two scoped guards; existing true midplane point coefficients lie in the whole-cell enclosures.
- [x] **FACTORED-NEAR-MIDPLANE-SOURCE:** introduce zeta=(Pstar^11*Cstar^10)*Z and a typed finite zeta cell. Keep the original source family, Cstar/Pstar/delta, physical-Z map and L=1-delta*Z^2. Collect the carrier with Z before any interval conversion; do not turn tiny Z into0 and then multiply it by an astronomical factor. Acceptance: actual original E/V/p1/p2 C0/Z functions on a whole zeta interval crossing0, with source factors and a positive L budget.
- [x] **ODD-p2/PRESSURE-REMAINDER-CARRIERS:** factor the full p2 expression only after binding the real pressure parity. Its baseline pressure has analytic alpha jets, while the actual even remainder remains. From original beta2/flatten first-derivative factors, retain the bound |R_P_Z_late|<=10*|Z|*exp(3/5)/(Pstar*(1+Z^2)^3); retain the finite alpha coefficient enclosure with its own |Z| factor too. Acceptance: true p2=Z*even carrier and P0_Z=Z*even carrier, all ordinary Z rows evaluated from original templates, no saved envelope differentiated as a field.
- [x] **REGULAR-SMALL-u-WHOLE-NEIGHBORHOOD:** use actual u=p2*q/dstar after source-factor collection and prove |u|<=1/4 over an entire finite zeta cell. Use the accepted Fourier primitives, inverse denominator and directed tails through zero; retain p2_Z*q/dstar with its large original factor. If the guard fails, refine the carrier interval or zeta domain. Acceptance: whole-neighborhood C0/Z primitives and inverse enclosures; exact Z0 alone is insufficient.
- [ ] **MIDPLANE-yZ/C1-PHASE-AVERAGING:** derive actual radial/mixed source rows, retaining radius powers and pressure parity. At exact Z0, u_y=0 but u_yZ is nonzero, so differentiate the original regular formulas with those identities before forming complete products. Use actual G_Z/G_yZ, finite-N Q_Z and all density product rules; retain endpoints and pressure memory. Acceptance: source-backed C0/Z averaged contributions with formal Lambda, not relabelled direct records or an r^-1 bound at zero.
- [ ] **SMALL/SIGNED-OVERLAP/WHOLE-Z-REFERENCE:** connect the regular neighborhood to source-signed intervals and cover[-1,1] as functions. Existing strict signed inequalities require a nonzero sign interval; compare the same actual values/derivatives in overlap and preserve positive rho/s/hinv. Acceptance: complete reference-window C1 oracle with functional Z coverage, including the central carrier.
- [ ] **SYMBOLIC-FREQUENCY/ORACLE-ERROR-CONTRACT:** retain the source Lambda and other large factors in all-N coefficient and integral bounds. Derive admissible frequency conditions from genuine whole-route C1 constants; bounded integers160/320/16384 do not certify global N. Keep phase-independent envelopes and actual endpoint provenance, numerical error and Picard tail separate. Acceptance: an executable typed frequency/integral contract compatible with the accepted centered bridge, without materializing an impossible native exponential.
- [ ] **SIGNED-NONLINEAR-MEANS/DRIFT:** integrate genuine nonlinear means and controlled y/Z drift where signed terminal closure requires them. Leading reflection does not zero quadratic/exponential means. Keep original phase measure, Jacobian and all errors.
- [ ] **O2/AXIAL/11-UNIT-BUFFER:** build actual C0,y,Z,yZ source rows with inherited pressure/moment histories, nonconstant parameters and correct offsets. Integrate true densities with explicit seam/initial-data contracts.
- [ ] **O3/RESTORE/ALL17-24-SOURCE-ORACLE:** follow the predecessor detailed unit queue, implement remaining definitions and histories, and produce complete same-N C1 function/integral outputs. Fail closed on unsupported units.
- [ ] **ACTUAL-FIVE-CONTROLS/TERMINAL-CLOSURE/PICARD:** attach the full actual oracle to the accepted centered all-N/unit-C1-ball bridge, evaluate five terminal identities as Z functions, select admissible N and replay corrections with separate oracle/Picard errors.
- [ ] **JOINS/EXACT-HEAT/STRESS/FLAT:** assemble consistent global background joins, exact heat exterior, finite energy, admissible stress cone and genuine flat remainder.
- [ ] **ACTUAL-n-RECURSION/TWO-PULSES/CORRECTED-UVW:** execute separate n=1/n>=2 recovery equations and moment repairs, controlled summation, both oscillatory families and quadratic stress cancellation; recover the corrected Cartesian field and physical contraction/elongation/winding diagnostics.

Whole-Z terminal closure, all17/24 oracle, actual five controls, global N, true recursion, pulses and corrected Cartesian NS remain open.
