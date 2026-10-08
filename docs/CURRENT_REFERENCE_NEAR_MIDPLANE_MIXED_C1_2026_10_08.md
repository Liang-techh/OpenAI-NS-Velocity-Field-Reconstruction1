# Original near-midplane mixed source and C1 phase integration

Checked source [989880da](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/989880da5d32dce3f5141820b9a7a1643b6e0c88). Predecessor: [CURRENT_REFERENCE_NEAR_MIDPLANE_2026_10_08.md](CURRENT_REFERENCE_NEAR_MIDPLANE_2026_10_08.md). The existing GPT-5.6 Luna/max worker reviewed regular Fourier tails and the implicit mixed equations read-only. Root implemented, computed and published; no child was spawned. The persistent full reconstruction objective remains active.

## Installed original source and mixed primitive service

**OriginalReferenceNearMidplaneMixed** supplies actual original C0,y,Z,yZ source rows over the whole Rh_reference radial window y[-5,0] and signed zeta[-1e-6,1e-6], with physical Z=zeta/Lambda0 and Lambda0=Pstar^11*Cstar^10. It extends the accepted original whole-neighborhood source; the pressure datum, finite alpha enclosure, even remainder jets, odd first-jet carrier and positive L/Q hulls remain. Separate jet/hull enlargements retain original-function provenance but do not preserve joint arithmetic correlation.

For every original term R^r times a finite coefficient, the radial derivative is r times that coefficient plus the original(f,H,D,P) rates(1/10,1/10,1/5,1/5). This is applied to both the original C0 and original ordinary-Z rows. The compiler preserves original Z powers and one Pstar^-1 pressure-remainder factor before amplification. It does not differentiate an envelope, apply L_Z twice, or interpret Z derivatives as zeta derivatives.

This service is specifically the reference interval: a=4/5,b=t0=0 are constant. The actual original q definition has Delta=a-2<0 and exact cutoff sigma=1; its eta and other datum values are fixed original-family parameters. Thus q_y,q_Z,q_yZ and nu_y,nu_Z,nu_yZ are zero, where nu=1+2q^2. The checker binds this to the actual q implementation. O2 has varying a/q and needs the additional nu derivative terms; this implementation does not claim O2 coverage.

The original regular Fourier branch now has smooth mixed jets through r=0:

~~~
h=1/sqrt(1+u^2), r=u*h, u=p2*q/dstar
r_yZ=h^3*u_yZ-3*u*h^5*u_y*u_Z
h_yZ=-u*h^3*u_yZ+(3*u^2*h^5-h^3)*u_y*u_Z

W1=sum r^(k-1)*sin(k*psi)/k
T1=2*q*h*W1
T2/q^2=2*psi+sum c_k*sin(k*psi)/k
c_1=4*r
c_k=2*(k-1)*r^(k-2)+(6-2*k)*r^k, k>=2
~~~

The M48 tails bound the actual analytic series and both first and second r partials with closed geometric sums. Tail yZ is assembled from tail_r*r_yZ+tail_rr*r_y*r_Z; the radius supremum is never differentiated as a field. All powers evaluated at zero are nonnegative. The original t rational denominator stays positive, and the inverse denominator1+t^2 has lower bound1.

At fixed true phi, both inverse cross terms and curvature remain:

~~~
D=1+t^2
psi_yZ=-T2_yZ/D
       +2*t*(t_y*T2_Z+t_Z*T2_y)/D^2
       -2*t*t_psi*T2_y*T2_Z/D^3
~~~

The complete total T1_yZ is combined before forming B_yZ. The E-weighted product includes E_yZ*T1, E_y*T1_total_Z, E_Z*T1_total_y and E*T1_total_yZ. All five leading density products retain their actual E/V jets. Exact period/halfperiod traces are applied only at singleton phase coordinates; the whole zeta interval crossing zero is never declared exact-midplane.

Before repeated Fourier jet products, each finite radial offset is enclosed once into its directed coefficient, with offset anchor0. Every original Pstar/Cstar/L power remains formal. This prevents repeated subtraction of independent copies of the same radial-offset hull from artificially inflating mixed bounds. The checker verifies source containment of all three anchored derivative rows. The original source functions and large derivative factors are unchanged.

The accepted actual-N coefficient, nonlinear remainder and own-rate C0/Z integration programs are reused. In particular Q_Z=E_Z*A^2*R2(A/N)+E*A*A_Z*exprel(A/N), both coefficient orders and all product rules remain. Actual original endpoint phases, global endpoint terms, same-source internal telescoping, positive own-rate masses, radius Jacobian and pressure zero-rate memory remain. Incoming histories and P0 are not reset. Direct and phase-averaged bounds are both retained, with the tighter valid bound selected.

Four files with stem **lei_ren_part1_paper_compliant_current_original_reference_near_midplane_mixed_C1_integrals** are published under experiments/root_st073: producer, deterministic compressed manifest, checker and receipt.

## Result and limitations

Three actual levels were computed: (4 cells,N160),(16,N160),(16,N16384). At N16384/16 cells, outward-rounded whole-neighborhood contribution bounds are:

| Row | C0 averaged bound | Z/Lambda0 averaged bound | Z/Lambda0 direct bound | Z tightening, approximately |
| --- | ---: | ---: | ---: | ---: |
| m |5.685e-10|0.001532318|6.380501|4164 times|
| h |1.170e-9|0.001055892|1.293145|1225 times|
| k |6.209e-10|0.001854987|4.983273|2686 times|
| e |1.351e-9|0.001355490|1.588825|1172 times|
| p |2.087e-9|0.001817067|3.260120|1794 times|

**Lambda0=Pstar^11*Cstar^10**. Z bounds must be multiplied by this actual astronomical factor. All variable-L powers are retained in the bounded normalization; this differs from the predecessor's Lambda0*L^-2 reporting convention. The table describes source-window contributions, not full terminal defects or NS momentum residuals. N16384 is a candidate local frequency, not a global admission.

Producer36.500s, checker63.141s. Evidence:19 exact original radial derivative identities;60 independent original C0 y/yZ derivatives with a consistent nonzero even pressure remainder;96 independent defining T1/T2 integral and A/B jet comparisons at negative/zero/positive u;9 exact period/halfperiod traces;actual constant q/nu source binding;2 smooth u-chain identities,2 implicit mixed identities and2 geometric tail sum identities;3 finite-offset source containment checks and3 frame/domain rejections;36 whole regular mixed cells and15 own-rate mass comparisons, pressure mass5;actual N-specific nonlinear coefficients and direct/averaged selection. **1224 Git-index dependency hashes PASS**. No ancestor producer or unrelated test suite was rebuilt.

This closes the current near-midplane mixed-source/regular-primitive/C1-averaging dependency. Full-Z reference coverage still needs regular/signed overlap and formal axial interval handling. All17/24 units, five actual controls, global N, global joins/heat, stress/flat remainder, n-dependent recursion, both pulses and corrected Cartesian NS remain open.

## Ordered executable queue

Read this before historical handoffs. Implement the first unfinished dependency and mark DONE only with code, source-scoped receipt and commit. Preserve the full objective and predecessor queues. Repeat accepted checks only when relevant inputs or behavior change.

- [x] **REGULAR-y/yZ-SOURCE:** original templates, radial factors, source powers and pressure remainder jets on the whole signed zeta window.
- [x] **REGULAR-FOURIER-MIXED-PRIMITIVES:** M48 analytic first/second tails, both implicit inverse cross terms, no r^-1 through zero, exact singleton traces.
- [x] **COMPLETE-E/V-MIXED-PRODUCTS:** full E/E_y/E_Z/E_yZ and V products before density and integral bounds.
- [x] **WHOLE-NEIGHBORHOOD-C1-AVERAGING:** actual-N C0/Z IBP, endpoints, nonlinear terms, own rates and direct/averaged selector on the whole near-midplane window.
- [x] **FINITE-OFFSET-ANCHORING:** source containment with original parameter powers retained; repeated radial-hull inflation removed from Fourier jet products.

- [ ] **NEXT TYPED-AXIAL-SOURCE-CELLS:** generalize the accepted exact carrier collector to signed axial cells with an explicit physical-Z map, original family and defining L/Q. Preserve ordinary Z/yZ derivatives. Allow source intervals beyond the current central window while retaining their formal positive lower factors; do not turn an astronomical physical endpoint into native0. Acceptance: actual original C0,y,Z,yZ source enclosures with sign/domain certificates and positive L/Q, no sampled sign or independently selected field caps.
- [ ] **REGULAR/SIGNED-OVERLAP:** adapt radial/axial partitioning to the actual original u guard. Prove a nonempty same-source overlap where regular |u|<=1/4 and signed |u|>=1/8 both hold. In that overlap compare actual C0/Z/yZ roots and defining primitive/implicit derivatives, not just plotted values. Keep positive rho,s,hinv and exact source phase.
- [ ] **WHOLE-Z-REFERENCE-FUNCTION-ORACLE:** connect the central factored neighborhood to positive/negative physical Z all the way to[-1,1]. Native rational endpoint samples cannot cover the astronomical gap between zeta/Lambda0 and ordinary Z; use typed formal endpoints or source inequalities/interval maps with complete coverage. Combine regular and signed primitive services under one C1 API, fail closed on unsupported cells, and retain source-backed seam contracts. Acceptance: functional full-Z reference values, derivatives and original-N integrals.
- [ ] **SOURCE-CORRELATION-WHERE-REQUIRED:** when overlap or nonlinear mean cancellation needs joint pressure/L information, retain the defined shared remainder/jet relationships and exact L map in a stronger representation. Current independent hulls cannot be relabelled exact correlated functions.
- [ ] **WHOLE-ROUTE-FREQUENCY/ORACLE-ERROR:** include Lambda0 and all other original source factors in all-N C1 constants, endpoint and numerical/oracle errors, compatible with the accepted centered all-N/unit-C1-ball bridge. Derive admissible symbolic N conditions from actual complete route bounds, not160/16384 samples. Do not materialize impossible native exponentials.
- [ ] **SIGNED-NONLINEAR-MEANS/DRIFT:** integrate actual quadratic/exponential means with original phase measure and controlled y/Z drift where terminal closure requires signs. Leading reflection does not cancel nonlinear means. Keep all finite-N terms and error provenance.
- [ ] **O2-NONCONSTANT-PARAMETERS:** implement actual a_y,q_y,nu_y and mixed source jets. Use F=psi+T2-2*pi*phi*nu for the implicit equation, retaining -2*pi*phi*nu_i and nu_yZ. Include a derivatives in A/B. Acceptance: independent defining source/inverse comparisons and complete density products; the constant-reference routine must keep rejecting this case.
- [ ] **AXIAL/11-UNIT-BUFFER:** recover original source rows, nonconstant parameters, inherited moment/pressure histories and every offset. Integrate genuine densities with explicit seams/initial data; do not reset incoming corrections.
- [ ] **O3/RESTORE/ALL17-24-ORACLE:** execute the predecessor unit queue, recover all definitions/histories, and publish complete same-N C1 function/integral outputs. Unsupported units fail closed.
- [ ] **ACTUAL-FIVE-CONTROLS/TERMINAL-PICARD:** connect the full actual oracle to the centered bridge, solve all five terminal identities as Z functions, admit a global N and replay controls with oracle/Picard errors separated.
- [ ] **JOINS/EXACT-HEAT/FINITE-ENERGY:** assemble the original core/annulus/flatten/collar/exterior with analytic pressure compatibility, required smooth joins and genuine heat/radial tails.
- [ ] **ADMISSIBLE-STRESS/FLAT:** separate the actual background residual into -div(T_B)+E_B, establish cone margins in all original regions and genuine flat decay independently.
- [ ] **ACTUAL-n-RECURSION:** execute distinct n=1/n>=2 recovery equations, common core interval, independent moment repairs, divergence-preserving cutoffs and controlled smooth summation.
- [ ] **TWO-FAMILIES/CORRECTED-UVW:** realize both oscillatory pulse families and averaged quadratic stress cancellation, recover corrected Cartesian u/v/w and independently validate forced NS residual plus physical contraction, elongation and material winding.

The full original reconstruction goal remains active.
