# Original reference cumulative m15 and pressure to Rm - 2026-10-16

`OriginalReferencePressureRows` continues genuine U/Q ordinary Z2, mean ordinary Z3 and Uy/Qy ordinary Z1 across reference, restoration and postrestore to Rm. `OriginalReferenceM15` appends the reference interval to the unchanged106-term long endpoint, producing124 cumulative sectors. **Restoration/postrestore cumulative m15 is not yet integrated.** Pressure data reaching Rm does not complete the pressure moment there.

The full goal remains **ACTIVE / INCOMPLETE**; tasks14/15 remain IN PROGRESS. Selected-N signed centers/useful relative cancellation, the independent Rp seam, all-five global histories/common N, exact heat/energy/stress/flat remainder, genuine temporal n>=2 and corrected NS remain open.

## Genuine source and changing baseline

The source uses the same original long object and exact Taylor contexts. Genuine five-jets are continued before any older public velocity/six-history truncation. Generic derivative owners are `inner_reference`, `axial_restore`, and `restore_buffer`; the last has source offsets[-7,-6] while its wrapper coordinate p runs0..1. The complete upstream prefix is retained with reference measure Lref*dq and later dt/dp once.

Write U=U0+epsilon*Psi and M=M0+epsilon*Mpsi. The j baseline is restored separately from the small epsilon sector:

```
reference:
  U0=M0=4Z+j; Psi=Psi_R110
  Mpsi=Psi_R110+(Mpsi_Rsh-Psi_R110)*exp(-gap)
  Uy=0, gap=Lref*q

restoration, alpha=1-sigma(t), K=integral_0^t exp(-(t-s))*alpha(s) ds:
  U0=4Z+j*alpha; M0=4Z+j*(exp(-t)+K)
  Psi=Psi_R110*alpha; Mpsi=Mpsi_Rz*exp(-t)+Psi_R110*K
  U0_y=-j*sigma_prime; Psi_y=-Psi_R110*sigma_prime

postrestore:
  U0=4Z; Psi=Uy=0
  M0=4Z+(M0_restore1-4Z)*exp(-p)
  Mpsi=Mpsi_restore1*exp(-p)
```

The independent numerical quadrature is a point proxy. Its truth enclosure comes from the original positive restoration kernel cells; a long decimal query retains its exact interval coordinate rather than substituting the rounded proxy coordinate. Exact chart0 inlets preserve the actual incoming jets.

Each sector is recovered with A=1-delta*Z^2, D=1-Z^2:

```
Q   = (2Z*U-(1-delta)Z*M-D*M_Z)/A
Q_y = (2Z*Uy-(1-delta)Z*(U-M)-D*(U-M)_Z)/A
```

Mean Z3 is required for QZZ. The retained five-jets also provide MZZ for QyZ; no missing row is zero-padded. `physical_pressure_rows` are ordinary derivatives; the extra `baseline_pressure_rows` are explicitly Taylor coefficients D_Z^k/k!. The RQy export separates R*Q0_y, R*epsilon*Qpsi_y and the finite-N R*deltaQy error. The latter has no extra epsilon, N denominator or hb multiplier on its fast terms.

## Finite reference geometry and coefficient correlation

```
logref=10*(logC+logP)
Lref=logref-T-8 > 0
y=T+Lref*q, R=110exp(y)
restoration: y=logref-8+t
postrestore: y=logref-7+p; p1 is Rm
```

The unit offsets -8,-7,-6 remain separate log atoms. The reference fraction uses T*(1-q),10logC*q,10logP*q,-8q; native `g:gap` uses T and gap. No enormous physical radius, exp(T), or exp(Lref) is materialized.

For H_dn=(2-n)F_dn+G_dn and k=2-n>0, the reference increment is the shared original long coefficient H_dn/k times

```
exp(k*(T+gap))-exp(k*T),
```

with common110^2*epsilon^d*exp(-n*(Y+log(11/10))). The new lower coefficient and log atom are copied exactly from the accepted long upper endpoint with opposite coefficient sign and the same correlation ID. The new upper starts from the long lower common atom and appends the current absolute offset atoms. This avoids counting T twice and preserves internal endpoint cancellation. For n2, the coefficient is copied unchanged and log(gap) is a separate positive measure atom.

The constant-U source error envelope remains valid under the exact positive mean relaxation: each Z derivative obeys eM'=eU-eM, and the existing long M bound covers its constant U bound. Only the known integrated mass1/2 is rescaled to min(gap,1/2); no original width or source parameter is selected from a cap. Source/inlet errors retain hbB and hbS^2 factors and all quadratic products. Finite-N endpoint, volume and viscous errors retain both physical endpoints and the full-prefix uniform background. Only the unmodified constant-U Delta(RQy) cancels; finite-N errors do not.

## Acceptance receipt

The focused run exited0 in137.512 seconds; saved-report inspection exited0. The report contains19 reference124-term m15 packets and61 reference/restoration/postrestore pressure packets, including independent Z=.17 and a100-digit restoration coordinate.

- 19 unchanged106-term inlet replays;854 ordinary physical pressure rows and61 nonzero meanZ3 finite-N error checks.
- 684 serialized increment-error conditions;108 exact empty rows;320 exact shared long coefficient rows;95 internal telescope atoms;16 positive resonance atoms.
- 64 retained lower finite-N endpoint conditions;126 source-seam overlap conditions;55 separate unit-offset conditions.
- 63 actual restoration baseline conditions;21 independent numerical mean-kernel inclusions;61 complete-prefix conditions; eight rejected invalid queries.
- Three identities executing the original restoration/postrestore source,18 independent shifted physical pressure integrals,90 signed kernel comparisons,90 internal endpoint cancellations,12 mean-error Z0..3 envelope conditions and two exact mean/viscous identities.

The report is27,355,240 bytes; SHA256 `2531c3bcabe67f47727175ba0fd08364ac8f2d45a78569a7b90366da687f8064`. It pins1767 input sources and70 implementations. Inspect with `python scripts/run_original_reference_m15.py --inspect-report`; do not rerun accepted ancestors. These checks establish the stated adapters and error enclosures, not a global five-moment closure or a blow-up theorem.

Existing read-only `/root/n1_atlas_map`: **GPT-5.6 Luna / max**, reused without spawning another agent. The worker confirmed source sectors, derivative orders, chart measures, the signed telescope and both finite-N endpoints; it identified the constant-U mean-error invariance proof obligation, covered by the explicit source binding and independent propagation checks. Its next-step advice is reflected below.

Two terminal development attempts are retained in the external recovery checkpoint: a wrong prerequisite filename, then an invalid nonzero assertion for the odd baseline Q0_y at Z=0. The first corrected the source filename to the accepted reference_velocity report; the second changed only the driver assertion. No assertion was weakened around a physical mismatch.

## Next implementation: changing-source cumulative pressure to Rm

1. Start with all124 accepted reference terms at q1. Preserve exact selected N and the genuine reference-endpoint pressure/mean jets; do not reset the moment history.
2. For restoration use R=Rz*exp(t), t0..1. Keep epsilon degrees0,1,2 separate. The degree0 coefficient contains the changing j*alpha and j*gamma baseline, including j and j^2; it cannot use the earlier constant-U frozen coefficients.
3. Use the exact identity integral R^2*(2F+F_y)dt=[R^2 F], leaving the actual integral R^2*G dt and both RQy endpoints. Keep all signed endpoints with source correlations; do not extend the constant-U viscous cancellation.
4. Build positive scalar cell enclosures for alpha, exp(-t), and gamma=exp(-t)+K. A recurrence with alpha(a),alpha(b) and exp(-(b-a)) encloses gamma at endpoints; the monotone whole-cell enclosure avoids nested kernel quadrature. A point quadrature is only a proxy and must retain its cell/truncation error.
5. Compute C1 density coefficients using U/Q through Z2 and M through Z3. Keep the physical Rref^2 base atom, offsets -16/-14 and local growth atoms separate, then apply epsilon^d once. Preserve all source/inlet uncertainty and quadratic cross terms; add actual finite-N upper/lower endpoint, uniform-prefix volume and viscous errors.
6. Continue postrestore p0..1 from the actual restoration endpoint. U=4Z is constant, but Qy/mean history is nonzero and the moment inlet is the full restoration packet. Verify the endpoint primitive and exact native dp measure through Rm.
7. Continue active/compact pressure toward Rp and independently evaluate the source-correlated value/Z1 seam. Complete all-five histories/global joins/common N, original flatten/preheat/exact heat and finite-energy tail, admissible stress/flat remainder, genuine temporal n>=2/smooth summation, oscillatory correction and independent corrected Cartesian NS.
