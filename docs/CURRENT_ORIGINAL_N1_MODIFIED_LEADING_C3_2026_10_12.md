# Finite-N corrected leading velocity and own-mean C3 interface

The pressure-tail route now has a source adapter for the **actual finite-N corrected leading velocities**, through ordinary Z3, and the corrected own axial mean increment at chart endpoints. It reuses the accepted original C3 inverse, finite-N density and own-rate transport graphs. It keeps the same N, source family, common velocity unit and global phase. The full reconstruction goal remains **ACTIVE / INCOMPLETE**.

Read this before [the n1 radial-cutoff checkpoint](CURRENT_ORIGINAL_N1_RADIAL_CUTOFF_2026_10_11.md). Three n1 velocity moments already have their full support values. Completing the other two requires the fixed leading pressure source outside Rin. This adapter closes the high-Z velocity/endpoint-history input layer; mixed radial/phase derivatives and the complete pressure-tail source are still open.

## Why the outer leading source needs this adapter

The selected original Rin has first phase Sin=sc/4, while the finite-N modification begins at r_minus, phase sc/2. Before r_minus the leading field is unchanged. To its right, the prescribed-shear background alone is insufficient: the actual leading field uses the original A/B loop and modified own histories.

The draft `current_original_n1_leading_forcing_tail_prefix.py` now refuses prescribed-shear queries beyond r_minus. It is included as an **unaccepted integration prototype**, without a numerical acceptance receipt. It must not be used as a completed finite-N tail provider, and phase 1/2 cannot be substituted for the unchanged-left endpoint.

## Signed corrected functions and their original units

The accepted original C3 graph defines F_N=N*(E_N-E), the original angular coefficient increment, and B, the axial coefficient increment. For the same original N,

```text
E_N = E + F_N/N = E*exp(A/N)
V_N = V + B/N
delta_m_at_chart_right = outgoing_N_scaled_m/N
delta_m_at_chart_left  = incoming_N_scaled_m/N.
```

Every row is an ordinary Z derivative through order three, with the original derivative/factorial conversion applied once upstream. E/V remain in their existing common velocity unit S. The physical axial mean increment is S*delta_m. The generic normalized history is m=Mz/(R*S); no additional radius Jacobian or Pstar conversion is applied to the already transported history.

The m transport has original rate one:

```text
partial_y delta_m = delta_V-delta_m, y=logR relative to the original origin
delta_m=0 at r_minus
delta_m(y)=exp(-(y-y_left))*delta_m(y_left)
           + integral_y_left^y exp(-(y-s))*delta_V(s) ds.
```

The same predecessor history and global phase are retained across every chart. A quiet chart has zero local velocity increment but can retain a nonzero incoming mean increment. It must not reset the latter to zero.

The bridge installs signed function handles for all 17 original charts. Its numerical receipt exercises the first- and second-micro right endpoints on the original whole-Z cell [0,1/2]. Directed C3 magnitude bounds use the accepted source ranges and outgoing histories, uniformly for N>=N0. These are enclosures of signed functions; magnitude caps do not become point values or derivatives of caps.

## Scoped evidence and entry

The independent exact fixture checks 24 finite-N velocity and mean Z0…Z3 rows, the physical radial-mean/rate-one identity, recovery of the radial divergence relation, and nonzero incoming-memory decay on a quiet chart. Source/function/scope mutations and unsupported queries are rejected. The existing read-only GPT-5.6 Luna / max review confirmed the finite-N density identities, rate-one m transport and common-S units.

```python
owner = CurrentOriginalN1ModifiedLeadingC3(accepted_C3_range_owner)
packet = owner.endpoint(('0', '.5'), 'first_micro')
view = owner.report(packet)
```

Module: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_modified_leading_C3.py`.

Scoped gate: `current_original_modified_finite_N_velocity_Z3_and_own_mean_endpoint_increment_installed`.

This supplies endpoint mean increments and corrected velocity function handles. It does not supply a whole-cell modified mean oracle, physical Cartesian velocity dispatcher, mixed radial/phase derivatives, complete Q_N/N1p_N, the two complete pressure moments, or the n1 repair.

## Detailed next work

- [x] Identify the original signed finite-N F_N/B C3 functions and N-scaled own mean transport instead of reusing unmodified roots alone.
- [x] Expose corrected E_N/V_N and delta_E/delta_V through Z3 using the existing N handle.
- [x] Divide incoming/outgoing N-scaled m by N once, preserving rate-one memory and common-S units.
- [x] Connect accepted C3 magnitude ranges to the signed corrected functions and current micro endpoint queries.
- [x] Prevent the draft prescribed-shear tail path from continuing beyond the actual modification origin.
- [ ] Build the same signed partial Duhamel mean function inside each chart: retain its incoming history and integrate B/N with the original rate-one kernel up to the queried native endpoint. Use whole source-cell enclosures; the current endpoint history is not a whole-cell oracle.
- [ ] Add the actual slow-y and phase derivatives of the finite-N B primitive through the mixed Z grades consumed by Q_y. Use the original global phase phi=frac(N*y), including its chart-specific Jacobian. Differentiate original defining functions, not range caps.
- [ ] Supply delta_V_y=B_slow_y/N+B_phi in ordinary logR units, and its Z derivatives. Verify the original periodic phase/seam identity before passing the phase wrap.
- [ ] Recover the own mean radial row from delta_m_y=delta_V-delta_m, and mixed Z rows from this same identity. Preserve all incoming histories on micro/macro and quiet joins.
- [ ] Recover Q_N from the actual corrected velocity/own mean: Q=(2Z*V-(1-delta)Z*m-(1-Z^2)*m_Z)/(1-delta*Z^2), in common units. Q_Z2 needs m_Z3; Q_y/Z1 also consumes the corrected radial velocity rows.
- [ ] Bind baseline and correction to the same original physical source frame; restore common-S/radius prefactors exactly once. P0 stays the independent original datum.
- [ ] Assemble the actual finite-N S=T_-2(q)+q*(q/2+R*q_R)+Uz0*Z_-2(q), then integrate S and R*S and use the exact finite-endpoint viscous identities. Do not drop upper endpoint terms before the actual leading source-support proof.
- [ ] Append those contributions to the n1 cutoff m13/m15 prefixes. Finish the remaining leading charts and actual zero-source support, then admit the two complete defects.
- [ ] Bind the positive power interval/Astar and actual normalized 2+3 bump matrices; solve the n1 repair with the complete signed defects, reconstruct M/V/P, and verify all five terminal identities and pressure compatibility.
- [ ] Close derivative joins/axis remainder, genuine n2, stress admissibility and oscillatory/flat-remainder/corrected-NS layers only at their full stated scope.

Keep all full-inner/solved/repair/completed-leading-support/temporal-recursion/global-stress/flat-remainder/full-corrected-NS flags false until their actual contracts are proved.
