# Original pressure and cumulative m15 through R110 - 2026-10-16

`OriginalR110PressureRows` and `OriginalR110M15` now continue the actual pressure source and cumulative m15 across the second switch and post-power interval to exact R110. The complete51-sector first-switch inlet is retained once; the second-switch packet has69 terms and the post-power packet has88. Pressure fields retain U/Q ordinary Z2, genuine meanZ3, Uy/Qy ordinary Z1, separate normalized epsilon sectors and original selected-N logarithmic errors.

Full goal: **ACTIVE / INCOMPLETE**. Tasks14/15 remain IN PROGRESS. This is a local cumulative pressure construction with signed proxies and complete factored absolute errors. Useful total relative cancellation, resolved selected-N signed centers, independent Rp value/Z1 cancellation, global five-moment closure/common N, heat/stress and genuine temporal n>=2 remain open.

## Genuine inlet and exact source ownership

The pressure adapter composes the actual first-switch five-jet source, then the actual second-switch source. It does not call the older public r110 mean adapter as a replacement for those data. Isolated chart views share immutable source objects and exact Taylor contexts; their public numerical caches remain separate. This avoids explicit base-method dispatch in the older r110 velocity implementation from truncating/resetting the pressure inlet.

The first-switch cumulative51 terms come from the checked, hash-pinned first-switch m15 report. Genuine incoming Psi/Mpsi five-jets come separately from `first.rows.inputs(Z,iv.mpf(1))`. The second/post moment inlets are not substituted by comparison histories or by the older public six-history output.

Generic derivative owners are **switch_second** (source domain1..2) and **switch_power** (0..1). The transport chart labels second_switch/post_power are different names and cannot select those bounds by themselves. All earlier micro/macro/switch histories and every axial frame survive. The common selected N, Pstar and fast B_phi/B_Z_phi units remain unchanged.

## Native geometries

For the second switch, local s in[0,1] has absolute switch phase1+s:

```
R=100*exp(hbS*(1+s)), dy=hbS ds, R*dR=R^2*hbS ds.
```

The source combines angular shear1-sigma(s) and angular power sigma(s). The original U is constant and Uy=0. Its local source-width interval is an error cover, not a selected physical parameter. hbS retains the original checked upstream K1/logh source identity and remains a distinct graph node from hbB.

For post power, q in[0,1]:

```
L=log(11/10), Lp=L-2hbS>0,
R=100*exp(2hbS+Lp*q), dy=Lp dq,
Phi_y=-2Phi/5, Uy=0, R(q=1)=110 exactly.
```

There is no additional hbS on the post-power measure. The physical radius is handled in logarithms, and R110 is evaluated directly at q1 so omitted tiny-width terms cannot spoil the exact endpoint.

## Mean transport, Qy and viscous cancellation

On both original charts, U is constant and My=U-M. The genuine mean therefore obeys

```
M(y)=U+(M_in-U)*exp(-y).
Q=(2ZU-(1-delta)ZM-(1-Z^2)M_Z)/(1-delta Z^2).
```

The pressure source retains the nonzero incoming mean and its derivative jets. Qy is recovered from this same mean ODE; it is not set to zero because Uy=0. A separate uniform-prefix background encloses every partial mean kernel when bounding finite-N volume products. Endpoint background values alone would not provide that uniform cover on the post-power interval.

Since Qy is the decaying mean mode and R grows by exp(y), **Delta(R Qy)=0 exactly for the unmodified prescribed source**, including any correlated incoming mean. The code uses this exact background identity. It still retains both finite-N viscous endpoints, both finite-N half-S1 endpoints, all volume products and quadratic cross errors.

## Signed half-S1 and native-width errors

For the original boundary/volume polynomials F/G and inherited mean-mode degree n,

```
H_dn=(2-n)*F_dn+G_dn,
half DeltaS1 = integral R^2*(2F+G+F_y) dy.
```

Each physical atom is anchored at R100 and retains `100^2*epsilon^d*exp(-nY)`. The signed inherited mode is not reconstructed by subtracting rounded R110 mean values. There are six modes d0..2,n0..d.

On second switch, the signed primitive retains hbS*s*H. The exp((2-n)*hbS*(1+s)) correction is bounded with a separate hbS^2 atom; n2 is exactly resonant. The actual field discrepancy is normalized before epsilon:

```
DeltaPsi <= hbB*bp + hbS^2*P,
DeltaM   <= hbB*bm + hbS^3*(1+s)*P.
```

Here P is the actual first-switch before-width source primitive and bm includes the full source-correlated macro inlet. Native integration adds one hbS: the source errors retain hbS^3 and the inlet errors retain hbS*hbB. Z2 jets are kept until C1 density derivatives are formed.

For post power, the proxy uses positive kernels J_k=expm1(k*L*q)/k, J_0=L*q. The omitted absolute offset2hbS and the actual Lp measure are covered together by a separate hbS geometry sector, with interval-directed 2/L. The mean source error retains hbS^2*P*(2hcap+1-exp(-Lp*q)); no cap is selected as the field definition. Both inlet and source error products keep their quadratic terms.

At each local coordinate0 every new increment/error coefficient is exactly zero. The complete upstream packet is returned unchanged.

## Executed evidence

The full run exited0 in111.881 seconds; saved-report inspection exited0:

- 38 actual second-switch/post-power pressure/m15 packets and38 exact complete inlet replays;
- 532 ordinary physical pressure rows,1406 serialized increment-error conditions and222 exact empty rows;
- 128 retained lower finite-N endpoint conditions,114 zero-Uy/native-measure conditions;
- 18 exact source-seam/R110 endpoint conditions and10 rejected invalid queries;
- eight value/Z1 identities executing the actual second-switch source body and three original post-power source identities;
- 48 independent physical half-S1 value/Z1 integrals,288 per-mode native geometry-error conditions;
- 48 nonzero-viscous-endpoint cancellation comparisons and48 source-mean Z0..3 error conditions.

The30,773,448-byte report SHA256 is `7226a0546b338c6834d035f6daa86707091c4c1e2901f631ea8a2351231bada1`; it pins1762 inputs and57 implementations. Inspect with `python scripts/run_original_r110_m15.py --inspect-report`. It includes the actual pressure row packets as well as cumulative moment packets. Do not rerun accepted ancestors.

Existing read-only `/root/n1_atlas_map`: **GPT-5.6 Luna / max**; no new agent. Its review found no material underbound and required preservation of the existing hbS/hbB source-identity guard. The signed kernels are anchored at R100; any future endpoint-anchored refactor must remove the offset from those kernels to avoid counting it twice.

## Next actual construction

1. Continue the original long reshape from the complete88-term post-power endpoint and genuine `post.rows.inputs(Z,iv.mpf(1))` five-jets. Keep the pressure inlet distinct from older public velocity/six-history adapters.
2. Use `OriginalLongNumeric` and its checked finite original T singleton, long profile and actual controls. Keep physical radius/axial amplitude/phase in their native normalized or logarithmic units; do not materialize the enormous radius or select an error cap as T.
3. Recover genuine meanZ3/QZZ and Uy/QyZ1 from the long source and its own positive-kernel mean equation, with the correct generic derivative owner and complete R110/long source prefix. The constant-U cancellation proved here must not be extended to a chart with varying U.
4. Integrate only the new long pressure interval with its actual native length, epsilon/mean/amplitude sectors, both leading/finite-N pressure endpoints and all quadratic products. Continue the subsequent physical transition/compact pressure toward Rm and Rp.
5. Measure independent source-correlated Rp value/Z1 cancellation; then complete global all-five histories/mixed joins/common N, original flatten/preheat/exact heat/finite-energy tail, admissible stress/flat remainder, genuine temporal n>=2/smooth summation, oscillatory stress correction and independent corrected Cartesian NS validation.
