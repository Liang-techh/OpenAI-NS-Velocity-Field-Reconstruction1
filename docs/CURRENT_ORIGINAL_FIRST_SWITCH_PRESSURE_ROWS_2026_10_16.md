# Original first-switch pressure source rows - 2026-10-16

`OriginalFirstSwitchPressureRows.evaluate(Z,s)` now provides the missing pressure rows on the original first-switch chart, s in[0,1]. It retains physical U/Q ordinary Z2, mean M ordinary Z3, Uy/Qy ordinary Z1, normalized epsilon sectors, the original factored Uy source and separate selected-N logarithmic errors. The actual R100 inlet uses the genuine five-jet macro source before the older velocity interface truncates the mean.

The full goal is **ACTIVE / INCOMPLETE**; tasks14/15 remain IN PROGRESS. This completes the first-switch pressure source/endpoint interface, **not its cumulative m15 integral**. The next integral must preserve the complete checked **28-sector** `OriginalMacroM15` q=1 inlet once. No new temporal coefficient order or global five-moment closure is claimed.

## Native switch width and controls

The checked source calls the switch parameter `hbS`. This is a native width parameter, not hb raised to a power. The source graph keeps bridge `hbB` and switch `hbS` as distinct nodes. Its parameter geometry guard proves the identical upstream K1 ledger/logh source and exact switch reboxing; the adapter does not identify parameters merely by their names or caps.

```
R(s) = 100*exp(hbS*s)
dy = hbS ds
R*dR = R^2*hbS ds
R1 = 100*exp(hbS)
```

The first-switch angular weight is1. Its axial weight is1-sigma(s). The actual prescribed source is

```
Phi_y = -(hbS/2)*R*(D/R)*Phi
U_y = -hbS*(1-sigma(s))*(Phi/barPhi)
      *(R*hydro + R*Pstar^2*pressure + R^2*F0^2*swirl).
```

The phase integral for U applies the separate `dy=hbS ds` measure. Angular and axial changes therefore have two width factors, whereas the first radial derivative Uy has one. Neither this source nor its control is taken from a frozen-macro continuation outside that chart.

## Genuine incoming jets and directed source errors

The retained macro chart view shares immutable source objects and the exact Taylor contexts with the switch owner, but has the real `OriginalMacroNumeric` method dispatch. Its q=1 evaluation cannot accidentally use switch geometry or coordinates. This supplies the actual nonzero incoming Psi/mean functions and genuine jets through order5.

The current-radius comparison modes are continued from their original R0 source into the switch. The source includes the macro angular exponential error, then the additional local switch bound

```
norm(raw_source_truth)*ell_switch*exp(ell_switch).
```

The full result is multiplied by the original directed1-sigma(s). The raw physical derivative uses a zero numerical proxy plus a complete width error; the signed coefficient is separately retained before `epsilon*hbS`. No numerical cap is selected as the physical switch width or the field definition.

The actual mean retains its full R100 inlet. Its change is bounded from the defining positive kernel by

```
hbS*s*(norm(Psi100-Mpsi100)+uniform switch Psi-change error).
```

The incoming large-Y mean relaxation is already part of this inlet. It is not repeated or reset on the short switch. The genuine M_ZZZ row is retained in Q_ZZ recovery, with its possibly large finite-N logarithmic uncertainty. A corrected small-error M_ZZZ/Q_ZZ value is not asserted.

`inputs()['geometry']['length']`, theta and mass are explicitly labeled **error covers only**. They must not be used as exact physical parameters in a subsequent integrator. Read the original hbS log atom and source identity; factor the unmaterializable width before evaluating a short-window kernel.

## Pressure derivatives and endpoint units

The physical pressure row recovery remains the exact common U/M mean-ODE construction. Its normalized epsilon sectors preserve small terms before addition to the O(1) affine field. The background finite endpoint is R*epsilon*Qpsi_y; the finite-N endpoint error is R*deltaQ_y. The latter receives no additional epsilon.

The derivative source is **switch_first**, with the complete first/second/macro/switch history and all axial frames. Ordinary y derivatives of the selected-N correction use

```
deltaU_y  = Pstar*(B_y/N + B_phi)
deltaU_yZ = Pstar*(B_yZ/N + B_Z_phi).
```

The fast phase term has no extra hbS and no extra N denominator. Slow y derivative normalization already accounts for the native phase-to-y conversion.

At s=1, the unmodified leading Uy source is exactly zero because1-sigma(1)=0. This does not zero Qy, its incoming mean contribution, or the finite-N Uy/Qy error. The source and corrected-field statements remain separate.

## Focused evidence

The focused run exited0 in83.88 seconds; saved-report inspection exited0. It produced19 actual switch packets,266 ordinary physical rows,570 serialized-error conditions,38 fast-phase conditions,15 exact macro-inlet/factored-derivative conditions, six exact terminal unmodified-Uy rows,114 width/control identity conditions,38 finite-RQy unit conditions and six invalid-query rejections. Repeated cached queries also passed.

Six value/Z1 identities execute the actual original first-switch source bodies. The actual phase U equation retains both width factors. Sixteen independent ordinary Z0..3 mean quadratures use a changing source and nonzero inlet. Existing read-only worker `/root/n1_atlas_map`: **GPT-5.6 Luna / max**; no new agent. Its review found no material underbound and required the explicit geometry-cover labels now present.

The1,487,379-byte report SHA256 is `6bcbd35931ddc35a3d743ae0a46c879f9479d4c6d77aba43ea3412d16881ae33`; it pins1754 inputs and42 implementations. Inspect with `python scripts/run_original_first_switch_pressure_rows.py --inspect-report`. Do not rerun checked ancestors.

## Next actual construction

1. Load the checked original macro q=1 **28-term** m15 inlet with its source/implementation pins. Keep it unchanged and integrate only first-switch0..s.
2. Use the original reduced half-S1 minus both RQy endpoints. Keep epsilon0/1/2 sectors separate, the actual hbS phase measure, nonzero mean relaxation and the changed1-sigma source.
3. Factor hbS before evaluating short exponential kernels. Retain higher width errors, rather than rounding exp(-hbS*s) to1 and dropping its effect. Geometry-cover arrays are not physical lengths.
4. Retain the actual leading Uy contributions at both endpoints, with the distinct source hbB/hbS atoms and their guarded identity. The frozen mean's constant RQy proxy identity cannot zero those source terms or the finite-N endpoints.
5. Continue second_switch/post_power to R110 and remaining active/compact pressure. Then measure source-correlated independent Rp value/Z1 cancellation; interval overlap or a terminal theorem is insufficient.

Global all-five joins/common N, original flatten/preheat/exact heat/energy, admissible stress/flat remainder, genuine temporal n>=2/smooth sum, oscillatory correction and independent corrected Cartesian NS remain open.
