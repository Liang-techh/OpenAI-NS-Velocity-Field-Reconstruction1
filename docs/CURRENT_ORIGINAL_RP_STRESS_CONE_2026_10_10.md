# Original strict admissible stress cone at the supported point — 2026-10-10

The actual original pulse_exit background now satisfies the paper's two-component admissible stress cone **at the supported point**. Directed source bounds prove a>0, Utheta>0, v_s>2 and both strict signed cone margins positive. This is the first accepted current-Rp cone point, following the full stress/remainder decomposition. Full goal **ACTIVE / INCOMPLETE**; regional cone, flat remainder and genuine coefficient recursion remain open.

Source commit: [e4a6bd56](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/e4a6bd560a9e86fe4cca53d3ddcddc1390d93cfd). Source quartet: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_stress_cone.py`, `.json.gz`, `_check.py`, `_check.json`.

## Actual original cone and positive factors

The checked delivery is pulse_exit, xi=21/2, Z=371/1000, finite log-time offset=1/5, theta=7/10. The accepted original stress owner supplies Ttheta=Itheta+Stheta and Tz=Iz+Sz, both with the same positive physical prefactor lambda**(-2-delta). Original viscosity is 1. The actual incoming histories, five cumulative moments, independent P0 and absolute pressure remain unchanged.

OpenAI equations (4.20)-(4.23) use the following cone coordinates:

```text
F = Utheta / sqrt(2R) > 0
a = 1 - 2*Utheta_y/Utheta = 2 + 2*mu
b_s = 2*Uz_y/Utheta
t_s = -b_s/a
v_s = a*(1+t_s**2)
G1 = Ttheta + t_s*Tz
G2 = 2*G1**2 - (v_s-2)*(Tz-t_s*Ttheta)**2
```

Admission requires a>0, v_s>2, G1>0 and G2>0. The completed diagonal T_thetatheta=r*partial_z(Tz) ensures the tensor-divergence identity and is **excluded from these two-component cone margins**.

The adapter avoids interval division by Utheta. It collects the following exact shared-source polynomials before enclosure:

```text
B = 2*Uz_y
H0 = (a-2)*a*Utheta**2 + B**2
   = a*Utheta**2*(v_s-2)
H1 = a*Utheta*Ttheta - B*Tz
   = a*Utheta*G1
Hperp = a*Utheta*Tz + B*Ttheta
H2 = 2*a*Utheta**2*H1**2 - H0*Hperp**2
   = a**3*Utheta**4*G2
```

All multipliers are proved positive using the actual original mu and swirl sources. The original theta derivative law is inherited from the accepted axial owner. Ttheta and Tz retain their actual source-scale factors and nonlinear histories. No stress target, independent rounded ratio or conservative cone cap replaces these values.

## Actual checked result

- a, Utheta, H0, H1 and H2 all have strictly positive source enclosures; no scale ratio is unresolved.
- The independent checker derives the positive-factor equivalences, checks the actual original margin definitions and common physical stress prefactor, and includes 347 collected source-product coefficients and 5 signed output bounds with directed 800-digit arithmetic.
- Total 0.1% width flags are reported separately: `{'a': True, 'Utheta_positive': False, 'H0_v_minus_2': False, 'H1_first_margin': False, 'H2_quadratic_margin': False}`. Cone admission uses strict signed margins; it does not promote the preceding 25 nonzero physical stress/remainder width flags.
- Only a=2+2mu, a finite dimensionless parameter, is ordinarily materialized. Utheta and the three nonzero cone margins remain source-scale values. No huge physical velocity/stress number is expanded.
- Copied packets, invalid targets, changed cone source DAGs and a substituted quadratic margin are rejected. The accepted constructor and public evaluate/report calls pass.
- One existing read-only reviewer is reused: **GPT-5.6 Luna / max**. No additional or Astra worker is spawned.

The point-cone flag is distinct from all inherited regional pulse/end cone, uniform/global/axis, flat-remainder, recursion and corrected-NS flags, which remain false. A strict point result does not establish a cone throughout a region or a blow-up theorem.

## Executable use

Create the live accepted background and field using [the background checkpoint](CURRENT_ORIGINAL_RP_BACKGROUND_STRESS_2026_10_10.md#executable-use). Append:

```python
from lei_ren_part1_paper_compliant_current_original_Rp_stress_cone import CurrentOriginalRpStressCone as StressCone

cone = StressCone(background)
point = cone.evaluate(field, '1/1000')
report = cone.report(point)
assert report['current_original_Rp_admissible_stress_cone_at_supported_point_certified']
assert report['cone_margins']['H2_quadratic_margin']['signed_log_value']['sign'] == 1
assert not report['regional_cone_certified']
```

Use the original live typed chain; JSON reports cannot hydrate new owners. Prepare new original queries before freezing source views, and do not reset source guards.

## Completed bounded tasks

- [x] **CURRENT-RP-SIGNED-CONE-MARGINS-POINT** Bind the original a/b_s and Ttheta/Tz sources, remove only strictly positive denominators, retain original pressure/incoming histories and prove all strict margins positive at the supported point.
- [x] **CURRENT-RP-ORIGINAL-POINT-CONE-EQUIVALENCE** Independently derive H0/H1/H2 positive-factor equivalences and check actual shared-source product/output bounds without adding the completed diagonal to the cone.

## Next executable tasks

- [ ] **CURRENT-RP-FLAT-REMAINDER-SCALE-CANCELLATION** Use the accepted separate E_r/E_theta/E_z from the background packet. Partition exact lambda/R/source scales and form ratios to the original reference decay scales before interval expansion. Retain original time offsets and Jacobians. Acceptance: controlled relative decay at actual deliveries; distinguish point, finite-order and uniform flatness.
- [ ] **CURRENT-RP-OFFCENTER-AND-TIME-CONE** Prepare independent nonzero Z/pulse sections and exact log-time offsets before freezing owners. Evaluate original derivative, axial, background and cone packets. Independently check each query and compatible cross-query scale relations. Acceptance: actual signed margins and measured velocity/vorticity/stress/remainder scale changes at multiple deliveries.
- [ ] **CURRENT-RP-BACKGROUND-RELATIVE-WIDTH** Separate coefficient and centered scale width errors for the 25 nonzero physical background outputs. Refine only dominant original dependencies and retain exact mu/delta/logC, incoming histories and pressure. Report included narrower bounds and actual 0.1% pass counts; no source replacement or flag reset.
- [ ] **CURRENT-RP-PRESSURE-END-AND-ALL-CHART-FIELD** Resolve the xi=13 remaining-pressure lower-bound guard, then implement original postpulse, flatten, steep, waiting, collar and exact-heat derivative/stress/cone adapters. Preserve joins, pressure compatibility and exact outside-heat definitions.
- [ ] **CURRENT-SIGNED-STRESS-CONE-FLAT-REGIONS** Extend controlled cone margins and separate remainder bounds across inner exit, matching annuli, pulse/end, flatten and collar. Use regional error control and original physical volume Jacobians for max/L2 norms. Sampled positive points do not close a regional task.
- [ ] **CURRENT-PHYSICAL-MATERIAL-WINDING-ENERGY** Integrate actual trajectories and cumulative material winding at multiple times separately from instantaneous streamlines. Bound total energy including actual radial tails and record physical integration/truncation error.
- [ ] **CURRENT-GENERAL-GLOBAL-AXIS-FIELD** Admit independent xyz/t with the original inverse, preserve support joins and pressure, and implement axis parity and Cartesian limits. Off-axis pulse point admission is insufficient.
- [ ] **CURRENT-N-DEPENDENT-RECOVERY** After leading functional moment/stress conditions hold, implement distinct n=1 and n>=2 recovery equations on a common core domain, independent moment repair, finite-order remainder and smooth sum. Truncate streamfunction/vector potential before curl. Acceptance requires actual order-dependent coefficients; coordinate rescaling alone is insufficient.
- [ ] **CURRENT-OSCILLATORY-CORRECTION-FULL-NS** Implement both pulse families and mean corrections, check realizable averaged quadratic stress cancellation, then validate fixed-forcing corrected Cartesian residual max/physical L2 against 1e-3. Do not define forcing as the residual.

Next priority: relative remainder scale cancellation and controlled time/query cone coverage. Keep the current point-cone progress distinct from regional admissible-stress completion and actual scale recursion.
