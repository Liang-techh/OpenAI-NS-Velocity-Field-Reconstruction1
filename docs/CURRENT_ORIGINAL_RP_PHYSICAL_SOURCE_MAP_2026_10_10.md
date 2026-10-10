# Current original Rp physical Cartesian/time source map (2026-10-10)

Full reconstruction **ACTIVE / INCOMPLETE**. [Source 5f4f9d0f](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/5f4f9d0f0d7cc2d83a15c190868ff19994dc05e1) installs the direct original physical-coordinate and fixed-physical-x time map on the accepted current fifteen-chart Rp-to-heat source chain. Each chart delivers 140 signed scaled Cartesian spatial rows through total order four and four first-time rows. The exact absolute radius, original delta, source amplitudes, analytic P0, complete Gamma future and original finite correction integer remain unchanged. Arbitrary numerical [u(x,y,z,t),v(x,y,z,t),w(x,y,z,t)], axis limits, global norms, stress/remainder, genuine coefficient recursion and oscillatory cancellation remain open.

## Current callable interface

```python
owner = CurrentOriginalRpPhysicalSourceMap(accepted_fourteen_interface_owner)
view = owner.evaluate('heat_exterior', '.521', 3,
                      log_tau='-10', theta='7/10')
view['physical_velocity_pressure']['ux'].report()
view['Cartesian_spatial_rows']['uz']['x1_y1_z2'].report()
view['fixed_x_time_rows']['p'].report()
view['coordinates']['x']       # exact FunctionRef on the current graph
view['coordinates']['log_r']   # lazy log-radius; no enormous exp needed
```

Use `evaluate` to obtain actual source rows from the accepted current owner. `map_source_view` is the trusted live-observation reuse interface used by the focused run/checker; it validates typed row/schema/graph/context/unit data but is not a provider-token-authenticated import API. Do not feed receipt scalar rows or independently fabricated jets into it. The accepted path reuses the existing current runtime without replaying unchanged repair, selection or future solves.

The domain is **r>0, -1<Z<1, tau=exp(log_tau)>0**. Native coordinates and theta are explicit exact rationals. `log_tau` is an exact rational or an exact FunctionRef belonging to the same current graph. The original critical physical time is 1, so physical t=1-tau; the radius parameter named T is the reshape length and is not substituted for this critical time. Z=+/-1 and an axis point are not admitted as finite points in this restricted Rp interface.

Returned rows are finite sums of exact positive source factors times signed directed derivative enclosures. No midpoint, cap, finite small replacement frequency or practical radius is chosen. `physical_velocity_pressure` aliases the order-zero Cartesian rows, not a second field. For an actual physical-coordinate query, the inverse and point/error source contracts below still have to be implemented.

## Original map and operators

```
lambda^2 * (1-Z^2) = tau
loglambda = (log_tau - log(1-Z^2))/2
r = lambda * sqrt(2R)
z = Z * lambda^(1-delta)
x = r*cos(theta), y = r*sin(theta)
u_r = lambda^-1 * Ur
u_theta, u_z = lambda^(-1-delta) * Utheta, Uz
p = lambda^(-2-2delta) * P
u_x = cos(theta)*u_r - sin(theta)*u_theta
u_y = sin(theta)*u_r + cos(theta)*u_theta
```

The exact current delta is exp(-4 logPstar - 30). Its native directed enclosure bounds operator coefficients and never defines delta. All rows consume **full ordinary similarity derivatives** D_logR^k D_Z^n of the actual velocity and pressure. Radial scale derivatives and the one Taylor-to-ordinary factorial are already present; neither is applied again. In particular current Ur already includes sqrt(R/2), so there is no second radial sqrt(2) conversion.

For Cartesian index (i,j,b), each template term (label,a,q), a+q=i+j, appends

```
R^(-(i+j)/2) * 2^((a-q)/2) * lambda^gamma
gamma = beta_label - (i+j) + b*(delta-1)
beta_Ur = -1
beta_Utheta = beta_Uz = -1-delta
beta_P = -2-2delta
```

Moving cylindrical basis derivatives use the original Cartesian templates. Native J^k factors are not consumed or appended; the physical operators act on ordinary logR rows. Fixed-physical-x time uses

```
lambda^(beta-2) * [ -beta*G/(2L)
                   + (1-delta)*Z*G_Z/(2L) + G_logR/L ]
L = 1-delta*Z^2 > 0
```

This derivative is physical partial_t at fixed x/y/z, with partial_t=-partial_tau. It is not differentiation along a source stage coordinate. Original finite correction N=2^(2^J), J=1358356628656378313, is retained separately from transverse derivative order i+j.

Signed terms merge only when their exact source-unit, power and log-scale nodes coincide. Pressure value/axial-only rows retain Pstar^2; positive radial pressure rows retain the actual theta-squared scale, Fpulse^2 Pstar^2 in pulse charts or Ev0^2 after the pulse. Independent P0 remains separate and its positive radial derivatives are zero. A single pressure amplitude is not imposed on these distinct source rows.

The graph can accept log_tau=-logR from the same exact current radius. Then log_r = (log2-log(1-Z^2))/2 by exact cancellation before numeric arithmetic. This preserves the physically relevant correlation of very small time and very large similarity radius without materializing either quantity.

## Focused evidence

- Fifteen actual current typed chart observations: 2100 Cartesian spatial rows and 60 fixed-x time rows.
- 31560 exact source/operator term identities, 3315 signed grouped consistency diagnostics against the unchanged original linear pullback, 135 exact physical-coordinate identities, and 225 pressure/Mp/P0 consistency rows.
- Original moving-basis/cylindrical operators and fixed-basepoint/divergence identities are rebound to the same current graph. The independently accepted original finite implicit-root fixture (140 spatial and four time derivatives) is reused by immutable source hashes; it is not substituted for current paper coefficients.
- Reject axial endpoints, a changed axial point, float/interval time input, foreign graph time, native-Jacobian rows, wrong pressure radial units, and a zero/cap substituted for defining delta.
- Exact correlated radius/time cancellation passes. All current source terms retain signs; zero-containing enclosures are not dropped or replaced with absolute upper bounds.
- An accepted constructor and a fresh public `evaluate('heat_exterior',...)` call pass on the current runtime. No old global physical provider is instantiated.
- The existing sole read-only worker **GPT-5.6 Luna / max** reports no material algebra, sign, unit or scope defect. Current dependency bytes are audited against the Git index before upload.

Files: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_physical_source_map.py`, `.json.gz`, `_check.py`, `_check.json`. The compressed candidate contains the current exact expression graph, signed grouped row reports and operator definitions. `exact_operator_integer_power` nodes retain the sign of a denominator written by the original symbolic operator as +/- L^b; the focused PhysicalInterpreter checks them together with the new sin/cos nodes.

## Completed bounded tasks

- [x] **CURRENT-PHYSICAL-TIME-MAPPING** Original physical/time operators installed on all fifteen actual current Rp-to-heat charts, at signed scaled row scope.
- [x] **CURRENT-RP-CARTESIAN-SPATIAL4-ROWS** All four Cartesian velocity/pressure components, all 35 spatial indices through total order four.
- [x] **CURRENT-RP-FIXED-X-TIME1-ROWS** Four original first-time rows with correct partial_t sign and fixed angular basis.
- [x] **CURRENT-RP-PHYSICAL-SOURCE-UNIT-SPLIT** Preserve distinct Ur, swirl/axial lambda exponents, baseline/radial pressure units, P0 and exact native/ordinary separation.
- [x] **CURRENT-RP-CORRELATED-RADIUS-TIME-GRAPH** Same-graph exact log-time input and demonstrated radius/time cancellation.

## Next implementation tasks and acceptance

Priority is point/error delivery and global source coverage. Reuse accepted current owners and the pure operators; do not replay unchanged selection/repair solves merely to refresh a receipt.

- [ ] **CURRENT-RP-PHYSICAL-INVERSE-CONTRACT** Define inputs for numerical physical x/y/z/t and a log-radius/log-tau form that preserves near-critical scales. Separate representable doubles from exact/log-scaled inputs. Solve lambda^2-lambda^(2delta)z^2=tau with a monotonicity/bracketing proof, then compute Z, logR and theta with directed error. Require source-domain admission and report inverse errors; do not use an arbitrary lambda or replace the source radius with an endpoint.
- [ ] **CURRENT-RP-NATIVE-CHART-LOCATION** Locate the inverse logR in the actual fifteen-chart registry. Keep logRp, inverse-mu length and short offsets separate and compare differences before arithmetic. Recover rational/interval native coordinates with true (Lrel-4), Ts or integral-defined W Jacobians. Handle exact boundaries by the accepted same-function seam callers, with left/right unit conversions explicit.
- [ ] **CURRENT-RP-RAW-POINT-ERROR-UNITS** Give each actual source evaluation an error budget in its true amplitude units: phase inversion, signed quadrature, narrow beta supports, coefficient recovery and Gamma tail. Propagate interval width through each physical scale group. An existing directed jet is an enclosure, not a selected value or automatically a sufficiently accurate physical point.
- [ ] **CURRENT-RP-PHYSICAL-ROW-MATERIALIZER** Consume exact FunctionRefs plus an admitted source oracle to evaluate finite grouped physical rows. Combine correlated log factors before exponentiation, preserve signs/cancellation, include the operator-angle errors and expose absolute/relative directed error. Avoid midpoint selection and practical smaller N. Support the new sin/cos and signed integer-power graph nodes.
- [ ] **CURRENT-RP-PHYSICAL-POINT-API** Combine inverse, chart admission, point-error source and materializer into a query returning u/v/w/p plus error/provenance at requested (x,y,z,t). Exercise pulse, active flatten/angular, waiting/collar and heat exterior points; compare at inherited interfaces after exact common-unit conversion. Leave unsupported core/axis inputs explicit until admitted.
- [ ] **CURRENT-RP-PHYSICAL-POINT-INDEPENDENT-CHECK** Check the completed query against an independent implicit-coordinate differentiation path at bounded representative points. Cover fixed-x time sign, moving basis, pressure baseline/radial split and a correlated near-critical log-time. State current source accuracy rather than reporting only a synthetic fixture.
- [ ] **CURRENT-GLOBAL-MIXED-DERIVATIVES** Add same-current core, transition, Rh/O2/O3, compact/quiet and leading-to-Rp rows and prior seams. Map their actual full derivative schema onto this same physical operator graph. The fifteen outer charts do not supply global field coverage or norms.
- [ ] **CURRENT-PHYSICAL-AXIS-LIMIT** Recover u_r and u_theta axis parity/limits from the same actual core definition, and give Cartesian derivatives at r=0 without division by r. Check axis pressure and divergence independently; do not transfer an old core-only axis exception to Rp.
- [ ] **CURRENT-CARTESIAN-VELOCITY** Deliver the complete callable [u(x,y,z,t),v(x,y,z,t),w(x,y,z,t)] after the point/error and global/axis contracts exist. Document source-family, units, accuracy and supported physical domain.
- [ ] **CURRENT-RP-POSTPULSE-QUANTITATIVE-SEAMS** Bound mixed derivatives uniformly over full source domains. Existing arbitrary-function joins and finite signed diagnostics are complete; they do not yield global quantitative admission.
- [ ] **CURRENT-RP-PULSE-UNIFORM-BOUNDS** Bound active pulse mixed rows and narrow supports with the exact selected amplitudes and full original integrals. Give real supremum and error estimates in physical units.
- [ ] **CURRENT-RP-STRESS-DERIVATIVE-ORDER-LEDGER** State actual row orders needed for stress, divergence and flat remainder. Compute missing mixed space/time derivatives only where those equations require them; no zero padding of missing Ur orders.
- [ ] **CURRENT-STRESS-CONE-AND-REMAINDER** Construct the signed original stress and independent flat remainder; compute full-domain cone margins and max/physical-volume L2 bounds. Pressure diagnostics cannot be patched into a new tail or datum.
- [ ] **CURRENT-PHYSICAL-ENERGY** Integrate actual physical-volume kinetic energy with the correct lambda-dependent Jacobian and full radial tail. Prove finite energy over the required domain/time sector; formal point bounds are insufficient.
- [ ] **CURRENT-TEMPORAL-RECURSION-N1** Implement the real n=1 recovery equations and independent five-moment repair on the shared inner interval, then map actual coefficients physically. Coordinate scaling is not recursion.
- [ ] **CURRENT-TEMPORAL-RECURSION-HIGHER** Implement n>=2 recovery/repairs, curl-preserving truncation, finite-order error and controlled smooth summation.
- [ ] **CURRENT-OSCILLATORY-FAMILIES** Construct both original oscillatory families and mean corrections from admitted stress and coefficients. Demonstrate averaged quadratic stress cancellation and bound the independent remainder.
- [ ] **CURRENT-DYNAMICS-DIAGNOSTICS** From the admitted physical field, measure radial/axial widths, aspect ratio, growth, true material winding, coefficient recurrence and energy. Geometric scaling alone does not prove recurrence or accumulated particle winding.
- [ ] **CURRENT-FULL-CARTESIAN-RESIDUAL** Independently validate complete corrected forced Navier-Stokes max and physical-volume L2 residual after stress/recursion/oscillatory correction is implemented.

Predecessor: [complete raw postpulse function chain](CURRENT_ORIGINAL_RP_POSTPULSE_MIXED_SEAMS_2026_10_10.md). Next agent starts with the actual physical inverse and point/error source contract, using this accepted adapter and current source objects. Keep the full reconstruction ACTIVE.
