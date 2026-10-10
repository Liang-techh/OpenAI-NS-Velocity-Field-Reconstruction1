# Current handoff update (2026-10-10): actual Rh C3 continuation

Read [the latest handoff](CURRENT_ORIGINAL_RH_C3_CONTINUATION_2026_10_10.md) first. Source [a40d1b63](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/a40d1b6355a118fa067c682dc47836c5a0a42316) installs actual Rh reference partial correction/complete five histories through Z3, first-y FTC and pressure, with original nonzero incoming, P0, phase and selected repair N retained. Exact Rh geometry/source registry, direct/full partitioned integral endpoint identities and four-cell whole-reference bounds pass. Next: outer partial-history chain and absolute heat matching; corrected high-y Rh rows remain open. Full physical time/Cartesian field, signed cone and temporal scale recursion stay ACTIVE / INCOMPLETE.

---

# Actual repaired mixed radial recovery (2026-10-10)

Full reconstruction **ACTIVE / INCOMPLETE**. Source [fe535e00](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/fe535e003f37a32051ecb05fa41f9a76f0225af1) extends the same accepted repaired band to ordinary y0..4/Z0..3 profiles and complete five histories, then recovers cylindrical velocity/pressure through y0..4/Z0..2 and all four inertial plus two shear sectors through y0..3/Z0..2. All four axial cells have directed whole-band magnitude bounds. Existing source family, selected repair N, original pressure, lower handles and exact radius are retained. Next: C2 Rh correction continuation and absolute heat matching. The global Cartesian/time field, signed stress cone and temporal coefficient recursion remain open.

## What is now installed

The same reserved band uses x in [1,2], y=log(x), native source offset t=2+log(x), and R=Rc*x. Actual signed high rows extend the existing C3 normalized profiles E/V, with their original nonzero incoming moment histories. Accepted y0/y1/y2 profile handles and y0/y1 history handles remain exact aliases.

For each original compact bump, arg=(y-center)/ell and

```
D_y^j bump = exp(-y)/J0 * sum(k=0..j,
    binomial(j,k)*(-1)^(j-k)*beta^(k)(arg)/ell^(k+1))
beta^(n)(arg) = exp(-1/(1-arg^2))*P_n(arg)/(1-arg^2)^(2*n)
```

The new ordinary beta third/fourth functions are bound to the original `flat_pulse_derivatives.py::beta_polynomials` and `beta_jets` source. They include the exp(-y), width and normalization factors. Support exterior and endpoint jets are zero; polynomial/tail bounds are magnitude evidence, not field definitions. C3 amplitude and all five actual repair-limit controls contribute to every mixed derivative. The same original power has y-j row (-alpha)^j*x^-alpha, alpha=1/2+mu.

Repeated complete-history rows use the signed own-rate FTC, including binomial ordinary-y products:

```
m_(j+1) = V_j-m_j
h_(j+1) = E_j-3h_j/2
k_(j+1) = (E*V)_j-3k_j/2
e_(j+1) = (V*V)_j-(E*E)_j/2-e_j
p_(j+1) = (E*E)_j/2
```

P0 is the same independent source function. Its positive-y rows are exact zero. The pressure is P0+p only in row zero, and p_j in every positive-y row. No incoming history or pressure constant is reset.

## Recovery orders and physical source factors

The unchanged generic recovery formulas are reimplemented in ordinary C2 graph algebra and independently compared with the original `GenericMomentRecovery.field_rows` method body. The actual C3 history shift supplies (q_Z,q_ZZ,q_ZZZ), consuming one Z derivative with no extra factorial. Explicit five-row/C3 guards reject inadequate source orders.

| Source/output | Ordinary y orders | Ordinary Z orders |
|---|---|---|
| Corrected E/V and all five complete histories | 0..4 | 0..3 |
| Q and original Utheta/Uz/Ur/Pi | 0..4 | 0..2 |
| Four signed inertial sectors, C/B and generic It/Iz | 0..3 | 0..2 |
| Four physical inertial and two physical shear sectors | 0..3 | 0..2 |

The original source-coordinate cylindrical fields are

```
Utheta=Pstar*E; Uz=Pstar*V
Ur=Pstar*sqrt(R/2)*Q
Pi=Pstar^2*(P0+p)
```

Physical inertial linear sectors multiply Pstar*sqrt(R/2), quadratic sectors additionally multiply Pstar. Physical shear sectors multiply Pstar/sqrt(2R). Their ordinary-y rows therefore include rates +1/2 and -1/2 respectively. Generic It/Iz=R*(linear+Pstar*quadratic) use rate +1. Every prefactor is differentiated once. These are physical **source** factors before the final time/similarity-to-Cartesian mapping.

Directed high-row bounds use the same selected repair integer, original ell/J0, genuine C3 amplitude/control bounds, analytic global beta caps and the existing complete history covers. Binomial y products and Z product/quotient recurrences remain separate. The inverse-square-root shear bound uses the strict **lower** actual Rc radius; using an upper radius cap in its denominator would be invalid. Recovery uses the accepted same-source positive L margin.

## Dedicated evidence

- 5 independent beta derivative identities, 6 y3/y4 bump checks, and 56 new ordinary Z0/Z1/Z2/Z3 profile rows.
- 60 independently differentiated repeated signed FTC history rows, plus exact zero positive-y P0.
- 273 independent mixed recovery and physical-prefactor ordinary y/Z rows.
- 273 exact directed polynomial-fixture comparisons, including pressure, meridional terms, quotient derivatives and physical shifts.
- 44 physical velocity/stress rows compared with the unchanged generic method body under symbolic scalar adapters.
- 124 accepted C3 profile/history alias rows, 51 accepted C2 recovered alias rows, and 150 direct C3-to-C2 history-consumption aliases.
- Whole reserved band on all 4 actual Z cells, retaining source Rc, Pstar, delta, independent P0 and original incoming histories.

The reused static reviewer is **GPT-5.6 Luna / max**, read-only. Producer admission flags stay false. Its separate checker admits only the local mixed profile/history functions, recovered mixed source functions and their magnitude ranges. These checks do not admit absolute exterior closure at the new orders, stress-cone signs, a global point oracle or temporal scale recursion.

```python
from lei_ren_part1_paper_compliant_current_original_mixed_repaired_recovery import CurrentMixedRepairedRecovery

mixed = CurrentMixedRepairedRecovery()  # pass recovery_field=live_accepted_recovery when available
velocity_rows = mixed.velocity_functions()  # Utheta, Uz, Ur, Pi: y0..4, ordinary Z0..2
stress_rows = mixed.stress_functions()      # six source sectors: y0..3, ordinary Z0..2
```

Installed quartet: `lei_ren_part1_paper_compliant_current_original_mixed_repaired_recovery.py`, `.json.gz`, `_check.py`, `_check.json`.

## Completed bounded tasks

- [x] **CURRENT-PATCHED-RADIAL-Y4-CONTRACT** Same corrected C3 profiles y0..4, original compact beta y3/y4, exact lower aliases, complete signed FTC histories y0..4 and P0 positive-y zeros.
- [x] **CURRENT-PATCHED-MIXED-RECOVERY** Actual Q/velocity y0..4 and inertial/shear y0..3 through Z2, correct physical R and half-power shifts, original generic equivalence and four-cell directed bounds.

## Next tasks in dependency order

- [x] **CURRENT-C2-RH-CORRECTION-JOIN** Completed through C3 by [a40d1b63](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/a40d1b6355a118fa067c682dc47836c5a0a42316); see [actual Rh continuation evidence](CURRENT_ORIGINAL_RH_C3_CONTINUATION_2026_10_10.md).
- [ ] **CURRENT-C2-RH-INTERFACE-MIXED-ROWS** Use the same leading patched Rh radius Rm*exp(1), reference Rref*exp(-5), and phase identity. Prove corrected value/Z/ZZ interface equality and required ordinary-y rows as functions; directed overlap alone is not equality. Preserve leading/correction separation and add correction once.
- [ ] **CURRENT-C2-EXTERIOR-SEGMENT-MEMORY** Extend the accepted source-owned quiet-power, selected-pulse, selected-postpulse and collar cumulative histories through Z2. Use `current_limit_absolute_future_integrals.py` as the dependency map. Each segment takes the previous complete incoming function; keep exact units, source radius, current phase and the original independent P0.
- [ ] **CURRENT-C2-ABSOLUTE-FUTURE-INTEGRALS** Recover the five absolute exterior terminal identities through Z2 with true final heat amplitude, native future-integral FTC and full Gamma tails. Keep final reference power distinct from the current mu power. Prove original incoming cancellations; relative repaired-band terminal zeros alone do not close the absolute exterior.
- [ ] **CURRENT-C2-PRESSURE-HEAT-ASSEMBLY** Assemble signed pressure restoration with source-owned axis datum and exact/controlled heat exterior. Prove matching pressure/derivatives, axis regularity and finite-energy radial tails before global field admission. No fitted pressure tail or reset datum.
- [ ] **CURRENT-HEAT-MIXED-ORDER-CONTRACT** Inventory the actual highest y/Z orders consumed by each velocity, Cartesian residual and stress recovery. Attach an explicit contract per chart/interface. Use the minimum necessary high-row extension; do not rebuild already accepted profiles or target bounds.
- [ ] **CURRENT-PATCHED-RECOVERY-Z3-PREREQUISITE** If a downstream stress formula needs recovered Q/inertial Z3, first extend genuine source/phase/target/limit/history/P0 data to Z4 on the same branch, with quantitative control tails and exact lower aliases. Current C3 inputs admit recovered Z2 only.
- [ ] **CURRENT-GLOBAL-MIXED-DERIVATIVES** Connect all core, transition, repair, pulse/end, flatten and collar charts with the required mixed rows. Preserve one source family, pressure registry and actual phase. Install interface function identities before advertising a global field.
- [ ] **CURRENT-PHYSICAL-TIME-MAP** Implement the original similarity variables and anisotropic space/time map, including all chain-rule factors, critical-time domain and source axis limits. Preserve radial/axial scales separately; do not interpret source-coordinate derivatives as physical time derivatives.
- [ ] **CURRENT-CARTESIAN-VELOCITY-API** Assemble the shared cylindrical dispatcher into u(x,y,z,t), v(x,y,z,t), w(x,y,z,t), with stable axis limits, domain guards and ordinary numeric outputs. First attach explicit uncertainty/tail estimates from the same source.
- [ ] **CURRENT-FREQUENCY-AND-NUMERIC-ORACLE** Select one finite N satisfying repaired-limit, physical/matching, stress and numerical constraints. Current N is repair-only. Implement source-based phase/Picard/integral evaluation with certified errors; all-N logarithmic caps are not a point field.
- [ ] **CURRENT-STRESS-CONE-REPAIR-BAND** Consume the complete pressure/history and newly installed mixed rows to assemble the local divergence-form stress. Prove signed cone margins, rather than treating magnitude bounds or the discriminant function as positivity. Keep admissible stress and flat remainder distinct.
- [ ] **CURRENT-STRESS-CONE-OTHER-REGIONS** Admit corresponding signed margins at core exit, transition, pulse/end, flatten and heat collar. Report smallest margin and controlling region, with common units and source parameters.
- [ ] **CURRENT-STRESS-REMAINDER-NORMS** Compute stress/remainder maxima, volume L2 norms, decay and scale dependence separately on admitted regions. Do not force the leading background total momentum residual below the final corrected threshold by pressure fitting.
- [ ] **CURRENT-TEMPORAL-FIRST-COEFFICIENT** Implement the original n=1 recovery equations on the shared core domain, using the completed background and independent first-order moment repair. Source spatial y-derivative rows do not themselves implement coefficient recursion.
- [ ] **CURRENT-TEMPORAL-HIGHER-RECURSION** Implement the original n-dependent n>=2 equations, each order's own moment repair, common core domain and curl-preserving truncation. Prove finite-order remainder and the final smooth sum, retaining source-order contracts.
- [ ] **CURRENT-OSCILLATORY-CANCELLATION** After accepted background stress and recursion, construct both original oscillatory families, common frequency hierarchy and averaged quadratic momentum-flux cancellation. Recover the flat forcing remainder from the same corrected field.
- [ ] **CURRENT-DYNAMICS-MEASUREMENTS** On the completed time-dependent field, measure radial contraction, axial aspect ratio, swirl/vorticity amplification, true material-line winding, interscale recurrence and finite energy. Keep instantaneous geometry and accumulated trajectory winding separate.
- [ ] **CURRENT-INDEPENDENT-CARTESIAN-RESIDUAL** Independently differentiate final Cartesian velocity/pressure and compute divergence and full forced NS residual after oscillatory correction. Only this completed stage targets Linfinity/L2 residuals below 1e-3.

Read this handoff first. Claim one bounded dependency, implement it, publish source/checker evidence and mark its task [x] only after completion. Keep accepted tasks complete. Reuse live owners and the existing read-only GPT-5.6 Luna / max reviewer; do not spawn Astra children. Full reconstruction remains ACTIVE / INCOMPLETE.
