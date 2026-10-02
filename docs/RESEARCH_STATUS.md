# Current research status and evidence

**Snapshot:** 2026-10-02 UTC. **Research:** ST073, `codex/st073-transition-next`, pinned at [`e0bfdad0`](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/e0bfdad05fc64f1f0fa4c662d13bd36b8b340f4c).

This assessment reviews saved implementation, reports and receipts. It does not independently rerun the numerical pipeline or audit the full mathematical proof. The research implementation has not been promoted into the main ST006 compatibility API.

## What changed since the previous homepage

The old main summary stopped at `9a6c1098` with two unresolved global offsets. The later [absolute leading moment report][closure] records a source-composition identity for both offsets, completing its five leading terminal moment conditions on `Z in [-1,1]`. The [angular high-jet report][jets] then supplies axial Taylor coefficients through order four for the two actual implicit angular correction functions.

The source-bound absolute identities are recorded as:

```text
Ctheta(Z) = 0
P0(Z) + Mp(Rv,Z) + Ev0^2*Prv(Z) = 0
```

These preserve the original analytic pressure datum, forward primitives and exact positive formal scales; they are not pressure fits or selections of zero from overlapping intervals. Old conservative receipts remain historical evidence, not the latest closure state.

**This closes a leading-moment subproblem, not the complete NS field.** Full-field C4, physical energy, stress/cone, temporal recursion, oscillatory cancellation and independent Cartesian residual validation remain incomplete in the [reviewed driver][driver].

## Assessment against the eight project goals

| Goal | Recorded evidence at this snapshot | Still required |
|---|---|---|
| Incompressibility | Structural radial recovery and local primitive identities are implemented | Coherent full physical field, required derivatives and independent Cartesian divergence validation |
| Finite energy | Corrected future radial swirl energy and selected positive axial amplitude are recorded | Full physical-domain/time energy with all components and radial/axial Jacobians |
| Anisotropic shrinking core | Core, exit and similarity-profile construction are developed | Verified physical time evolution and measured radial/axial scale laws |
| Swirl/axial growth | Source-bound angular field and selected axial pulse are available in partial charts | Genuine temporal coefficient recovery and multi-time growth diagnostics |
| Stress and remainder | Local identities and bounded construction components are recorded | Whole-outer cone, global admissible stress and a controlled flat remainder |
| Functional moment cancellation | Inner functional repair plus all five absolute leading terminal identities are recorded | Higher-regularity and later-order requirements; no automatic global NS acceptance |
| Inner/transition/outer matching | Corrected pulse-to-Gamma assembly and absolute leading pressure/angular matching are recorded | Full mixed derivatives, C4 interfaces and global stress-compatible matching |
| Full NS residual <= 1e-3 | No independently validated complete ST073 Cartesian residual is established | Both full-vector maximum and physical spatial volume L2 under the registered problem |

No completion percentages are assigned: these are dependent acceptance gates, not equal-sized tasks.

## Source identity and parameter provenance

The reviewed compliant source is `epsilon_source=.001*delta`:

```text
source SHA: 5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae
inner five-moment family SHA: 3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894
```

It is distinct from legacy `epsilon=.01*delta`; those old receipts must not be relabeled as evidence for the new source. The source epsilon is also distinct from the core parameter `1/Lambda`. Meeting one epsilon bound does not establish every paper hypothesis. The 144-order view uses sensitivity/envelope inclusion and is not a newly rerun point-coefficient recurrence.

The [absolute closure report][closure] identifies the supplied Lei–Ren `2609.35406v2` source and its text hash. That is the version used by this implementation; this documentation review does not independently establish arXiv publication history.

## What the five terminal identities mean

The closure report identifies the five moments as `Mz`, `Mtheta_z`, `Mztheta`, `Mtheta` and `Mp`. Their targets and normalizations are not five interchangeable scalar zero tests.

The selected axial equations give the zero `Mz` and `Mtheta_z` histories after `Rv`. The energy equation leaves `Mztheta(Rv,Z)` equal to **half the positive remaining corrected swirl energy**, tending to its zero target at infinity. The angular condition uses its renormalized asymptotic target; the pressure condition uses the absolute analytic datum. Resetting these finite-radius histories to zero would change the constructed source.

The earlier inner repair records approximately `||d||_C1 <= 4.90e-23` and `||h||_C1 <= 3.37e-20`; the selected amplitude was enclosed near `[1.0086895652, 1.0114075202]`. These are source/repair/amplitude bounds, **not** momentum residuals. See the [source rebuild][source] and [energy/amplitude report][energy].

## Higher derivatives: the precise new boundary

The two angular correction functions now have ordinary axial Taylor coefficients through order four, meaning derivative divided by `n!`. Their actual implicit equations and Jacobian are differentiated; the selected source branch and nonzero Gamma corrections are retained. The recorded checks include 15 symbolic identities and 50 independent branch/derivative/integral fixtures. These are saved checker results, not a new run in this review.

This does not supply fourth-order regularity for every velocity component, full pulse `Ur_Z`, mixed derivatives or all region interfaces. It also does not supply a fifth-derivative remainder for treating a fourth-degree Taylor jet as a finite-cell approximation. [Source and limits][jets].

## Next dependency sequence

1. Recover C4 corrected future energy and the actual incoming moment/energy functions; differentiate the selected `ap/c1/c2` equations on the same branch and with the same units.
2. Propagate complete jets through pulse and outer charts, recover radial/axial mixed derivatives and certify every required C4 interface, including support endpoints.
3. Establish full physical energy and whole-outer cone margins; construct global admissible stress and a separately controlled flat remainder.
4. Implement genuine order-dependent temporal recovery and required oscillatory corrections, then independently evaluate the complete Cartesian residual and time-dependent geometry.

The [pinned task handoff][tasks] gives the executable subtasks. A coordinate rescaling alone is not temporal recursion; a local or radial energy integral is not the full energy certificate.

## Historical numerical results are separate

The original ST061/ST063 tables, including unfavorable controls and validation limits, are preserved byte-for-byte in [historical numerical results](LEGACY_NUMERICAL_RESULTS.md). ST063-G2R's two paired sampled maxima are about `0.02439` and `0.02250`; spatial volume L2 is about `0.03426` and `0.03421`. Both original momentum gates fail. These values must not be assigned to the new ST073 construction.

The [old asset catalog](research_catalog.json), ST054 viewer and ST006 loader retain their original identities. The [structured reconstruction snapshot](reconstruction_status.json) is a separate navigation record, not a scientific acceptance receipt.

[closure]: https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/blob/e0bfdad05fc64f1f0fa4c662d13bd36b8b340f4c/docs/ABSOLUTE_LEADING_MOMENT_CLOSURE_2026_10_01.md
[jets]: https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/blob/e0bfdad05fc64f1f0fa4c662d13bd36b8b340f4c/docs/ACTUAL_ANGULAR_C4_JETS_2026_10_01.md
[driver]: https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/blob/e0bfdad05fc64f1f0fa4c662d13bd36b8b340f4c/experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py
[source]: https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/blob/e0bfdad05fc64f1f0fa4c662d13bd36b8b340f4c/docs/COMPLIANT_INNER_OUTER_HEAT_FAMILY_2026_10_01.md
[energy]: https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/blob/e0bfdad05fc64f1f0fa4c662d13bd36b8b340f4c/docs/COMPLETE_FUTURE_ENERGY_AND_AXIAL_AMPLITUDE_2026_10_01.md
[tasks]: https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/blob/e0bfdad05fc64f1f0fa4c662d13bd36b8b340f4c/docs/AGENT_TASKS.md
