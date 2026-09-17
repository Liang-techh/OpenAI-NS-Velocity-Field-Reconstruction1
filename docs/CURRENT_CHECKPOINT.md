## Energy-normalized bipolar child (2026-09-17)

The new energy quadrature compares orders 24, 48, 96 and 192 on the full
axisymmetric support with the cylindrical volume Jacobian. Both original
candidates fail the preregistered E(0.25)=1 condition; the bipolar seed has
E(0.25)=0.3053128140696095. A separately serialized child scales all existing
Phi and swirl coefficients by 1.8097870818686452 within the original bounds.
Its evaluated energies are 1.0000000000000004 at t=0.25 and 0.5842436855426052
at t=0.75, with all six validation times in the prescribed energy range.
Positive common scaling preserves flow signs and parity; it changes nonlinear
momentum balance, so parent residual results must not be transferred to it.

Reproduce: `PYTHONPATH=src python -m openai_ns_reconstruction.constrained_bipolar_energy`.
Report and explicit child: `artifacts/bipolar_energy/`. This remains experimental;
the public default and all PDE/source-correspondence truth flags are unchanged.
Two newly fetched branches add MATLAB export and visual-evidence governance;
they were not wholesale merged, as this continuation prioritizes field feasibility.

### Quartic dependency review

Read-only review found PR168 source headbcf49a2 diverges from live at1a98cff:
38 source-branch commits versus7 live commits. Direct candidate/test commits
1055a2f/bcf49a2 depend on absent capacity helpers. A minimal future replay can
extract d3=early_delta/-3*tau*(tau-.5)*(tau-1), g=(tau+1)*tau*(tau-.5)*(tau-1),
and a=-dot(gprime,d3prime)/dot(gprime,gprime) on tau={0,.5,1}; for early_delta=-1.4,
a=-.2831460674157303. Audit-only imports are unnecessary. No stacked PR was
merged in this root continuation. This eta-even temporal mode cannot alone
correct the observed axial parity mismatch; keep it a separate optional layer.

## Root continuation after GitHub synchronization (2026-09-17)

Fetched and fast-forwarded7dd05bb, preserving all integrated agent work. Added
portable baseline and explicit bipolar-seed bundles with Python/MATLAB grids.
The new bipolar seed fixes a directly measured central-motion parity mismatch:
existing default w has the same sign at both axial sides, whereas the source
describes opposite axial outflow. Existing Phi(0,1)=1 with other poloidal
coefficients zero gives inward swirl and opposite axial outflow at96 fresh
central probes. All swirl/support settings and bounds are unchanged.

This is an opt-in qualitative candidate, not canonical or final acceptance.
Entry and exact hash: docs/EQ45_DELIVERY_API.md. Full momentum, energy,
whole-domain geometry and source correspondence remain unverified for it.
Prior temporal quartic/compact/taper evidence applies to its original candidates,
not to this changed poloidal seed. Before further shape-only optimization,
compare central flow direction/parity through the same public interface.
Archived local compact-core failures are in docs/LOCAL_COMPACT_EXPERIMENT_ARCHIVE.md.

# Current checkpoint — support-connected Eq45 delivery

Snapshot date: **2026-09-17**. The machine-readable authority is [`project_status.json`](../project_status.json). This checkpoint describes the live constrained-integration path; older `coupled_joint`, paper-core-series, FUN/SCH, and exact-reconstruction notes are historical evidence rather than the active routing source.

At the start of this snapshot, the active integration branch was `codex/cr001-constraints@cd556a20d4cdd7d63c7626300932119dd4c3fdb4`. GitHub Actions run `35226585877` completed successfully on that head.

## Immediate objective

Deliver a directly callable, saveable/loadable, MATLAB/Python-usable time-varying 3D field

`velocity(x, y, z, t) -> [u, v, w]`

whose public observable geometry, streamlines/vorticity structure, and time evolution are made progressively closer to the public OpenAI velocity-field visualization. This is **not** a paper-exact reconstruction target and is **not** a complete Navier–Stokes blow-up proof target. Full PDE acceptance remains an independent scientific gate rather than a blocker for exporting a clearly labeled visualization candidate.

## Integrated delivery path

The active delivered family is `eq45_supported_velocity_candidate_v1`.

- supported child SHA-256: `2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d`
- parent SHA-256: `48f1845fcfd71ec95495c98bb7aac3fca4653e748a8856a5421234e9601525a7`
- public evaluator: `openai_ns_reconstruction.eq45_supported_delivery:velocity`
- grid evaluator: `openai_ns_reconstruction.eq45_supported_delivery:default_field().grid`
- candidate load path: `openai_ns_reconstruction.eq45_supported_delivery:Eq45SupportedDeliveryField.load_candidate`
- delivery capsule: `artifacts/constrained/eq45_supported_delivery_capsule.json`
- registered working box: `[-2,2]^3`
- registered time interval: `[0.25, 0.75]`
- physical support connection: `r < 2`, `|z| < 2`

The support-connected child can therefore be evaluated directly as Cartesian `[u,v,w]`, serialized/reloaded, sampled at arbitrary points, and exported on a `(time,x,y,z,component)` grid without downstream code having to reconstruct the Eq45 parent plus taper manually.

## Independent state flags

These states must remain independent:

- `velocity_export_ready = true`
- `visualization_ready = false`
- `pde_validated = false`
- `physical_support_connection_implemented = true`
- `physical_support_validated = false`
- `visual_correspondence_verified = false`
- `paper_exact = false`
- `openai_field_identified = false`
- `blowup_proved = false`

A green CI run, successful serialization, training convergence, a visually appealing render, or a target-free morphology improvement does not promote any of the scientific states above.

## Integrated visualization-facing blocker

The integrated whole-domain vorticity-envelope audit shows that the static support transform is **not** visually neutral at early time. At `t=0.25`, the parent radial q99 is about `1.3086`, while the supported child radial q99 is about `1.8581`, approximately a **+42%** increase. The refinement audit found this early widen/flatten branch to be resolution-stable rather than a one-grid artifact. Accordingly the current static supported child remains export-ready but not visualization-ready.

This blocker is about morphology. It is not evidence that the field is or is not the hidden OpenAI numerical field, and it is not a PDE acceptance result.

## Open candidate-development evidence — not yet integrated

The current open stack explores existing `Phi(1,0)` temporal freedom rather than growing the spatial basis indiscriminately. The derivative-balanced quartic candidate is the lower-collateral baseline currently worth replaying first. Existing open evidence reports strong late-time morphology preservation, while fresh-seed PDE ordering versus the cubic candidate remains unresolved/seed-sensitive.

A slope-capped compact-C2 schedule improves the early off-keyframe morphology further but carries a measurable PDE-side diagnostic cost. At the audited `t=0.3125` / `81^3` comparison, open evidence reports radial q99 of about `1.151` for compact-C2, `1.193` for quartic, and `1.254` for static supported; support-collar vorticity-squared fractions are about `1.59%`, `2.20%`, and `3.85%` respectively, while axial q99 stays unchanged. A separate pressure-free vorticity diagnostic reports the compact candidate at roughly `8.9%` worse than quartic. PR #184 then screens a single compact–quartic interpolation degree; it is capacity evidence only and does not select a canonical blend.

These open results are useful for routing, but none of #168–#184 is automatically part of the live integration branch merely because its CI is green. Many are stacked on one another.

## Agent-8 integration discipline

Do not merge a stacked PR wholesale merely because GitHub reports it mergeable. Before integration, inspect ancestry and changed files, identify the smallest owned delta, and replay/rebase that delta onto the current integration head when necessary. Verify the resulting exact head with CI and record the real run ID. Supersede or close the old stacked PR when a clean replay replaces it so another agent does not harvest the same work twice.

Optimizer convergence, candidate identity, public velocity evaluation, independent validator output, visualization diagnostics, and MATLAB/Python export must remain separable layers. In particular, `velocity_export_ready=true` is allowed while `visualization_ready=false` and `pde_validated=false`.

## Shortest next delivery chain

1. Replay the smallest callable/serializable derivative-balanced quartic candidate delta onto the live integration branch after dependency review; do not import unrelated stacked history.
2. Rebind downstream evidence to that replayed candidate identity: public evaluator/save-load/grid first, then whole-domain morphology and independent stability/PDE diagnostics.
3. Harvest visualization smoke only with correct semantics. A meridional `(u,w)` line overlay is a projected 2D streamline diagnostic, not a true 3D streamline when swirl is omitted.
4. Compare quartic, compact-C2, and any bounded blend through the same public `[u,v,w]` interface and the same whole-domain fingerprint before choosing a visualization candidate. PDE failure does not prevent a clearly labeled visualization-candidate artifact.
5. Produce the final MATLAB/Python-facing report from the frozen candidate identity: one-command load/evaluate/export, sample grid, vorticity/true-3D-streamline diagnostics, hashes, time/domain mapping, and explicit truth-boundary states.

## Historical checkpoint material

The earlier function-first paper-core-series work, `coupled_joint` optimizer lineage, FUN/VIS/SCH queues, and exact-reconstruction-era experiments remain useful provenance and failure evidence, but they no longer define the active merge blockers. Consult `docs/MIGRATION.md`, `reports/CONSTRAINED_PROGRESS.md`, `docs/SCHEDULED_AGENT_TASKS.md`, and the relevant artifacts when a current task specifically needs that history. Do not route new integration work back to those queues by default.
