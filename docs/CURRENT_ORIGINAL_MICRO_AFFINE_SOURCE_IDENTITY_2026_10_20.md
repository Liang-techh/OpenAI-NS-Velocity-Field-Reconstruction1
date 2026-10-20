# Original micro initial and affine source identity - 2026-10-20

Tasks 14/15 and the long-term goal remain **ACTIVE / INCOMPLETE**. This
publication records the saved source-only micro initial projection and the
saved native affine offset result. No accepted ancestor computation is
replayed. Source-only validation uses the saved component packets; the
separate worker audit records an execution deviation from the read-only brief,
including uncheckpointed helper calls.

## Saved micro projection

The passed micro result covers the background inlet at `s=0`, `rho=4`, with
finite coefficients plus infinite-tail enclosures projected onto the same
unique analytic core. It binds the six nonzero inlet moments `H`, `M`, `K`,
`A`, `B`, and `C`, retaining ordinary `Z` orders `0`, `1`, and `2`,
factorials, radial units, and the physical `Pstar` scale. It does not claim
the finite-`N` exit correction at `sc/2`.

The definitive saved micro result has `passed = True` and
`execution_seconds = 0.6020156999584287`. Its flags are
`common_core_initial_projection_proved =
True`,
`six_inlet_moment_projection_proved = True`,
and `all_tail_enclosures_retained = True`.
Full micro evolution, `F0`, `hb`, phase cutoff, density, full fields, and
global closure remain false in the saved result.

## Saved affine offset result

The corrected access component covers the native affine chart source domains,
offset derivative/J, kernel forms, and global phase/source rows according to
the saved result. Its definitive saved status is:

```text
passed = True
execution_seconds = 1.5423346000025049
cross_archive_kernel_identity_proved = True
global_offset_identity_proved = True
offset_derivative_equals_jacobian = True
```

The formal affine source is pinned by SHA256
`5c74e45263f97bd5a2c40cd5e35c19375f5f925f3f4ed17deebc9c239949a4b4`. The access runner and access
kernel are pinned by SHA256 `958bb1a6adf057d3b8bcc0d8fe916cdb062bf9b0265e53e1a7c026486d049ba9` and
`a76b4554c79cac5e40975b800ff3e20fe0cd6e5f02c22bec968a2c2179fe1702`.

The affine result is a source identity only. The publication scope fixes
`same_density_certificate = False`, `finite_N_history_identity = False`,
and `physical_Q_y_Z1_seam = False`; density, history, and `Q_y` `Z1` remain
open.

The original runner's negative affine attempt is preserved separately under
`preserved_initial_native_affine_kernel_identity`, with its exact packet under
`initial_offset_component_input`. It remains a negative, one-gate source result
with `source_formula_token_missing` at the lower-collar token `name=='hbB'`;
the corrected access result is the publication's `native_affine_kernel_identity`.

## Worker execution audit

The exact audit record `affine_worker_uncheckpointed_helper_audit_2026_10_20.json` is embedded by SHA256
and base64 in the report and execution receipt. It records
`worker_execution_deviated_from_read_only_brief = True`, helper counts
`C2descriptor = 4`, `C5operands = 1`, `source_evidence = 1`, and
`read_record = 37`, with toy-algebra call counts unknown. During those
worker checks, no `derive`, `_expr`, `_row_result`, or ancestor-operator
calls ran. The worker's source
code 893 preimage is preserved, and the audit records `formalderive` as
never-entered before the first checkpoint. The independent reviewer
deviation, including its two-mode actual source evaluator, is preserved in the
same audit and is not adopted as certification. The top-level checkpointed
source-component count is `3`; this audit is recorded separately.

## Next micro evolution gate

The implementation-ready C2/C6 source link requires four premises:

1. Bind the exact comparison `Phi/barPhi` and `U/barU` functions on both
   `[0,1]` and `[1,2]`, including the `s=1` seam. `WidthPolynomial` and
   `value_enclosure` rows remain enclosures, not function values.
2. Bind C2 `bridge.direction` to C6 `comparison_directions` with the same
   `D_over_R`, hydro/pressure/swirl drives, `F0^2` dressing, `Pstar` scale,
   and argument mapping.
3. Bind the same Z-independent `hb`, first/second `chi` controls,
   `R=Ra*exp(hb*s)`, and physical measure `dy=hb*ds` in both owners.
4. Apply positivity and exact ODE uniqueness to the same initial `Phi/U` and
   six nonzero moments, then identify the two field and six-history solutions.

The common source equations are `Phi_y=-chi*R*D_over_R*Phi/2` and
`U_y=-chi*(Phi/barPhi)*(R*drive_hydro + R*Pstar^2*drive_pressure +
R^2*F0^2*drive_swirl)`. History rates are `(2,1,2,1,2,1)` with densities
`(2*Phi,U,2*Phi*U,U^2,Phi^2,Phi^2)`. These formulas are implementation
premises, not a completed full-evolution admission.

All `212` cumulative terms, raw `213`, reduced `183`, six signed channels,
all old errors, and one reserved contribution remain retained. No forbidden
continuation archive is decoded.

`P1` differs from `m15`, physical `Q` differs from loop `q`, and the accepted
physical baseline remains commit `8f36f19e3c524b424a9286c46cc8e0b855f75307`.
The stopped attempt remains stopped at session `17273`; no ancestor is
replayed.

Parent report `lei_ren_part1_paper_compliant_current_original_Rm_physical_Qy_Z_calculus.json` has SHA256 `8903b59267bae01d5182bf2dcc84ed3c02c44eb8c7d3dd9d7b14af67ac092746`. The accepted
source stage `lei_ren_part1_paper_compliant_current_original_Rm_native_transport_pressure_sources.json` has SHA256 `7c88052348a78ab7e26728f492296fb57858cead1794b9d8cb6e57ad873a7b9b`. All nine
source, entry, result checkpoints and every preimage are preserved in the
execution receipt with exact bytes and base64.

Current report: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_micro_affine_source_identity.json`, SHA256 `3ba7c5829ebc1243ef63d04f092c17613cfda7e76913980e7e16b71c810fe37c`.
