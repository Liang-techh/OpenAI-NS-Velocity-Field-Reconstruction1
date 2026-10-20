# Original macro C2/C6 evolution source projection - 2026-10-20

The long-term reconstruction goal remains **ACTIVE / INCOMPLETE**.
This publication records one source-only saved macro component after the accepted
micro source projection. The saved macro result is preserved verbatim with
`passed=True`; its unresolved gates and every actual boolean
flag remain in the report. A normal-negative result is not promoted to a pass.

The actual saved background projection flag is `True` and the actual
saved own-history projection flag is `True`.

The exact macro inlet is separated from finite enclosures. The C6 exact function
ports are `actual_micro_end` and `actual_moment_end`. The C2 core-inlet
`phi0,V0`, second-exit `ell_in,deltaV_in`, and six incoming moment memories
remain saved enclosure rows. No enclosure cap or covering
row is treated as a function value, and no cover-equals-rowvalue identity is
claimed.

The saved result includes `exact_macro_inlet_function_ports`, with authenticated
function-graph references for both velocity functions and all six incoming
histories. Its `graph_ports_valid` flag is
`True`.
The source graph keeps bridge `hbB` and switch `hbS` as separate handles;
their numerical equality or inequality is not asserted.

The exact geometry is `R=Ra*exp(2*hbB+length*m)`,
`length=log(100/Ra)-2*hbB`, with `dy=length*dm` on `m in [0,1]`.
The C2 physical log-radius history derivative is converted by
`d/dm=length*d/dt`. Its full angular prefix is `ell_in+ell_macro`,
whereas the C6 macro starts from `phi_micro_end` and integrates only
`ell_macro`. The endpoint port must preserve this factorization once.

The source stage uses ordinary Z orders `[0, 1, 2]`; these are spatial ordinary
Z derivatives and are not temporal `n` recursion. Density, finite-N history,
physical `Q_y Z1`, global five-moment closure, and the full-goal flag remain
false. The actual macro projection flags are copied from the saved result into
the report scope without filling omitted flags.

If the saved background projection is true, the next source gate propagates it
through R100, the switch, and the remaining native charts before admitting actual
finite-N density and differentiating nonreset histories with the physical measure
once. If it is false, the exact unresolved gates remain the next gate.

All `212` cumulative terms, raw `213`, reduced `183`, six signed channels, old
errors, and one reserved contribution remain retained. The physical baseline is
`8f36f19e3c524b424a9286c46cc8e0b855f75307`. The full parent report SHA256 is `8903b59267bae01d5182bf2dcc84ed3c02c44eb8c7d3dd9d7b14af67ac092746` and the
accepted micro source-stage SHA256 is `9916b2aae49a624b5b505ba0c4f7b3cfbdf92f9a7a30d7852fb247efbd085dbd`. Three source, entry, and
result checkpoints are embedded byte-for-byte in the execution receipt.

Source component calls are `1`; accepted-ancestor, mathematical/history operator,
numeric, quadrature, recurrence, and norm calls are all `0`. Runtime model
metadata is unavailable (`null`).
