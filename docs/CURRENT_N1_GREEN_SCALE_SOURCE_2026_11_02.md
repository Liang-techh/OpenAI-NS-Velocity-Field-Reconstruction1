# Current n1 Green scale source record — 2026-11-02

The canonical source `green_apply` is the preconditioned operator

\[
G=(I+J)G_{\mathrm{bare}}.
\]

Its six output components use source rows `(5, 6, 6, 4, 5, 6)`.  The corresponding diagonal
source scale is

\[
S=\operatorname{diag}(F_0/\epsilon,\epsilon,\epsilon,
\epsilon P_\star^2,F_0/\epsilon,\epsilon).
\]

Each output row and its incoming source column has the same scale.  Therefore,
for the common source function, radially independent scale factors, and invertible
`S`,

\[
S^{-1}GS=G.
\]

The statement `G_Z=0` is only the fixed-`x`, fixed-`h_b` kernel statement
with `Z`-independent endpoints.  It does not say that an output derivative is
zero.

Since `S_1=S_5=F_0/\epsilon`, the supplied conditional support statement
`A=-2\chi E_{51}` gives `S^{-1}AS=A`.  The support of `A` is not reverified
by this source leaf.

For ordinary `Z` differentiation with `epsilon` and `Pstar` independent of `Z`,

\[
\Gamma=S^{-1}S_Z=
\operatorname{diag}(F_{0,Z}/F_0,0,0,0,F_{0,Z}/F_0,0).
\]

For a same-source transformation `\widehat Q=S^{-1}QS`, retaining the
`B_1S_Z` term,

\[
\widehat d=S^{-1}d_{\rm raw},\qquad
\widehat d_Z=S^{-1}d_{{\rm raw},Z}-\Gamma\widehat d.
\]

This record does not establish a complete `Kchi` resolvent, a same-source
residual function, the signed raw residual or its tail terms, or any closure.
Those require the real complete resolver and same-source residual mapping.

## Provenance state

The v3 source proof was called once.  Its structural checks passed, but its
saved status was `failed`: the only three false checks were the old first-raw
parent-admission checks, which looked for the first-raw SHA in
`parent.input_hashes` even though the authoritative direct pin is in
`parent.implementation_hashes`.

The provenance reconciliation authenticated the unchanged v3 result and
parent mappings:

- v3 result raw SHA: `87bc982493ca1f3dcce0183e5eeb368c439b1ddf57ec295296895fbe5363db95`
- parent report SHA: `4c36c2dc2e709d7ca740a689af9c9e16802569f300dfe9997851b9bb72142856`
- complete receipt: `experiments/root_st073/n1_green_scale_provenance_bridge_v1_receipt.json.gz`
- receipt SHA: `680f974a0d318b2566a398b27eb4247a3205b8fd05e519749f292b9c6a9a5243`
- manifest SHA: `345dc361263ab18d726f42b17d5b08fe2a5185a582d9bc02f81fba12ac50c653`
- entry SHA: `eb2285277198d3a3e5d882581f2f385bf478cadc40ce2ab6b71c413bda02b62a`
- result SHA: `6ffcf3e8a0f07670d2f750252246f2f5cf2e5e96ad107b2def504fa24bb64cef`

The reconciliation reported 69 actual checks and 7 new synthetic mapping
checks passed.  It did not rerun AST analysis and did not alter the original
v3 failure.

The v1 fixture-key-prefix failure is retained under raw SHA
`d0ccd6e1ffd646eaa702568e9cae3d930acc440fa8b977fde4eaba48c068edf0`.
The v2 preflight omission (not entered), its Store/Load AST-target fixture
bug, raw failure-archive SHA
`719a57707fa5b478bcd4e3ff61905ff50f1fdb9b084983f388444db8be449bc4`, and
the separate traceback text remain part of the repair record.  The v3
fixture had 6/6 checks and one actual leaf entry; its structural checks passed
while the parent provenance checks were false.  Its raw failure-archive SHA
is `3b07192d99b21d8127775145f444e2d64767f0b2987a39cbbd5a32b47253393a`.

Retain the 212 / six-item / errors record, and retain `P1 != m15`.

The next required work is the complete same-source `A`/`A_Z` and genuine
`Kchi` resolvent, followed by the signed raw residual mapping needed before
calling the earlier bridge.  A finite `D` or a high-order norm must not be
substituted for the complete function.
