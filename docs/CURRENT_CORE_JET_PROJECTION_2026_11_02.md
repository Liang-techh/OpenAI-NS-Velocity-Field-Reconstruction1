# Current Core Jet Projection Handoff — 2026-11-02

## v2 result established

The v2 exact symbolic rational algebra establishes the first-order projection

\[
\pi(a_0+a_1h)=(a_0,a_1),
\]

for all six operations: addition, negation, multiplication, squaring,
reciprocal, and division with a nonzero denominator. It also establishes the
\(d/dh\) projection of second-order coefficients.

The exponential-coefficient recurrence is

\[
r_0=1,\qquad r_1=\ell_0,\qquad
r_2=\frac{\ell_1+\ell_0^2}{2}.
\]

For a \(\Phi\) jet written as \((\Phi_0,\Phi_1,\Phi_2/2)\), the first two
ordinary derivatives of the normalized product \(\Phi r\) are

\[
D_1=\Phi_1+\ell_0\Phi_0,
\qquad
D_2=\Phi_2+2\ell_0\Phi_1+(\ell_1+\ell_0^2)\Phi_0.
\]

Here \(D_1,D_2\) denote derivatives, not the perturbation-state fields. All six operations, the recurrence, and the \(F/f_z\) calculation matched the
exact scalar calculation. Eight partial AST anchors against the frozen source
code matched. Nine synthetic checks passed, including changed
signature, denominator, factor, omission, and syntax cases. The actual new
proof used one call and exited cleanly with process exit 0.

## Validity and limits

This establishes an exact field-algebra result under the same source atoms,
proper ordinary derivatives, a nonzero denominator, and fixed
\(\rho,\delta,\Lambda,j\). The AST anchors are partial shape checks; they are
not full verification of constructors, operations, or the definition of
\(\ell\), and there is no live owner.

This does not claim interval endpoint equality, a tighter error bound, a
physical field, \(K_\chi\), or a five-moment closure.

## Source and experiment notes

The docstring at line 12 of the previous frozen file
`n1_core_low_matrix_v3.py` is misleading: \(\Phi_Z+\ell\Phi\) is
\((F_0\Phi)_Z/F_0\), not the derivative of \(\Phi/F_0\). The original source
remains unchanged for source-hash purposes.

The v1 synthetic attempt failed on a high relative loop index and a low unary
AST shape; no actual source was entered. Its failure archive was
`experiments/root_st073/core_jet_projection_v1_failure.json.gz`, SHA-256
`c12fe1f1a232f9cc88f25115694945fd614b0620d5ba724dfa078f9b2ce5d718`.

The v2 provenance is:

- source SHA-256: `cbed28bc2fc107e1785d00c8cbcf75d6c0c8d791310f91bcacf74460da76c5e0`
- entry SHA-256: `c2373dfa87a98af9a8b8f31e74699ffce55b68c0b7e3555c7f4e3da1986cbdd4`
- result SHA-256: `928ad38e8a5e6661cf8ccc42404abfacb2536a990049795ce0936aed8247ba5a`
- [receipt](../experiments/root_st073/core_jet_projection_v2_receipt.json.gz)
- receipt SHA-256: `a3d6bd3bcfbb896be9c738d0023198ad64b6a86a10a74fbebf4d103627b6e780`

No source modules or ancestor computations were invoked. The v3 operator
receipt is recorded as a raw hash only.
Both archives were byte-checked against saved source, entry, result and executable bytes. All staged inputs match the frozen manifests, and the staged whitespace check passed.

## Next actionable step

Bind \(F_0,\Phi,U_z,M,\ell\) from the same live source to \(E_{00}/K_0\) and
the signed \(Q/P\) endpoint. Require full source paths and one measure, and
preserve the 212 ledger together with all errors. Do not merge until this is
proven.
