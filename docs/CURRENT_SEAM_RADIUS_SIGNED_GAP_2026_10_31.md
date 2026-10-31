# Conditional old/current endpoint gap — 2026-10-31

Status: ACTIVE / INCOMPLETE. A new symbolic packet keeps the signed difference between an old endpoint and the current F endpoint separately at `main.upper` and `gap.lower`. It combines the previous saved F endpoint candidates and opaque C3 endpoint references with a conditional radius/field gap formula. These are **proposed bindings**, not a source-function identity or a numerical evaluation.

For a paired endpoint, suppose the radii have compatible units, \(R_o=R_f e^D\), and define independently

\[
K_o=R_o^2 H_o-R_o q_{y,o},\qquad K_f=R_f^2 H_f-R_f q_{y,f}.
\]

The candidate \(H\) is \(Q/(2L)+Q^2/4-Z^2U^2/L^2\), with each owner's own fields and \(L=1-\delta Z^2\). None of these proposed source-field or derivative bindings is proved by the packet. Define \(e=\operatorname{expm1}(D)\), \(e_2=e(2+e)\), \(G=H_o-H_f+H_o e_2\), and \(T=q_{y,o}-q_{y,f}+q_{y,o}e\). The exact conditional signed difference is

\[
K_o-K_f=R_f^2G-R_fT.
\]

Its ordinary Z1 is retained as

\[
2R_f(R_f)_Z G+R_f^2G_Z-(R_f)_Z T-R_fT_Z,
\]

where \(G_Z=(H_o)_Z-(H_f)_Z+(H_o)_Ze_2+2H_o e^{2D}D_Z\) and \(T_Z=(q_{y,o})_Z-(q_{y,f})_Z+(q_{y,o})_Ze+q_{y,o}e^D D_Z\). No derivative is defaulted to zero. `D` is the full logarithmic radius ratio, including any old/current coordinate, \(\mu\), reference-scale or radius difference; the packet does not assume the two \(y\) coordinates agree. These old-minus-current differences are per-side quantities, not the separate main-minus-gap seam orientation.

The prior saved current-F \(R,U,Q,Q_y\) and Z1 candidates and all four old C3 composite endpoint reference pairs remain adjacent to each side's expression. The ports are symbolic labels; the packet does not evaluate them or assert that old composites equal products of the proposed old fields. Both `H` owners and both `q_y` owners remain independent. Ordinary Z semantics, actual endpoint substitution, coordinate pairing, positive compatible radii and units, physical scales, the source-level relation for \(D\), and source identity remain open. Existing errors are untouched: this algebra does not certify a smaller error or justify cancellation.

The independent synthetic fixture compares the expression against direct old-minus-current evaluation using dual-number differentiation. It passed all 23 checks, including nonzero \(D,D_Z\), independent fields and derivatives, and \(D=0\) with unequal fields. The numbers are fixture-local floating-point values, not source values or certified physical bounds. The first frozen v1 attempt failed before evaluating a scenario because its fixture parser required an `arguments` list from an exact literal. Its complete failure receipt is retained; v2 corrects that parser in a separate frozen component. No ancestor helper, source provider, graph traversal, numerical quadrature, or recurrence was run. The previous receipt was raw-hash authenticated before decompression.

## Provenance

| Artifact | SHA-256 |
| --- | --- |
| v1 failure receipt, raw gzip | `73e5a8e06601fde2cc01fe5beb329c94829cdf8a8598ec9b4197303208cc3b9f` |
| v2 frozen source manifest | `a656bf98cd84a8c52c54e54e6307dcf24b9e1f3675e7efad83f9e55b3f044cda` |
| v2 exclusive entry | `5f313db1d80b6747711ba64bd47163259a1b94506459de7f2f561ebf80eae6d8` |
| v2 complete result | `bf10558427cb1de2ed9c69a685e1caecefd1c4f4b3e0d39fdeff4ce8e031a8f5` |
| v2 receipt, raw gzip | `d15fc07138161bcb788cde8226d6c232b356b309967e3c89e7ea2de27a2deebe` |

## Next executable step

Bind the separate old and F \(H,q_y\) ports to actual source functions and ordinary Z1, including capture-safe endpoint substitution and the old C3 composite-product identities. Supply a correlated source definition and derivative for the full \(D\), incorporating the existing cross-owner Rp gap without assuming common coordinates. Preserve the old/C3 and current/F namespace distinction until proved. Only a same-source endpoint and error-packet theorem can justify removing an inherited error term.

The forward 212 terms, six signed inner/lower channels, reserved own-swirl once, independent negative complete-pulse target, four geometry corrections and all inherited errors remain unchanged; `P1 != m15`. Global five-moment closure, exact heat exterior, admissible stress/flat remainder, order-dependent coefficient recursion, oscillatory correction and independent Cartesian NS validation remain open.
