# Compact M Z0 recipe attempt and physical pulse gap — 2026-10-28

Goal ACTIVE / INCOMPLETE. A new, source-only saved-recipe linkage was frozen and entered exactly once. It failed in its amplitude-certificate reader with `KeyError: 'rebuilt_source_leaf_ports'`; the full traceback and the frozen source, capsule, exclusive entry and result bytes are embedded in `experiments/root_st073/lei_ren_part1_paper_compliant_current_compact_M_z0_recipe_linkage_failure.json`. Result SHA-256: `2cb4f418867d5cdb44e2c3fa1cf302ffcf4d08a67bb19747487d02bfe175ac3a`. The saved `source_binding` has `rebuilt_source_leaf_A3` and `primary_source_targets_A3`, not the attempted plural keys. This failed leaf and the earlier generic M/Z0 failure must not be re-entered; no new binding from this invocation is accepted. The previously accepted 20 D/I/S/Cp source-port correspondences remain intact. Runtime model metadata is unavailable (`null`).

Independent static inspection located two already-saved base-source recipes in the previously published joint numerator report: C3 `M` endpoint 4938 is a positive quotient of m source handle 4920 by amplitude 4934; C5 endpoint 258913 uses m 258883 and amplitude 251009. These are authenticated source-role recipes at ordinary Z0, not values, Z1–Z3 endpoint identities, common physical points, same-unit assertions or cross-archive function equality. The saved scoped amplitude identity is not connected by a completed joint-vector binding, so retain `a=A_N-A_C` rather than setting it to zero. `M` and physical `m15` remain separate data.

## Physical pulse correction remains explicit

For each physical (`P`) and frame-only (`F`) source, with `y=log R`, `h=1-delta`, `d=1-Z^2`, `L=1-delta Z^2`, the source defines

```
S_X = (Q_X + h*Z*Q_X,Z/2 + Q_X,y)/L
    + Q_X*(Q_X/2 + Q_X,y)
    + U_X*(-2*Z*Q_X + d*Q_X,Z - 2*Z*Q_X,y)/L
N_X = -S_X/2 + (Q_X,yy + Q_X,y)/R_X
m15_X = -integral R_X*N_X dR_X.
```

The integral is the sum over main `xi=[0,11]`, gap `xi=[11,13-4mu]`, and end `s=[-4,0]`, with respective signed `m15` weights `R^2/mu`, `R^2/mu`, and `R^2`. Under the already stated source-specific regularity and Z-independent radii and endpoints, its ordinary Z1 is the same three integrals with `N_X` replaced by `N_X,Z = -S_X,Z/2 + (Q_X,yyZ+Q_X,yZ)/R_X`. This is an analytic consequence of the source definition, not a newly evaluated function or an independent source certificate. Without those assumptions, retain the radius and Leibniz endpoint derivatives.

In a common `y` chart *if subsequently certified*, write `rho=R_P/R_F`; then for either power `k=1,2`,

```
I_P^(k)[f_P]-I_F^(k)[f_F]
 = integral R_F^k*(rho^k*f_P-f_F) dy
 = rho^k*integral R_F^k*(f_P-f_F) dy
   + (rho^k-1)*I_F^(k)[f_F].
```

The last line assumes `rho` is constant in the chart variable; comparison also requires the same endpoints, source normalization and the correct orientation. For a general chart change the measure must include its Jacobian. This separates the geometry and field mismatches; it does not establish either is zero. Existing source joins establish telescoping *within* the respective pulse construction, but do not establish cross-source equality of `U`, `Q`, their required Z jets, or the outer boundary terms. Therefore preserve the full physical terminal difference `D_ac + (physical_pulse-frame_only_pulse)` and its ordinary Z1; do not cancel the correction or old errors.

Preserved: 212 cumulative / 213 raw / 183 reduced terms, six signed channels, 16 J, 17 K, two outer B terms, the once-only reserved contribution, the independent unresolved reference offset, all historical errors, and the original dirty worktree files. No new physical value, finite-N error reduction, global five-moment closure, exact heat exterior, stress correction or NS validation is claimed.

Next: independently bind the physical-pulse/current-frame chart fields, measures, endpoints, Pstar and ordinary Z1 to the saved source owners, with the separately accepted radius/log-radius witnesses as explicit dependencies. Use a new frozen leaf only; do not repeat either entered M source-check family, the entered amplitude reader, an ancestor solver, graph comparison or numerical quadrature. Derive the signed correction and same-source internal endpoint relation before any error elimination. Continue the higher-order and NS stages only when their stated dependencies are reliable.
