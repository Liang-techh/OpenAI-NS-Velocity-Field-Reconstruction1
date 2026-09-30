# Incremental analytic five-bump inverse

The actual common-source five-bump correction now supports incremental Picard
inversion of `A h + Q(h,h) = -d`, replacing a fixed degree-three response when
more finite accuracy is required. This is the paper's local analytic moment
inverse. It is not temporal n-dependent coefficient recursion.

## Implementation

- `lei_ren_part1_paper_five_bump_inverse.py`: generic scalar/PWJet/AxialDual
  inverse. All increments and their first-Z tangents remain separately exposed.
- `lei_ren_part1_paper_five_bump_inverse_field.py`: callable same-source field
  adapter, obtaining each Z's defects and background from the same reference.
  It caches only source evaluations actually obtained, never extrapolating a
  cached point to another Z.
- `lei_ren_part1_paper_five_bump_inverse_check.py/json`: actual Z=.3, x=1.25
  value and first-Z receipt from the internally generated source snapshot.
- `lei_ren_part1_paper_five_bump_inverse_fixture.py/json`: resolved analytic
  surrogate checked against an independent nonlinear root solver and a
  four-point axial stencil of that root solver.

Let `H0 = -A^-1 d`. Its residual is `Q(H0,H0)`. If `s_k` is the last increment,
the next residual is computed by direct polarization:

```
R(H_k) = 2 sum_{j<k} Q(s_j,s_k) + Q(s_k,s_k),  s_0 = H0
s_(k+1) = -A^-1 R(H_k).
```

This keeps the very small tail visible without subtracting large physical
moment totals. Materialized coefficient sums can still round away tiny parts;
the increment list is authoritative. A supplied contraction majorant yields
only a conditional bound `radius * lipschitz^(steps+1)`, dependent on its
declared uniform C1 source norm and finite quadrature constants. No actual
uniform bound is inferred from a sample.

## Results

The actual P9/W2 common-source run used 10 nonlinear updates at 473 working
digits. Maximum nominal centered residual decreased as follows:

| Update | Maximum absolute nominal residual |
| --- | ---: |
| 1 | 4.39225e-50 |
| 2 | 4.49310e-73 |
| 4 | 1.71076e-107 |
| 6 | 6.51375e-142 |
| 8 | 2.48012e-176 |
| 10 | 9.44313e-211 |

At update 10 the maximum nominal first-Z residual is `6.23766e-210`.
Aggregate replay agrees with the polarized residual to a maximum retained-atom
error of `1.41e-488` across values and first-Z tangents.

The field was joined at x=1.25, retaining the common P0 and P0_Z and using the
same stress evaluator. The JSON retains every increment and residual's value
and first-Z PW atoms. An aggregate replay verifies the retained-jet algebra;
it is not independent source quadrature or source-error enclosure.

On the resolved analytic surrogate, 16 updates agree with the independent
root to value error `1.86e-73` and first-Z error `1.85e-72`. The conditional
majorant is checked only for that declared surrogate family.

## Remaining construction gates

The actual fifth-row source is far flatter than even the 10-update residual.
Small absolute finite residuals therefore do not prove its relative hierarchy,
the exact infinite inverse, or the five terminal identities as functions of Z.
Uniform actual C1/C2 source bounds, quadrature and finite-ring errors, and the
infinite inverse tail remain open. The previous 38-point relaxed-cone sample
was for the degree-three field; it is not transferred as a certificate of this
new inverse field. Exact heat/energy matching, global admissible stress,
temporal recursion, and oscillatory correction are still incomplete.

Reproduce after generating the actual local reference snapshot:

```powershell
python experiments/root_st073/lei_ren_part1_paper_five_bump_inverse_fixture.py
python experiments/root_st073/lei_ren_part1_paper_five_bump_inverse_check.py
```
