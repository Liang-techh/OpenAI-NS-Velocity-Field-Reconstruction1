# KOKUNO-LP-NATIVE-COORDINATES-002

## Scope

Kokuno Agent 1 minimal executable increment: implement the native similarity-coordinate map needed before assembling the public leading profiles into a 3-D velocity candidate.

This increment does **not** construct the missing radial leading profile and does not identify Kokuno's convention with the repository's separate user-Eq45 convention.

## Source and version

Formula evidence is pinned to:

- repository: `KokunoYumeto/yang-mills-interacting-workbench`
- commit: `143f6773feb424ad9ed3a8d116653200f20346b7`
- path: `navier-stokes/navier_stokes_workbench.tex`
- corrected release provenance: Zenodo record `22678406`, dated 2026-09-09, file `released_ns_reader_corrected_20260909.pdf`

The workbench's leading-field coordinate section states, with `tau=1-t`, `D=1/2-h`,

```text
z = q^D eta
tau = q(1-eta^2)
X = (x^2+y^2)/(2q)
```

and therefore

```text
q - z^2 q^(2h) = tau.
```

For fixed `tau>0,z`, the source lower endpoint is `q_star=|z|^(1/D)` (zero for `z=0`), and above it

```text
g'(q) = 1 - 2h z^2 q^(2h-1) >= 1-2h > 0.
```

Thus the source relation has exactly one `q>q_star`, with `|eta|<1`.

The completed source construction fixes the stricter parameter range `0<h<1/100`; the executable class keeps that construction bound.

## Executable mapping

`KokunoNativeSimilarityCoordinates` provides vectorized:

- `solve_q(z,t)`;
- `eta(z,t)`;
- `evaluate(x,y,z,t)` returning `q, eta, X, d=1-eta^2, L=1-2h eta^2`;
- the source exact fixed-physical-coordinate Jacobian
  `q_t, eta_t, X_t, q_z, eta_z, X_z`;
- deterministic SHA and fail-closed JSON save/load.

The root solver uses

```text
q_{n+1} = tau + z^2 q_n^(2h),
q_0 = q_star + tau.
```

On the source domain its derivative satisfies

```text
0 <= 2h z^2 q^(2h-1) < 2h < 0.02,
```

so the iteration is strongly contractive while remaining inside the unique positive source branch. The implementation independently checks the defining-equation residual before returning.

## Exact Jacobian implemented

At fixed physical coordinates, the pinned reconstruction gives

```text
q_t   = -1/L
eta_t = D eta/(q L)
X_t   = X/(q L)
q_z   = 2 eta q^(1-D)/L
eta_z = (1-eta^2)/(q^D L)
X_z   = -2 eta X/(q^D L)
```

The regression compares each expression against centered finite differences of the public coordinate evaluator.

## Convention boundary

The native Kokuno relation implemented here is

```text
q - z^2 q^(2h) = 1-t.
```

Elsewhere in the repository/user Eq45 path the relation has been written

```text
q - z^2 q^(-2h) = 1-t.
```

This increment does not silently equate them. Serialization records `cross_convention_equivalence_claimed=false`, and a regression confirms that a nontrivial native solution does not in general satisfy the negative-exponent equation.

## Direct contribution to the final velocity candidate

This closes the coordinate blocker identified by Agent 1 #213 and Agent 5 #217: Kokuno leading profiles can now be evaluated on their own public `q/X/eta` backbone with analytic coordinate derivatives. The next Agent-1 increment may assemble source radial/base profile data on this native coordinate object without forcing it through the incompatible convention.

## Truth boundary

`native_q_eta_coordinates_executable=true`, but `full_leading_profile_reconstructed=false`, `full_3d_velocity_candidate=false`, `pde_validated=false`, `paper_exact=false`, `openai_field_identified=false`, and `blowup_proved=false`.
