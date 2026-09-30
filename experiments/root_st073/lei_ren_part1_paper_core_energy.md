# Finite regular-core energy receipt

This bounded computation uses the actual core in
`build_joined_field().inner.core`. It integrates over

```text
0 <= R <= inner.r = 4e-36 = R_a
-0.5 <= Z <= 0.5
nu = 0.01
delta = 1e-200
```

The upper radius is the regular-core exit `inner.r`; the polynomial itself
continues to `4.1e-36`, so the integration stays inside the actual core. The
radial quadrature uses normalized `R / inner.r` coordinates. Every node calls
`core.evaluate(R, Z)` and uses its returned `Ur`, `F`, and `Uz`.

The physical cylindrical Jacobian is

```text
2 pi nu^(3/2) q^((3-delta)/2) (1-delta Z^2)/(1-Z^2),
q = tau/(1-Z^2).
```

With the half kinetic-energy factor, the radial part is evaluated from
`nu q^(-1) Ur^2`, and the swirl-plus-axial part from
`nu q^(-(1+delta)) ((2 R) F^2 + Uz^2)`. Factoring out `tau` gives the two
reported scales

```text
radial:          tau^((1-delta)/2)
swirl plus axial: tau^((1-3delta)/2)
```

The receipt evaluates `log(tau) = 0, -2e200, -4e200, -6e200`. Extreme values
are retained as sign, `log_abs`, and arbitrary-exponent MP values. The direct
Jacobian quadrature agrees with the separated expressions to roughly
`e^-562` to `e^-561` in the extreme rows at 443 digit working precision. The
measured radial log gains follow the requested radial scale, while the
swirl-plus-axial gains differ by the small `delta` exponent as recorded in the
JSON.

The bounded quadrature uses Gauss-Legendre orders radial `16/32` and axial
`12/24`. The coarse and fine shape integrals are close, but this is sample
refinement and does not establish a uniform quadrature bound. At a moderate
point `(R, Z, tau) = (0.37 inner.r, 0.2, 1)`, the actual
`inner.physical_field()` callable is round-tripped and its physical integrand
agrees with the separated expression with relative error around `e^-1021`.

This receipt covers the finite regular core only. It does not include the outer
or heat tail, alter the mean tail, certify global finite energy, establish a
recursive scale transition, provide uniform-in-time estimates, or prove NS
closure.

Run command:

```text
python experiments/root_st073/lei_ren_part1_paper_core_energy.py
```

Machine-readable results are in
`lei_ren_part1_paper_core_energy.json`.
