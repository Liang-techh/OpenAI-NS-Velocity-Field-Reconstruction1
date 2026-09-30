# Same-profile future-pressure Z derivative receipt

`lei_ren_part1_paper_pressure_tail_Z_bounds.py` evaluates a pointwise
envelope for the uncorrected pressure tail beginning at `R_v`, using the
same `AxisPressureJets` profile and Decimal stage coordinates.  On the
finite flattening interval it uses

```text
|∂Z(Uθ²/2)| <= 4 |Z|/(1+Z²) Uθ² <= 2 Uθ².
```

The nominal interval `[y_f, y_tail]` is marked exactly Z-independent.  The
heat connection and exact heat exterior use `|H'(xi)| <= h(1+h)` and retain
separate bounds.  Values are reported in normalized `P/Pstar²` and physical
pressure units.

The angular bump correction is not merged into the baseline bound.  Its
value at `R_v` is included, while the derivative remains explicitly
unresolved because the current signed-log coefficient receipt supplies
`d1(Z), d2(Z)` but no certified `Z` derivatives.  The module exposes the
formal relation

```text
dZ I_b = d1_Z (L1 + d1 Q1) + d2_Z (L2 + d2 Q2).
```

No float finite-difference derivative is used.  Gauss order 64 versus 128
refinement is recorded as numerical quadrature evidence, not as a global
quadrature proof.

The executable receipt uses the explicit `Lambda = 10^36` source schedule,
`log Pstar = 14`, `delta = 10^-32`, and samples `Z = 0, 0.3, 0.7, 1`.
The baseline derivative bound is exactly zero at `Z = 0`; for the other
samples it is represented by arbitrary-exponent signed logs because the
source `y_v` is astronomical.  The corrected-tail C1 bound remains open.
