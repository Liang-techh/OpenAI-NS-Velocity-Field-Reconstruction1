# Live continuous axial atom runtime

`SharedContinuousAxialRuntime(seeded, receipt)` owns one
`ContinuousAxialBump`, one `ContinuousAxialPulse`, the live complete matrix,
the two live complete pulse rows, their numeric `p_i` values, the bump energy
atoms, and the solved numeric `a` and `c`. Its `.component` is a compatible
`ContinuousAxialCorrection` adapter that shares the basis, pulse, coefficients,
and incoming base list by object identity.

The complete end primitive is exposed through
`runtime.full_end(coefficients=None, row=1)`. It evaluates the stored matrix
row directly. `component.terminal_balance` uses the same stored `p_i` and
matrix row, so the residual is evaluated from one atom owner even when the
caller requests a higher display precision. The pulse object's complete-row
method is cached; partial rows and value jets still delegate to the original
provider and preserve omitted-piece metadata.

The installed incoming outer profile uses this runtime for point coefficients,
complete pulse values, cumulative means and coefficient Z tangents. Fully
passed or wholly remaining end bumps use the stored matrix atoms; partial
bump integrals use the same live basis. The regenerated terminal receipt
identifies this path as `complete_atom_binding=shared_runtime`. Actual mean
and independent full-atom replay now agree at the reported precision, including
the Z derivative. Their common terminal residual remains nonzero.

The runtime receipt is JSON safe and retains the inherited energy atom
precision and its uncertified-accuracy flag. The exact continuous statement

\[
M c + b + a p = 0
\]

is reported separately from the materialized residual. The runtime never
assigns the terminal balance or its derivative to zero.

The focused CLI fixture reads
`lei_ren_part1_paper_continuous_incoming_outer.json` and the seeded candidate,
without constructing the joined field. It checks object identity, matrix-dot
end primitives against the shared bump primitive, complete-row balance
consistency, and nonzero residual retention. Quadrature enclosure, incoming
uncertainty, global mean closure, and finite-energy certification remain false.

Run from the repository root:

```powershell
python experiments/root_st073/lei_ren_part1_paper_continuous_axial_runtime.py
```

Evidence is written to
`lei_ren_part1_paper_continuous_axial_runtime.json`.
