# Original n1 signed pre-Rp pressure functions, 2026-10-13

The original `r_minus -> Rc -> 2Rc -> Rp` pressure contributions now have
genuine signed C1 function handles in the accepted n1 linear operator graph.
This advances function binding after the complete numerical pressure-C1
milestone. The full original goal remains **ACTIVE / INCOMPLETE**.

The checked source quartet is
`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions`
(`.py`, `.json.gz`, `_check.py`, `_check.json`). Its local receipt gate is
`current_original_n1_rminus_to_Rp_signed_pressure_C1_functions_installed`.

## What is installed

- All 17 active source charts retain the actual finite-N fields, their own
  nonlocal mean, and all ten signed physical reduced-S integrals. Original
  Jacobians, positive microscopic radii and the chosen integer N are retained.
- The lower finite endpoint uses actual background `V_y` and `V_yZ`, plus
  the checked zero incoming/local finite-N driver. The upper O3 endpoint
  uses `q_y=-q`; its nonlocal q history is retained.
- Compact `Rc -> 2Rc` uses the accepted C3 velocity/own-mean graph and
  actual `q_y`, with `dR=Rc dx` and `R dR=Rc^2 x dx` exactly once.
- Quiet `2Rc -> Rp` uses the actual native-frame mean derivative as kappa,
  the exact identity `m_common=kappa Z`, and the signed pressure primitives
  `K^2(1/(2Rc)-1/Rp)/4`, `-K^2 log(Rp/(2Rc))/4`.
- The three contributions are summed as real expressions. All quotient and
  product derivatives are ordinary Z derivatives of the same functions.
  Neither interval centers nor cover slopes are selected as coefficients.

The restored graph has **31,222 nodes**. Aggregate handles are
`m13=(31218,31217)` and `m15=(31221,31220)` for `(value,Z)`.
These are expression identifiers, not numerical moment values. The quiet
derivatives are exactly zero; the active and compact derivatives are retained.

## How the next agent continues

Run from the repository root:

```powershell
python scripts/run_original_n1_preRp_moment_functions.py
```

Optional focused independent fixtures:

```powershell
python scripts/run_original_n1_preRp_moment_functions.py --fixtures
```

For implementation, import `SavedCurrentOriginalN1PreRpMomentFunctions` from
the source module. Its `.functions` contains the two signed C1 pressure
functions; `.pieces` contains active/compact/quiet separately; `.operator`
is the checked physical linear repair operator on the **same graph**.
Append subsequent genuine expressions to `.g`. Preserve its accepted prefix.
The archive-only constructor verifies source receipts and exact bytes and
does not reconstruct the leading-field ancestors.

The graph importer relocates explicit expression operands only. Numbers in
source recipes and derivative orders retain their meaning. The metadata
`original_native_source_partitions` and `free_query_endpoint` retains the
original graph node space; the recorded active import mapping maps those
provenance identifiers to the imported expression graph when needed.
They must not be interpreted directly as new-graph node identifiers.

## Evidence and limits

The independent fixture directly integrates the original physical N1p
forcing over `R=3.7 -> 7.4`, with physical unit `1.5` and delta `0.2`.
It checks both signed moments, both true first-Z derivatives, and both
finite viscous endpoints symbolically. A separate relocation fixture checks
literal flat values, nullable collar operands, exact rational metadata and
foreign-graph rejection. The live-source check compares **22,811 source
nodes** after relocation and verifies all 17 charts and the actual-patch
inner Duhamel variable. Archive restoration also passed.

These checks establish signed function composition on the stated radial
scope. They do not select complete defect coefficients, evaluate the full
repaired positive-order velocity family, close the five terminal moments,
or establish temporal recursion, admissible stress, oscillatory cancellation
or the final corrected Cartesian NS residual. Existing full numerical
pressure-C1 accuracy remains a separate accepted result.

## Detailed next tasks

- [x] Claim and implement `N1-FUNCTION-PRE-RP`: combine genuine active,
  compact and quiet pressure functions in the n1 linear operator graph.
- [x] Preserve the actual common Pstar, selected N, physical radii and both
  finite pressure endpoints; apply the external m15 minus exactly once.
- [x] Supply a checked archive-only restore path for the next GitHub agent.
- [ ] `N1-FUNCTION-INNER`: bind actual leading F0, Uz0, Q0 and P0 recipes
  and their radial/Z derivatives. Translate the existing `known_forcing`,
  regular matrices and regular Green operator into the common graph.
- [ ] Define the positive-order inner W as the actual unique regular
  Volterra/Neumann solution, with its source, axis normalization and domain.
  Retain the accepted source error when enclosing it. A finite prefix center
  alone is not the solution function.
- [ ] Bind the true first-Z solution derivative using the same defining
  equation and regular data. Do not differentiate a remainder cap.
- [ ] Form the existing smooth radial cutoff on F1/Uz1, then recover its
  own M/V/P from the original equations. Preserve the common P0 datum.
- [ ] Integrate its five signed cutoff densities on the full inner interval,
  with the exact radial measure and first-Z integrands. Keep the m15 pressure
  boundary subtraction and all tiny nonzero sectors.
- [ ] `N1-FUNCTION-LOWER`: rebuild the actual `Rin -> r_minus` leading N1p
  forcing as expressions and integrate `N1p dR`, `-R N1p dR` with the
  original physical radius/Jacobian. Use the same selected source family.
- [ ] `N1-FUNCTION-PULSE`: translate the actual refined original main/gap/end
  pulse coefficients and fields, including their own mean and finite-mu
  saddle dependence. Build both pressure integrals and their Z derivatives
  from these defining bodies, preserving the signed weighted endpoint.
- [ ] Bind the checked exact meridional support theorem beyond Rv. Do not
  zero angular/heat fields or append a second isolated end-pulse contribution.
- [ ] `N1-FULL-FUNCTIONS`: combine inner, lower collar, this pre-Rp adapter,
  and the complete pulse exactly once. The result must have all five genuine
  same-graph C1 handles m11/m12/m13/m14/m15. Compare enclosure queries with
  the accepted complete numerical value/C1 archives on full axial cells.
- [ ] Call `.operator.coefficient_functions(full_moments)` with those genuine
  handles. Persist actual alpha/beta function definitions and their Z rows.
- [ ] `N1-OWN-RECOVERY`: install the five fresh bump corrections, integrate
  repaired own M/P, recover V from the n1 equation, and verify the five
  terminal functional identities and axis/support/joins.
- [ ] Only after n1 functional closure, proceed to the original n>=2
  recovery/repair and temporal recursion. Keep global heat/stress, flat
  remainder, oscillatory cancellation and corrected NS validation open.

Ownership for this checkpoint: root implemented/integrated the new module;
the existing GPT-5.6 Luna/max worker ran it against retained checked owners,
the focused checker and the archive restore. No checked ancestors were edited
or reconstructed, and no additional subagent was spawned.
