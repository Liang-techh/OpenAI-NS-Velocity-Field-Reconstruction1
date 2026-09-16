# Agent 9 external-method migration log

## CR-A9-001 — Sobol off-grid held-out sampling

- **Purpose:** reduce structured-grid blind spots in independent CR006/CR009 diagnostics without changing the optimizer, candidate, forcing model, or preregistered acceptance thresholds.
- **Classification:** **可直接迁移（public API only）**. The repository calls SciPy's public `scipy.stats.qmc.Sobol` / `random_base2` API. No SciPy source code is copied.
- **Source repository:** `scipy/scipy`
- **Source commit inspected:** `eff78058ed0d4cb8e1f0b5c5585d62d712e948d2`
- **Source file inspected:** `scipy/stats/_qmc.py` (`Sobol`, `random_base2`)
- **License:** SciPy BSD-3-Clause (`LICENSE.txt` at the source commit).
- **Migration scope:** one generic rectangular-box sampler with explicit seed, `2**m` sample sizing, SciPy-version-compatible RNG construction, and a fail-closed training/validation seed separation helper.
- **Difference from source:** no QMC engine internals or direction numbers are vendored or reimplemented. This repository only wraps the public SciPy API already present in `pyproject.toml`.
- **Regression evidence:** tests cover deterministic replay, bounds/count invariants, fail-closed invalid inputs/seed reuse, and a narrow synthetic hotspot that a coarse structured grid misses but the held-out Sobol sample detects.
- **Truth boundary:** this is a supplemental independent sampling primitive. It does not evaluate the PDE, does not define a pass threshold, does not replace refinement studies, and cannot turn sampled validation or green CI into a Navier–Stokes proof. Training loss must not be passed off as this independent validation result.
- **Integration limitation:** this increment is branched from `main` as required. The production constrained candidate currently lives in open work, so candidate-specific wiring is intentionally deferred to avoid overlapping those PRs.
