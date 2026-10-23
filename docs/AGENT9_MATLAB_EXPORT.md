# Agent 9 external method note: MATLAB velocity export

- task_id: `CR-A9-003`
- classification: **directly reusable public API**
- source repository: `scipy/scipy`
- screened commit: `eff78058ed0d4cb8e1f0b5c5585d62d712e948d2`
- source API: `scipy.io.savemat` / `scipy.io.loadmat`
- source file at screened commit: `scipy/io/matlab/_mio.py`
- license: SciPy BSD-3-Clause
- migration scope: call the already-declared SciPy public API to write/read MATLAB v5 MAT files. No SciPy writer/reader implementation is copied.
- local implementation difference: this repository adds a fail-closed velocity-grid schema with explicit MATLAB `meshgrid` axis order, candidate identity and truth-boundary fields; it does not expose arbitrary MAT dictionaries as a reconstruction artifact contract.

## Direct contribution to final velocity delivery

The final project target is a callable and serializable `velocity(x,y,z,t)->[u,v,w]` usable from Python and MATLAB. This increment samples any frozen velocity callable on a structured 3D/time grid and writes `x,y,z,times,U,V,W` in an axis order that directly matches MATLAB `[X,Y,Z]=meshgrid(x,y,z)`. The paired Python loader validates the same schema so the exported field can be round-tripped before plotting or sharing.

## Truth boundary

Successful export/load is only serialization and visualization-readiness evidence. It does not test Navier--Stokes residuals, divergence, support, forcing, independent-validation thresholds, similarity to OpenAI's public visualizations, paper exactness, or blow-up. The MAT file therefore writes `claim_scope="visualization_export_only"` and false flags for `pde_validated`, `paper_exact`, and `openai_field_identified`.
