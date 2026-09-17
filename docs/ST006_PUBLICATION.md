# ST006 publication pointer — no scientific/default change

The user requested publishing the retained root-research result and organizing the repository. Publication is tracked by [PR #245](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/pull/245); read its current GitHub merge/CI status rather than assuming it from this note.

Frozen publication source: [`749999cb093906db32b859eb983cf9e8387da67f`](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/tree/749999cb093906db32b859eb983cf9e8387da67f).

The exact retained ST006 JSON is under `artifacts/research/ST006/candidate.json` in that publication tree. SHA256:
`6b4d84b48ab9dbcd2ee1a1858d3e56ef81523f5864369d7e96c6431fccf107a3`.
It includes 1946 stored parameters; three original independent evidence files and the original runtime are included. The repository-local API is `from research_baseline import load_best`, with `scripts/ns_candidate.py` for integrity, evaluation and independent validation. This API exists in the publication tree/main after merging, not automatically in this integration checkout.

On original held-out seed9172801,4096points,six times: full momentum max0.1082289305 and volumeL20.1075843288, both ABOVE0.001; the original finite-difference divergence-max gate also fails. The original reports and false flags are retained. "Retained best" means this root research lane's current baseline, not a globally best candidate or a validated NS solution.

## Integration boundary

This note does not replace Eq45/Kokuno or any live integration candidate, change `project_status.json`, alter the force/support/core/normalization/threshold contract, modify agent instructions or schedules, or merge stacked scientific PRs. Do not apply an Eq45 or core-only result to ST006, or an ST006 result to a different artifact.

Use the published candidate as a reproducible reference only. Active construction stays on `codex/cr001-constraints` under its latest instructions and task claims. Root experiments #210/#240 remain separately reviewable. Repository navigation and archived main-entry docs are part of the publication tree.

`pde_validated=false`, `paper_exact=false`, `openai_field_identified=false`, `blowup_proved=false`.
