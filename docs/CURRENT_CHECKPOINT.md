# Current checkpoint and candidate selection

Updated 2026-09-21. This page separates latest documented research from what the `main` checkout can currently execute. It is not the integration branch's task queue.

## Four distinct entry points

| Role | Candidate | Where to start | State |
|---|---|---|---|
| Latest geometry experiment | ST063-G2R | [ST063 record](research_snapshots/ST063.md), [comparison guide](../visualization/README.md) | Local study and complete offline bundle; full array/runtime import to `main` not performed |
| Residual-oriented controls | ST061-D and ST061-P | [ST061 record](research_snapshots/ST061.md) | D lower L2, P lower peak within the ST061 paired samples; no universal winner |
| Viewer bundled on `main` | ST054-Q2 and ST054-M3 | [MATLAB viewer](../visualization/matlab/README.md) | Existing frozen export and original native-test receipts; not the latest geometry field |
| Backward-compatible Python baseline | ST006 | [Frozen manifest](../artifacts/research/ST006/manifest.json), `research_baseline.load_best()` | Historical retained API default; not a claim that subsequent studies do not exist |

ST063-G1R remains an explicitly reported control, including its second-sample residual regressions and one off-grid radial-pressure sign miss. It is not hidden, but is not the preferred geometry example.

## What changed most recently

G2R has a longer moderate-strength axial rotation plateau and a 3.5%–7.3% fixed-window aspect increase at t=0.25/0.5/0.75. Its paired sampled peaks decrease, but paired volume L2 increases by 0.65%/1.84%. Strong-axis continuity at thresholds 0.25 and 0.35 is not established. The new comparison UI is supplied, not natively MATLAB-tested in that round.

Every candidate above remains unvalidated for the original full NS target. No geometry result, source-code upload or software test changes that status.

## Resume the right experiment

A working continuation needs the actual frozen candidate, compatible runtime, registration, source identity and evidence. A branch README or hash alone is not a recoverable field. [The catalog](research_catalog.json) binds the available complete bundle names and hashes, source commits, exact raw identities and test limits. Do not run an ST063 command from a checkout that lacks its complete bundle.

The main ST006 API and its original failure reports remain unchanged. See [results](RESEARCH_STATUS.md) for current paired comparisons and [repository guide](REPOSITORY_GUIDE.md) for runnable commands.

For multi-agent routing use [the live integration task file](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/blob/codex/cr001-constraints/docs/AGENT_TASKS.md). This organization change does not close, duplicate, unblock or reassign those tasks.
