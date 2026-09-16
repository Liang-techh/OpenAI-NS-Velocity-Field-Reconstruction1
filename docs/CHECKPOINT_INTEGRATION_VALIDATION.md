# Checkpoint integration validation

Task: NS021  
Owner: Agent 8  
Selected integration baseline: `codex/stage1-complete-slow2`  
Baseline commit: `77ed17c1e81b7abce99fa4c7ffc5b3d96ba389d0`

This report records executable validation evidence for the exact checkpoint selected by NS020. It does not accept any other task and does not promote mathematical status. `paper_exact_velocity_available=false` and `full_reconstruction=false` remain mandatory.

## Baseline CI evidence

GitHub Actions run `35045650971` (`tests`, run 2532) was associated with checkpoint head `77ed17c1e81b7abce99fa4c7ffc5b3d96ba389d0` through PR #367. Actions checked out PR merge ref `03dc7654b1ee2c9d0f7a3fc176925c25ff8dbd59`. GitHub comparison shows that merge commit is exactly one merge commit ahead of `77ed17c1...` with **no file differences**, so the executed tree content is identical to the selected checkpoint tree. This distinction is retained here rather than calling a PR merge-ref run a raw branch-head run.

The workflow configuration at the selected checkpoint runs four Python 3.11 slice jobs and two full jobs. Each full job executes:

```text
python -m pip check
python -m pytest -q -W error
```

then, only after the full test step succeeds, builds a wheel, installs it into `/tmp/ns-wheel`, runs `ns-reconstruct audit`, requires `ns-reconstruct audit --require-paper-exact` to exit 2, and runs `ns-reconstruct demo --output artifacts`.

Observed baseline results:

| job | result |
| --- | --- |
| `slice-coordinates` | success |
| `slice-velocity` | success |
| `slice-forcing` | success |
| `slice-provenance` | success |
| `full-3.10-numpy<2` | failure in full pytest |
| `full-3.13-numpy>=2,<3` | failure in full pytest |

Both full environments passed `python -m pip check` with `No broken requirements found.` Both then failed on the same two tests in `tests/test_axis_coefficient_wide_first_picard_slow2.py`:

```text
test_production_aggregate_sums_injected_nonzero_channels_and_preserves_scales
test_production_aggregate_cancellation_and_bad_amplitude_metadata_fail_closed
```

Python 3.10 / NumPy 1.26.4 reported `2 failed, 1991 passed in 3634.89s (1:00:34)`. Python 3.13 / NumPy 2.5.3 reported `2 failed, 1991 passed in 1390.96s (0:23:10)`.

The common exception was raised at `axis_coefficient_wide_first_picard_slow2.py:408`:

```text
ValueError: slow2 branch jet missing amplitude log source metadata
```

Because pytest failed, both full jobs skipped wheel build, installed-wheel CLI audit, required-paper-exact refusal, demo, and diagnostics upload. Therefore run `35045650971` supplies **no** successful wheel/demo/audit evidence and must not be reported as such.

## Failure diagnosis and minimal fix

The production `ActualScheduleWideFirstPicardSlow2State.jet(...)` now fail-closes on an `AmplitudeLogSource` for every branch jet and requires every source to match the actual schedule source exactly. That guard is intentional provenance/truth-boundary behavior and is not weakened here.

The two failing regressions replace the four already-constructed production branch `jet` callables with test-only finite polynomial rows. Their helper `_injected_branch(...)` still populated only the legacy scalar `amplitude_log` field, so the synthetic jets no longer satisfied the production metadata contract. The underlying aggregate arithmetic assertions were not reached.

Commit `2aadbd57cd92fd972887b323b3c4e0c93f591a0f` updates only that regression fixture:

- each injected jet obtains the actual schedule `AmplitudeLogSource` at the requested eta;
- the source midpoint remains the scalar `amplitude_log` supplied to the jet dataclass;
- the malformed-amplitude branch constructs a deliberately shifted source midpoint, so the existing negative test still exercises fail-closed amplitude identity rejection;
- production code and mathematical status flags are unchanged.

This is a test-contract repair, not a new mathematical construction.

## Current validation status

The candidate branch is `agent8/ns021-checkpoint-integration-validation`, based directly on `77ed17c1...`. The final NS021 acceptance evidence must come from the exact final candidate tree after this report/task ledger are committed. Until the final workflow run completes, the following remain **not yet established for the candidate head**:

- both full pytest matrices green;
- wheel build and outside-checkout installation;
- ordinary `ns-reconstruct audit` execution;
- `ns-reconstruct audit --require-paper-exact` exit code 2;
- demo generation.

No old PR CI result is substituted for those missing checks. Lean build is not part of the current GitHub workflow and is not claimed here.

## Truth boundary

Even if every repository test and CLI gate becomes green, NS021 is only integration/runtime validation of the selected checkpoint tree. It does not by itself establish the missing decay/geometry/pressure/fixed-point/Section 5-10 construction chain. In particular:

```text
paper_exact_velocity_available=false
full_reconstruction=false
```

remain unchanged.