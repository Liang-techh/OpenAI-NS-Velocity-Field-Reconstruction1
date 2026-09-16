# Migration Smoke Results

Date: 2026-09-16

The local migration smoke batch was run once from the repository root with:

```text
python -m pytest -q -W error tests/test_coordinates.py tests/test_velocity.py tests/test_forcing.py
```

Actual output and status:

```text
........                                                                 [100%]
8 passed in 0.36s
```

The command exited with status `0`. This focused batch covers the coordinate, velocity, and forcing tests only. The default push and pull request workflow also checks that the package imports and that `ns-reconstruct --help` succeeds. The legacy full suite and wheel checks run only when a manually dispatched workflow sets the `run_legacy_full` boolean input to `true`.

Limits: the full suite was not run as part of this local smoke check. These results do not establish full-suite status, exact paper reproduction, or candidate validation.
