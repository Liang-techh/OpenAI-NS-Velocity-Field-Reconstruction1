# Current checkpoint and version selection

Reviewed 2026-10-02 UTC against research commit [`e0bfdad05fc64f1f0fa4c662d13bd36b8b340f4c`](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/e0bfdad05fc64f1f0fa4c662d13bd36b8b340f4c). This is a pinned snapshot, not a live scheduler report.

| Layer | Identity | What to use it for |
|---|---|---|
| Current reconstruction | ST073 / `codex/st073-transition-next` | Source-bound Part I construction and staged checks |
| Earlier numerical candidates | ST061-D/P; ST063-G1R/G2R | Frozen residual/geometry comparisons, with original failures retained |
| Main MATLAB viewer | ST054-Q2/M3 | Existing checked-in visualizer and its own native-test history |
| Main Python API | ST006 | Backward-compatible `research_baseline.load_best()` and CLI |

The latest source records admit absolute leading five-moment closure and angular coefficient axial jets through order four. They do not establish a complete C4 physical field, energy/stress/temporal closure or independent NS residual acceptance. Read [the evidence and limits](RESEARCH_STATUS.md) before continuing.

## ST073 research checkout

From an existing clone, use a separate worktree so local main work is not overwritten:

```bash
git fetch origin main codex/st073-transition-next
git worktree add --detach ../ns-st073-e0bfdad e0bfdad05fc64f1f0fa4c662d13bd36b8b340f4c
cd ../ns-st073-e0bfdad
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage all --list
```

Choose an unused worktree directory. The immutable commit above is the reviewed source; a later branch head may have additional work. `--list` prints module names without executing numerical construction or validating dependencies.

The source driver orders these stages:

```text
source -> inner -> outer -> angular -> energy -> pulse
       -> postpulse -> closure -> angularjets
```

At this snapshot, `--stage all` selects 45 modules. On a complete checkout with its required scientific dependencies, full reproduction is:

```bash
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage all
```

This command regenerates source-bound receipts in the research worktree. A selected stage such as `--stage angularjets` does not automatically run all earlier stages; use it only with matching prerequisites. The driver stops on a failed stage. These commands were inspected in source, **not numerically executed for this documentation update**.

## Historical runtime and assets

The [visualization hub](../visualization/README.md) distinguishes the bundled ST054 viewer from the separately delivered ST063 comparison ZIP. The [historical catalog](research_catalog.json) records candidate/bundle hashes and missing branch-only dependencies. Do not substitute a different field when an array is unavailable.

The main ST006 loader and original failed validation remain unchanged. See [historical numerical results](LEGACY_NUMERICAL_RESULTS.md) and [experiment navigation](../experiments/README.md).

## Continuation boundary

The next dependency is complete same-source higher derivatives: corrected future energy, incoming moment/energy functions and selected amplitude/end coefficients, then full pulse/outer mixed derivatives and C4 interfaces. Energy, outer cone/stress, flat remainder, temporal recursion and full residual validation follow; see [task routing](AGENT_TASKS.md).

The September 22 pause snapshot and older supervision headings are historical. Neither their presence nor recent branch activity establishes that an external task is currently running.
