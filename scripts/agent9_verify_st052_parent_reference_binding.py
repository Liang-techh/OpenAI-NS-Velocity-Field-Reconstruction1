"""Replay the package parent-binding guard against exact ST052-M source."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

from openai_ns_reconstruction.st052_parent_reference_binding import (
    DEFAULT_PARENT_REFERENCE,
    St052ParentBindingError,
    St052ReferenceBoundParent,
    TRUTH_BOUNDARY,
)

TASK_ID = "CR-A9-067"
SOURCE_REPO = "Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1"
SOURCE_PR = 508
SOURCE_HEAD = "b3b8bfdbe1077f9ec967d158602951997d81e17d"
SOURCE_RECIPE_BLOB_SHA1 = "e30c769052379f72afeee46ca264482884cc5ac7"


def _git_blob_sha1(data: bytes) -> str:
    return hashlib.sha1(f"blob {len(data)}\0".encode() + data).hexdigest()


def run(source_root: Path, out: Path) -> dict:
    recipe_path = source_root / "experiments" / "root_st052" / "recipe.json"
    data = recipe_path.read_bytes()
    if _git_blob_sha1(data) != SOURCE_RECIPE_BLOB_SHA1:
        raise RuntimeError("exact ST052-M recipe blob identity drift")
    recipe = json.loads(data)
    reference = recipe["reference"]
    if reference["points"] != [list(row) for row in DEFAULT_PARENT_REFERENCE.points]:
        raise RuntimeError("migrated parent reference points drifted from exact source")
    if reference["times"] != list(DEFAULT_PARENT_REFERENCE.times):
        raise RuntimeError("migrated parent reference times drifted from exact source")
    if reference["velocity"] != [list(row) for row in DEFAULT_PARENT_REFERENCE.velocity]:
        raise RuntimeError("migrated parent reference velocity drifted from exact source")
    if recipe.get("pde_validated") is not False:
        raise RuntimeError("exact source truth boundary unexpectedly changed")

    exp052 = source_root / "experiments" / "root_st052"
    exp051 = source_root / "experiments" / "root_st051"
    sys.path.insert(0, str(exp052))
    sys.path.insert(0, str(exp051))
    import replay_st052

    field, raw = replay_st052.reconstruct()

    def exact_parent(points, time):
        points = np.asarray(points, dtype=float)
        velocity, _ = field.fields(raw, points, np.full(len(points), float(time)))
        return np.asarray(velocity, dtype=float)

    bound = St052ReferenceBoundParent(exact_parent)
    receipt = dict(bound.receipt)

    rejected_perturbation = False

    def perturbed_parent(points, time):
        value = exact_parent(points, time).copy()
        value[:, 1] += 2.0e-6
        return value

    try:
        St052ReferenceBoundParent(perturbed_parent)
    except St052ParentBindingError:
        rejected_perturbation = True
    if not rejected_perturbation:
        raise RuntimeError("behavioral parent binding failed to reject perturbation")

    report = {
        "task_id": TASK_ID,
        "source": {
            "repo": SOURCE_REPO,
            "source_pr": SOURCE_PR,
            "source_head": SOURCE_HEAD,
            "source_recipe_path": DEFAULT_PARENT_REFERENCE.source_recipe_path,
            "source_recipe_git_blob_sha1": SOURCE_RECIPE_BLOB_SHA1,
            "migration_class": "direct_internal_migration_of_source_native_reference",
            "license": "repository-native internal migration; no third-party code copied",
        },
        "binding": {
            "reference_spec_sha256": DEFAULT_PARENT_REFERENCE.sha256(),
            "behavioral_binding_id": bound.binding_id,
            "exact_source_parent_passed": receipt["passed"],
            "probe_count": receipt["probe_count"],
            "max_abs": receipt["max_abs"],
            "max_scaled_error": receipt["max_scaled_error"],
            "perturbed_parent_rejected": rejected_perturbation,
        },
        "direct_contribution": (
            "adds a runtime source-native behavioral guard for externally injected ST052-M "
            "parent callables, closing the metadata-string-only binding seam at the frozen probes"
        ),
        "remaining_limitations": [
            "ten source-native probes do not prove whole-domain parent function equality",
            "the complete ST052-M parent evaluator is still not serialized by this increment",
            "no immutable full-parent artifact SHA is assigned by this increment",
            "the static/temporal child is not yet a whole-candidate save/load capsule",
            "pressure and restricted forcing are unchanged and fresh CR001 PDE validation was not run",
        ],
        **TRUTH_BOUNDARY,
    }
    if not report["binding"]["exact_source_parent_passed"]:
        raise RuntimeError("exact-source parent failed migrated runtime binding")

    out.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    out.write_text(text, encoding="utf-8")
    print("exact_source_parent_passed=True")
    print(f"probe_count={receipt['probe_count']}")
    print(f"max_abs={receipt['max_abs']:.17g}")
    print(f"max_scaled_error={receipt['max_scaled_error']:.17g}")
    print("perturbed_parent_rejected=True")
    print("report_sha256=" + hashlib.sha256(text.encode()).hexdigest())
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    run(args.source_root, args.out)


if __name__ == "__main__":
    main()
