"""Enter the frozen source-factorized P-minus-F leaf exactly once.

The root agent must create and freeze ``tmp/pulse_gap_source_factorized_delta_v2/
source.json`` and the matching pending brief before invoking this file.  This
runner does not create a checkpoint, does not execute a preflight expression,
and never re-enters an existing ``entered.json`` or ``result.json``.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path


R = Path(__file__).resolve().parents[1]
B = R.parent.parent / "RESUME_BRIEF_2026_10_16.json"
L = R / "tmp/pulse_gap_source_factorized_delta_v2"
SOURCE = L / "source.json"
ENTRY = L / "entered.json"
RESULT = L / "result.json"
CODE = R / "src/openai_ns_reconstruction/pulse_gap_source_factorized_delta_v2.py"
COMPONENT = "pulse_gap_source_factorized_delta_v2"
PHASE = "pulse_gap_source_factorized_delta_v2_pending"

EXPECTED_INPUTS = {
    "tmp/pulse_gap_actual_source_program_v1/result.json": "e611cdc7e049ce40e6667275b4e9706c3478044880b92c93baebf8ab8d06ebf2",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_native_physical_pulse_source.json": "8b21f8e868025d8cc95149ffab80b54603833f5f299ec9d100602af355cea82b",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz": "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_pulse_kernel_functions.json.gz": "9584201b517fe3058ea6d1a169ddb7d319776eb97479b36024af4bffd4973b75",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_pulse_kernel_functions_check.json": "47343a839590c3cb08c8d726ef225415e15f4fb12a5f7cf888049ea1f20cffc0",
}


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def save_new(path: Path, value: object) -> str:
    raw = (
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")
    with path.open("xb") as stream:
        stream.write(raw)
    return hashlib.sha256(raw).hexdigest()


def _path_from_entry(raw_path: object) -> Path:
    if not isinstance(raw_path, str) or not raw_path:
        raise RuntimeError("source file entry has no path")
    candidate = Path(raw_path)
    return candidate if candidate.is_absolute() else R / candidate


def _checkpoint_sha(brief: dict[str, object]) -> str:
    for key in ("source_checkpoint", f"{COMPONENT}_source_checkpoint"):
        value = brief.get(key)
        if isinstance(value, dict) and isinstance(value.get("sha256"), str):
            if not isinstance(value.get("path"), str) or Path(value["path"]).resolve() != SOURCE.resolve():
                raise RuntimeError("pending brief source_checkpoint path mismatch")
            return value["sha256"]
    raise RuntimeError("pending brief source_checkpoint missing")


def _verify_source_checkpoint(frozen: dict[str, object]) -> str:
    if frozen.get("schema") != "pulse_gap_source_factorized_delta_checkpoint_v2":
        raise RuntimeError("source checkpoint schema mismatch")
    if frozen.get("component") != COMPONENT:
        raise RuntimeError("source checkpoint component mismatch")
    if frozen.get("one_entry_only") is not True:
        raise RuntimeError("source checkpoint is not one-entry-only")
    entries = frozen.get("files")
    if not isinstance(entries, list) or not entries:
        raise RuntimeError("source checkpoint files missing")
    seen: set[str] = set()
    listed: dict[str, str] = {}
    for item in entries:
        if not isinstance(item, dict):
            raise RuntimeError("source checkpoint file entry is not an object")
        raw_path = item.get("path")
        expected = item.get("sha256")
        path = _path_from_entry(raw_path)
        canonical = path.as_posix().lower()
        if canonical in seen:
            raise RuntimeError("duplicate source checkpoint path: " + str(raw_path))
        seen.add(canonical)
        if not isinstance(expected, str) or len(expected) != 64:
            raise RuntimeError("source checkpoint file SHA missing: " + str(raw_path))
        if not path.is_file() or sha(path) != expected:
            raise RuntimeError("frozen source file changed: " + str(raw_path))
        listed[str(raw_path).replace("\\", "/")] = expected

    # The checkpoint must pin both new code files and every authenticated
    # source artifact.  The code/runner bytes are checked from ``files``;
    # these five entries additionally pin the exact raw hashes used by the
    # leaf itself.
    required_paths = {
        CODE.relative_to(R).as_posix().lower(),
        Path(__file__).resolve().relative_to(R).as_posix().lower(),
        *(path.lower() for path in EXPECTED_INPUTS),
    }
    normalized_listed = {key.lower(): value for key, value in listed.items()}
    for required in required_paths:
        if not any(key.endswith(required) for key in normalized_listed):
            raise RuntimeError("source checkpoint omits required file: " + required)
    for required, expected in EXPECTED_INPUTS.items():
        matches = [
            value
            for key, value in normalized_listed.items()
            if key.endswith(required.lower())
        ]
        if matches != [expected]:
            raise RuntimeError("source checkpoint raw SHA pin mismatch: " + required)
    return sha(SOURCE)


def main() -> int:
    if ENTRY.exists() or RESULT.exists():
        raise RuntimeError("already entered; never reenter")
    if not SOURCE.is_file():
        raise RuntimeError("frozen source checkpoint is missing")
    if not B.is_file():
        raise RuntimeError("pending brief is missing")

    frozen = json.loads(SOURCE.read_bytes())
    brief = json.loads(B.read_bytes())
    if not isinstance(frozen, dict) or not isinstance(brief, dict):
        raise RuntimeError("checkpoint or pending brief is not an object")
    if (
        brief.get("component") != COMPONENT
        or brief.get("phase") != PHASE
        or brief.get("goal_active") is not True
        or brief.get("full_goal_complete") is not False
        or brief.get("do_not_reenter") is not True
    ):
        raise RuntimeError("pending brief mismatch")

    source_sha = _verify_source_checkpoint(frozen)
    if _checkpoint_sha(brief) != source_sha:
        raise RuntimeError("pending brief source_checkpoint SHA mismatch")

    # Persist the exclusivity marker before importing the leaf.  The target
    # directory is part of the frozen checkpoint and is not created here.
    entry_sha = save_new(
        ENTRY,
        {
            "component": COMPONENT,
            "phase": PHASE,
            "source_sha256": source_sha,
            "entered": True,
            "consumer_calls_started": 0,
            "ancestor_math_calls": 0,
            "runtime_model_metadata": None,
        },
    )

    try:
        spec = importlib.util.spec_from_file_location(
            "_pulse_gap_source_factorized_delta", CODE
        )
        if spec is None or spec.loader is None:
            raise RuntimeError("source program load failed")
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        builder = getattr(module, "build_pulse_gap_source_factorized_delta")
        value = builder(R)
        output = {
            "component": COMPONENT,
            "phase": PHASE,
            "status": "returned",
            "source_sha256": source_sha,
            "entry_sha256": entry_sha,
            "consumer_calls": 1,
            "ancestor_math_calls": 0,
            "runtime_model_metadata": None,
            "result": value,
        }
    except Exception as exc:  # preserve a complete failure result
        output = {
            "component": COMPONENT,
            "phase": PHASE,
            "status": "failed",
            "source_sha256": source_sha,
            "entry_sha256": entry_sha,
            "consumer_calls": 1,
            "ancestor_math_calls": 0,
            "runtime_model_metadata": None,
            "exception": type(exc).__name__ + ": " + str(exc),
            "traceback": traceback.format_exc(),
        }

    result_sha = save_new(RESULT, output)
    print(
        json.dumps(
            {
                "component": COMPONENT,
                "phase": PHASE,
                "status": output["status"],
                "result_sha256": result_sha,
                "consumer_calls": 1,
                "ancestor_math_calls": 0,
                "runtime_model_metadata": None,
            },
            sort_keys=True,
        )
    )
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
