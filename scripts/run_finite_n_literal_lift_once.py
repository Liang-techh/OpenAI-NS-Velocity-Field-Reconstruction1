"""Run the finite-N literal-lift fixture exactly once.

The pending brief and frozen source manifest must already exist.  This runner
only authenticates them, records an exclusive entry, invokes the new direct
guard fixture once, and saves its complete result.  It never creates or edits
the checkpoint or pending brief and it does not import an old consumer.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import traceback
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
COMPONENT = "finite_n_literal_lift_v2"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
DIRECTORY = ROOT / "tmp/finite_n_literal_lift_v2"
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"
FIXTURE = ROOT / "src/openai_ns_reconstruction/finite_n_literal_lift_fixture.py"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def _write_new(path: Path, value: Any) -> str:
    _require(not path.exists(), "artifact already exists: " + str(path))
    raw = (
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")
    with path.open("xb") as stream:
        stream.write(raw)
    return hashlib.sha256(raw).hexdigest()


def _source_path(name: str) -> Path:
    path = Path(name)
    return path if path.is_absolute() else ROOT / path


def _verify_source_manifest(source: dict[str, Any]) -> dict[str, str]:
    pins = source.get("sha256")
    _require(isinstance(pins, dict) and pins, "frozen source sha256 manifest is missing")
    normalized: dict[str, str] = {}
    for name, expected in pins.items():
        _require(isinstance(name, str) and isinstance(expected, str), "invalid source pin")
        path = _source_path(name)
        _require(path.is_file(), "frozen source is missing: " + name)
        _require(_sha(path) == expected, "frozen source changed: " + name)
        normalized[name] = expected
    return normalized


def _pending_source_sha(state: dict[str, Any]) -> str | None:
    direct = state.get("source_sha256")
    if isinstance(direct, str):
        return direct
    checkpoint = state.get("source_checkpoint")
    if isinstance(checkpoint, dict) and isinstance(checkpoint.get("sha256"), str):
        return checkpoint["sha256"]
    return None


def main() -> int:
    _require(BRIEF.is_file(), "pending brief is missing")
    _require(SOURCE.is_file(), "frozen source manifest is missing")
    _require(FIXTURE.is_file(), "adapter guard fixture is missing")
    _require(not ENTRY.exists() and not RESULT.exists(), "adapter guard entry already exists")

    state = json.loads(BRIEF.read_bytes())
    _require(
        state.get("phase") == "finite_n_literal_lift_pending"
        and state.get("component") == COMPONENT,
        "pending brief is not in the finite-N literal-lift phase",
    )
    source = json.loads(SOURCE.read_bytes())
    _require(isinstance(source, dict), "frozen source manifest is not an object")
    source_pins = _verify_source_manifest(source)
    source_sha256 = _sha(SOURCE)
    expected_source_sha256 = _pending_source_sha(state)
    _require(
        expected_source_sha256 is not None
        and expected_source_sha256 == source_sha256,
        "pending brief source hash does not match frozen source manifest",
    )

    entry = {
        "component": COMPONENT,
        "entered": True,
        "source_sha256": source_sha256,
        "source_pins": source_pins,
        "fixture_sha256": _sha(FIXTURE),
        "kernel_sha256": _sha(
            ROOT / "src/openai_ns_reconstruction/original_finite_N_density_kernel.py"
        ),
        "fixture_calls_started": 0,
        "old_math_calls": 0,
        "runtime_model_metadata": None,
    }
    entry_sha256 = _write_new(ENTRY, entry)

    try:
        sys.path.insert(0, str(ROOT / "src"))
        spec = importlib.util.spec_from_file_location(
            "finite_n_literal_lift_fixture_once", FIXTURE
        )
        _require(spec is not None and spec.loader is not None, "fixture import failed")
        fixture = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = fixture
        spec.loader.exec_module(fixture)
        _require(callable(getattr(fixture, "run_fixture", None)), "fixture API missing")
        fixture_result = fixture.run_fixture()
        _require(isinstance(fixture_result, dict), "fixture result is not an object")
        outcome: dict[str, Any] = {
            "component": COMPONENT,
            "status": "returned",
            "source_sha256": source_sha256,
            "entry_sha256": entry_sha256,
            "fixture_sha256": _sha(FIXTURE),
            "kernel_sha256": _sha(
                ROOT / "src/openai_ns_reconstruction/original_finite_N_density_kernel.py"
            ),
            "fixture_calls": 1,
            "old_math_calls": 0,
            "runtime_model_metadata": None,
            "result": fixture_result,
        }
    except Exception:
        outcome = {
            "component": COMPONENT,
            "status": "failed",
            "source_sha256": source_sha256,
            "entry_sha256": entry_sha256,
            "fixture_sha256": _sha(FIXTURE),
            "kernel_sha256": _sha(
                ROOT / "src/openai_ns_reconstruction/original_finite_N_density_kernel.py"
            ),
            "fixture_calls": 1,
            "old_math_calls": 0,
            "runtime_model_metadata": None,
            "traceback": traceback.format_exc(),
        }

    result_sha256 = _write_new(RESULT, outcome)
    print(
        json.dumps(
            {
                "component": COMPONENT,
                "status": outcome["status"],
                "fixture_passed": outcome.get("result", {}).get("passed"),
                "fixture_calls": 1,
                "old_math_calls": 0,
                "result_sha256": result_sha256,
                "runtime_model_metadata": None,
            },
            sort_keys=True,
        )
    )
    return 0 if outcome["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
