"""Run the fixed-physical-radius gap transport leaf exactly once.

The runner is deliberately narrow.  It authenticates the pending brief,
manifest, new leaf, and this runner before consuming the one-entry slot.  Only
the new fixed-radius transport leaf is imported after the entry receipt is
written.  The checks are exact rational fixture checks; they do not inspect a
source graph, read an ancestor result, or certify a physical/source identity.

The primary fixture is ``m(Z, y) = Z*exp(-y)`` at ``Z = 0`` with
``delta = 0``, ``r = (log Rp)_Z = 1``, and ``r_Z = r_ZZ = 0``.  The expected
physical-radius rows are ``Q_R = -1``, ``Q_RZ = -2``, and ``Q_RZZ = -3``.
If the ``s*r**2*m1`` contribution is omitted, the last row is ``-1``; that
omission is recorded as an explicit negative control.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import traceback
from collections.abc import Mapping
from fractions import Fraction
from pathlib import Path
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[1]
COMPONENT = "gap_fixed_radius_transport_v2"
PHASE = COMPONENT + "_pending"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"

RUNNER = Path("scripts/run_gap_fixed_radius_transport_v2_once.py")
LEAF = Path("src/openai_ns_reconstruction/gap_fixed_radius_transport_v2.py")

F = Fraction
SHA256_LENGTH = 64


class RunnerInputError(ValueError):
    """Raised when the pending one-shot contract is not satisfied."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RunnerInputError(message)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            value,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
        + b"\n"
    )


def _json_safe(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, (set, frozenset)):
        return sorted((_json_safe(item) for item in value), key=repr)
    if isinstance(value, Path):
        return value.as_posix()
    return repr(value)


def _write_once(path: Path, value: Any) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = _json_bytes(_json_safe(value))
    with path.open("xb") as stream:
        stream.write(raw)
    return hashlib.sha256(raw).hexdigest()


def _load_json_object(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise RunnerInputError(label + " is not valid JSON") from exc
    _require(isinstance(value, dict), label + " is not a JSON object")
    return value


def _normal_path(value: Any, label: str) -> str:
    _require(isinstance(value, str) and value, label + " path is missing")
    normalized = value.replace("\\", "/")
    path = Path(normalized)
    _require(not path.is_absolute() and ".." not in path.parts, label + " must be repository-relative")
    return path.as_posix()


def _manifest_hashes(source: Mapping[str, Any]) -> dict[str, str]:
    files = source.get("files")
    _require(isinstance(files, list) and files, "source manifest files are missing")
    hashes: dict[str, str] = {}
    for index, item in enumerate(files):
        _require(isinstance(item, Mapping), f"source.files[{index}] is not an object")
        relative = _normal_path(item.get("path"), f"source.files[{index}]")
        expected = item.get("sha256")
        _require(
            isinstance(expected, str) and len(expected) == SHA256_LENGTH,
            f"source.files[{index}] SHA-256 is invalid",
        )
        path = (ROOT / relative).resolve()
        _require(path.is_file(), "manifest input is missing: " + relative)
        actual = _sha256_file(path)
        _require(actual.lower() == expected.lower(), "manifest SHA-256 mismatch: " + relative)
        hashes[relative] = actual
    return hashes


def _brief_scope(raw: Mapping[str, Any]) -> dict[str, Any]:
    """Select the component-specific pending section from the real brief."""

    if raw.get("component") == COMPONENT:
        return dict(raw)
    for key in ("pending_components", "components", "pending"):
        table = raw.get(key)
        if isinstance(table, Mapping):
            candidate = table.get(COMPONENT)
            if isinstance(candidate, Mapping):
                return dict(candidate)
    raise RunnerInputError("real pending brief has no component section: " + COMPONENT)


def _preflight() -> tuple[dict[str, Any], dict[str, str], str, str]:
    _require(not ENTRY.exists() and not RESULT.exists(), "one-entry component already entered")
    _require(BRIEF.is_file(), "real pending brief is missing")
    _require(SOURCE.is_file(), "source manifest is missing")

    brief_raw = _load_json_object(BRIEF, "real pending brief")
    brief = _brief_scope(brief_raw)
    _require(brief.get("component") == COMPONENT, "pending brief component changed")
    _require(brief.get("phase") == PHASE, "pending brief phase changed")
    _require(brief.get("goal_active") is True, "pending brief goal is not active")
    _require(brief.get("full_goal_complete") is False, "pending brief goal is already complete")
    _require(brief.get("do_not_reenter") is True, "pending brief re-entry guard is missing")

    source_path = SOURCE.resolve()
    source_sha = _sha256_file(source_path)
    checkpoint = brief.get("source_checkpoint")
    _require(isinstance(checkpoint, Mapping), "pending brief source checkpoint is missing")
    checkpoint_path_raw = checkpoint.get("path")
    _require(
        isinstance(checkpoint_path_raw, str) and checkpoint_path_raw,
        "pending brief source checkpoint path is missing",
    )
    checkpoint_path = Path(checkpoint_path_raw.replace("\\", "/"))
    if not checkpoint_path.is_absolute():
        checkpoint_path = ROOT / checkpoint_path
    _require(
        checkpoint_path.resolve() == source_path,
        "pending brief source checkpoint path does not equal source manifest",
    )
    _require(
        isinstance(checkpoint.get("sha256"), str)
        and checkpoint.get("sha256", "").lower() == source_sha,
        "pending brief source checkpoint SHA-256 does not match source manifest",
    )

    source = _load_json_object(SOURCE, "source manifest")
    _require(source.get("schema") == COMPONENT + "_checkpoint_v1", "source manifest schema changed")
    _require(source.get("component") == COMPONENT, "source manifest component changed")
    _require(source.get("one_entry_only") is True, "source manifest one-entry guard is missing")
    _require(source.get("synthetic_only") is True, "source manifest synthetic-only guard is missing")
    _require(source.get("source_identity_proved") is False, "source identity must remain false")
    _require(source.get("physical_scale_proved") is False, "physical scale proof must remain false")
    _require(source.get("error_certified") is False, "error certification must remain false")

    files = source.get("files")
    _require(isinstance(files, list), "source manifest files are missing")
    expected_paths = {LEAF.as_posix(), RUNNER.as_posix()}
    manifest_paths = [
        _normal_path(item.get("path"), f"source.files[{index}]")
        for index, item in enumerate(files)
        if isinstance(item, Mapping)
    ]
    _require(len(files) == 2, "source manifest must contain exactly two file entries")
    _require(len(manifest_paths) == len(files), "source manifest contains a non-object file entry")
    _require(len(set(manifest_paths)) == 2, "source manifest contains duplicate file entries")
    _require(set(manifest_paths) == expected_paths, "source manifest contains an extra or missing file")

    hashes = _manifest_hashes(source)
    _require(set(hashes) == expected_paths, "authenticated manifest paths changed")
    return dict(source), hashes, source_sha, _sha256_file(BRIEF)


def _import_leaf() -> Any:
    path = (ROOT / LEAF).resolve()
    _require(path.is_file(), "new fixed-radius transport leaf is missing")
    spec = importlib.util.spec_from_file_location("gap_fixed_radius_transport_v2_new_leaf", path)
    _require(spec is not None and spec.loader is not None, "new fixed-radius leaf load failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _callable(module: Any, name: str) -> Callable[..., Any]:
    function = getattr(module, name, None)
    _require(callable(function), "new fixed-radius leaf callable is missing: " + name)
    return function


def _field(value: Any, name: str) -> Any:
    if isinstance(value, Mapping):
        return value.get(name)
    return getattr(value, name)


def _rows(value: Any, name: str) -> tuple[Any, ...]:
    rows = getattr(value, "rows", None)
    if rows is None and isinstance(value, Mapping):
        rows = value.get("rows")
    _require(rows is not None, name + " has no rows")
    return tuple(rows)


def _packet(value: Any) -> Any:
    as_dict = getattr(value, "as_dict", None)
    if callable(as_dict):
        return _json_safe(as_dict())
    return _json_safe(value)


def _gate_state() -> tuple[Callable[[Any], bool], dict[str, Any]]:
    state: dict[str, Any] = {"calls": 0, "values": []}

    def gate(value: Any) -> bool:
        state["calls"] += 1
        state["values"].append(value)
        return True

    return gate, state


def _assert_false_metadata(value: Any) -> dict[str, bool]:
    names = (
        "source_identity_proved",
        "physical_scale_proved",
        "error_certified",
        "certified",
        "endpoint_terms_cancelled",
        "errors_removed",
        "common_scale_materialized",
    )
    flags = {name: _field(value, name) for name in names}
    _require(all(item is False for item in flags.values()), "certification/source flag became true")
    return flags


def _run_case(
    builder: Callable[..., Any],
    *,
    label: str,
    scaled_m: tuple[Fraction, Fraction, Fraction, Fraction],
    Z: Fraction,
    r: Fraction,
    r_Z: Fraction,
    r_ZZ: Fraction,
    expected_scaled_Q_R: tuple[Fraction, Fraction, Fraction],
    expected_Q_R: Fraction,
    expected_Q_RZ: Fraction,
    expected_Q_RZZ: Fraction,
) -> dict[str, Any]:
    gate, gate_state = _gate_state()
    result = builder(
        scaled_m,
        Z,
        F(0),
        F(0),
        F(0),
        r,
        r_Z,
        r_ZZ,
        denominator_nonzero_gate=gate,
    )
    _require(gate_state["calls"] == 1, label + " denominator gate was not called once")
    _require(
        _rows(_field(result, "fixed_y_scaled_Q_R"), label + " scaled Q_R")
        == expected_scaled_Q_R,
        label + " fixed-y Q_R rows changed",
    )
    _require(_field(result, "Q_R") == expected_Q_R, label + " Q_R changed")
    _require(_field(result, "Q_RZ") == expected_Q_RZ, label + " Q_RZ changed")
    _require(_field(result, "Q_RZZ") == expected_Q_RZZ, label + " Q_RZZ changed")
    flags = _assert_false_metadata(result)
    return {
        "label": label,
        "gate_calls": gate_state["calls"],
        "gate_values": gate_state["values"],
        "expected": {
            "fixed_y_scaled_Q_R_rows": expected_scaled_Q_R,
            "Q_R": expected_Q_R,
            "Q_RZ": expected_Q_RZ,
            "Q_RZZ": expected_Q_RZZ,
        },
        "returned_scalars": {
            "Q_R": _field(result, "Q_R"),
            "Q_RZ": _field(result, "Q_RZ"),
            "Q_RZZ": _field(result, "Q_RZZ"),
        },
        "returned": _packet(result),
        "flags_checked_false": flags,
    }


def main() -> int:
    source: dict[str, Any] | None = None
    hashes: dict[str, str] = {}
    source_sha: str | None = None
    brief_sha: str | None = None
    entry_sha: str | None = None
    steps: list[dict[str, Any]] = []
    cases: list[dict[str, Any]] = []

    source, hashes, source_sha, brief_sha = _preflight()
    entry_sha = _write_once(
        ENTRY,
        {
            "schema": COMPONENT + "_entry_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "entered": True,
            "source_sha256": source_sha,
            "brief_sha256": brief_sha,
            "input_hashes": hashes,
            "leaf_sha256": hashes.get(LEAF.as_posix()),
            "runner_sha256": hashes.get(RUNNER.as_posix()),
            "new_leaf_imported": False,
            "consumer_calls_started": 0,
            "graph_traversal": False,
            "old_module_imported": False,
            "source_identity_proved": False,
            "physical_scale_proved": False,
            "error_certified": False,
            "runtime_model_metadata": None,
        },
    )

    output: dict[str, Any]
    try:
        module = _import_leaf()
        _require(getattr(module, "SCHEMA", None) == COMPONENT, "new leaf schema changed")
        builder = _callable(module, "build_gap_fixed_radius_transport")
        steps.append({"step": "import_new_leaf", "invoked": True, "returned": True})

        base_m = (F(0), F(1), F(0), F(0))  # scaled m for m(Z,y) = Z*exp(-y)
        steps.append({"step": "primary_nonzero_r_counterexample", "invoked": True})
        primary = _run_case(
            builder,
            label="primary_nonzero_r_counterexample",
            scaled_m=base_m,
            Z=F(0),
            r=F(1),
            r_Z=F(0),
            r_ZZ=F(0),
            expected_scaled_Q_R=(F(-1), F(-1), F(0)),
            expected_Q_R=F(-1),
            expected_Q_RZ=F(-2),
            expected_Q_RZZ=F(-3),
        )
        primary["omitted_sr2m1_would_return"] = F(-1)
        _require(primary["returned_scalars"]["Q_RZZ"] == F(-3), "primary Q_RZZ is not -3")
        _require(primary["omitted_sr2m1_would_return"] != F(-3), "omission negative control collapsed")
        cases.append(primary)

        steps.append({"step": "combined_nonzero_r_Z_r_ZZ_hand_packet", "invoked": True})
        cases.append(
            _run_case(
                builder,
                label="combined_nonzero_r_Z_r_ZZ_hand_packet",
                scaled_m=(F(1), F(1), F(1), F(1)),
                Z=F(0),
                r=F(1),
                r_Z=F(2),
                r_ZZ=F(3),
                expected_scaled_Q_R=(F(-2), F(-5), F(-7)),
                expected_Q_R=F(-2),
                expected_Q_RZ=F(-7),
                expected_Q_RZZ=F(-23),
            )
        )

        steps.append({"step": "r_zero_fixed_y_reduction", "invoked": True})
        cases.append(
            _run_case(
                builder,
                label="r_zero_fixed_y_reduction",
                scaled_m=base_m,
                Z=F(0),
                r=F(0),
                r_Z=F(0),
                r_ZZ=F(0),
                expected_scaled_Q_R=(F(-1), F(0), F(0)),
                expected_Q_R=F(-1),
                expected_Q_RZ=F(0),
                expected_Q_RZZ=F(0),
            )
        )

        steps.append({"step": "nonzero_r_Z_coefficient", "invoked": True})
        cases.append(
            _run_case(
                builder,
                label="nonzero_r_Z_coefficient",
                scaled_m=base_m,
                Z=F(0),
                r=F(0),
                r_Z=F(2),
                r_ZZ=F(0),
                expected_scaled_Q_R=(F(-1), F(0), F(-4)),
                expected_Q_R=F(-1),
                expected_Q_RZ=F(0),
                expected_Q_RZZ=F(-6),
            )
        )

        steps.append({"step": "nonzero_r_ZZ_coefficient", "invoked": True})
        cases.append(
            _run_case(
                builder,
                label="nonzero_r_ZZ_coefficient",
                scaled_m=(F(1), F(0), F(0), F(0)),
                Z=F(1, 2),
                r=F(0),
                r_Z=F(0),
                r_ZZ=F(2),
                expected_scaled_Q_R=(F(-1, 2), F(-1), F(-3, 2)),
                expected_Q_R=F(-1, 2),
                expected_Q_RZ=F(-1),
                expected_Q_RZZ=F(-3, 2),
            )
        )

        output = {
            "schema": COMPONENT + "_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "returned",
            "source_sha256": source_sha,
            "brief_sha256": brief_sha,
            "entry_sha256": entry_sha,
            "input_hashes": hashes,
            "leaf_sha256": hashes.get(LEAF.as_posix()),
            "runner_sha256": hashes.get(RUNNER.as_posix()),
            "source_manifest": _json_safe(source),
            "consumer_calls": [item["step"] for item in steps],
            "consumer_call_count": len(cases),
            "graph_traversal": False,
            "old_module_imported": False,
            "source_identity_proved": False,
            "physical_scale_proved": False,
            "error_certified": False,
            "certified": False,
            "runtime_model_metadata": None,
            "negative_control": {
                "omitted_sr2m1_expected_Q_RZZ": F(-1),
                "correct_Q_RZZ": F(-3),
            },
            "cases": cases,
            "steps": steps,
        }
    except Exception as exc:  # noqa: BLE001 - complete one-shot failure receipt
        steps.append(
            {
                "step": "exception",
                "invoked": True,
                "exception_type": type(exc).__name__,
                "exception": str(exc),
                "traceback": traceback.format_exc(),
            }
        )
        output = {
            "schema": COMPONENT + "_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "failed",
            "source_sha256": source_sha,
            "brief_sha256": brief_sha,
            "entry_sha256": entry_sha,
            "input_hashes": hashes,
            "leaf_sha256": hashes.get(LEAF.as_posix()),
            "runner_sha256": hashes.get(RUNNER.as_posix()),
            "source_manifest": _json_safe(source),
            "consumer_calls": [item["step"] for item in steps],
            "consumer_call_count": len(cases),
            "graph_traversal": False,
            "old_module_imported": False,
            "source_identity_proved": False,
            "physical_scale_proved": False,
            "error_certified": False,
            "certified": False,
            "runtime_model_metadata": None,
            "negative_control": {
                "omitted_sr2m1_expected_Q_RZZ": F(-1),
                "correct_Q_RZZ": F(-3),
            },
            "cases": cases,
            "steps": steps,
        }

    result_sha = _write_once(RESULT, output)
    print(
        json.dumps(
            {
                "component": COMPONENT,
                "status": output["status"],
                "source_sha256": source_sha,
                "brief_sha256": brief_sha,
                "entry_sha256": entry_sha,
                "result_sha256": result_sha,
                "consumer_call_count": len(cases),
            },
            sort_keys=True,
        )
    )
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())


