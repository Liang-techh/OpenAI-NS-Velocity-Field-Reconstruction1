"""Authenticated one-shot runner for the complete Gamma heat-delta leaf.

The root agent supplies ``tmp/gamma_weighted_energy_delta_v1/source.json`` and
the matching pending brief.  This runner authenticates every declared raw
input before the new leaf is imported or any gzip is decoded.  It never
imports a historical source, calls an older Gamma kernel/quadrature/provider,
walks a source graph, or computes a source value.  After the exclusive entry
receipt is written it invokes the new leaf exactly once and writes a complete
success/failure result before printing.

Expected checkpoint shape (additional direct metadata is allowed)::

    {
      "schema": "gamma_weighted_energy_delta_checkpoint_v1",
      "component": "gamma_weighted_energy_delta_v1",
      "one_entry_only": true,
      "no_source_graph_traversal": true,
      "all_raw_inputs_authenticated_before_decompression": true,
      "goal_flags": {"goal_active": true, "full_goal_complete": false,
                      "do_not_reenter": true},
      "files": [{"path": "...", "sha256": "..."}],
      "source_pointers": {
        "gamma": {
          "path": "...Gamma_future_functions.json.gz",
          "encoding": "gzip_json",
          "pointers": {
            "source_S": "/function_bindings/schedule/S",
            "Dhat": "/function_bindings/functions/Dhat",
            "collar_energy": "/function_bindings/pieces/collar/energy",
            "outer_energy": "/function_bindings/pieces/exterior/energy",
            "actual_scaled_deficit": "/function_bindings/pieces/actual_scaled_deficit",
            "preheat_bracket": "/function_bindings/pieces/preheat_bracket",
            "a": "/function_bindings/pieces/a"
          },
          "receipt": {
            "path": "...Gamma_future_functions_check.json",
            "pointers": {"archive_sha256": "/input_hashes/...json.gz"}
          }
        },
        "complete_future": {
          "path": "...complete_future_functions.json.gz",
          "encoding": "gzip_json",
          "pointers": {
            "energy_complete": "/function_bindings/energy/complete",
            "energy_weighted": "/function_bindings/energy/weighted",
            "Nrel": "/function_bindings/pieces/Nrel",
            "tail_multiplier": "/function_bindings/pieces/tail_multiplier",
            "Gamma_heat_difference": "/function_bindings/pieces/Gamma_heat_difference"
          },
          "receipt": {
            "path": "...complete_future_functions_check.json",
            "pointers": {"archive_sha256": "/input_hashes/...json.gz"}
          }
        }
      }
    }

Pointers are direct RFC 6901 lookups into the authenticated saved JSON.  A
missing saved role remains a false binding gate in the result; it is never
replaced by a nearby node or interpreted as a numerical source value.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import re
import sys
import traceback
from pathlib import Path
from typing import Any, Mapping


ROOT = Path(__file__).resolve().parents[1]
DIRECTORY = ROOT / "tmp" / "gamma_weighted_energy_delta_v1"
SOURCE = DIRECTORY / "source.json"
ENTERED = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
LEAF = ROOT / "src" / "openai_ns_reconstruction" / "gamma_weighted_energy_delta.py"
COMPONENT = "gamma_weighted_energy_delta_v1"
PHASE = "gamma_weighted_energy_delta_pending"
SCHEMA = "gamma_weighted_energy_delta_checkpoint_v1"
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")

GAMMA_STEM = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions"
)
COMPLETE_STEM = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_complete_future_functions"
)
PREVIOUS_RESULT = "tmp/gamma_correlated_scale_delta_v1/result.json"
REQUIRED_INPUTS = {
    GAMMA_STEM + ".py",
    GAMMA_STEM + "_check.json",
    GAMMA_STEM + ".json.gz",
    COMPLETE_STEM + ".py",
    COMPLETE_STEM + "_check.json",
    COMPLETE_STEM + ".json.gz",
    PREVIOUS_RESULT,
    "src/openai_ns_reconstruction/gamma_weighted_energy_delta.py",
    "scripts/run_gamma_weighted_energy_delta_once.py",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


def _json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if hasattr(value, "_mpf_"):
        return str(value)
    if type(value).__module__.split(".", 1)[0] == "sympy":
        return str(value)
    raise TypeError(f"not JSON serializable: {type(value).__name__}")


def _write_exclusive(path: Path, value: Mapping[str, Any]) -> str:
    raw = (
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            default=_json_default,
        )
        + "\n"
    ).encode("utf-8")
    with path.open("xb") as stream:
        stream.write(raw)
    return hashlib.sha256(raw).hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def _relative_path(raw: Any, label: str) -> str:
    _require(isinstance(raw, str) and bool(raw), f"{label} path is not text")
    path = Path(raw)
    if path.is_absolute():
        try:
            return path.resolve().relative_to(ROOT.resolve()).as_posix()
        except ValueError as exc:
            raise RuntimeError(f"{label} path is outside repository") from exc
    _require(".." not in path.parts, f"{label} path escapes repository")
    return path.as_posix()


def _manifest_hashes(source: Mapping[str, Any]) -> dict[str, str]:
    """Read direct path/hash records only; never walk saved source metadata."""

    result: dict[str, str] = {}

    def add(raw_path: Any, raw_digest: Any, label: str) -> None:
        relative = _relative_path(raw_path, label)
        _require(
            isinstance(raw_digest, str) and SHA256_RE.fullmatch(raw_digest) is not None,
            f"{label} has no valid SHA-256",
        )
        digest = raw_digest.lower()
        prior = result.get(relative)
        _require(prior is None or prior == digest, f"conflicting SHA-256 for {relative}")
        result[relative] = digest

    for key in ("sha256", "input_hashes", "input_sha256", "raw_input_sha256", "code_sha256"):
        mapping = source.get(key)
        if isinstance(mapping, Mapping):
            for raw_path, raw_digest in mapping.items():
                add(raw_path, raw_digest, f"source.{key}[{raw_path!r}]")
    files = source.get("files")
    _require(isinstance(files, list) and bool(files), "source files manifest missing")
    for index, item in enumerate(files):
        _require(isinstance(item, Mapping), f"source.files[{index}] is not an object")
        add(item.get("path"), item.get("sha256"), f"source.files[{index}]")
    return result


def _brief_checkpoint(brief: Mapping[str, Any]) -> tuple[str, str]:
    for key in (f"{COMPONENT}_source_checkpoint", "source_checkpoint"):
        value = brief.get(key)
        if isinstance(value, Mapping):
            path, digest = value.get("path"), value.get("sha256")
            if isinstance(path, str) and isinstance(digest, str):
                return path, digest
    path = brief.get("source_checkpoint_path")
    digest = brief.get("source_checkpoint_sha256")
    _require(
        isinstance(path, str) and isinstance(digest, str),
        "pending source checkpoint pointer missing",
    )
    return path, digest


def _validate_previous_result(hashes: Mapping[str, str]) -> None:
    expected = hashes.get(PREVIOUS_RESULT)
    _require(expected is not None, "previous Gamma delta result is not pinned")
    path = ROOT / PREVIOUS_RESULT
    _require(path.is_file(), "previous Gamma delta result is missing")
    _require(sha256(path) == expected, "previous Gamma delta result changed")
    previous = json.loads(path.read_bytes())
    _require(
        previous.get("component") == "gamma_correlated_scale_delta_v1"
        and previous.get("status") == "returned",
        "previous Gamma delta receipt is not returned",
    )
    body = previous.get("result")
    _require(isinstance(body, Mapping) and body.get("checks_passed") is True,
             "previous Gamma delta synthetic receipt did not pass")
    _require(
        previous.get("source_values_computed") is False
        and previous.get("source_value_computation") is False,
        "previous Gamma delta source-value flags are not accepted",
    )


def _preflight() -> tuple[dict[str, Any], dict[str, Any], str, dict[str, str]]:
    _require(not ENTERED.exists() and not RESULT.exists(), "already entered; never re-enter")
    _require(BRIEF.is_file(), f"pending brief is missing: {BRIEF}")
    _require(SOURCE.is_file(), f"source checkpoint is missing: {SOURCE}")
    brief = json.loads(BRIEF.read_bytes())
    source = json.loads(SOURCE.read_bytes())
    _require(isinstance(brief, Mapping), "pending brief is not a JSON object")
    _require(isinstance(source, Mapping), "source checkpoint is not a JSON object")
    _require(brief.get("component") == COMPONENT, "pending brief component mismatch")
    _require(brief.get("phase") == PHASE, "pending brief phase mismatch")
    _require(
        brief.get("goal_active") is True
        and brief.get("full_goal_complete") is False
        and brief.get("do_not_reenter") is True,
        "pending goal flags mismatch",
    )
    _require(source.get("schema") == SCHEMA, "source checkpoint schema mismatch")
    _require(source.get("component") == COMPONENT, "source checkpoint component mismatch")
    _require(source.get("one_entry_only") is True, "source checkpoint is not one-entry-only")
    _require(source.get("no_source_graph_traversal") is True,
             "source graph traversal guard missing")
    _require(
        source.get("all_raw_inputs_authenticated_before_decompression") is True,
        "raw-before-decompression guard missing",
    )
    goal_flags = source.get("goal_flags")
    _require(
        isinstance(goal_flags, Mapping)
        and goal_flags.get("goal_active") is True
        and goal_flags.get("full_goal_complete") is False
        and goal_flags.get("do_not_reenter") is True,
        "source checkpoint goal flags mismatch",
    )

    raw_path, expected_checkpoint_sha = _brief_checkpoint(brief)
    checkpoint_relative = _relative_path(raw_path, "brief source checkpoint")
    _require((ROOT / checkpoint_relative).resolve() == SOURCE.resolve(),
             "brief points to another checkpoint")
    checkpoint_sha = sha256(SOURCE)
    _require(
        SHA256_RE.fullmatch(expected_checkpoint_sha) is not None
        and checkpoint_sha == expected_checkpoint_sha.lower(),
        "source checkpoint SHA mismatch",
    )

    hashes = _manifest_hashes(source)
    _require(REQUIRED_INPUTS.issubset(hashes),
             "source checkpoint omits required Gamma/complete-future/code input")
    for relative, expected in hashes.items():
        path = ROOT / relative
        _require(path.is_file(), f"frozen input is missing: {relative}")
        _require(sha256(path) == expected, f"frozen raw input changed: {relative}")
    return dict(brief), dict(source), checkpoint_sha, hashes


def main() -> int:
    _, source, checkpoint_sha, input_hashes = _preflight()
    entry_sha = _write_exclusive(
        ENTERED,
        {
            "schema": f"{COMPONENT}_entered_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "entered": True,
            "source_checkpoint_sha256": checkpoint_sha,
            "input_hashes": dict(input_hashes),
            "consumer_calls_started": 0,
            "ancestor_math_calls": 0,
            "source_values_computed": False,
            "source_value_computation": False,
            "runtime_model_metadata": None,
        },
    )

    consumer_calls = 0
    try:
        spec = importlib.util.spec_from_file_location(
            "gamma_weighted_energy_delta_new_leaf", LEAF
        )
        if spec is None or spec.loader is None:
            raise RuntimeError("new Gamma weighted-energy leaf load failed")
        leaf = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = leaf
        spec.loader.exec_module(leaf)
        consumer_calls = 1
        value = leaf.verify_synthetic_fixture(ROOT, source, input_hashes)
        _require(isinstance(value, Mapping), "synthetic fixture did not return an object")
        _require(value.get("checks_passed") is True, "synthetic fixture checks did not pass")
        _require(value.get("certified") is False, "synthetic fixture must not claim certification")
        output: dict[str, Any] = {
            "schema": f"{COMPONENT}_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "returned",
            "source_checkpoint_sha256": checkpoint_sha,
            "entry_sha256": entry_sha,
            "input_hashes": dict(input_hashes),
            "consumer_calls": consumer_calls,
            "ancestor_math_calls": 0,
            "numeric_evaluation": False,
            "exact_symbolic_fixture": True,
            "certified": False,
            "source_values_computed": False,
            "source_value_computation": False,
            "runtime_model_metadata": None,
            "error": None,
            "result": value,
        }
    except Exception as exc:  # complete failure artifact is required
        output = {
            "schema": f"{COMPONENT}_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "failed",
            "source_checkpoint_sha256": checkpoint_sha,
            "entry_sha256": entry_sha,
            "input_hashes": dict(input_hashes),
            "consumer_calls": consumer_calls,
            "ancestor_math_calls": 0,
            "numeric_evaluation": False,
            "exact_symbolic_fixture": True,
            "certified": False,
            "source_values_computed": False,
            "source_value_computation": False,
            "runtime_model_metadata": None,
            "error": {"type": type(exc).__name__, "message": str(exc)},
            "traceback": traceback.format_exc(),
        }

    result_sha = _write_exclusive(RESULT, output)
    print(
        json.dumps(
            {
                "component": COMPONENT,
                "phase": PHASE,
                "status": output["status"],
                "result_sha256": result_sha,
                "consumer_calls": consumer_calls,
                "numeric_evaluation": False,
                "exact_symbolic_fixture": True,
                "certified": False,
            },
            sort_keys=True,
        )
    )
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
