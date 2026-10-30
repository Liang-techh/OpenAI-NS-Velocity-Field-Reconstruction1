"""One-entry runner for the bounded F source-jet binding consumer.

The pending checkpoint must pin exactly three raw files: the new leaf, this
runner, and the one authenticated prior result
``tmp/pulse_gap_actual_source_program_v1/result.json``.  The prior result is
hashed as raw bytes only; this runner never parses or imports it.  No archive,
graph, child row, ancestor, numerical evaluator, integral, or quadrature is an
input to this component.

The leaf API is deliberately small and explicit:

``verify_synthetic_fixture()``
    One in-memory guard fixture.  Its packet must be marked ``synthetic_only``
    and must not make any source, derivative, radius, coordinate, numerical,
    closure, or error-removal claim.

``bind_f_gap_source_jets(root, prior_result_path, prior_result_sha256,
                         source_manifest, input_hashes)``
    The one new source consumer.  It may bind source-role/jet records from the
    authenticated prior packet, but must return a non-empty
     ``adapter_input_order`` and ``adapter_inputs`` packet.  The packet must
     contain exactly the sixteen fixed Q-adapter operands; every operand may
     remain unresolved, and no opaque handle may be promoted to a value or
     authentication claim.  It must not evaluate values or walk a graph.

The runner writes ``entered.json`` before importing the leaf and writes a
complete JSON-safe success or failure receipt before printing.  There are no
retries: at most one synthetic call and one source-consumer call are recorded.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import re
import sys
import traceback
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[1]
COMPONENT = "f_gap_source_jet_bindings_v1"
PHASE = COMPONENT + "_pending"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"

RUNNER = Path("scripts/run_f_gap_source_jet_bindings_once.py")
LEAF = Path("src/openai_ns_reconstruction/f_gap_source_jet_bindings.py")
PRIOR = Path("tmp/pulse_gap_actual_source_program_v1/result.json")
PRIOR_SHA256 = "e611cdc7e049ce40e6667275b4e9706c3478044880b92c93baebf8ab8d06ebf2"
REQUIRED_INPUTS = (LEAF.as_posix(), RUNNER.as_posix(), PRIOR.as_posix())
FIXED_SHA256 = {PRIOR.as_posix(): PRIOR_SHA256}
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")
EXPECTED_ADAPTER_INPUT_ORDER = (
    "U0",
    "U1",
    "U2",
    "m0",
    "m1",
    "m2",
    "m3",
    "Uy0",
    "Uy1",
    "my0",
    "my1",
    "my2",
    "delta0",
    "delta1",
    "delta2",
    "Z",
)


# These names are rejected when a returned packet sets them to true.  The
# consumer is structural/source-role binding only; a packet cannot promote
# this component into a numerical or closure result by adding metadata.
FORBIDDEN_TRUE_KEYS = frozenset(
    {
        "archive_read",
        "archive_decoded",
        "archive_whole_json_decoded",
        "graph_access",
        "graph_traversed",
        "source_graph_traversed",
        "children_access",
        "children_indexed",
        "child_rows_indexed",
        "child_dereferenced",
        "graph_rows_indexed",
        "graph_edges_followed",
        "ancestor_replayed",
        "ancestor_called",
        "old_leaf_imported",
        "old_leaf_called",
        "old_producer_imported",
        "old_producer_called",
        "source_constructor_called",
        "source_checker_called",
        "source_values_selected",
        "source_value_selected",
        "source_values_computed",
        "actual_source_values",
        "numeric_values_computed",
        "numeric_evaluation",
        "numeric_use_admissible",
        "integrals_evaluated",
        "integral_evaluated",
        "quadrature_evaluated",
        "source_identity",
        "source_identity_proved",
        "same_source",
        "same_source_proved",
        "same_source_function_proved",
        "function_identity_proved",
        "m3_ordinary_derivative_proved",
        "ordinary_Z_derivative_proved",
        "ordinary_Z_derivative_semantics_proved",
        "ordinary_z_derivative_proved",
        "ordinary_Z_semantics_proved",
        "ordinary_Z3_semantics_proved",
        "ordinary_Z_transport_proved",
        "ordinary_z_transport_proved",
        "m3_ordinary_Z3_semantics_proved",
        "m3_extra_row_promoted_to_ordinary_Z3",
        "delta1_binding_proved",
        "delta2_binding_proved",
        "Z_binding_proved",
        "source_capture_proved",
        "physical_scale_proved",
        "namespace_identity_proved",
        "namespace_candidates_verified",
        "owner_candidates_verified",
        "namespace_verified",
        "owner_verified",
        "field_radius_owner_proved",
        "radius_metadata_promoted_to_field_owner",
        "coordinate_literal_is_value",
        "node_ids_are_values",
        "admissible_for_numeric_use",
        "physical_radius_owner_bound",
        "physical_radius_bound",
        "radius_identity_proved",
        "coordinate_identity_proved",
        "coordinates_identity_proved",
        "closure",
        "closure_proved",
        "m15_closure_proved",
        "error_deletion",
        "error_removed",
        "errors_removed",
        "old_error_removed",
        "old_errors_removed",
        "quantitative_error_reduced",
        "ledger_212_213_183_changed",
    }
)

SEMANTIC_FALSE_KEYS = frozenset(FORBIDDEN_TRUE_KEYS)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_load(path: Path, label: str) -> Mapping[str, Any]:
    try:
        value = json.loads(path.read_bytes().decode("utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise RuntimeError(label + " is not valid JSON") from exc
    _require(isinstance(value, Mapping), label + " is not an object")
    return value


def _json_safe(value: Any) -> Any:
    """Serialize a receipt without evaluating or narrowing backend values."""

    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_json_safe(item) for item in value]
    # A source-role object is receipt data, not a value to coerce.  repr keeps
    # the runner from importing a backend or extracting an interval midpoint.
    return repr(value)


def _write_once(path: Path, value: Any) -> str:
    raw = (
        json.dumps(
            _json_safe(value),
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


def _repo_path(raw_path: Any, label: str) -> tuple[str, Path]:
    _require(isinstance(raw_path, str) and raw_path, label + " path missing")
    candidate = Path(raw_path)
    path = candidate if candidate.is_absolute() else ROOT / candidate
    path = path.resolve()
    _require(path.is_relative_to(ROOT.resolve()), label + " escapes repository")
    return path.relative_to(ROOT.resolve()).as_posix(), path


def _manifest_hashes(source: Mapping[str, Any]) -> dict[str, str]:
    files = source.get("files")
    _require(isinstance(files, list), "source checkpoint files missing")
    seen: dict[str, str] = {}
    for index, item in enumerate(files):
        _require(isinstance(item, Mapping), f"source.files[{index}] is not an object")
        key, path = _repo_path(item.get("path"), f"source.files[{index}]")
        _require(key not in seen, "duplicate source input: " + key)
        expected = item.get("sha256")
        _require(
            isinstance(expected, str) and SHA256_RE.fullmatch(expected) is not None,
            "invalid SHA-256 for " + key,
        )
        _require(path.is_file(), "source input is missing: " + key)
        actual = _sha256_file(path)
        _require(actual.lower() == expected.lower(), "raw SHA changed: " + key)
        fixed = FIXED_SHA256.get(key)
        if fixed is not None:
            _require(expected.lower() == fixed, "fixed prior SHA changed: " + key)
        seen[key] = actual

    expected_paths = set(REQUIRED_INPUTS)
    _require(
        set(seen) == expected_paths,
        "source checkpoint must pin exactly leaf, runner, and prior result; "
        "missing/extra=" + ",".join(sorted(expected_paths ^ set(seen))),
    )
    return seen


def _checkpoint_sha(brief: Mapping[str, Any]) -> str:
    checkpoint = brief.get(f"{COMPONENT}_source_checkpoint")
    if not isinstance(checkpoint, Mapping):
        checkpoint = brief.get("source_checkpoint")
    _require(isinstance(checkpoint, Mapping), "pending brief source checkpoint missing")
    raw_path = checkpoint.get("path")
    _, path = _repo_path(raw_path, "pending source checkpoint")
    _require(path == SOURCE.resolve(), "pending brief points to another checkpoint")
    expected = checkpoint.get("sha256")
    _require(
        isinstance(expected, str) and SHA256_RE.fullmatch(expected) is not None,
        "pending source checkpoint SHA missing",
    )
    actual = _sha256_file(SOURCE)
    _require(actual.lower() == expected.lower(), "source checkpoint SHA changed")
    return actual


def _walk_true_keys(value: Any, path: str = "$") -> list[str]:
    found: list[str] = []
    if isinstance(value, Mapping):
        for key, item in value.items():
            key_text = str(key)
            child = path + "." + key_text
            if key_text in FORBIDDEN_TRUE_KEYS and item is True:
                found.append(child)
            found.extend(_walk_true_keys(item, child))
    elif isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            found.extend(_walk_true_keys(item, path + "[" + str(index) + "]"))
    return found


def _walk_false_keys(value: Any, path: str = "$") -> list[tuple[str, Any]]:
    found: list[tuple[str, Any]] = []
    if isinstance(value, Mapping):
        for key, item in value.items():
            key_text = str(key)
            child = path + "." + key_text
            if key_text in SEMANTIC_FALSE_KEYS:
                found.append((child, item))
            found.extend(_walk_false_keys(item, child))
    elif isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            found.extend(_walk_false_keys(item, path + "[" + str(index) + "]"))
    return found


def _validate_scope(value: Mapping[str, Any], label: str) -> None:
    true_claims = _walk_true_keys(value)
    _require(not true_claims, label + " claims forbidden true scope: " + ",".join(true_claims))
    bad_false = [path for path, marker in _walk_false_keys(value) if marker is not False]
    _require(not bad_false, label + " has non-false restricted flags: " + ",".join(bad_false))


def _validate_opaque_handles(value: Any, label: str, path: str = "$") -> None:
    """Require every advertised source handle to remain explicitly opaque."""

    if isinstance(value, Mapping):
        if "source_node_id" in value:
            _require(
                value.get("node_id_is_not_a_numeric_value") is True,
                label + " exposes a source_node_id without an opaque-value guard at " + path,
            )
            if "source_value_selected" in value:
                _require(
                    value.get("source_value_selected") is False,
                    label + " promotes a source handle to a selected value at " + path,
                )
        for key, item in value.items():
            _validate_opaque_handles(item, label, path + "." + str(key))
    elif isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            _validate_opaque_handles(item, label, path + "[" + str(index) + "]")


def _preflight() -> tuple[dict[str, Any], dict[str, str], str]:
    # Preflight is read-only.  In particular, the prior result is hashed as
    # raw bytes and is never decoded or imported here.
    _require(not ENTRY.exists() and not RESULT.exists(), "component already entered")
    _require(DIRECTORY.is_dir(), "component directory is missing")
    _require(BRIEF.is_file(), "pending brief is missing")
    _require(SOURCE.is_file(), "source checkpoint is missing")
    brief = _json_load(BRIEF, "pending brief")
    source = _json_load(SOURCE, "source checkpoint")
    _require(brief.get("component") == COMPONENT, "pending brief component mismatch")
    _require(brief.get("phase") == PHASE, "pending brief phase mismatch")
    _require(
        brief.get("goal_active") is True
        and brief.get("full_goal_complete") is False
        and brief.get("do_not_reenter") is True,
        "pending brief goal flags mismatch",
    )
    source_sha = _checkpoint_sha(brief)
    _require(source.get("schema") == COMPONENT + "_checkpoint_v1", "source schema mismatch")
    _require(source.get("component") == COMPONENT, "source component mismatch")
    _require(source.get("one_entry_only") is True, "one-entry guard missing")
    _require(source.get("raw_sha256_authenticated_before_entry") is True,
             "raw SHA pre-entry guard missing")
    _require(source.get("allow_missing_m3") is True, "missing-m3 guard missing")
    _require(source.get("allow_unresolved_radius_metadata") is True,
             "unresolved-radius guard missing")
    for name in (
        "archive_access", "graph_access", "children_access", "ancestor_replay",
        "numeric_evaluation", "integrals_evaluated", "source_identity_proved",
        "m3_ordinary_derivative_proved", "physical_radius_owner_bound",
        "coordinate_identity_proved", "closure_proved", "error_deletion",
        "ledger_212_213_183_changed",
    ):
        _require(source.get(name) is False, "source checkpoint guard must be false: " + name)
    _validate_scope(source, "source checkpoint")
    hashes = _manifest_hashes(source)
    return dict(source), hashes, source_sha


def _import_leaf() -> Any:
    spec = importlib.util.spec_from_file_location("f_gap_source_jet_bindings_v1_new_leaf", ROOT / LEAF)
    _require(spec is not None and spec.loader is not None, "new source-jet leaf load failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _callable(module: Any, name: str) -> Callable[..., Any]:
    value = getattr(module, name, None)
    _require(callable(value), "new leaf callable is missing: " + name)
    return value


def _require_synthetic(value: Any) -> Mapping[str, Any]:
    _require(isinstance(value, Mapping), "synthetic fixture did not return an object")
    _require(
        value.get("passed") is True or value.get("checks_passed") is True,
        "synthetic fixture did not pass",
    )
    _require(value.get("synthetic_only") is True, "synthetic fixture is not synthetic-only")
    _validate_scope(value, "synthetic fixture")
    return value


def _require_source_packet(value: Any) -> Mapping[str, Any]:
    _require(isinstance(value, Mapping), "source consumer did not return an object")
    _require(value.get("schema") == COMPONENT, "source consumer schema mismatch")
    _require(value.get("source_only") is True, "source consumer is not source-only")
    _require(value.get("status") not in {None, "failed"}, "source consumer did not return")
    for flag in (
        "admissible_for_numeric_use",
        "source_identity_proved",
        "source_values_selected",
        "source_values_computed",
        "source_graph_traversed",
        "graph_edges_followed",
        "node_ids_are_values",
        "ordinary_Z_semantics_proved",
        "ordinary_Z3_semantics_proved",
        "source_capture_proved",
        "physical_scale_proved",
        "field_radius_owner_proved",
        "namespace_identity_proved",
    ):
        _require(value.get(flag) is False, "source consumer flag must be false: " + flag)
    order = value.get("adapter_input_order")
    _require(
        isinstance(order, list)
        and tuple(order) == EXPECTED_ADAPTER_INPUT_ORDER,
        "source consumer adapter_input_order is not the fixed sixteen-field order",
    )
    inputs = value.get("adapter_inputs")
    _require(isinstance(inputs, list) and len(inputs) == 16,
             "source consumer adapter_inputs must contain exactly sixteen entries")
    _require(
        value.get("adapter_input_count") == 16,
        "source consumer adapter_input_count must be exactly sixteen",
    )
    names: list[Any] = []
    for index, item in enumerate(inputs):
        _require(isinstance(item, Mapping), f"adapter_inputs[{index}] is not an object")
        name = item.get("input_name")
        names.append(name)
        _require(
            name == EXPECTED_ADAPTER_INPUT_ORDER[index],
            f"adapter_inputs[{index}] input_name/order mismatch",
        )
        binding = item.get("binding")
        _require(isinstance(binding, Mapping),
                 f"adapter_inputs[{index}] binding is not an object")
        _require(
            binding.get("status") in {"bound_opaque_source_handle", "unresolved"},
            f"adapter_inputs[{index}] has an unrecognized binding status",
        )
        _require(
            binding.get("node_id_is_not_a_numeric_value") is True,
            f"adapter_inputs[{index}] lacks the opaque-node guard",
        )
        _require(
            binding.get("source_value_selected") is False,
            f"adapter_inputs[{index}] promotes a handle to a selected value",
        )
        if binding.get("status") == "bound_opaque_source_handle":
            _require(
                "source_node_id" in binding,
                f"adapter_inputs[{index}] bound handle has no source_node_id",
            )
        _validate_opaque_handles(binding, "source adapter binding", f"$.adapter_inputs[{index}]")
    _require(tuple(names) == EXPECTED_ADAPTER_INPUT_ORDER,
             "source consumer adapter input names are duplicated or reordered")
    # This recursive check covers F rows, geometry, coordinate metadata, and
    # the adapter packet, so opaque node ids cannot become numeric values or
    # implicit source-authentication claims in a nested field.
    _validate_opaque_handles(value, "source consumer packet")
    _validate_scope(value, "source consumer")
    return value


def main() -> int:
    source, hashes, source_sha = _preflight()
    entry_sha = _write_once(
        ENTRY,
        {
            "schema": COMPONENT + "_entry_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "entered": True,
            "source_checkpoint_sha256": source_sha,
            "input_hashes": hashes,
            "leaf_imported": False,
            "synthetic_fixture_invocations": 0,
            "source_consumer_invocations": 0,
            "archive_access": False,
            "graph_access": False,
            "children_access": False,
            "ancestor_replay": False,
            "numeric_evaluation": False,
            "integrals_evaluated": False,
            "source_identity_proved": False,
            "m3_ordinary_derivative_proved": False,
            "physical_radius_owner_bound": False,
            "coordinate_identity_proved": False,
            "closure_proved": False,
            "error_deletion": False,
            "ledger_212_213_183_changed": False,
        },
    )

    steps: list[dict[str, Any]] = []
    fixture_result: Any = None
    source_result: Any = None
    leaf_imported = False
    try:
        module = _import_leaf()
        leaf_imported = True

        fixture = _callable(module, "verify_synthetic_fixture")
        steps.append({"step": "verify_synthetic_fixture", "invoked": True})
        fixture_result = fixture()
        steps[-1]["returned_packet"] = _json_safe(fixture_result)
        _require_synthetic(fixture_result)

        consumer = _callable(module, "bind_f_gap_source_jets")
        steps.append({"step": "bind_f_gap_source_jets", "invoked": True})
        source_result = consumer(
            root=ROOT,
            prior_result_path=ROOT / PRIOR,
            prior_result_sha256=PRIOR_SHA256,
            source_manifest=source,
            input_hashes=hashes,
        )
        steps[-1]["returned_packet"] = _json_safe(source_result)
        _require_source_packet(source_result)

        output: dict[str, Any] = {
            "schema": COMPONENT + "_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "returned",
            "source_checkpoint_sha256": source_sha,
            "entry_sha256": entry_sha,
            "input_hashes": hashes,
            "leaf_imported": leaf_imported,
            "consumer_calls": [item["step"] for item in steps],
            "synthetic_fixture_invocations": 1,
            "source_consumer_invocations": 1,
            "archive_access": False,
            "graph_access": False,
            "children_access": False,
            "ancestor_replay": False,
            "numeric_evaluation": False,
            "integrals_evaluated": False,
            "source_identity_proved": False,
            "m3_ordinary_derivative_proved": False,
            "physical_radius_owner_bound": False,
            "coordinate_identity_proved": False,
            "closure_proved": False,
            "error_deletion": False,
            "ledger_212_213_183_changed": False,
            "old_errors_removed": False,
            "reference_offset_removed": False,
            "reserved_contribution_removed": False,
            "m3_missing_allowed": True,
            "radius_metadata_unresolved_allowed": True,
            "result": {
                "synthetic_fixture": _json_safe(fixture_result),
                "adapter_source_packet": _json_safe(source_result),
            },
            "source_manifest": _json_safe(source),
            "steps": _json_safe(steps),
            "error": None,
        }
        _validate_scope(output, "runner success result")
    except BaseException as exc:
        output = {
            "schema": COMPONENT + "_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "failed",
            "source_checkpoint_sha256": source_sha,
            "entry_sha256": entry_sha,
            "input_hashes": hashes,
            "leaf_imported": leaf_imported,
            "consumer_calls": [item["step"] for item in steps],
            "synthetic_fixture_invocations": sum(
                item.get("step") == "verify_synthetic_fixture" for item in steps
            ),
            "source_consumer_invocations": sum(
                item.get("step") == "bind_f_gap_source_jets" for item in steps
            ),
            "archive_access": False,
            "graph_access": False,
            "children_access": False,
            "ancestor_replay": False,
            "numeric_evaluation": False,
            "integrals_evaluated": False,
            "source_identity_proved": False,
            "m3_ordinary_derivative_proved": False,
            "physical_radius_owner_bound": False,
            "coordinate_identity_proved": False,
            "closure_proved": False,
            "error_deletion": False,
            "ledger_212_213_183_changed": False,
            "old_errors_removed": False,
            "reference_offset_removed": False,
            "reserved_contribution_removed": False,
            "m3_missing_allowed": True,
            "radius_metadata_unresolved_allowed": True,
            "result": {
                "synthetic_fixture": _json_safe(fixture_result),
                "adapter_source_packet": _json_safe(source_result),
            },
            "source_manifest": _json_safe(source),
            "steps": _json_safe(steps),
            "error": {
                "type": type(exc).__name__,
                "message": str(exc),
                "traceback": traceback.format_exc(),
            },
        }

    # The complete receipt is persisted before the one-line summary is shown.
    result_sha = _write_once(RESULT, output)
    print(
        json.dumps(
            {
                "component": COMPONENT,
                "phase": PHASE,
                "status": output["status"],
                "source_checkpoint_sha256": source_sha,
                "entry_sha256": entry_sha,
                "result_sha256": result_sha,
                "consumer_calls": output["consumer_calls"],
                "synthetic_fixture_invocations": output["synthetic_fixture_invocations"],
                "source_consumer_invocations": output["source_consumer_invocations"],
                "m3_missing_allowed": True,
                "radius_metadata_unresolved_allowed": True,
            },
            sort_keys=True,
        )
    )
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
