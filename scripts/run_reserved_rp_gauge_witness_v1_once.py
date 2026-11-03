"""Run the reserved-Rp/gauge witness consumer exactly once.

The runner authenticates every raw byte listed by ``source.json`` before it
creates the exclusive entry receipt.  Only after entry does it import the new
pure consumer modules or decode the authenticated inner-Volterra archive.
Historical producer modules are never imported or called.  Historical source
files and check receipts are read as static AST/JSON evidence only.
"""

from __future__ import annotations

import ast
import dataclasses
import gzip
import hashlib
import importlib.util
import json
import re
import sys
import traceback
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
NAME = "reserved_rp_gauge_witness_v1"
DIRECTORY = ROOT / "tmp" / NAME
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entry.json"
RESULT = DIRECTORY / "result.json"

GAUGE = "src/openai_ns_reconstruction/n1_core_gauge_binding_witness_v2.py"
BRIDGE = "src/openai_ns_reconstruction/reserved_rp_support_bridge_v1.py"
VARIABLE_Z = "src/openai_ns_reconstruction/reserved_m15_variable_z.py"
RUNNER = "scripts/run_reserved_rp_gauge_witness_v1_once.py"
INNER_VOLTERRA = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.json.gz"
)
FIXED_POINT_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py"
)
RESERVED_LINEAR_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_reserved_linear_operator.py"
)
RESERVED_LINEAR_CHECK = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_reserved_linear_operator_check.json"
)
C3_POWER_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_C3_power_to_Rp.py"
)
C3_POWER_CHECK = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_C3_power_to_Rp_check.json"
)

GRAPH_SHA256 = "1b44be8c6cbfffb692c7a0430ecbc8a590f46001a233c99c59ca8c5e2c44a6fe"
INTEGRAL_NODE_ID = 35646
BINDER_NODE_ID = 35644
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")


class RunnerInputError(ValueError):
    """Raised when the frozen one-entry input contract is not satisfied."""


def _require(condition: object, message: str) -> None:
    if not condition:
        raise RunnerInputError(message)


def _sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_constant(token: str) -> Any:
    raise ValueError("non-finite JSON value is not admitted: " + token)


def _load_json(raw: bytes, label: str) -> Any:
    try:
        return json.loads(raw.decode("utf-8"), parse_constant=_json_constant)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise RunnerInputError(label + " is not finite UTF-8 JSON") from exc


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


def _repository_path(raw_path: Any, label: str) -> tuple[str, Path]:
    _require(isinstance(raw_path, str) and raw_path, label + " path is missing")
    relative = Path(raw_path.replace("\\", "/"))
    _require(not relative.is_absolute() and ".." not in relative.parts, label + " must be repository-relative")
    key = relative.as_posix()
    path = (ROOT / relative).resolve()
    _require(path.is_relative_to(ROOT.resolve()) and path.is_file(), label + " is missing")
    return key, path


def _manifest_hashes(source: Mapping[str, Any]) -> dict[str, str]:
    rows = source.get("files")
    _require(isinstance(rows, list), "source.files is missing")
    hashes: dict[str, str] = {}
    for index, row in enumerate(rows):
        label = f"source.files[{index}]"
        _require(isinstance(row, Mapping), label + " is not an object")
        key, path = _repository_path(row.get("path"), label)
        _require(key not in hashes, "duplicate source file: " + key)
        expected = row.get("raw_sha256", row.get("sha256"))
        if row.get("raw_sha256") is not None and row.get("sha256") is not None:
            _require(str(row["raw_sha256"]).lower() == str(row["sha256"]).lower(), label + " SHA fields disagree")
        _require(isinstance(expected, str) and SHA256_RE.fullmatch(expected) is not None, label + " has invalid SHA-256")
        actual = _sha_file(path)
        _require(actual.lower() == expected.lower(), "raw SHA-256 mismatch: " + key)
        if key.endswith(".gz"):
            _require(
                row.get("raw_sha256_before_decode") is True
                or row.get("raw_bytes_authenticated_before_decode") is True,
                "gzip input is not marked raw-authenticated: " + key,
            )
        hashes[key] = actual
    required = {
        GAUGE,
        BRIDGE,
        VARIABLE_Z,
        RUNNER,
        INNER_VOLTERRA,
        FIXED_POINT_SOURCE,
        RESERVED_LINEAR_SOURCE,
        RESERVED_LINEAR_CHECK,
        C3_POWER_SOURCE,
        C3_POWER_CHECK,
    }
    _require(required.issubset(hashes), "source manifest is missing required frozen inputs")
    _require(hashes[INNER_VOLTERRA].lower() == GRAPH_SHA256, "inner_volterra graph raw SHA pin mismatch")
    return hashes


def _preflight() -> tuple[dict[str, Any], dict[str, str], str]:
    _require(not ENTRY.exists() and not RESULT.exists(), "component already entered")
    _require(SOURCE.is_file(), "source manifest is missing")
    source_raw = SOURCE.read_bytes()
    source = _load_json(source_raw, "source manifest")
    _require(isinstance(source, Mapping), "source manifest must be an object")
    _require(source.get("component") in (None, NAME), "source manifest component mismatch")
    _require(source.get("one_entry_only") in (None, True), "source manifest one-entry guard")
    hashes = _manifest_hashes(source)
    return dict(source), hashes, _sha(source_raw)


def _write_exclusive(path: Path, value: Any) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = _json_bytes(value)
    try:
        with path.open("xb") as stream:
            stream.write(raw)
    except FileExistsError as exc:
        raise RunnerInputError("one-entry file already exists: " + str(path)) from exc
    return _sha(raw)


def _recheck_raw(hashes: Mapping[str, str], key: str) -> bytes:
    _require(key in hashes, "unhashed input: " + key)
    raw = (ROOT / key).read_bytes()
    _require(_sha(raw).lower() == hashes[key].lower(), "input changed after entry: " + key)
    return raw


def _import_new(path: str, name: str) -> Any:
    _require(path in {GAUGE, BRIDGE, VARIABLE_Z}, "unexpected import: " + path)
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    _require(spec is not None and spec.loader is not None, "module load failed: " + path)
    loaded = importlib.util.module_from_spec(spec)
    sys.modules[name] = loaded
    spec.loader.exec_module(loaded)
    return loaded


def _line_hits(text: str, needle: str) -> list[int]:
    return [line_no for line_no, line in enumerate(text.splitlines(), 1) if needle in line]


def _normalized_line_hits(text: str, needle: str) -> list[int]:
    return [
        line_no
        for line_no, line in enumerate(text.splitlines(), 1)
        if needle in "".join(line.split())
    ]


def _literal_assignment(tree: ast.AST, name: str) -> tuple[Any, list[int]]:
    hits: list[int] = []
    value: Any = None
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        target = node.targets[0]
        if not isinstance(target, ast.Name) or target.id != name:
            continue
        hits.append(int(node.lineno))
        try:
            value = ast.literal_eval(node.value)
        except (ValueError, TypeError, SyntaxError):
            value = None
    return value, hits


def _static_source_audit(linear_raw: bytes, c3_raw: bytes, linear_check_raw: bytes, c3_check_raw: bytes) -> dict[str, Any]:
    """Perform source/check evidence collection without importing history."""

    linear_text = linear_raw.decode("utf-8")
    c3_text = c3_raw.decode("utf-8")
    try:
        linear_tree = ast.parse(linear_text, filename=RESERVED_LINEAR_SOURCE, mode="exec")
        ast.parse(c3_text, filename=C3_POWER_SOURCE, mode="exec")
    except SyntaxError as exc:
        raise RunnerInputError("static source AST parse failed") from exc
    linear_check = _load_json(linear_check_raw, RESERVED_LINEAR_CHECK)
    c3_check = _load_json(c3_check_raw, C3_POWER_CHECK)
    _require(isinstance(linear_check, Mapping), "reserved linear check is not an object")
    _require(isinstance(c3_check, Mapping), "C3 power check is not an object")
    centers, centers_lines = _literal_assignment(linear_tree, "CENTERS")
    angular, angular_lines = _literal_assignment(linear_tree, "ANGULAR")
    ell_hits = _normalized_line_hits(linear_text, "ell=g.mul(L,g.constant('1/400'))")
    tw_hits = _normalized_line_hits(c3_text, "Tw=g.mul(g.constant(-60),g.unary('log',mu))")
    rp_hits = _line_hits(c3_text, "Rp=Rc*exp(Tw-2)")
    quiet_constructor_hits = _normalized_line_hits(
        c3_text,
        "length=g.sub(g.sub(Tw,two),g.unary('log',two))",
    )
    quiet_coordinate_hits = _line_hits(c3_text, "Tw-2-log2")
    _require(centers == (12, 13, 14, 15, 16), "reserved source CENTERS literal mismatch")
    _require(angular == (1, 2, 3), "reserved source ANGULAR literal mismatch")
    _require(ell_hits, "reserved source ell=log(2)/400 construction not found")
    _require(tw_hits, "C3 source Tw construction not found")
    _require(rp_hits, "C3 source Rp=Rc*exp(Tw-2) record not found")
    _require(quiet_constructor_hits, "C3 source quiet-length constructor not found")
    _require(quiet_coordinate_hits, "C3 source Tw-2-log2 coordinate record not found")
    positive_keys = (
        "actual_mu_positive_and_retained",
        "original_amplitude_positive_on_all_four_Z_cells",
        "new_supports_strictly_disjoint",
    )
    positive_records = {key: linear_check.get(key) for key in positive_keys}
    _require(linear_check.get("all_passed") is True, "reserved linear check all_passed is not True")
    _require(all(linear_check.get(key) is True for key in positive_keys), "reserved linear positive checks incomplete")
    geometry = None
    if isinstance(c3_check.get("independent_native_frame_and_geometry"), Mapping):
        geometry = c3_check["independent_native_frame_and_geometry"].get("exact_Rc_2Rc_Rp_and_quiet_length_identity")
    _require(c3_check.get("all_passed") is True, "C3 power check all_passed is not True")
    _require(geometry is True, "C3 positive-quiet identity check is not True")

    def bool_records(value: Any, prefix: str = "") -> dict[str, bool]:
        result: dict[str, bool] = {}
        if isinstance(value, Mapping):
            for key, item in value.items():
                path = prefix + "." + str(key) if prefix else str(key)
                if isinstance(item, bool) and any(token in str(key).lower() for token in ("positive", "quiet", "rp", "support")):
                    result[path] = item
                elif isinstance(item, Mapping):
                    result.update(bool_records(item, path))
        return result

    return {
        "static_only": True,
        "imports_executed": False,
        "cross_owner_binding_claimed": False,
        "reserved_linear_operator": {
            "source_ast_parsed": True,
            "CENTERS": {"value": list(centers), "lines": centers_lines},
            "ANGULAR": {"value": list(angular), "lines": angular_lines},
            "ell_log2_over_400_lines": ell_hits,
            "positive_quiet_records": positive_records,
            "all_passed": linear_check.get("all_passed"),
        },
        "C3_power_to_Rp": {
            "source_ast_parsed": True,
            "Tw_lines": tw_hits,
            "Rp_formula_lines": rp_hits,
            "quiet_length_formula_lines": quiet_constructor_hits,
            "quiet_coordinate_Tw_minus_2_minus_log2_lines": quiet_coordinate_hits,
            "positive_quiet_identity": geometry,
            "all_passed": c3_check.get("all_passed"),
        },
        "linear_check_positive_quiet_boolean_records": bool_records(linear_check),
        "c3_check_positive_quiet_boolean_records": bool_records(c3_check),
    }


def _safe(value: Any, *, context: str = "value") -> Any:
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, (bytes, bytearray)):
        return {"kind": "bytes", "sha256": _sha(bytes(value)), "length": len(value)}
    if hasattr(value, "_mpi_"):
        lower, upper = value._mpi_
        return {
            "kind": "interval",
            "lower_exact_mpf_tuple": list(lower),
            "upper_exact_mpf_tuple": list(upper),
        }
    if hasattr(value, "_mpf_"):
        return {"kind": "mpf", "exact_mpf_tuple": list(value._mpf_)}
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return {
            field.name: _safe(getattr(value, field.name), context=context + "." + field.name)
            for field in dataclasses.fields(value)
        }
    if isinstance(value, Mapping):
        return {str(key): _safe(item, context=context + "." + str(key)) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_safe(item, context=context) for item in value]
    if isinstance(value, (set, frozenset)):
        return sorted((_safe(item, context=context) for item in value), key=repr)
    if isinstance(value, Path):
        return value.as_posix()
    raise TypeError("cannot serialize " + context + " of type " + type(value).__name__)


def _case(name: str, function: Any, *, expect_error: bool = False) -> dict[str, Any]:
    try:
        output = function()
    except BaseException as exc:
        if not expect_error:
            return {
                "name": name,
                "status": "failed",
                "expected_error": False,
                "error": {
                    "type": type(exc).__name__,
                    "message": str(exc),
                    "traceback": traceback.format_exc(),
                },
            }
        return {
            "name": name,
            "status": "passed",
            "expected_error": True,
            "rejection": {
                "type": type(exc).__name__,
                "message": str(exc),
                "traceback": traceback.format_exc(),
            },
        }
    if expect_error:
        return {
            "name": name,
            "status": "failed",
            "expected_error": True,
            "error": {
                "type": "ExpectedRejectionMissing",
                "message": "case returned without the required rejection",
                "traceback": "",
            },
        }
    try:
        safe = _safe(output, context=name)
    except BaseException as exc:
        return {
            "name": name,
            "status": "failed",
            "expected_error": False,
            "error": {
                "type": type(exc).__name__,
                "message": str(exc),
                "traceback": traceback.format_exc(),
            },
        }
    return {"name": name, "status": "passed", "expected_error": False, "output": safe}


def _bridge_inputs(bridge: Any, token: object, rc_owner: object, parameter_owner: object) -> tuple[Any, Any, Any]:
    ports = (
        bridge.LiveSourceOwnerPort("reserved_bump", rc_owner, parameter_owner, token),
        bridge.LiveSourceOwnerPort("quiet_rp", rc_owner, parameter_owner, token),
    )
    mapping = bridge.PhysicalRpMapping(
        quiet_length_positive=bridge.PositiveSourceWitness("Tw-2-log2>0", token),
        rc_positive=bridge.PositiveSourceWitness("Rc>0", token),
        exponential_order=bridge.PositiveSourceWitness("exp_strictly_increasing_on_real", token),
    )
    matrix = bridge.FullSupportMatrixPremise(
        complete_m15_moment_row=object(),
        angular_matrix_third_row=object(),
        row_identity_witness=bridge.SourcePremiseWitness(
            "same_three_bump_complete_m15_moment_row_equals_angular_matrix_third_row", token
        ),
        invertibility_witness=bridge.SourcePremiseWitness(
            "same_three_bump_angular_matrix_invertible", token
        ),
        row_C1_identity_witness=bridge.SourcePremiseWitness(
            "same_three_bump_complete_m15_moment_row_equals_angular_matrix_third_row_as_full_Z_C1_function", token
        ),
    )
    return ports, mapping, matrix


def _interval_zero(value: Any) -> bool:
    raw = getattr(value, "_mpi_", None)
    if raw is None:
        return False
    lower, upper = raw
    return tuple(lower) == tuple(upper) and tuple(lower)[1] == 0


def _term_map(value: Any) -> dict[str, Any]:
    return {str(name): term for name, term in value}


def _conditional_variable_z_inputs(witness: Any, iv: Any) -> dict[str, Any]:
    """Request conditional response inputs with explicit acknowledgement."""
    method = getattr(witness, "conditional_variable_z_response_inputs", None)
    _require(callable(method), "support bridge lacks conditional_variable_z_response_inputs")
    return method(iv, acknowledge_unverified_conditions=True)


def _consume(hashes: Mapping[str, str], progress: dict[str, Any] | None = None) -> dict[str, Any]:
    # Every re-read is checked against the pre-entry manifest before use.
    raw = {
        key: _recheck_raw(hashes, key)
        for key in (
            GAUGE,
            BRIDGE,
            VARIABLE_Z,
            INNER_VOLTERRA,
            FIXED_POINT_SOURCE,
            RESERVED_LINEAR_SOURCE,
            RESERVED_LINEAR_CHECK,
            C3_POWER_SOURCE,
            C3_POWER_CHECK,
        )
    }
    if progress is not None:
        progress["raw_recheck_complete"] = True
    gauge = _import_new(GAUGE, NAME + "_gauge")
    bridge = _import_new(BRIDGE, NAME + "_bridge")
    variable_z = _import_new(VARIABLE_Z, NAME + "_variable_z")
    if progress is not None:
        progress["existing_pure_consumer_imported"] = True
    static_evidence = _static_source_audit(
        raw[RESERVED_LINEAR_SOURCE],
        raw[C3_POWER_SOURCE],
        raw[RESERVED_LINEAR_CHECK],
        raw[C3_POWER_CHECK],
    )
    graph = _load_json(gzip.decompress(raw[INNER_VOLTERRA]), INNER_VOLTERRA)
    _require(isinstance(graph, Mapping), "inner_volterra archive is not an object")
    exact_graph_nodes = graph.get("exact_graph_nodes")
    _require(isinstance(exact_graph_nodes, list), "inner_volterra exact_graph_nodes is missing")
    source_text = raw[FIXED_POINT_SOURCE].decode("utf-8")
    validator_result = gauge.build_core_gauge_binding_witness(
        exact_graph_nodes,
        source_text,
        integral_node_id=INTEGRAL_NODE_ID,
        binder_node_id=BINDER_NODE_ID,
        source_sha256=_sha(raw[FIXED_POINT_SOURCE]),
        integral_target="G",
    )
    if progress is not None:
        progress["validator"] = validator_result
        progress["static_source_evidence"] = static_evidence
        progress["graph_node_count"] = len(exact_graph_nodes)

    # Check two new structural counterexamples without evaluating the graph.
    altered = list(exact_graph_nodes)
    altered[35645] = dict(altered[35645], value=20)
    rejected_value = gauge.validate_core_gauge_binding_witness(
        altered, source_text, integral_node_id=INTEGRAL_NODE_ID,
        binder_node_id=BINDER_NODE_ID, integral_target="G")
    _require(rejected_value.get("passed") is False, "wrong substitution binder accepted")
    altered = list(exact_graph_nodes)
    altered[INTEGRAL_NODE_ID] = dict(altered[INTEGRAL_NODE_ID], upper=BINDER_NODE_ID)
    rejected_upper = gauge.validate_core_gauge_binding_witness(
        altered, source_text, integral_node_id=INTEGRAL_NODE_ID,
        binder_node_id=BINDER_NODE_ID, integral_target="G")
    _require(rejected_upper.get("passed") is False, "wrong upper endpoint accepted")
    if progress is not None:
        progress["graph_negative_cases"] = {"wrong_substitution_rejected": True, "wrong_upper_rejected": True}

    # The same interval context is passed to the bridge and the new probe;
    # no historical fixture or producer is called.
    iv_module = __import__("mpmath.ctx_iv", fromlist=["MPIntervalContext"])
    iv = iv_module.MPIntervalContext()
    token = object()
    rc_owner = object()
    parameter_owner = object()
    ports, mapping, matrix = _bridge_inputs(bridge, token, rc_owner, parameter_owner)
    shared: dict[str, Any] = {"iv": iv, "ports": ports, "mapping": mapping, "matrix": matrix}
    cases: list[dict[str, Any]] = []
    if progress is not None:
        progress["synthetic_cases"] = cases
    cases.append(
        _case(
            "support_bridge_valid_input",
            lambda: (
                shared.update(
                    witness=bridge.prove_reserved_rp_full_support(
                        owner_ports=ports, mapping=mapping, matrix_premise=matrix
                    )
                )
                or shared["witness"].report()
            ),
        )
    )
    cases.append(
        _case(
            "support_bridge_wrong_owner_rejected",
            lambda: bridge.prove_reserved_rp_full_support(
                owner_ports=(
                    bridge.LiveSourceOwnerPort("reserved_bump", rc_owner, parameter_owner, token),
                    bridge.LiveSourceOwnerPort("quiet_rp", object(), parameter_owner, token),
                ),
                mapping=mapping,
                matrix_premise=matrix,
            ),
            expect_error=True,
        )
    )
    fake_matrix = bridge.FullSupportMatrixPremise(
        complete_m15_moment_row=object(),
        angular_matrix_third_row=object(),
        row_identity_witness=bridge.SourcePremiseWitness("fake matrix premise", token),
        invertibility_witness=matrix.invertibility_witness,
        row_C1_identity_witness=matrix.row_C1_identity_witness,
    )
    cases.append(
        _case(
            "support_bridge_fake_matrix_premise_rejected",
            lambda: bridge.prove_reserved_rp_full_support(
                owner_ports=ports, mapping=mapping, matrix_premise=fake_matrix
            ),
            expect_error=True,
        )
    )
    cases.append(
        _case(
            "support_bridge_certified_row_rejected",
            lambda: shared["witness"].certified_exact_response_row(),
            expect_error=True,
        )
    )

    def response_inputs() -> dict[str, Any]:
        supplied = _conditional_variable_z_inputs(shared["witness"], iv)
        shared["response_inputs"] = supplied
        return {
            "same_iv_context_used": True,
            "r_i": supplied["r_i"],
            "r_i_Z": supplied["r_i_Z"],
        }

    cases.append(_case("support_bridge_variable_z_inputs_shared_iv", response_inputs))

    def exact_row_variable_z_probe() -> dict[str, Any]:
        supplied = shared["response_inputs"]
        zero = iv.mpf("0")

        def mpf(value: int) -> Any:
            return iv.mpf(str(value))

        result = variable_z.evaluate_partial_reserved_minus_pulse(
            iv,
            supplied["r_i"],
            supplied["r_i_Z"],
            zero,
            mpf(29),
            (mpf(3), mpf(5)),
            (mpf(7), mpf(11)),
            mpf(13),
            mpf(17),
            mpf(19),
            mpf(23),
            zero,
        )
        c0 = _term_map(result.C0_terms)
        z1 = _term_map(result.Z1_baseline_terms)
        sensitivity = _term_map(result.sensitivity_coefficients)
        _require("plus_r2_J" in c0 and "plus_one_minus_r2_P" in c0, "exact row C0 terms missing")
        _require("plus_r2_J_Z" in z1 and "plus_one_minus_r2_P_Z" in z1, "exact row Z1 terms missing")
        _require(_interval_zero(c0["plus_one_minus_r2_P"]), "exact row mixed a pulse P term")
        _require(_interval_zero(z1["plus_one_minus_r2_P_Z"]), "exact row mixed a pulse P_Z term")
        _require(not _interval_zero(sensitivity["r2_Z"]), "exact row erased J-P sensitivity")
        _require(_interval_zero(result.C0 - mpf(13)), "full-support C0 does not equal J")
        _require(result.Z1 is not None and _interval_zero(result.Z1 - mpf(17)), "full-support Z1 does not equal J_Z")
        return {
            "same_iv_context_used": True,
            "nonzero_P_and_P_Z_supplied": True,
            "nonzero_theta_and_mu_Z_supplied": True,
            "pulse_terms_zero_under_(0,0,1)": True,
            "J_minus_P_sensitivity_retained": True,
            "C0_terms": result.C0_terms,
            "Z1_baseline_terms": result.Z1_baseline_terms,
            "sensitivity_coefficients": result.sensitivity_coefficients,
        }

    cases.append(_case("variable_z_exact_response_row_probe", exact_row_variable_z_probe))
    synthetic_checks_passed = all(case.get("status") == "passed" for case in cases)
    return {
        "validator": validator_result,
        "graph": {
            "archive_path": INNER_VOLTERRA,
            "raw_sha256": hashes[INNER_VOLTERRA],
            "raw_sha256_pin": GRAPH_SHA256,
            "raw_sha256_verified_before_gzip_decode": True,
            "top_level_keys": sorted(str(key) for key in graph.keys()),
            "exact_graph_nodes_count": len(exact_graph_nodes),
            "selected_rows": {
                str(INTEGRAL_NODE_ID): exact_graph_nodes[INTEGRAL_NODE_ID],
                str(BINDER_NODE_ID): exact_graph_nodes[BINDER_NODE_ID],
                "35645": exact_graph_nodes[35645],
            },
        },
        "static_source_evidence": static_evidence,
        "graph_negative_cases": {"wrong_substitution_rejected": True, "wrong_upper_rejected": True},
        "synthetic_cases": cases,
        "synthetic_checks_passed": synthetic_checks_passed,
        "partial_progress_preserved": True,
        "shared_interval_context": {
            "same_iv_object_used_for_bridge_and_variable_z": True,
            "ancestor_producer_called": False,
            "historical_producer_imported": False,
            "historical_fixture_called": False,
        },
        "claims": {
            "candidate_witness_supported": True,
            "live_runtime_owner_proved": False,
            "archived_owner_same_proved": False,
            "cross_owner_binding_claimed": False,
            "certified_matrix_row": False,
            "numeric_or_physical_solution": False,
        },
    }


def _error_record(exc: BaseException) -> dict[str, str]:
    return {
        "type": type(exc).__name__,
        "message": str(exc),
        "traceback": traceback.format_exc(),
    }


def main() -> int:
    DIRECTORY.mkdir(parents=True, exist_ok=True)
    if ENTRY.exists() or RESULT.exists():
        print(json.dumps({"component": NAME, "status": "rejected", "error": "component already entered"}))
        return 1
    source_sha: str | None = None
    hashes: dict[str, str] = {}
    preflight_error: dict[str, str] | None = None
    try:
        _, hashes, source_sha = _preflight()
    except BaseException as exc:
        preflight_error = _error_record(exc)
    entry_value: dict[str, Any] = {
        "component": NAME,
        "source_sha256": source_sha,
        "input_hashes": hashes,
        "entered": True,
        "preflight_passed": preflight_error is None,
        "raw_sha256_verified_before_import_or_decode": preflight_error is None,
        "ancestor_producer_called": False,
        "existing_pure_consumer_imported": False,
        "historical_producer_imported": False,
        "historical_fixture_called": False,
        "manifest_created_by_runner": False,
    }
    if preflight_error is not None:
        entry_value["error"] = preflight_error
    try:
        entry_sha = _write_exclusive(ENTRY, entry_value)
    except BaseException as exc:
        print(json.dumps({"component": NAME, "status": "failed", "error": _error_record(exc)}))
        return 1
    progress: dict[str, Any] = {
        "raw_recheck_complete": False,
        "existing_pure_consumer_imported": False,
    }
    outcome: dict[str, Any] = {
        "component": NAME,
        "status": "failed",
        "entry_sha256": entry_sha,
        "source_sha256": source_sha,
        "input_hashes": hashes,
        "ancestor_producer_called": False,
        "existing_pure_consumer_imported": False,
        "historical_producer_imported": False,
        "historical_fixture_called": False,
        "raw_sha256_verified_before_import_or_decode": preflight_error is None,
        "progress": progress,
        "output": None,
    }
    if preflight_error is not None:
        outcome["error"] = preflight_error
    else:
        try:
            outcome["output"] = _consume(hashes, progress)
            outcome["status"] = "passed" if outcome["output"].get("synthetic_checks_passed") is True else "failed"
            if outcome["status"] != "passed":
                failed_cases = [
                    case for case in outcome["output"].get("synthetic_cases", [])
                    if case.get("status") != "passed"
                ]
                outcome["error"] = {
                    "type": "SyntheticValidationFailed",
                    "message": "one or more new synthetic cases failed",
                    "traceback": "\n\n".join(
                        case.get("error", {}).get("traceback", "")
                        for case in failed_cases
                        if isinstance(case.get("error"), Mapping)
                    ),
                }
        except BaseException as exc:
            outcome["error"] = _error_record(exc)
        finally:
            outcome["progress"] = progress
            outcome["existing_pure_consumer_imported"] = progress.get("existing_pure_consumer_imported") is True
    try:
        result_sha = _write_exclusive(RESULT, outcome)
    except BaseException as exc:
        print(json.dumps({"component": NAME, "status": "failed", "entry_sha256": entry_sha, "error": _error_record(exc)}))
        return 1
    print(
        json.dumps(
            {
                "component": NAME,
                "status": outcome["status"],
                "entry_sha256": entry_sha,
                "result_sha256": result_sha,
                "error": (outcome.get("error") or {}).get("message"),
            },
            sort_keys=True,
        )
    )
    return 0 if outcome["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
