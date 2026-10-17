from __future__ import annotations

"""Bounded source interface for the original m 15 full-minus-pulse graph.

This reads pinned receipts and source text only.  It does not import an
owner, execute an integral, replay a fixture, or turn an enclosure into a
function value.  The returned relation is a typed C 1 identity before
boxing; existing consumers still need a bridge to consume it.
"""

import ast
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

_PARENT_RELATIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_macro_coupled_m15_numeric.json"
)
_PARENT_SHA256 = "b63476443d5de5c1fca24ed5979ccb8a71273c42ff7a0b9b9094d8f80782f177"
_PREFIX = "lei_ren_part1_paper_compliant_current_original_"
_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}
_REPORTS = {
    "full_defect": _PREFIX + "n1_full_defect_linear_functions.json.gz",
    "full_bridge": _PREFIX + "n1_full_five_numeric_bridge.json.gz",
    "pulse_bridge": _PREFIX + "n1_pulse_numeric_bridge.json.gz",
    "inner": _PREFIX + "n1_inner_cutoff_functions.json.gz",
    "lower": _PREFIX + "n1_lower_collar_functions.json.gz",
    "preRp": _PREFIX + "n1_preRp_moment_functions.json.gz",
    "pulse_functions": _PREFIX + "n1_pulse_pressure_functions.json.gz",
}
_SOURCES = {
    "compose": _PREFIX + "n1_full_defect_linear_functions.py",
    "full_bridge": _PREFIX + "n1_full_five_numeric_bridge.py",
    "pulse": _PREFIX + "n1_pulse_pressure_functions.py",
    "preRp": _PREFIX + "n1_preRp_moment_functions.py",
    "lower": _PREFIX + "n1_lower_collar_functions.py",
    "inner": _PREFIX + "n1_inner_cutoff_functions.py",
}
_IMPLEMENTATIONS = {
    "reserved": "src/openai_ns_reconstruction/original_reserved_m15_response.py",
    "reserved_n1": "src/openai_ns_reconstruction/original_reserved_n1_numeric.py",
    "forward_native": "src/openai_ns_reconstruction/original_native_pulse_m15.py",
}
_FORWARD_PIN_REPORT = _PREFIX + "native_pulse_m15.json"
_FORWARD_PIN_REPORT_SHA256 = "311c61f1bc38143498ba443acac266b7160001037ef4bc022e52d9faf52940f0"


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _sha(path: Path) -> str | None:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError:
        return None


def _repo(root: Any, directory: Any) -> Path:
    for value in (root, directory):
        if value is None:
            continue
        try:
            path = Path(value)
        except TypeError:
            continue
        if path.is_file():
            path = path.parent
        if (path / "src").is_dir() and (path / "experiments").is_dir():
            return path
    return Path(__file__).resolve().parents[2]


def _directory(repo: Path, directory: Any) -> Path:
    if directory is None:
        return repo / "experiments" / "root_st073"
    path = Path(directory)
    if path.is_file():
        return path.parent
    if (path / "experiments" / "root_st073").is_dir():
        return path / "experiments" / "root_st073"
    return path


def _relative(path: Path, repo: Path) -> str:
    try:
        return path.relative_to(repo).as_posix()
    except ValueError:
        return path.as_posix()


def _load(path: Path) -> Mapping[str, Any]:
    value = json.loads(gzip.decompress(path.read_bytes()))
    if not isinstance(value, Mapping):
        raise ValueError(f"saved mapping required: {path}")
    return value


def _load_plain(path: Path) -> Mapping[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, Mapping):
        raise ValueError(f"saved mapping required: {path}")
    return value


def _digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _function(tree: ast.AST, name: str, cls: str | None = None) -> ast.AST:
    scope: ast.AST = tree
    if cls is not None:
        found = [n for n in ast.walk(tree) if isinstance(n, ast.ClassDef) and n.name == cls]
        if len(found) != 1:
            raise ValueError(f"unique class required: {cls}")
        scope = found[0]
    found = [n for n in ast.walk(scope) if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)) and n.name == name]
    if len(found) != 1:
        raise ValueError(f"unique function required: {name}")
    return found[0]


def _statement(text: str) -> ast.AST:
    tree = ast.parse(text)
    if len(tree.body) != 1:
        raise ValueError(text)
    return tree.body[0]


def _contains(node: ast.AST, wanted: ast.AST) -> bool:
    target = ast.dump(wanted, include_attributes=False)
    return any(ast.dump(part, include_attributes=False) == target for part in ast.walk(node))


def _ast_guard(
    path: Path,
    repo: Path,
    name: str,
    *,
    cls: str | None = None,
    statements: tuple[str, ...] = (),
    text: tuple[str, ...] = (),
) -> dict[str, Any]:
    row: dict[str, Any] = {
        "source_path": _relative(path, repo), "function": name, "class": cls,
        "required_AST": {s: False for s in statements}, "required_text": {s: False for s in text},
        "function_AST_sha256": None, "source_executed": False, "passed": False,
    }
    try:
        source = path.read_text(encoding="utf-8")
        fn = _function(ast.parse(source, filename=str(path)), name, cls)
        row["required_AST"] = {s: _contains(fn, _statement(s)) for s in statements}
        segment = ast.get_source_segment(source, fn) or ""
        compact = "".join(segment.split())
        row["required_text"] = {s: "".join(s.split()) in compact for s in text}
        row["function_AST_sha256"] = hashlib.sha256(ast.dump(fn, include_attributes=False).encode()).hexdigest()
        row["passed"] = all(row["required_AST"].values()) and all(row["required_text"].values())
    except (OSError, SyntaxError, UnicodeError, ValueError) as exc:
        row["error"] = str(exc)
    return row


def _pin(path: Path, repo: Path, expected: Any, name: str, kind: str) -> dict[str, Any]:
    actual = _sha(path)
    return {
        "name": name, "kind": kind,
        "key": path.name if kind == "input" else _relative(path, repo),
        "path": _relative(path, repo), "sha256": actual, "expected_sha256": expected,
        "parent_pin_present": isinstance(expected, str), "passed": isinstance(expected, str) and actual == expected,
        "source_executed": False,
    }


def _node(nodes: Any, index: int) -> Mapping[str, Any]:
    return nodes[index] if isinstance(nodes, list) and 0 <= index < len(nodes) and isinstance(nodes[index], Mapping) else {}


def _node_proofs(nodes: Any, expected: Mapping[int, Mapping[str, Any]]) -> dict[str, Any]:
    out = {}
    for index, wanted in expected.items():
        actual = dict(_node(nodes, index))
        out[str(index)] = {"expected": dict(wanted), "actual": actual, "passed": actual == dict(wanted)}
    return out


def _normalization_node(report: Mapping[str, Any], index: int) -> Mapping[str, Any]:
    span, rows = report.get("normalization_graph_span"), report.get("normalization_graph_nodes")
    if not isinstance(span, list) or len(span) != 2 or not isinstance(rows, list):
        return {}
    if not isinstance(span[0], int) or not span[0] <= index < span[1]:
        return {}
    row = index - span[0]
    return rows[row] if row < len(rows) and isinstance(rows[row], Mapping) else {}


def _pulse_endpoints(report: Mapping[str, Any]) -> dict[str, Any]:
    wanted = {
        "main": ("original_own_pulse_main_xi", "(R^2/mu)*dxi"),
        "gap": ("original_own_pulse_gap_xi", "(R^2/mu)*dxi"),
        "end": ("original_own_pulse_end_s", "R^2*ds"),
    }
    result: dict[str, Any] = {}
    for index, row in enumerate(report.get("exact_graph_nodes", ())):
        if not isinstance(row, Mapping) or row.get("operation") != "definite_integral":
            continue
        chart = row.get("original_complete_pulse_chart")
        if chart not in wanted or chart in result or row.get("ordinary_slow_Z_derivative_order") != 0:
            continue
        variable, measure = wanted[chart]
        if row.get("variable") == variable and row.get("measure") == measure:
            result[chart] = {"graph_node": index, "variable": variable, "lower": row.get("lower"),
                             "upper": row.get("upper"), "measure": measure, "ordinary_Z_derivative_order": 0}
    return result


def _empty(gaps: list[str], repo: Path, directory: Path, pin: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "passed": False, "compact_proofs_passed": False, "correlation_proved": False,
        "global_closure_claimed": False, "parent_report": _PARENT_RELATIVE,
        "parent_report_sha256": _PARENT_SHA256, "parent_report_pin": dict(pin), "source_family": dict(_FAMILY),
        "input_hashes": {}, "implementation_hashes": {}, "source_pins": {}, "AST": {},
        "ast_function_hashes": {}, "graph_identity": {}, "partition": {}, "normalized_expression": {},
        "source_function_bindings": {}, "source_generator_functions": {},
        "cached_bounds_refs": {}, "numeric_value_available": False,
        "numeric_joint_response_complete": False, "full_forward_partition_mapping_complete": False,
        "forward_partition_consumer_bound": False, "physical_scale_guards": {},
        "numeric_value_unavailable_reason": "producer bounds are not function values; no interval arithmetic performed",
        "consumer_dependencies": {}, "unresolved_consumer_center_angle_coupling": {}, "source_executed": False,
        "root": str(repo), "directory": str(directory), "gaps": list(gaps),
    }


def _parent_path(repo: Path, directory: Path, value: Any) -> Path:
    value = Path(value)
    if value.is_absolute():
        return value
    candidate = repo / value
    return candidate if candidate.exists() else directory / value


def _load_parent(value: Any, repo: Path, directory: Path) -> tuple[Mapping[str, Any], dict[str, Any], list[str]]:
    gaps: list[str] = []
    latest_path = repo / _PARENT_RELATIVE
    latest_sha = _sha(latest_path)
    pin: dict[str, Any] = {"path": _relative(latest_path, repo), "sha256": latest_sha,
        "expected_sha256": _PARENT_SHA256, "passed": latest_sha == _PARENT_SHA256, "source_executed": False}
    if not pin["passed"]:
        gaps.append("latest parent report pin mismatch")
    try:
        latest = _load_plain(latest_path)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        return {}, pin, gaps + [f"cannot load latest parent: {exc}"]
    try:
        supplied = value if isinstance(value, Mapping) else _load_plain(_parent_path(repo, directory, value))
    except (OSError, ValueError, TypeError, json.JSONDecodeError) as exc:
        return {}, pin, gaps + [f"cannot load supplied parent: {exc}"]
    supplied_digest, latest_digest = _digest(supplied), _digest(latest)
    pin.update({"supplied_mapping_sha256": supplied_digest, "latest_mapping_sha256": latest_digest,
                "mapping_match": supplied_digest == latest_digest})
    if supplied_digest != latest_digest:
        gaps.append("supplied parent mapping does not match pinned latest report")
    return supplied, pin, gaps


def build_joint_m15_source(parent_report: Any, *, root: Any, directory: Any) -> dict[str, Any]:
    """Export the exact full-minus-pulse m 15 source relation."""
    repo, work = _repo(root, directory), _directory(_repo(root, directory), directory)
    parent, parent_pin, gaps = _load_parent(parent_report, repo, work)
    if not parent:
        return _empty(gaps, repo, work, parent_pin)
    if _mapping(parent.get("source_family")) != _FAMILY:
        gaps.append("parent source_family mismatch")
    if _mapping(parent.get("checks")).get("passed") is not True:
        gaps.append("parent checks.passed is not true")
    parent_inputs, parent_impl = _mapping(parent.get("input_hashes")), _mapping(parent.get("implementation_hashes"))
    source_pins: dict[str, Any] = {}
    input_hashes: dict[str, str] = {}
    implementation_hashes: dict[str, str] = {}
    paths: dict[str, Path] = {}
    forward_pin_path = work / _FORWARD_PIN_REPORT
    forward_pin: Mapping[str, Any] = {}
    try:
        forward_pin = _load_plain(forward_pin_path)
        forward_report_sha = _sha(forward_pin_path)
        forward_report_ok = forward_report_sha == _FORWARD_PIN_REPORT_SHA256
        if not forward_report_ok:
            gaps.append("forward native pin report SHA mismatch")
        if _mapping(forward_pin.get("source_family")) != _FAMILY or _mapping(forward_pin.get("checks")).get("passed") is not True:
            gaps.append("forward native pin report family/check failed")
        source_pins["forward_native_pin_report"] = {
            "key": forward_pin_path.name, "path": _relative(forward_pin_path, repo),
            "sha256": forward_report_sha, "expected_sha256": _FORWARD_PIN_REPORT_SHA256,
            "sha_passed": forward_report_ok,
            "source_family": _mapping(forward_pin.get("source_family")),
            "checks_passed": _mapping(forward_pin.get("checks")).get("passed") is True,
            "passed": forward_report_ok and _mapping(forward_pin.get("source_family")) == _FAMILY and _mapping(forward_pin.get("checks")).get("passed") is True,
            "source_executed": False,
        }
        if source_pins["forward_native_pin_report"]["passed"]:
            input_hashes[forward_pin_path.name] = forward_report_sha
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        gaps.append(f"cannot load forward native pin report: {exc}")
    for key, basename in _REPORTS.items():
        path = work / basename
        row = _pin(path, repo, parent_inputs.get(path.name), key + "_report", "input")
        source_pins[key + "_report"], paths[key] = row, path
        if row["passed"]:
            input_hashes[path.name] = row["sha256"]
        else:
            gaps.append("report pin failed: " + path.name)
    for key, basename in _SOURCES.items():
        path = work / basename
        row = _pin(path, repo, parent_inputs.get(path.name), key + "_source", "input")
        source_pins[key + "_source"], paths[key + "_source"] = row, path
        if row["passed"]:
            input_hashes[path.name] = row["sha256"]
        else:
            gaps.append("source pin failed: " + path.name)
    for key, relative in _IMPLEMENTATIONS.items():
        path = repo / relative
        expected = parent_impl.get(relative)
        if expected is None and key == "forward_native":
            expected = _mapping(forward_pin.get("implementation_hashes")).get(relative)
        row = _pin(path, repo, expected, key + "_consumer", "implementation")
        if key == "forward_native":
            row["pin_source"] = "native_pulse_m15.json implementation_hashes" if relative not in parent_impl else "parent implementation_hashes"
        source_pins[key + "_consumer"], paths[key + "_consumer"] = row, path
        if row["passed"]:
            implementation_hashes[relative] = row["sha256"]
        else:
            gaps.append("consumer pin failed: " + relative)

    reports: dict[str, Mapping[str, Any]] = {}
    for key in _REPORTS:
        try:
            reports[key] = _load(paths[key])
        except (OSError, ValueError, json.JSONDecodeError, gzip.BadGzipFile) as exc:
            reports[key] = {}
            gaps.append(f"cannot load {key}: {exc}")
    family_rows = {key: _mapping(value.get("source_family")) == _FAMILY for key, value in reports.items()}
    if not all(family_rows.values()):
        gaps.append("producer source_family mismatch")

    full_nodes = reports["full_defect"].get("exact_graph_nodes", [])
    expected = {
        53161: {"operation": "function_substitution", "expression": 53141, "variable": 35848, "value": 35824, "Z_independent_substitution": True},
        53162: {"operation": "function_substitution", "expression": 53142, "variable": 35848, "value": 35824, "Z_independent_substitution": True},
        53171: {"operation": "sum", "arguments": [53161, 36592]},
        53172: {"operation": "sum", "arguments": [53162, 36593]},
        53173: {"operation": "sum", "arguments": [53171, 31221]},
        53174: {"operation": "sum", "arguments": [53172, 31220]},
        53175: {"operation": "sum", "arguments": [53173, 35556]},
        53176: {"operation": "sum", "arguments": [53174, 35557]},
    }
    node_proofs = _node_proofs(full_nodes, expected)
    graph_passed = all(row["passed"] for row in node_proofs.values())
    full_bind = _mapping(reports["full_defect"].get("full_defect_linear_bindings"))
    segments, moments = _mapping(full_bind.get("segments")), _mapping(full_bind.get("moments"))
    binding_checks = {
        "inner_m15": _mapping(segments.get("inner")).get("m15") == [53161, 53162],
        "lower_m15": _mapping(segments.get("lower")).get("m15") == [36592, 36593],
        "preRp_m15": _mapping(segments.get("preRp")).get("m15") == [31221, 31220],
        "pulse_m15": _mapping(segments.get("pulse")).get("m15") == [35556, 35557],
        "full_m15": moments.get("m15") == [53175, 53176],
    }
    inner_bind = _mapping(reports["inner"].get("inner_cutoff_bindings"))
    lower_bind = _mapping(reports["lower"].get("lower_collar_bindings"))
    pulse_bind = _mapping(reports["pulse_functions"].get("pulse_pressure_bindings"))
    component_checks = {
        "inner_endpoint": _mapping(inner_bind.get("endpoint")).get("m15", [])[:2] == [53161, 53162],
        "lower_function": _mapping(lower_bind.get("functions")).get("m15") == [36592, 36593],
        "preRp_function": _mapping(reports["preRp"].get("signed_pressure_C1_functions")).get("m15") == [31221, 31220],
        "pulse_function": _mapping(pulse_bind.get("functions")).get("m15") == [35556, 35557],
        "pulse_pieces": (
            _mapping(pulse_bind.get("pieces")).get("main", {}).get("m15") == [34708, 34709]
            and _mapping(pulse_bind.get("pieces")).get("gap", {}).get("m15") == [35007, 35008]
            and _mapping(pulse_bind.get("pieces")).get("end", {}).get("m15") == [35352, 35353]
        ),
    }
    graph_passed = graph_passed and all(binding_checks.values()) and all(component_checks.values())
    endpoints = _pulse_endpoints(reports["pulse_functions"])
    endpoints_passed = set(endpoints) == {"main", "gap", "end"}
    if not endpoints_passed:
        gaps.append("pulse chart endpoint graph incomplete")
    lower_domain = _mapping(lower_bind.get("fields")).get("domain")
    partition = {
        "inner": {"scope": "0 to Rin", "C1_m15": [53161, 53162], "cutoff_endpoint_node": 35824},
        "lower": {"scope": "Rin to r_minus", "C1_m15": [36592, 36593], "domain_nodes": lower_domain},
        "preRp": {"scope": "r_minus to Rp", "C1_m15": [31221, 31220], "scope_receipt": reports["preRp"].get("radial_scope")},
        "pulse": {"scope": "Rp to Rv", "C1_m15": [35556, 35557], "charts": endpoints},
    }

    ast_rows = {
        "compose_full_moments": _ast_guard(
            paths["compose_source"], repo, "compose_full_moments",
            statements=("moments[key]=g.c1add(moments[key],segment[key])",),
            text=("for segment in (lower,preRp,pulse)",),
        ),
        "pulse_pressure_functions": _ast_guard(
            paths["pulse_source"], repo, "pulse_pressure_functions",
            statements=(
                "pieces[chart]=parent.pressure_functions(g,delta,integrals[chart],*ends)",
                "functions=parent.pressure_functions(g,delta,total,*outer)",
            ),
            text=("for chart in ('main','gap','end')", "controls.integral", "physical_Pstar_not_reapplied=True"),
        ),
        "preRp_pressure_functions": _ast_guard(
            paths["preRp_source"], repo, "pressure_functions",
            statements=("return dict(m13=a.add(a.neg(half(S0)),boundary('q_plus_q_y')),m15=a.sub(half(S1),boundary('R_q_y')))",),
        ),
        "lower_direct_pressure_functions": _ast_guard(
            paths["lower_source"], repo, "direct_pressure_functions",
            text=("weighted_physical_measure", "original_m15_minus_applied_once=True", "signed_integrand_not_a_cap=True"),
        ),
        "inner_own_cutoff_fields": _ast_guard(
            paths["inner_source"], repo, "own_cutoff_fields",
            text=("actual_signed_n1_inner_moment=True", "full_core_to_first_history_retained=True"),
        ),
        "full_bridge_attachment": _ast_guard(
            paths["full_bridge_source"], repo, "__init__", cls="CurrentOriginalN1FullFiveNumericBridge",
            statements=("self.functions=self.full.moments",),
            text=("same positive reference*A*Rv power, Z independent",),
        ),
    }
    ast_passed = all(row.get("passed") is True for row in ast_rows.values())
    if not ast_passed:
        gaps.append("producer AST source guard failed")

    full_true, pulse_true = _mapping(reports["full_bridge"].get("true_function_bindings")), _mapping(reports["pulse_bridge"].get("true_function_bindings"))
    binding_checks.update({
        "full_true_physical": _mapping(full_true.get("physical_five")).get("m15") == [53175, 53176],
        "full_true_normalized": _mapping(full_true.get("normalized_pressure")).get("m15") == [53441, 53443],
        "pulse_true_physical": _mapping(pulse_true.get("physical")).get("m15") == [35556, 35557],
        "pulse_true_normalized": _mapping(pulse_true.get("normalized")).get("m15") == [53428, 53431],
    })
    graph_passed = graph_passed and all(binding_checks.values())

    pulse_norm = {str(i): dict(_normalization_node(reports["pulse_bridge"], i)) for i in (53399, 53421, 53427, 53428, 53429, 53430, 53431)}
    full_norm = {str(i): dict(_normalization_node(reports["full_bridge"], i)) for i in (53441, 53442, 53443)}
    unit_checks = {
        "u0_substitution": pulse_norm["53399"] == {"operation": "function_substitution", "expression": 31393, "variable": 20, "value": 0, "Z_independent_substitution": True},
        "unit_log_node": pulse_norm["53421"].get("operation") == "sum" and pulse_norm["53421"].get("arguments") == [53409, 5, 53401, 31320, 53418, 53420, 31343],
        "unit_exp_node": pulse_norm["53427"] == {"operation": "analytic_unary", "name": "exp", "argument": 53421},
        "unit_Z_exact_zero": pulse_norm["53399"].get("Z_independent_substitution") is True,
    }
    normalization_checks = {
        "shared_unit_log": pulse_norm["53427"].get("operation") == "analytic_unary" and pulse_norm["53427"].get("name") == "exp",
        "full_C0_quotient": full_norm["53441"] == {"operation": "positive_quotient", "numerator": 53175, "denominator": 53427, "source_positive_certificate": "same positive reference*A*Rv power, Z independent"},
        "full_C1_product": full_norm["53442"] == {"operation": "product", "arguments": [53176, 53427]},
        "full_C1_quotient": full_norm["53443"] == {"operation": "positive_quotient", "numerator": 53442, "denominator": 53430, "source_positive_certificate": "same positive reference*A*Rv power, Z independent squared"},
        "pulse_C0_quotient": pulse_norm["53428"] == {"operation": "positive_quotient", "numerator": 35556, "denominator": 53427, "source_positive_certificate": "positive reference*A*Rv power; all factors Z independent"},
        "pulse_C1_product": pulse_norm["53429"] == {"operation": "product", "arguments": [35557, 53427]},
        "shared_unit_square": pulse_norm["53430"] == {"operation": "product", "arguments": [53427, 53427]},
        "pulse_C1_quotient": pulse_norm["53431"] == {"operation": "positive_quotient", "numerator": 53429, "denominator": 53430, "source_positive_certificate": "positive reference*A*Rv power; all factors Z independent squared"},
        "unit_Z_exact_zero": all(unit_checks.values()),
    }
    normalization_passed = all(normalization_checks.values())
    graph_passed = graph_passed and normalization_passed

    full_identity, pulse_identity = _mapping(reports["full_bridge"].get("source_identity")), _mapping(reports["pulse_bridge"].get("source_identity"))
    scale_guards = {
        "full_source_identity_passed": full_identity.get("passed") is True,
        "full_composed_partition": full_identity.get("exact_five_signed_function_composition_replayed") is True,
        "full_same_positive_A_Rv_units": full_identity.get("preRp_pulse_same_original_positive_reference_A_Rv_units") is True,
        "full_Z1_unit_exact_zero": full_identity.get("pressure_unit_Z1_exact_zero") is True,
        "full_signed_m15_once": full_identity.get("m15_minus_already_in_original_signed_components") is True,
        "pulse_source_identity_passed": pulse_identity.get("passed") is True,
        "pulse_physical_measure": pulse_identity.get("physical_measure_main_gap_R_over_mu_end_R") is True,
        "pulse_Pstar_and_m15_sign_once": pulse_identity.get("Pstar_applied_once_and_signed_m15_once") is True,
        "pulse_reference_positive_once": pulse_identity.get("reference_is_positive_and_used_once") is True,
        "common_normalization_denominators": normalization_checks["shared_unit_log"] and normalization_checks["full_C1_quotient"],
        "shared_unit_Z_derivative_zero": unit_checks["unit_Z_exact_zero"],
    }
    source_identity_passed = all(scale_guards.values())
    if not source_identity_passed:
        gaps.append("physical source scale/sign/measure guard failed")

    relation = {
        "ordinary_Z_orders": [0, 1],
        "physical": {
            "C0": {"operation": "negative_full_plus_pulse", "full": 53175, "pulse": 35556, "reduced": [53161, 36592, 31221]},
            "C1": {"operation": "negative_full_plus_pulse", "full": 53176, "pulse": 35557, "reduced": [53162, 36593, 31220]},
        },
        "normalized": {
            "C0": {"full_node": 53441, "pulse_node": 53428, "shared_denominator": 53427, "reduced": "-(53161+36592+31221)/53427"},
            "C1": {"full_node": 53443, "pulse_node": 53431, "shared_denominator_square": 53430, "reduced": "-(53162+36593+31220)*53427/53430"},
        },
        "shared_unit": {
            "log_node": 53421, "unit_node": 53427, "square_node": 53430,
            "Z_variable_node": 20, "U0_substitution_node": 53399,
            "U_Z_exact_zero": unit_checks["unit_Z_exact_zero"],
            "ordinary_C1_is_scaled_physical_derivative": unit_checks["unit_Z_exact_zero"],
            "quotient_derivative_keeps_U_Z_term": False,
        },
        "symbolic": "formal joint pressure source = -full + pulse = -(inner + lower + preRp)",
        "relation_typed_exact_before_boxing": True,
        "same_physical_inner_lower_preRp_generator": bool(graph_passed and source_identity_passed),
        "correct_derivative_quotient": normalization_checks["full_C1_quotient"],
    }

    whole = _mapping(_mapping(reports["full_bridge"].get("samples")).get("whole_Z"))
    components = whole.get("pressure_component_enclosures", {}).get("m15") if isinstance(whole.get("pressure_component_enclosures"), Mapping) else None
    cached_refs = {
        "full_bridge": paths["full_bridge"].name + ":samples.whole_Z.pressure_component_enclosures.m15",
        "component_names": list(components) if isinstance(components, Mapping) else [],
        "component_C1_rows_present": {
            key: isinstance(value, Mapping) and isinstance(value.get("C1_enclosures"), list) and len(value["C1_enclosures"]) == 2
            for key, value in (components.items() if isinstance(components, Mapping) else ())
        },
        "pulse_box": paths["pulse_bridge"].name + ":samples.whole_Z.normalized_true_pulse_C1_enclosures.m15",
        "bounds_are_function_values": False,
        "source_enclosure_midpoint_used": False,
    }

    consumer_ast = {
        "reserved_response": _ast_guard(
            paths["reserved_consumer"], repo, "evaluate", cls="OriginalReservedM15Response",
            statements=("inputs=[None,None,[self.decode(v) for v in pressure]]",),
            text=("response[i]*inputs[i][k]", "normalized_full_pressure_moment_C1_enclosures"),
        ),
        "forward_inlet": _ast_guard(
            paths["forward_native_consumer"], repo, "inlet",
            cls="OriginalNativePulseM15Numeric",
            statements=("value=-self.J(self.iv,[self.native._decode(v) for v in rows])",),
            text=("normalized_true_pulse_C1_enclosures",),
        ),
    }
    ast_rows.update(consumer_ast)
    ast_passed = all(row.get("passed") is True for row in ast_rows.values())
    if not all(row.get("passed") is True for row in consumer_ast.values()):
        gaps.append("consumer AST source guard failed")
    consumer_flags = {
        "reserved_typed_joint_input_consumed": False,
        "reserved_partial_response_input_correlations_resolved": False,
        "reserved_signed_source_function_centers_resolved": False,
        "reserved_angle_center_coupling_resolved": False,
        "forward_typed_joint_input_consumed": False,
        "forward_relative_pressure_cancellation_accuracy": False,
        "forward_global_five_moment_numeric_join": False,
        "forward_angle_center_coupling_resolved": False,
    }
    consumer_dependencies = {
        "reserved_response": {
            "current_input": "full_bridge.samples.whole_Z.normalized_full_pressure_moment_C1_enclosures.m15",
            "typed_input_exported": True, "typed_input_consumed": False,
            "remaining_dependency": "replace the independent full-pressure box in evaluate with this relation while retaining response/angle source correlation",
            "AST": consumer_ast["reserved_response"],
        },
        "forward_m15": {
            "current_input": "pulse_bridge.samples.whole_Z.normalized_true_pulse_C1_enclosures.m15",
            "typed_input_exported": True, "typed_input_consumed": False,
            "remaining_dependency": "bind inlet and forward chart arithmetic to the same full/pulse C 1 source bundle before cancellation",
            "AST": consumer_ast["forward_inlet"],
        },
    }
    physical_units = _mapping(pulse_identity.get("original_physical_units"))
    physical_scale_guards = {
        "m15_source_units": physical_units.get("m15"),
        "source_A_formula": physical_units.get("source_A"),
        "Pstar_applied_once": pulse_identity.get("Pstar_applied_once_and_signed_m15_once") is True,
        "physical_measure": pulse_identity.get("physical_measure_main_gap_R_over_mu_end_R") is True,
        "signed_m15_once": full_identity.get("m15_minus_already_in_original_signed_components") is True,
    }
    source_function_bindings = {
        "inner": [53161, 53162], "lower": [36592, 36593], "preRp": [31221, 31220],
        "pulse": [35556, 35557], "full": [53175, 53176],
    }
    source_generator_functions = {
        "full_composition": {"function": "compose_full_moments", "source_pin": source_pins["compose_source"], "AST": ast_rows["compose_full_moments"]},
        "pulse": {"function": "pulse_pressure_functions", "source_pin": source_pins["pulse_source"], "AST": ast_rows["pulse_pressure_functions"]},
        "preRp": {"function": "pressure_functions", "source_pin": source_pins["preRp_source"], "AST": ast_rows["preRp_pressure_functions"]},
        "lower": {"function": "direct_pressure_functions", "source_pin": source_pins["lower_source"], "AST": ast_rows["lower_direct_pressure_functions"]},
        "inner": {"function": "own_cutoff_fields", "source_pin": source_pins["inner_source"], "AST": ast_rows["inner_own_cutoff_fields"]},
    }

    # This module proves source identity, not a global repair or consumer
    # closure.  Consumer AST failures remain explicit unresolved dependencies.
    compact = bool(
        parent_pin.get("passed") and not gaps and all(family_rows.values()) and graph_passed
        and endpoints_passed and ast_passed and source_identity_passed
    )
    return {
        "passed": compact, "compact_proofs_passed": compact,
        "correlation_proved": bool(graph_passed and source_identity_passed),
        "global_closure_claimed": False,
        "parent_report": _PARENT_RELATIVE, "parent_report_sha256": _PARENT_SHA256,
        "parent_report_pin": parent_pin, "source_family": dict(_FAMILY),
        "input_hashes": input_hashes, "implementation_hashes": implementation_hashes,
        "source_pins": source_pins, "AST": ast_rows,
        "ast_function_hashes": {key: row.get("function_AST_sha256") for key, row in ast_rows.items()},
        "graph_identity": {
            "passed": graph_passed, "node_proofs": node_proofs,
            "binding_checks": binding_checks, "component_checks": component_checks,
            "normalization_checks": normalization_checks,
            "normalization_nodes": {"pulse": pulse_norm, "full": full_norm},
            "source_identity_guards": scale_guards,
            "full_is_actual_composed_integral_recipe": bool(graph_passed and ast_rows["compose_full_moments"]["passed"]),
            "free_bounded_atom_used": False,
        },
        "partition": partition, "partition_endpoints_passed": endpoints_passed,
        "source_function_bindings": source_function_bindings,
        "source_generator_functions": source_generator_functions,
        "normalized_expression": relation, "cached_bounds_refs": cached_refs,
        "physical_scale_guards": physical_scale_guards,
        "numeric_value_available": False,
        "numeric_joint_response_complete": False,
        "full_forward_partition_mapping_complete": False,
        "forward_partition_consumer_bound": False,
        "numeric_value_unavailable_reason": "typed source relation is proven, but cached C 1 enclosures are not function values and no joint interval arithmetic is performed",
        "consumer_dependencies": consumer_dependencies,
        "unresolved_consumer_center_angle_coupling": consumer_flags,
        "source_executed": False, "root": str(repo), "directory": str(work), "gaps": gaps,
    }


__all__ = ["build_joint_m15_source"]
