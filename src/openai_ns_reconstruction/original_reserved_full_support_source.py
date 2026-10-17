from __future__ import annotations

"""Bounded source binding for the reserved full-support m 15 response.

The helper is deliberately receipt and AST only.  It validates the saved
first-switch parent, the saved repaired response and bridge graphs, and the
current source methods which own the response.  It never imports an owner,
replays a producer, or turns a cached enclosure into a function value.
"""

import ast
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

_PARENT_RELATIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_first_switch_terminal_m15_numeric.json"
)
_PARENT_SHA256 = "1a3c4759f475b9eb615352ead770baa2b874f2e03485d2f4828eb0744ff36535"
_PREFIX = "lei_ren_part1_paper_compliant_current_original_"
_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}
_REPAIRED_RELATIVE = "experiments/root_st073/" + _PREFIX + "repaired_Rp_m15_numeric.json"
_REPAIRED_SHA256 = "bd9fb2fc23072c743eeeadd300041ba187aff908a3bff4a3399f7222b72e6871"
_ARCHIVE_RELATIVE = "experiments/root_st073/" + _PREFIX + "macro_endpoint_m15_numeric.json"
_ARCHIVE_SHA256 = "16eb888da56fad8e5696b5fd8400b827ff7634e2c8f524110267cfcc6f525b80"
_FORWARD_RELATIVE = "experiments/root_st073/" + _PREFIX + "compact_quiet_m15_numeric.json"
_FORWARD_SHA256 = "7ccb30dbb3e440a3e598eb30deaf03509cf8624e429478aeeb364671baaac071"

_REPORTS = {
    "operator": _PREFIX + "n1_reserved_linear_operator.json.gz",
    "operator_check": _PREFIX + "n1_reserved_linear_operator_check.json",
    "coefficients": _PREFIX + "n1_reserved_coefficient_enclosures.json.gz",
    "coefficients_check": _PREFIX + "n1_reserved_coefficient_enclosures_check.json",
    "terminal": _PREFIX + "n1_terminal_integral_identities.json.gz",
    "terminal_check": _PREFIX + "n1_terminal_integral_identities_check.json",
    "full_bridge": _PREFIX + "n1_full_five_numeric_bridge.json.gz",
    "full_bridge_check": _PREFIX + "n1_full_five_numeric_bridge_check.json",
    "pulse_bridge": _PREFIX + "n1_pulse_numeric_bridge.json.gz",
    "pulse_bridge_check": _PREFIX + "n1_pulse_numeric_bridge_check.json",
    "repaired": _PREFIX + "repaired_Rp_m15_numeric.json",
    "archive": _PREFIX + "macro_endpoint_m15_numeric.json",
    "forward_checkpoint": _PREFIX + "compact_quiet_m15_numeric.json",
}
_IMPLEMENTATIONS = {
    "response": "src/openai_ns_reconstruction/original_reserved_m15_response.py",
    "repaired": "src/openai_ns_reconstruction/original_repaired_Rp_m15.py",
    "reserved_n1": "src/openai_ns_reconstruction/original_reserved_n1_numeric.py",
    "native_pulse": "src/openai_ns_reconstruction/original_native_pulse_m15.py",
}
_REQUIRED_SOURCE = {
    "response": "src/openai_ns_reconstruction/original_reserved_m15_response.py",
    "repaired": "src/openai_ns_reconstruction/original_repaired_Rp_m15.py",
}
_PHYSICAL = {
    "inner": [53161, 53162],
    "lower": [36592, 36593],
    "preRp": [31221, 31220],
    "pulse": [35556, 35557],
    "full": [53175, 53176],
}
_NORMALIZED = {"full": [53441, 53443], "pulse": [53428, 53431]}
_ORDERS = [0, 1]


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


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
    raw = gzip.decompress(path.read_bytes()) if path.suffix == ".gz" else path.read_bytes()
    value = json.loads(raw)
    if not isinstance(value, Mapping):
        raise ValueError(f"saved mapping required: {path}")
    return value


def _find(repo: Path, work: Path, key: str) -> Path:
    raw = Path(key)
    candidates = [raw] if raw.is_absolute() else [work / raw, repo / raw, repo / "experiments" / "root_st073" / raw]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def _function(tree: ast.AST, name: str, cls: str | None = None) -> ast.AST:
    scope: ast.AST = tree
    if cls is not None:
        classes = [node for node in ast.walk(tree) if isinstance(node, ast.ClassDef) and node.name == cls]
        if len(classes) != 1:
            raise ValueError(f"unique class required: {cls}")
        scope = classes[0]
    functions = [node for node in ast.walk(scope) if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name]
    if len(functions) != 1:
        raise ValueError(f"unique function required: {name}")
    return functions[0]


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
        "source_path": _relative(path, repo),
        "function": name,
        "class": cls,
        "required_AST": {value: False for value in statements},
        "required_text": {value: False for value in text},
        "function_AST_sha256": None,
        "source_executed": False,
        "passed": False,
    }
    try:
        source = path.read_text(encoding="utf-8")
        function = _function(ast.parse(source, filename=str(path)), name, cls)
        row["required_AST"] = {value: _contains(function, _statement(value)) for value in statements}
        segment = ast.get_source_segment(source, function) or ""
        compact = "".join(segment.split())
        row["required_text"] = {value: "".join(value.split()) in compact for value in text}
        row["function_AST_sha256"] = hashlib.sha256(ast.dump(function, include_attributes=False).encode()).hexdigest()
        row["passed"] = all(row["required_AST"].values()) and all(row["required_text"].values())
    except (OSError, SyntaxError, UnicodeError, ValueError) as exc:
        row["error"] = str(exc)
    return row


def _pin(path: Path, repo: Path, expected: Any, name: str, kind: str, *, anchor: bool = False) -> dict[str, Any]:
    actual = _sha(path)
    expected_ok = isinstance(expected, str) and actual == expected
    return {
        "name": name,
        "kind": kind,
        "key": path.name if kind == "input" else _relative(path, repo),
        "path": _relative(path, repo),
        "sha256": actual,
        "expected_sha256": expected,
        "parent_pin_present": isinstance(expected, str),
        "anchor_source": anchor,
        "passed": expected_ok,
        "source_executed": False,
    }


def _node(nodes: Any, index: int) -> Mapping[str, Any]:
    if isinstance(nodes, list) and 0 <= index < len(nodes) and isinstance(nodes[index], Mapping):
        return nodes[index]
    return {}


def _normalization_node(report: Mapping[str, Any], index: int) -> Mapping[str, Any]:
    span, nodes = report.get("normalization_graph_span"), report.get("normalization_graph_nodes")
    if not isinstance(span, list) or len(span) != 2 or not isinstance(nodes, list):
        return {}
    if not isinstance(span[0], int) or not span[0] <= index < span[1]:
        return {}
    position = index - span[0]
    return nodes[position] if 0 <= position < len(nodes) and isinstance(nodes[position], Mapping) else {}


def _empty(gaps: list[str], repo: Path, work: Path, parent_pin: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "passed": False,
        "compact_proofs_passed": False,
        "correlation_proved": False,
        "global_closure_claimed": False,
        "parent_report": _PARENT_RELATIVE,
        "parent_report_sha256": _PARENT_SHA256,
        "parent_report_pin": dict(parent_pin),
        "source_family": dict(_FAMILY),
        "input_hashes": {},
        "implementation_hashes": {},
        "added_pins": {},
        "transitive_hash_validation": {"passed": False},
        "source_pins": {},
        "AST": {},
        "ast_function_hashes": {},
        "graph_identity": {},
        "source_full_support_identity": {},
        "response_minus_target": {},
        "archive_provenance": {},
        "normalized_expression": {},
        "cached_bounds_refs": {},
        "physical_scale_guards": {},
        "numeric_value_available": False,
        "numeric_full_closure": False,
        "source_executed": False,
        "root": str(repo),
        "directory": str(work),
        "gaps": list(gaps),
    }


def _parent_path(repo: Path, work: Path, value: Any) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    candidate = repo / path
    return candidate if candidate.exists() else work / path


def _load_parent(value: Any, repo: Path, work: Path) -> tuple[Mapping[str, Any], dict[str, Any], list[str]]:
    gaps: list[str] = []
    latest_path = repo / _PARENT_RELATIVE
    latest_sha = _sha(latest_path)
    parent_pin: dict[str, Any] = {
        "path": _relative(latest_path, repo),
        "sha256": latest_sha,
        "expected_sha256": _PARENT_SHA256,
        "passed": latest_sha == _PARENT_SHA256,
        "source_executed": False,
    }
    if not parent_pin["passed"]:
        gaps.append("latest parent report pin mismatch")
    try:
        latest = _load(latest_path)
    except (OSError, ValueError, json.JSONDecodeError, gzip.BadGzipFile) as exc:
        return {}, parent_pin, gaps + [f"cannot load latest parent: {exc}"]
    try:
        supplied = value if isinstance(value, Mapping) else _load(_parent_path(repo, work, value))
    except (OSError, ValueError, TypeError, json.JSONDecodeError, gzip.BadGzipFile) as exc:
        return {}, parent_pin, gaps + [f"cannot load supplied parent: {exc}"]
    latest_digest, supplied_digest = _digest(latest), _digest(supplied)
    parent_pin.update({"latest_mapping_sha256": latest_digest, "supplied_mapping_sha256": supplied_digest,
                       "mapping_match": latest_digest == supplied_digest})
    if latest_digest != supplied_digest:
        gaps.append("supplied parent mapping does not match pinned latest report")
    return supplied, parent_pin, gaps


def _verify_hash_map(
    hashes: Mapping[str, Any], repo: Path, work: Path, *, label: str, source_pins: dict[str, Any],
    selected: dict[str, str], gaps: list[str], added: dict[str, str], anchor_keys: Mapping[str, str] = (),
) -> tuple[int, int]:
    passed = total = 0
    anchors = dict(anchor_keys)
    for key, expected in hashes.items():
        if not isinstance(key, str):
            continue
        total += 1
        path = _find(repo, work, key)
        actual = _sha(path)
        is_anchor = key in anchors
        expected_for_pin = anchors.get(key, expected)
        ok = isinstance(expected_for_pin, str) and actual == expected_for_pin
        row = {
            "name": label + ":" + key,
            "kind": "transitive_input" if label == "input" else "transitive_implementation",
            "key": key,
            "path": _relative(path, repo),
            "sha256": actual,
            "expected_sha256": expected_for_pin,
            "parent_expected_sha256": expected,
            "anchor_source": is_anchor,
            "passed": ok,
            "source_executed": False,
        }
        source_pins["transitive:" + label + ":" + key] = row
        if ok:
            passed += 1
            selected[key] = actual
            if is_anchor and expected != actual:
                added[key] = actual
        else:
            gaps.append(f"{label} hash failed: {key}")
    return passed, total


def _check_saved_report(report: Mapping[str, Any], name: str) -> bool:
    checks = _mapping(report.get("checks"))
    if name in {"repaired", "archive", "forward_checkpoint"}:
        return checks.get("passed") is True
    return checks.get("passed") is True or checks.get("all_passed") is True or report.get("all_passed") is True


def _archive_provenance(parent: Mapping[str, Any], archive: Mapping[str, Any], repaired: Mapping[str, Any], gaps: list[str]) -> dict[str, Any]:
    views = parent.get("signed_m15_query_views")
    samples = archive.get("samples")
    repaired_samples = repaired.get("samples")
    inherited_packet = _mapping(parent.get("inherited_m15_packet"))
    inherited_samples = inherited_packet.get("samples")
    views_ok = isinstance(views, list) and len(views) == 3
    sample_rows: list[dict[str, Any]] = []
    if not views_ok or not isinstance(samples, list) or len(samples) != 3 or not isinstance(repaired_samples, list) or len(repaired_samples) != 3 or not isinstance(inherited_samples, list) or len(inherited_samples) != 3:
        gaps.append("saved signed/archive/repaired sample count is not three")
        return {"passed": False, "signed_view_count": len(views) if isinstance(views, list) else None}
    for index, (view, sample, repaired_sample, inherited_sample) in enumerate(zip(views, samples, repaired_samples, inherited_samples)):
        view_map, sample_map, repaired_map = _mapping(view), _mapping(sample), _mapping(repaired_sample)
        inherited_map = _mapping(inherited_sample)
        reduced = sample_map.get("correlated_reduced_m15_C0_C1_terms")
        raw = sample_map.get("raw_m15_C0_C1_terms")
        reserved = reduced[182] if isinstance(reduced, list) and len(reduced) > 182 else {}
        # The accepted parent stores the canonical packet digests.  Validate
        # those pins directly; the archive's current terms remain provenance
        # and are never treated as newly recomputed numerical values.
        view_digest = view_map.get("parent_original_reduced183_digest")
        reduced_digest_ok = view_digest == inherited_map.get("reduced_sha256")
        raw_digest_ok = view_map.get("raw213_digest") == inherited_map.get("raw_sha256")
        target = _mapping(sample_map.get("independent_target"))
        target_binding = _mapping(target.get("actual_function_binding"))
        target_ok = (
            target_binding.get("physical_handles") == _PHYSICAL["pulse"]
            and target_binding.get("normalized_handles") == _NORMALIZED["pulse"]
            and target_binding.get("sign") == -1
            and target_binding.get("physical_measure") == "-R*pressure_density*dR"
            and target_binding.get("ordinary_Z1") is True
        )
        response_saved = _mapping(repaired_map.get("reserved_response"))
        angular = response_saved.get("angular_inverse_response")
        diagnostics = _mapping(response_saved.get("diagnostics"))
        response_exact = (
            response_saved.get("native_compact_x") in (2, "2")
            and isinstance(angular, list) and len(angular) == 3
            and angular[0].get("source_enclosure_exact_mpf_tuples") == [[0, 0, 0, 0], [0, 0, 0, 0]]
            and angular[1].get("source_enclosure_exact_mpf_tuples") == [[0, 0, 0, 0], [0, 0, 0, 0]]
            and angular[2].get("source_enclosure_exact_mpf_tuples") == [[0, 1, 0, 1], [0, 1, 0, 1]]
            and diagnostics.get("exact_full") is True
        )
        own = repaired_map.get("m15_C0_C1_terms", [])
        own = own[-1] if isinstance(own, list) and own else {}
        repaired_binding = _mapping(_mapping(own).get("actual_function_binding"))
        archived_own = reduced[-1] if isinstance(reduced, list) and reduced else {}
        own_equal = archived_own == own
        entries = view_map.get("entries")
        parent_last = entries[-1] if isinstance(entries, list) and entries else {}
        parent_last_ok = (
            _mapping(parent_last).get("parent_reduced_index") == 182
            and _mapping(parent_last).get("name") == "reserved_n1_own_swirl_once"
            and _mapping(parent_last).get("operation") == "unchanged_parent_reference"
            and _mapping(parent_last).get("parent_term_digest") == _digest(own)
        )
        parent_target = _mapping(view_map.get("independent_target"))
        target_equal = target == _mapping(repaired_map.get("independent_target")) == parent_target
        own_ok = (
            response_exact
            and repaired_binding.get("terminal_patch") == [53905, 53906]
            and repaired_binding.get("reduced_response") == [53921, 53922]
            and repaired_binding.get("physical_measure") == "-R*theta0*Utheta1 ds"
        )
        row_ok = (
            isinstance(raw, list) and len(raw) == 213
            and isinstance(reduced, list) and len(reduced) == 183
            and _mapping(reserved).get("name") == "reserved_n1_own_swirl_once"
            and reduced_digest_ok and raw_digest_ok and target_ok and own_equal and parent_last_ok and target_equal
            and view_map.get("independent_target_digest") == inherited_map.get("independent_target_sha256")
            and own_ok
        )
        if not row_ok:
            gaps.append(f"archive provenance sample {index} failed")
        sample_rows.append({
            "index": index,
            "Z": sample_map.get("Z"),
            "raw_count": len(raw) if isinstance(raw, list) else None,
            "reduced_count": len(reduced) if isinstance(reduced, list) else None,
            "reserved_parent_index": 182,
            "reserved_term_digest": _digest(reserved) if isinstance(reserved, Mapping) else None,
            "view_parent_reduced183_digest": view_digest,
            "stored_reduced183_digest": inherited_map.get("reduced_sha256"),
            "parent_view_digest_match": reduced_digest_ok and raw_digest_ok,
            "archive_digest_recomputed": False,
            "target_binding": target_binding,
            "reserved_response_binding": repaired_binding,
            "response_exact_angular_identity": response_exact,
            "archive_repaired_own_term_equal": own_equal,
            "parent_last_reserved_entry_equal": parent_last_ok,
            "archived_repaired_parent_target_equal": target_equal,
            "passed": row_ok,
        })
    packet = _mapping(parent.get("inherited_m15_packet"))
    packet_ok = packet.get("report_sha256") == _ARCHIVE_SHA256 and packet.get("raw_term_count") in (None, 213)
    if not packet_ok:
        gaps.append("inherited m15 packet archive pin/count failed")
    return {
        "passed": bool(views_ok and all(row["passed"] for row in sample_rows) and packet_ok),
        "archive_report": _ARCHIVE_RELATIVE,
        "archive_report_sha256": _ARCHIVE_SHA256,
        "repaired_report": _REPAIRED_RELATIVE,
        "repaired_report_sha256": _REPAIRED_SHA256,
        "signed_view_count": len(views),
        "sample_rows": sample_rows,
        "parent_inherited_packet": packet,
        "retains_parent_reserved_index_182": True,
        "retains_independent_negative_pulse_target": True,
    }


def build_reserved_full_support_source(parent_report: Any, *, root: Any, directory: Any) -> dict[str, Any]:
    """Build a fail-closed, source-only full-support identity descriptor."""
    repo = _repo(root, directory)
    work = _directory(repo, directory)
    parent, parent_pin, gaps = _load_parent(parent_report, repo, work)
    if not parent:
        return _empty(gaps, repo, work, parent_pin)

    if _mapping(parent.get("source_family")) != _FAMILY:
        gaps.append("parent source_family mismatch")
    if _mapping(parent.get("checks")).get("passed") is not True:
        gaps.append("parent checks.passed is not true")
    joint = _mapping(parent.get("joint_m15_pressure_source"))
    joint_relation = _mapping(joint.get("normalized_expression"))
    joint_ok = (
        joint.get("passed") is True
        and joint.get("correlation_proved") is True
        and joint.get("global_closure_claimed") is False
        and joint.get("joint_pressure_source_used_as_total_reserved_value") is not True
        and joint_relation.get("ordinary_Z_orders") == _ORDERS
        and joint_relation.get("relation_typed_exact_before_boxing") is True
        and joint_relation.get("symbolic") == "formal joint pressure source = -full + pulse = -(inner + lower + preRp)"
    )
    if not joint_ok:
        gaps.append("accepted joint_m15_pressure_source relation guard failed")

    parent_inputs, parent_impl = _mapping(parent.get("input_hashes")), _mapping(parent.get("implementation_hashes"))
    source_pins: dict[str, Any] = {}
    input_hashes: dict[str, str] = {}
    implementation_hashes: dict[str, str] = {}
    added_pins: dict[str, str] = {}
    reports: dict[str, Mapping[str, Any]] = {}
    paths: dict[str, Path] = {}
    for key, basename in _REPORTS.items():
        path = work / basename
        expected = parent_inputs.get(basename)
        anchor = key == "repaired" and expected is None
        if anchor:
            expected = _REPAIRED_SHA256
        if key == "archive" and expected is None:
            expected = _ARCHIVE_SHA256
            anchor = True
        row = _pin(path, repo, expected, key + "_report", "input", anchor=anchor)
        source_pins[key + "_report"] = row
        paths[key] = path
        if row["passed"]:
            input_hashes[basename] = row["sha256"]
            if anchor:
                added_pins[basename] = row["sha256"]
        else:
            gaps.append("report pin failed: " + basename)
        try:
            reports[key] = _load(path)
        except (OSError, ValueError, json.JSONDecodeError, gzip.BadGzipFile) as exc:
            reports[key] = {}
            gaps.append(f"cannot load {key}: {exc}")

    for key, relative in _IMPLEMENTATIONS.items():
        path = repo / relative
        expected = parent_impl.get(relative)
        row = _pin(path, repo, expected, key + "_implementation", "implementation")
        source_pins[key + "_implementation"] = row
        if row["passed"]:
            implementation_hashes[relative] = row["sha256"]
        else:
            gaps.append("implementation pin failed: " + relative)

    # The parent carries a complete input and implementation closure.  Check
    # every path here, including paths not used by this compact identity.
    transitive_input_passed, transitive_input_total = _verify_hash_map(
        parent_inputs, repo, work, label="input", source_pins=source_pins,
        selected=input_hashes, gaps=gaps, added=added_pins,
        anchor_keys={_REPAIRED_RELATIVE: _REPAIRED_SHA256, _ARCHIVE_RELATIVE: _ARCHIVE_SHA256},
    )
    transitive_impl_passed, transitive_impl_total = _verify_hash_map(
        parent_impl, repo, work, label="implementation", source_pins=source_pins,
        selected=implementation_hashes, gaps=gaps, added=added_pins,
    )
    transitive_ok = transitive_input_passed == transitive_input_total and transitive_impl_passed == transitive_impl_total

    family_rows = {key: _mapping(value.get("source_family")) == _FAMILY for key, value in reports.items()}
    for key, report in reports.items():
        if not family_rows[key]:
            gaps.append("saved report source_family mismatch: " + key)
        if key.endswith("_check") and report and not _check_saved_report(report, key):
            gaps.append("saved check report did not pass: " + key)
    if not _check_saved_report(reports["repaired"], "repaired"):
        gaps.append("repaired report checks.passed is not true")
    if not _check_saved_report(reports["archive"], "archive"):
        gaps.append("archive report checks.passed is not true")

    operator = reports["operator"]
    terminal = reports["terminal"]
    full = reports["full_bridge"]
    pulse = reports["pulse_bridge"]
    matrix = _mapping(operator.get("exact_normalized_basis_and_matrix_functions"))
    matrix_checks = {
        "pressure_row": matrix.get("angular_rows", [None, None, None])[2] == [25072, 25075, 25078],
        "inverse_pressure_column": [
            _mapping(matrix.get("inverse_angular", [None, None, None])[row]).get("2")
            if isinstance(matrix.get("inverse_angular", [None, None, None])[row], Mapping) else matrix.get("inverse_angular", [None, None, None])[row][2]
            for row in range(3)
        ] == [25117, 25134, 25150],
        "pressure_row_nodes_present": all(_node(operator.get("exact_graph_nodes"), node) for node in [25072, 25075, 25078]),
        "inverse_column_nodes_present": all(_node(operator.get("exact_graph_nodes"), node) for node in [25117, 25134, 25150]),
    }
    if not all(matrix_checks.values()):
        gaps.append("reserved operator pressure row/inverse column guard failed")
    operator_check = reports["operator_check"]
    determinant_identities = _mapping(operator_check.get("determinant_identities"))
    operator_source_checks = {
        "check_all_passed": operator_check.get("all_passed") is True,
        "angular_exact_Vandermonde_identity": determinant_identities.get("angular_exact_Vandermonde_identity") is True,
        "angular_determinant_positive": determinant_identities.get("angular_determinant_positive_for_mu_nonnegative") is True,
        "actual_mu_positive": operator_check.get("actual_mu_positive_and_retained") is True,
        "source_functions_and_bounds_distinct": operator_check.get("source_functions_and_bounds_distinct") is True,
        "new_supports_strictly_disjoint": operator_check.get("new_supports_strictly_disjoint") is True,
    }
    if not all(operator_source_checks.values()):
        gaps.append("operator_check exact inverse/determinant source proof failed")

    terminal_bindings = _mapping(terminal.get("terminal_identity_bindings"))
    # Terminal identities are keyed by binding first, then by moment.
    terminal_m15 = terminal_bindings
    terminal_checks = {
        "patch_endpoint": _mapping(terminal_m15.get("patch_endpoint")).get("m15") == [53905, 53906],
        "reduced_patch": _mapping(terminal_m15.get("reduced_patch")).get("m15") == [53921, 53922],
        "terminal_sources_m15": _mapping(terminal_m15.get("terminal_sources")).get("m15") == [53832, 53833],
        "pressure_terminal_radius": "Rv" in str(_mapping(terminal_m15.get("definition")).get("pressure_moment_terminal_radius", "")),
        "velocity_terminal_radius": "2Rc" in str(_mapping(terminal_m15.get("definition")).get("velocity_moment_terminal_radius", "")),
    }
    if not all(terminal_checks.values()):
        gaps.append("terminal m15 patch/tail guard failed")

    full_true, pulse_true = _mapping(full.get("true_function_bindings")), _mapping(pulse.get("true_function_bindings"))
    physical_checks = {
        "inner": _PHYSICAL["inner"] == _mapping(joint.get("source_function_bindings")).get("inner"),
        "lower": _PHYSICAL["lower"] == _mapping(joint.get("source_function_bindings")).get("lower"),
        "preRp": _PHYSICAL["preRp"] == _mapping(joint.get("source_function_bindings")).get("preRp"),
        "full": _mapping(full_true.get("physical_five")).get("m15") == _PHYSICAL["full"],
        "pulse": _mapping(pulse_true.get("physical")).get("m15") == _PHYSICAL["pulse"],
        "full_normalized": _mapping(full_true.get("normalized_pressure")).get("m15") == _NORMALIZED["full"],
        "pulse_normalized": _mapping(pulse_true.get("normalized")).get("m15") == _NORMALIZED["pulse"],
    }
    if not all(physical_checks.values()):
        gaps.append("saved physical/normalized m15 handle guard failed")
    coefficient_identity = _mapping(_mapping(reports["coefficients"].get("source_identity")).get("actual_arc_source_identity"))
    coefficient_arc_checks = {
        "coefficient_source_passed": _mapping(reports["coefficients"].get("source_identity")).get("passed") is True,
        "selected_five_graph_action": _mapping(reports["coefficients"].get("source_identity")).get("actual_selected_five_coefficient_graph_action_replayed") is True,
        "actual_arc_source_passed": coefficient_identity.get("passed") is True,
        "canonical_Rd_amplitude": coefficient_identity.get("actual_canonical_Rd_amplitude_body_bound") is True,
        "J_one_half": coefficient_identity.get("actual_J_one_half") is True,
        "finite_mu_factors": coefficient_identity.get("actual_nonzero_finite_mu_factor_retained") is True,
        "Astar_Z_retained": _mapping(reports["coefficients"].get("source_identity")).get("actual_Astar_Z_retained") is True,
        "quotient_formula": coefficient_identity.get("exact_finite_coefficient") == "exp(-17/10-5*mu/2)/(1+Z^2)",
    }
    if not all(coefficient_arc_checks.values()):
        gaps.append("coefficient actual A_rc/theta2 source identity guard failed")

    full_nodes = {str(index): dict(_normalization_node(full, index)) for index in [53441, 53442, 53443]}
    pulse_nodes = {str(index): dict(_normalization_node(pulse, index)) for index in [53399, 53421, 53427, 53428, 53429, 53430, 53431]}
    normalization_checks = {
        "full_C0": full_nodes["53441"] == {"operation": "positive_quotient", "numerator": 53175, "denominator": 53427, "source_positive_certificate": "same positive reference*A*Rv power, Z independent"},
        "full_C1_product": full_nodes["53442"] == {"operation": "product", "arguments": [53176, 53427]},
        "full_C1": full_nodes["53443"] == {"operation": "positive_quotient", "numerator": 53442, "denominator": 53430, "source_positive_certificate": "same positive reference*A*Rv power, Z independent squared"},
        "pulse_C0": pulse_nodes["53428"] == {"operation": "positive_quotient", "numerator": 35556, "denominator": 53427, "source_positive_certificate": "positive reference*A*Rv power; all factors Z independent"},
        "pulse_C1_product": pulse_nodes["53429"] == {"operation": "product", "arguments": [35557, 53427]},
        "shared_unit_square": pulse_nodes["53430"] == {"operation": "product", "arguments": [53427, 53427]},
        "pulse_C1": pulse_nodes["53431"] == {"operation": "positive_quotient", "numerator": 53429, "denominator": 53430, "source_positive_certificate": "positive reference*A*Rv power; all factors Z independent squared"},
        "U_Z_exact_zero": pulse_nodes["53399"] == {"operation": "function_substitution", "expression": 31393, "variable": 20, "value": 0, "Z_independent_substitution": True},
    }
    if not all(normalization_checks.values()):
        gaps.append("normalized full/pulse graph identity guard failed")

    full_identity, pulse_identity = _mapping(full.get("source_identity")), _mapping(pulse.get("source_identity"))
    units = _mapping(pulse_identity.get("original_physical_units"))
    scale_guards = {
        "full_source_identity": full_identity.get("passed") is True,
        "full_partition_complete": full_identity.get("inner_lower_preRp_pulse_pressure_partitions_complete") is True,
        "full_same_A_Rv_unit": full_identity.get("preRp_pulse_same_original_positive_reference_A_Rv_units") is True,
        "full_Z1_unit_zero": full_identity.get("pressure_unit_Z1_exact_zero") is True,
        "pulse_source_identity": bool(pulse_identity.get("passed") is True or pulse_identity.get("exact_native_numeric_source_identity")),
        "pulse_measure": pulse_identity.get("physical_measure_main_gap_R_over_mu_end_R") is True,
        "pulse_sign_once": pulse_identity.get("Pstar_applied_once_and_signed_m15_once") is True,
        "positive_reference_once": pulse_identity.get("reference_is_positive_and_used_once") is True,
        "m15_unit": units.get("m15") == {"A_power": 1, "Rv_power": 2, "inverse_mu_numerator": 37, "inverse_mu_denominator": 2},
        "pulse_radial_tail": units.get("radial_interval") == "Rp<=R<=Rv" and units.get("terminal_M_at_Rv_exact_zero") is True,
        "actual_arc_source_identity": all(coefficient_arc_checks.values()),
    }
    if not all(scale_guards.values()):
        gaps.append("physical measure/unit/tail guard failed")

    response_path, repaired_path = repo / _REQUIRED_SOURCE["response"], repo / _REQUIRED_SOURCE["repaired"]
    ast_rows = {
        "inverse3_cofactors": _ast_guard(
            response_path, repo, "inverse3",
            statements=("row.append((-1)**(i+j)*(matrix[a][c]*matrix[b][d]-matrix[a][d]*matrix[b][c])/determinant)",),
            text=("for j in range(3)", "/determinant"),
        ),
        "reserved_init_inverse_binding": _ast_guard(
            response_path, repo, "__init__", cls="OriginalReservedM15Response",
            statements=("self.inverse=inverse3(self.matrix,det)",),
            text=("self.matrix=[[self.decode(v) for v in row] for row in values['angular_rows']]", "_bounds(self.ctx,det)[0]<=0"),
        ),
        "reserved_partial_row": _ast_guard(
            response_path, repo, "partial_row", cls="OriginalReservedM15Response",
            statements=("response.append(self.intersect(f,b))",),
            text=("if exact_empty:response=[iv.mpf(0)]*3", "elif exact_full:response=[iv.mpf(0),iv.mpf(0),iv.mpf(1)]", "1<=Fraction(x)<=2"),
        ),
        "reserved_evaluate": _ast_guard(
            response_path, repo, "evaluate", cls="OriginalReservedM15Response",
            statements=("response,diagnostic=self.partial_row(x)", "inputs=[None,None,[self.decode(v) for v in pressure]]", "result=[-sum((response[i]*inputs[i][k] for i in range(3)),iv.mpf(0)) for k in (0,1)]"),
            text=("normalized_full_pressure_moment_C1_enclosures", "-1<c.mpf(Z)<1", "constant=iv.exp(-iv.mpf(17)/10-iv.mpf(5)*self.mu/2)", "arc=[constant/(1+z*z),-2*z*constant/(1+z*z)**2]"),
        ),
        "repaired_init_saved_selection": _ast_guard(
            repaired_path, repo, "__init__", cls="OriginalRepairedRpM15",
            statements=("self.saved={Fraction(p['Z']):p for p in data['samples'] if p['chart']=='quiet' and Fraction(p['native_coordinate'])==1}",),
            text=("self.saved", "p['chart']=='quiet'", "Fraction(p['native_coordinate'])==1"),
        ),
        "repaired_evaluate": _ast_guard(
            repaired_path, repo, "evaluate", cls="OriginalRepairedRpM15",
            statements=("crosswalk=self.crosswalk(packet)", "response = self.response.evaluate(Z,'2')"),
            text=("reserved_n1_own_swirl_once", "-R*theta0*Utheta1 ds", "terminal_patch=self.bindings['patch_endpoint']", "reduced_response=self.bindings['reduced_patch']"),
        ),
    }
    ast_passed = all(row.get("passed") is True for row in ast_rows.values())
    if not ast_passed:
        gaps.append("reserved response/repaired AST guard failed")

    repaired_samples = reports["repaired"].get("samples")
    forward_samples = reports["forward_checkpoint"].get("samples")
    forward_selected = [
        {"Z": _mapping(sample).get("Z"), "chart": _mapping(sample).get("chart"), "native_coordinate": _mapping(sample).get("native_coordinate")}
        for sample in (forward_samples if isinstance(forward_samples, list) else ())
        if _mapping(sample).get("chart") == "quiet" and str(_mapping(sample).get("native_coordinate")) in {"1", "1.0"}
    ]
    repaired_Z = {str(_mapping(sample).get("Z")) for sample in (repaired_samples if isinstance(repaired_samples, list) else ())}
    forward_selection_ok = len(forward_selected) == 3 and {str(row.get("Z")) for row in forward_selected} == repaired_Z
    if not forward_selection_ok:
        gaps.append("pinned quiet q1/Rp forward selection is incomplete")
    response_identity_rows: list[dict[str, Any]] = []
    transport_rows: list[dict[str, Any]] = []
    if isinstance(repaired_samples, list):
        for sample in repaired_samples:
            sample_map = _mapping(sample)
            response_saved = _mapping(sample_map.get("reserved_response"))
            response_binding = _mapping(_mapping(sample_map.get("m15_C0_C1_terms", [{}])[-1] if isinstance(sample_map.get("m15_C0_C1_terms"), list) else {}).get("actual_function_binding"))
            response_contract = _mapping(response_saved.get("contract"))
            angular = response_saved.get("angular_inverse_response")
            diagnostics = _mapping(response_saved.get("diagnostics"))
            response_exact = (
                response_saved.get("native_compact_x") in (2, "2")
                and isinstance(angular, list) and len(angular) == 3
                and angular[0].get("source_enclosure_exact_mpf_tuples") == [[0, 0, 0, 0], [0, 0, 0, 0]]
                and angular[1].get("source_enclosure_exact_mpf_tuples") == [[0, 0, 0, 0], [0, 0, 0, 0]]
                and angular[2].get("source_enclosure_exact_mpf_tuples") == [[0, 1, 0, 1], [0, 1, 0, 1]]
                and diagnostics.get("exact_full") is True
            )
            crosswalk = sample_map.get("crosswalk")
            transport_ok = (
                response_saved.get("native_compact_x") in (2, "2")
                and isinstance(crosswalk, list) and bool(crosswalk)
                and response_binding.get("terminal_patch") == [53905, 53906]
                and response_binding.get("reduced_response") == [53921, 53922]
            )
            response_ok = (
                response_binding.get("physical_measure") == "-R*theta0*Utheta1 ds"
                and response_binding.get("terminal_patch") == [53905, 53906]
                and response_binding.get("reduced_response") == [53921, 53922]
                and response_contract.get("full_support_response_is_exact_matrix_row_identity") is True
                and response_contract.get("independent_complete_pulse_not_used_to_define_this_response") is True
                and response_exact
            )
            response_identity_rows.append({"Z": sample_map.get("Z"), "binding": response_binding, "contract": response_contract, "exact_angular_identity": response_exact, "passed": response_ok})
            transport_rows.append({"Z": sample_map.get("Z"), "compact_x": response_saved.get("native_compact_x"), "passed": transport_ok,
                                   "description": "compact x 2 response transported to quiet q1/Rp; physical pressure remains through Rv"})
    response_source_ok = bool(response_identity_rows) and all(row["passed"] for row in response_identity_rows)
    transport_ok = forward_selection_ok and bool(transport_rows) and all(row["passed"] for row in transport_rows) and ast_rows["repaired_init_saved_selection"]["passed"] and ast_rows["repaired_evaluate"]["passed"]
    if not response_source_ok:
        gaps.append("saved repaired full-support response identity failed")
    if not transport_ok:
        gaps.append("compact-to-quiet/Rp transport or pressure-tail guard failed")

    archive = _archive_provenance(parent, reports["archive"], reports["repaired"], gaps)
    graph_passed = all(matrix_checks.values()) and all(operator_source_checks.values()) and all(terminal_checks.values()) and all(physical_checks.values()) and all(coefficient_arc_checks.values()) and all(normalization_checks.values())
    source_full_support_identity = {
        "passed": bool(parent_pin.get("passed") and joint_ok and transitive_ok and ast_passed and response_source_ok and transport_ok and graph_passed and all(scale_guards.values())),
        "statement": "reserved_n1_own_swirl_once response = -normalized(full) as ordinary C1 functions",
        "ordinary_Z_orders": list(_ORDERS),
        "response_source": {"module": _REQUIRED_SOURCE["response"], "class": "OriginalReservedM15Response", "methods": ["partial_row", "evaluate"], "native_compact_x": "2", "exact_angular_response": [0, 0, 1]},
        "physical_pairs": dict(_PHYSICAL),
        "normalized_full_pair": list(_NORMALIZED["full"]),
        "terminal_patch": [53905, 53906],
        "reduced_patch": [53921, 53922],
        "matrix_action": {"pressure_row": [25072, 25075, 25078], "inverse_pressure_column": [25117, 25134, 25150], "exact_common_matrix_inverse_identity": [0, 0, 1], "boxes_not_identity": True},
        "actual_arc_source_identity": coefficient_identity,
        "identity_before_boxing": True,
        "response_contract_rows": response_identity_rows,
        "transport": transport_rows,
        "forward_quiet_q1_selection": {"passed": forward_selection_ok, "rows": forward_selected, "source": _FORWARD_RELATIVE},
        "numeric_value_available": False,
        "numeric_full_closure": False,
    }
    response_minus_target = {
        "passed": bool(source_full_support_identity["passed"] and archive.get("passed")),
        "statement": "response - negative pulse target = -(inner + lower + preRp) through ordinary C1 before boxing",
        "ordinary_Z_orders": list(_ORDERS),
        "source_signs": {"reserved_response": "-normalized(full)", "independent_target": "-normalized(pulse)", "difference": "-full + pulse"},
        "physical_pairs": dict(_PHYSICAL),
        "normalized_pairs": dict(_NORMALIZED),
        "shared_unit": {"log_node": 53421, "unit_node": 53427, "square_node": 53430, "U_Z_exact_zero_node": 53399, "U_Z_exact_zero": normalization_checks["U_Z_exact_zero"]},
        "expression": "-normalize_C1(full)+normalize_C1(pulse)=-normalize_C1(inner+lower+preRp)",
        "before_boxing": True,
        "bounds_are_not_function_values": True,
        "joint_numeric_closure": False,
        "component_source_digests": {key: _digest(value) for key, value in _mapping(_mapping(_mapping(full.get("samples")).get("whole_Z")).get("pressure_component_enclosures")).get("m15", {}).items()},
        "archive_provenance_passed": archive.get("passed") is True,
    }
    whole = _mapping(_mapping(full.get("samples")).get("whole_Z"))
    components = _mapping(whole.get("pressure_component_enclosures")).get("m15")
    def _c1_row(value: Any) -> bool:
        if isinstance(value, Mapping):
            value = value.get("C1_enclosures")
        return isinstance(value, list) and len(value) == 2

    cached_bounds_refs = {
        "full_bridge": paths["full_bridge"].name + ":samples.whole_Z.pressure_component_enclosures.m15",
        "component_names": list(components) if isinstance(components, Mapping) else [],
        "component_C1_rows_present": {key: _c1_row(value) for key, value in (components.items() if isinstance(components, Mapping) else ())},
        "bounds_are_function_values": False,
        "source_enclosure_midpoint_used": False,
    }
    if not cached_bounds_refs["component_C1_rows_present"] or not all(cached_bounds_refs["component_C1_rows_present"].values()):
        gaps.append("full bridge component C1 enclosure reference incomplete")

    compact = bool(
        parent_pin.get("passed") and joint_ok and transitive_ok and all(family_rows.values()) and ast_passed
        and response_source_ok and transport_ok and archive.get("passed") and graph_passed and not gaps
    )
    return {
        "passed": compact,
        "compact_proofs_passed": compact,
        "correlation_proved": bool(source_full_support_identity["passed"] and response_minus_target["passed"]),
        "global_closure_claimed": False,
        "parent_report": _PARENT_RELATIVE,
        "parent_report_sha256": _PARENT_SHA256,
        "parent_report_pin": parent_pin,
        "source_family": dict(_FAMILY),
        "input_hashes": input_hashes,
        "implementation_hashes": implementation_hashes,
        "added_pins": added_pins,
        "transitive_hash_validation": {"passed": transitive_ok, "input_passed": transitive_input_passed, "input_total": transitive_input_total, "implementation_passed": transitive_impl_passed, "implementation_total": transitive_impl_total},
        "source_pins": source_pins,
        "AST": ast_rows,
        "ast_function_hashes": {key: row.get("function_AST_sha256") for key, row in ast_rows.items()},
        "graph_identity": {"passed": graph_passed, "matrix_checks": matrix_checks, "operator_source_checks": operator_source_checks, "terminal_checks": terminal_checks, "physical_checks": physical_checks, "coefficient_arc_checks": coefficient_arc_checks, "normalization_checks": normalization_checks, "full_normalization_nodes": full_nodes, "pulse_normalization_nodes": pulse_nodes, "source_executed": False},
        "source_full_support_identity": source_full_support_identity,
        "response_minus_target": response_minus_target,
        "archive_provenance": archive,
        "normalized_expression": response_minus_target,
        "cached_bounds_refs": cached_bounds_refs,
        "physical_scale_guards": scale_guards,
        "numeric_value_available": False,
        "numeric_full_closure": False,
        "numeric_value_unavailable_reason": "source identity and saved C1 enclosure references are exported; no numeric full closure or interval arithmetic is performed",
        "source_executed": False,
        "root": str(repo),
        "directory": str(work),
        "gaps": gaps,
    }


__all__ = ["build_reserved_full_support_source"]
