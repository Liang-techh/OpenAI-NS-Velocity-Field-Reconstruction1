"""Numeric contraction from the saved outer-chain state through ``Rp``.

This adapter consumes the already accepted 190-term outer-chain report and
one saved compact/quiet numerical report.  It binds the two new source seams
to their actual packet endpoints, then replaces the old outer remainder and
the compact viscous sector by one directed positive error enclosure.  Source
owners, constructors, producers, and report checkers are deliberately not
loaded or executed here.
"""

from __future__ import annotations

import copy
import ast
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .original_core_numeric import PREFIX, ROOT


SAVED_Z = ("-.3", "0", ".3")
RAW_COUNT = 213
PREVIOUS_REDUCED_COUNT = 190
NEW_REDUCED_COUNT = 189
UNCHANGED_COUNT = 188

PARENT = PREFIX + "compliant_current_original_outer_chain_m15_numeric.json"
PARENT_SHA = "cf185c36affc4183984aff5ef9caf119c68ea90d0b8e95535b1f8698e2161e53"
COMPACT_QUIET_REPORT = PREFIX + "compliant_current_original_compact_quiet_m15_numeric.json"
COMPACT_QUIET_SHA = "7ccb30dbb3e440a3e598eb30deaf03509cf8624e429478aeeb364671baaac071"
COMPACT_PRESSURE_ROWS_SOURCE = "src/openai_ns_reconstruction/original_compact_pressure_rows.py"

OLD_EXTERNAL_TERM = "long_outer_chain_finite_N_viscous_external_endpoints"
COMPACT_VISCOUS_TERM = "compact_finite_N_viscous_error"
NEW_EXTERNAL_TERM = "long_compact2Rc_finite_N_viscous_external_endpoints"
QUIET_TERM = "quiet_exact_original_mean_square_integral"
RESERVED_TERM = "reserved_n1_own_swirl_once"

COMPACT_LOWER = "1"
COMPACT_UPPER = "2"
QUIET_LOWER = "0"
QUIET_UPPER = "1"
QUIET_SOURCE_COORDINATE_NODE = 263344
COMPACT_PREFIX_COUNT = 207
COMPACT_TERM_COUNT = 211
QUIET_PREFIX_COUNT = 211
QUIET_TERM_COUNT = 212


def _canonical(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _digest(value: Any) -> str:
    return hashlib.sha256(_canonical(value).encode("utf-8")).hexdigest()


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{label} must be a mapping")
    return value


def _find(rows: Any, *, label: str, **criteria: Any) -> Mapping[str, Any]:
    if not isinstance(rows, (list, tuple)):
        raise ValueError(f"Saved {label} rows are missing")
    for row in rows:
        if isinstance(row, Mapping) and all(row.get(key) == value for key, value in criteria.items()):
            return row
    raise ValueError(f"Saved {label} not found: {criteria}")


def _coerce_coordinate(value: Any) -> Any:
    if isinstance(value, str):
        return value
    if isinstance(value, int):
        return str(value)
    return value


def _path(owner: Any, relative: str) -> Path:
    candidates: list[Path] = []
    directory = getattr(owner, "directory", None)
    if directory is not None:
        base = Path(directory)
        candidates.extend((base / relative, base / "experiments" / "root_st073" / relative))
    candidates.extend((ROOT / relative, ROOT / "experiments" / "root_st073" / relative))
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return candidates[0] if candidates else ROOT / relative


def _pin_file(owner: Any, relative: str, expected: str) -> dict[str, Any]:
    path = _path(owner, relative)
    actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
    hashes = _mapping(getattr(owner, "hashes", {}), "owner.hashes")
    code_hashes = _mapping(getattr(owner, "code_hashes", {}), "owner.code_hashes")
    basename = Path(relative).name
    for key, known in ((relative, hashes.get(relative)), (basename, hashes.get(basename)),
                       (relative, code_hashes.get(relative))):
        if known is not None and known != actual:
            raise ValueError(f"Changed accepted source: {key}")
    if actual != expected:
        raise ValueError(f"Unexpected accepted source: {relative}")
    return {"path": relative, "sha256": actual, "expected_sha256": expected, "passed": True}


def _load_json(owner: Any, relative: str, expected: str) -> Mapping[str, Any]:
    _pin_file(owner, relative, expected)
    value = json.loads(_path(owner, relative).read_text(encoding="utf-8"))
    return _mapping(value, relative)


def _geometry_source_binding(owner: Any) -> dict[str, Any]:
    """Pin the saved q-to-s front-end geometry by AST and owner source hash.

    The compact/quiet numeric packets use the fraction ``q``.  The accepted
    source graph uses the log-offset ``s=L*q``; this check records that
    distinction without importing or constructing the historical owner.
    """
    relative = COMPACT_PRESSURE_ROWS_SOURCE
    path = _path(owner, relative)
    code_hashes = _mapping(getattr(owner, "code_hashes", {}), "owner.code_hashes")
    expected = code_hashes.get(relative)
    if not isinstance(expected, str):
        raise ValueError(f"Pinned geometry source hash missing: {relative}")
    actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
    if actual != expected:
        raise ValueError(f"Changed pinned geometry source: {relative}")
    text = path.read_text(encoding="utf-8")
    try:
        tree = ast.parse(text, filename=str(path))
    except SyntaxError as exc:
        raise ValueError("Pinned compact/quiet geometry source is not valid Python") from exc

    methods: dict[tuple[str, str], str] = {}
    for node in ast.walk(tree):
        if not isinstance(node, ast.ClassDef):
            continue
        for child in node.body:
            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)) and child.name == "geometry":
                methods[(node.name, child.name)] = ast.get_source_segment(text, child) or ""
    compact_text = methods.get(("OriginalCompactPressureRows", "geometry"), "")
    quiet_text = methods.get(("OriginalQuietPressureRows", "geometry"), "")
    compact_required = (
        "q,p,t=self.coordinate(value)",
        "sp,st=c.log(p),iv.log(t)",
        "physical_jacobian_point=1/p",
    )
    quiet_required = (
        "q,p,t=self.coordinate(value)",
        "sp,st=self.length*p,self.lengthi*t",
        "physical_jacobian_point=self.length",
        "native=q",
    )
    compact_passed = bool(compact_text) and all(token in compact_text for token in compact_required)
    quiet_passed = bool(quiet_text) and all(token in quiet_text for token in quiet_required)
    if not compact_passed or not quiet_passed:
        raise ValueError("Pinned compact/quiet q-to-s geometry AST binding changed")
    return {
        "passed": True,
        "path": relative,
        "sha256": actual,
        "owner_code_hash": expected,
        "ast_parse_passed": True,
        "compact_geometry": {
            "passed": True,
            "saved_coordinate": "q",
            "source_mapping": "s=log(p) for compact x=p=R/Rc",
            "physical_jacobian": "1/p",
        },
        "quiet_geometry": {
            "passed": True,
            "saved_coordinate": "q",
            "source_coordinate_node": QUIET_SOURCE_COORDINATE_NODE,
            "source_mapping": "s=L*q",
            "source_length": "L=-60*log(mu)-2-log(2)",
            "physical_jacobian": "L",
            "Rp_fraction_is_not_native_s_one": True,
        },
        "old_owner_or_constructor_replayed": False,
        "source_values_evaluated": False,
    }


def _parent_report(owner: Any) -> Mapping[str, Any]:
    for candidate in (getattr(owner, "parent", None), getattr(owner, "parent_report", None)):
        if isinstance(candidate, Mapping) and isinstance(candidate.get("samples"), list):
            return candidate
        report = getattr(candidate, "report", None)
        if isinstance(report, Mapping) and isinstance(report.get("samples"), list):
            return report
    return _load_json(owner, PARENT, PARENT_SHA)


def _term_names(rows: Any, label: str) -> list[str]:
    if not isinstance(rows, list):
        raise ValueError(f"{label} must be a list")
    names = [row.get("name") if isinstance(row, Mapping) else None for row in rows]
    if any(name is None for name in names):
        raise ValueError(f"{label} has an unnamed term")
    return names


def _exact(label: str, left: Any, right: Any) -> None:
    if left != right:
        raise ValueError(f"{label} full dictionary equality changed")


def _positive_atom(atom: Any, label: str) -> Mapping[str, Any]:
    value = _mapping(atom, label)
    if value.get("positive_error_bound_atom_only") is not True or value.get("source_value_selected") is not False:
        raise ValueError(f"{label} must remain an error-only positive atom")
    if not isinstance(value.get("source_enclosure_exact_mpf_tuples"), list):
        raise ValueError(f"{label} source enclosure is missing")
    return value


def _parent_long_witness(parent_sample: Mapping[str, Any]) -> Mapping[str, Any]:
    budget = _mapping(parent_sample.get("external_endpoint_budget_witness"), "parent external endpoint budget")
    witness = _mapping(budget.get("own_longreshape0"), "parent own long native0 witness")
    rows = witness.get("radius_weighted_C0_C1_log_bounds")
    if not isinstance(rows, list) or len(rows) != 2 or any(value is None for value in rows):
        raise ValueError("Parent own long native0 C0/Z1 positive bounds required")
    for index, value in enumerate(rows):
        _positive_atom(value, f"parent long native0 atom {index}")
    return witness


def _proofs(seam_proofs: Mapping[str, Any]) -> tuple[dict[str, Any], Mapping[str, Any] | None, Mapping[str, Any] | None]:
    if not isinstance(seam_proofs, Mapping) or seam_proofs.get("passed") is not True:
        raise ValueError("Compact/quiet source seam proofs must pass")
    joins = seam_proofs.get("joins")
    if not isinstance(joins, Mapping):
        joins = seam_proofs.get("seams")
    if not isinstance(joins, Mapping):
        joins = seam_proofs
    expected = {
        "power_compact": ("O3_power", "compact", 2, 1, (250215, 260181), (48, 1), ("2", "1")),
        "compact_quiet": ("compact", "quiet", 2, 0, (260181, 263344), (48, 0), ("2", "0")),
    }
    result: dict[str, Any] = {}
    for name, (left_chart, right_chart, left_endpoint, right_endpoint, coordinates, queries, query_values) in expected.items():
        proof = joins.get(name)
        if not isinstance(proof, Mapping) or proof.get("passed") is not True:
            raise ValueError(f"Accepted {name} source proof must pass")
        for key, expected_value in (("left_chart", left_chart), ("right_chart", right_chart),
                                    ("left_endpoint", left_endpoint), ("right_endpoint", right_endpoint),
                                    ("left_native_query", left_endpoint), ("right_native_query", right_endpoint)):
            if key in proof and proof.get(key) != expected_value:
                raise ValueError(f"{name} source proof has the wrong {key}")
        # Endpoint values (2 Rc, 1 Rc, or 0 Rc) are distinct from the
        # archived graph-node IDs.  When a source helper emits the latter,
        # bind those IDs explicitly without treating them as coordinates.
        node_values = {
            "left_endpoint_node": 48,
            "right_endpoint_node": right_endpoint,
            "left_query_node": 48,
            "right_query_node": right_endpoint,
            "left_native_query_node": 48,
            "right_native_query_node": right_endpoint,
        }
        for key, expected_value in node_values.items():
            if key in proof and proof.get(key) != expected_value:
                raise ValueError(f"{name} source proof has the wrong {key}")
        direct_scope_values = {
            "left_coordinate_node": coordinates[0],
            "right_coordinate_node": coordinates[1],
            "left_native_coordinate_node": coordinates[0],
            "right_native_coordinate_node": coordinates[1],
        }
        for key, expected_value in direct_scope_values.items():
            if key in proof and proof.get(key) != expected_value:
                raise ValueError(f"{name} source proof has the wrong {key}")
        # Source correspondence helpers normally carry these as left_scope
        # and right_scope.  Validate the archived coordinate and query nodes
        # when emitted, while allowing a smaller source proof schema.
        for side, coordinate_value, query_value, point_value in (
            ("left", coordinates[0], queries[0], query_values[0]),
            ("right", coordinates[1], queries[1], query_values[1]),
        ):
            scopes = [proof[key] for key in (f"{side}_scope", f"{side}_endpoint_scope")
                      if key in proof and isinstance(proof[key], Mapping)]
            for scope in scopes:
                for key in ("coordinate_node", "native_coordinate_node"):
                    if key in scope and scope.get(key) != coordinate_value:
                        raise ValueError(f"{name} source proof has the wrong {side} {key}")
                for key in ("native_query_node", "query_node", "query"):
                    if key in scope and scope.get(key) != query_value:
                        raise ValueError(f"{name} source proof has the wrong {side} {key}")
                for key in ("query_value", "native_query_value"):
                    if key in scope and str(scope.get(key)) != point_value:
                        raise ValueError(f"{name} source proof has the wrong {side} {key}")
        if proof.get("source_values_selected") is True or proof.get("numeric_caps_used") is True:
            raise ValueError(f"{name} source proof selected a value or cap")
        result[name] = copy.deepcopy(dict(proof))

    def optional_certificate(keys: tuple[str, ...]) -> Mapping[str, Any] | None:
        present = [(key, seam_proofs[key]) for key in keys if key in seam_proofs and seam_proofs[key] is not None]
        if not present:
            return None
        _, value = present[0]
        certificate = _mapping(value, "optional source certificate")
        if certificate.get("passed") is not True:
            raise ValueError("Explicit optional source certificate failed")
        return copy.deepcopy(dict(certificate))

    terminal = optional_certificate(("terminal_D_zero_certificate",))
    if terminal is not None:
        for key in ("D_upper_zero", "ordinary_Z1", "physical_R_once"):
            if terminal.get(key) is not True:
                raise ValueError(f"Terminal D-zero certificate lacks {key}")
    quiet = optional_certificate(("quiet_transport_certificate", "quietconstantD", "quiet_constant_D"))
    if quiet is not None:
        if not any(quiet.get(key) is True for key in ("constant_D", "D_constant", "lower_upper_D_equal", "quiet_D_constant")):
            raise ValueError("Quiet transport certificate lacks an explicit constant-D flag")
        for key in ("ordinary_Z1", "physical_R_once"):
            if quiet.get(key) is not True:
                raise ValueError(f"Quiet transport certificate lacks {key}")
    return result, quiet, terminal


def _sample(report: Mapping[str, Any], Z: str, chart: str, coordinate: str) -> Mapping[str, Any]:
    return _find(report.get("samples"), label=f"{chart} compact/quiet sample",
                 Z=Z, chart=chart, native_coordinate=_coerce_coordinate(coordinate))


def _validate_sample(
    owner: Any,
    parent_sample: Mapping[str, Any],
    compact_report: Mapping[str, Any],
    family: Mapping[str, Any],
) -> dict[str, Any]:
    Z = str(parent_sample.get("Z"))
    raw = parent_sample.get("raw_m15_C0_C1_terms")
    previous = parent_sample.get("correlated_reduced_m15_C0_C1_terms")
    if not isinstance(raw, list) or len(raw) != RAW_COUNT:
        raise ValueError("Parent raw 213 m15 state must remain verbatim")
    if not isinstance(previous, list) or len(previous) != PREVIOUS_REDUCED_COUNT:
        raise ValueError("Parent outer-chain 190-term view required")
    if compact_report.get("source_family") != family:
        raise ValueError("Compact/quiet report source family changed")
    selected_parent = parent_sample.get("selected_N_evidence")
    if not isinstance(selected_parent, Mapping):
        raise ValueError("Parent selected-N evidence required")

    compact_lower = _sample(compact_report, Z, "compact", COMPACT_LOWER)
    compact_upper = _sample(compact_report, Z, "compact", COMPACT_UPPER)
    quiet_lower = _sample(compact_report, Z, "quiet", QUIET_LOWER)
    quiet_upper = _sample(compact_report, Z, "quiet", QUIET_UPPER)
    rows = {
        "compact_lower": compact_lower,
        "compact_upper": compact_upper,
        "quiet_lower": quiet_lower,
        "quiet_upper": quiet_upper,
    }
    for label, row in rows.items():
        if row.get("source_family") != family:
            raise ValueError(f"{label} source family changed")
        if row.get("actual_selected_integer") != selected_parent:
            raise ValueError(f"{label} selected N changed")

    c_lower_terms = compact_lower.get("m15_C0_C1_terms")
    c_lower_incoming = compact_lower.get("incoming_m15_terms")
    c_upper_terms = compact_upper.get("m15_C0_C1_terms")
    c_upper_incoming = compact_upper.get("incoming_m15_terms")
    q_lower_terms = quiet_lower.get("m15_C0_C1_terms")
    q_lower_incoming = quiet_lower.get("incoming_m15_terms")
    q_upper_terms = quiet_upper.get("m15_C0_C1_terms")
    q_upper_incoming = quiet_upper.get("incoming_m15_terms")
    for label, terms, incoming, term_count, prefix_count in (
        ("compact lower", c_lower_terms, c_lower_incoming, COMPACT_TERM_COUNT, COMPACT_PREFIX_COUNT),
        ("compact upper", c_upper_terms, c_upper_incoming, COMPACT_TERM_COUNT, COMPACT_PREFIX_COUNT),
        ("quiet lower", q_lower_terms, q_lower_incoming, QUIET_TERM_COUNT, QUIET_PREFIX_COUNT),
        ("quiet upper", q_upper_terms, q_upper_incoming, QUIET_TERM_COUNT, QUIET_PREFIX_COUNT),
    ):
        if not isinstance(terms, list) or len(terms) != term_count:
            raise ValueError(f"{label} complete packet count changed")
        if not isinstance(incoming, list) or len(incoming) != prefix_count:
            raise ValueError(f"{label} incoming prefix count changed")
        _term_names(terms, label + " terms")
        _term_names(incoming, label + " incoming")

    raw_prefix = raw[:COMPACT_PREFIX_COUNT]
    _exact("compact lower incoming versus parent raw prefix", c_lower_incoming, raw_prefix)
    _exact("compact lower terms prefix versus incoming", c_lower_terms[:COMPACT_PREFIX_COUNT], c_lower_incoming)
    _exact("compact upper incoming versus parent raw prefix", c_upper_incoming, raw_prefix)
    _exact("compact upper terms prefix versus incoming", c_upper_terms[:COMPACT_PREFIX_COUNT], c_upper_incoming)
    _exact("compact lower and upper cumulative prefix", c_lower_terms[:COMPACT_PREFIX_COUNT], c_upper_terms[:COMPACT_PREFIX_COUNT])
    _exact("compact upper complete packet versus quiet lower incoming", c_upper_terms, q_lower_incoming)
    _exact("quiet lower terms prefix versus incoming", q_lower_terms[:QUIET_PREFIX_COUNT], q_lower_incoming)
    _exact("quiet lower complete packet versus quiet upper incoming", q_lower_terms[:QUIET_PREFIX_COUNT], q_upper_incoming)
    _exact("quiet upper terms prefix versus incoming", q_upper_terms[:QUIET_PREFIX_COUNT], q_upper_incoming)

    if compact_lower.get("exact_empty_increment") is not True or quiet_lower.get("exact_empty_increment") is not True:
        raise ValueError("Compact and quiet lower endpoints must carry exact empty increments")
    if raw[COMPACT_PREFIX_COUNT]["name"] != "compact_exact_original_mean_square_integral":
        raise ValueError("Parent raw compact prefix route changed")
    if raw[COMPACT_PREFIX_COUNT + 1]["name"] != "compact_finite_N_endpoint_error":
        raise ValueError("Parent raw compact endpoint sector changed")
    if raw[COMPACT_PREFIX_COUNT + 2]["name"] != "compact_finite_N_volume_error":
        raise ValueError("Parent raw compact volume sector changed")
    if raw[COMPACT_PREFIX_COUNT + 3]["name"] != COMPACT_VISCOUS_TERM:
        raise ValueError("Parent raw compact viscous sector changed")
    if raw[QUIET_PREFIX_COUNT]["name"] != QUIET_TERM or raw[QUIET_PREFIX_COUNT + 1]["name"] != RESERVED_TERM:
        raise ValueError("Parent raw quiet/reserved route changed")
    expected_compact_tail = [
        "compact_exact_original_mean_square_integral",
        "compact_finite_N_endpoint_error",
        "compact_finite_N_volume_error",
        COMPACT_VISCOUS_TERM,
    ]
    if _term_names(c_upper_terms[COMPACT_PREFIX_COUNT:], "compact upper tail") != expected_compact_tail:
        raise ValueError("Compact upper tail route changed")
    if _term_names(q_upper_terms[COMPACT_PREFIX_COUNT:], "quiet upper tail") != expected_compact_tail + [QUIET_TERM]:
        raise ValueError("Quiet upper tail route changed")

    compact_aggregate = c_upper_terms[COMPACT_PREFIX_COUNT + 3]
    if compact_aggregate.get("name") != COMPACT_VISCOUS_TERM:
        raise ValueError("Compact upper viscous aggregate changed")
    _exact("parent raw compact viscous term versus saved compact aggregate",
           raw[COMPACT_PREFIX_COUNT + 3], compact_aggregate)
    binding = _mapping(compact_aggregate.get("actual_function_binding"), "compact viscous source binding")
    if (binding.get("chart") != "compact" or binding.get("incoming_terms_retained_once") != COMPACT_PREFIX_COUNT
            or binding.get("error_enclosure_only") is not True):
        raise ValueError("Compact viscous aggregate source binding changed")
    logs = compact_aggregate.get("log_scale_terms")
    if not isinstance(logs, list) or len(logs) != 1:
        raise ValueError("Complete compact viscous aggregate must have one positive log atom")
    compact_atom = _positive_atom(logs[0], "compact upper viscous aggregate")
    _exact("compact viscous aggregate retained through quiet lower", compact_aggregate,
           q_lower_terms[COMPACT_PREFIX_COUNT + 3])
    _exact("compact viscous aggregate retained through quiet upper", compact_aggregate,
           q_upper_terms[COMPACT_PREFIX_COUNT + 3])

    quiet_term_lower = q_lower_terms[QUIET_PREFIX_COUNT]
    quiet_term_upper = q_upper_terms[QUIET_PREFIX_COUNT]
    if quiet_term_lower.get("name") != QUIET_TERM or quiet_term_upper.get("name") != QUIET_TERM:
        raise ValueError("Quiet terminal term changed")
    quiet_binding = _mapping(quiet_term_upper.get("actual_function_binding"), "quiet source binding")
    if quiet_binding.get("chart") != "quiet" or quiet_binding.get("incoming_terms_retained_once") != QUIET_PREFIX_COUNT:
        raise ValueError("Quiet source binding changed")

    old_external = [term for term in previous if term.get("name") == OLD_EXTERNAL_TERM]
    old_compact = [term for term in previous if term.get("name") == COMPACT_VISCOUS_TERM]
    if len(old_external) != 1 or len(old_compact) != 1:
        raise ValueError("Parent 190-term view must contain one outer external and one compact viscous term")
    _exact("parent reduced compact viscous term versus saved compact aggregate",
           old_compact[0], compact_aggregate)
    if sum(term.get("name") == QUIET_TERM for term in previous) != 1:
        raise ValueError("Parent quiet mean-square term must remain once")
    if sum(term.get("name") == RESERVED_TERM for term in previous) != 1:
        raise ValueError("Parent reserved own-swirl term must remain once")

    witness = {
        "passed": True,
        "parent_raw_count": len(raw),
        "parent_raw_prefix_count": len(raw_prefix),
        "compact_lower_terms_count": len(c_lower_terms),
        "compact_lower_incoming_count": len(c_lower_incoming),
        "compact_upper_terms_count": len(c_upper_terms),
        "compact_upper_incoming_count": len(c_upper_incoming),
        "quiet_lower_terms_count": len(q_lower_terms),
        "quiet_lower_incoming_count": len(q_lower_incoming),
        "quiet_upper_terms_count": len(q_upper_terms),
        "quiet_upper_incoming_count": len(q_upper_incoming),
        "compact_lower_incoming_equals_parent_raw": True,
        "compact_lower_terms_prefix_equals_incoming": True,
        "compact_upper_incoming_equals_parent_raw": True,
        "compact_upper_terms_prefix_equals_incoming": True,
        "compact_upper_equals_quiet_lower_incoming": True,
        "quiet_lower_terms_prefix_equals_incoming": True,
        "quiet_lower_equals_quiet_upper_incoming": True,
        "quiet_upper_terms_prefix_equals_incoming": True,
        "exact_dictionary_equality": True,
        "parent_raw_prefix_sha256": _digest(raw_prefix),
        "compact_lower_terms_sha256": _digest(c_lower_terms),
        "compact_upper_terms_sha256": _digest(c_upper_terms),
        "quiet_lower_terms_sha256": _digest(q_lower_terms),
        "quiet_upper_terms_sha256": _digest(q_upper_terms),
        "compact_aggregate_sha256": _digest(compact_aggregate),
        "selected_N_exact": True,
        "selected_N": copy.deepcopy(selected_parent),
    }
    return {
        "Z": Z,
        "parent_sample": parent_sample,
        "raw": raw,
        "previous": previous,
        "selected_N": selected_parent,
        "target": parent_sample.get("independent_target"),
        "old_external": old_external[0],
        "old_compact": old_compact[0],
        "compact_aggregate": compact_aggregate,
        "compact_atom": compact_atom,
        "compact_packet_witness": witness,
        "quiet_term_lower": quiet_term_lower,
        "quiet_term_upper": quiet_term_upper,
    }


def _trace_chain(*, terminal_zero: bool, quiet_transport: bool) -> dict[str, Any]:
    if terminal_zero and quiet_transport:
        function = "D(long native0)"
        upper = "Rp (D(compact2Rc)=0 and D is constant through quiet)"
    elif terminal_zero:
        function = "D(long native0)-D(compact2Rc=0)"
        upper = "compact native2Rc"
    elif quiet_transport:
        function = "D(long native0)-D(Rp)"
        upper = "Rp (quiet constant-D transport)"
    else:
        function = "D(long native0)-D(compact2Rc)"
        upper = "compact native2Rc"
    return {
        "external_endpoints": {"lower": "long native0", "upper": upper},
        "cancelled_interfaces_explicit": [
            "O3 power native2Rc = compact native1Rc",
            "compact native2Rc = quiet native0 2Rc",
        ],
        "function": function,
        "delta_definition": "D=R*(Qy_modified-Qy_leading)",
        "ordinary_Z1": True,
        "physical_R_once": True,
        "leading_subtraction_retained": True,
        "quiet_coordinate_binding": {
            "saved_quiet_lower_fraction": QUIET_LOWER,
            "saved_quiet_Rp_fraction": QUIET_UPPER,
            "actual_source_coordinate_node": QUIET_SOURCE_COORDINATE_NODE,
            "actual_source_coordinate": "s=L*q",
            "quiet_lower_source_s": "0",
            "quiet_Rp_source_s": "L",
            "quiet_Rp_fraction_is_not_native_s_one": True,
        },
        "terminal_D_zero_used": terminal_zero,
        "quiet_transport_used": quiet_transport,
    }


def _evaluate_sample(
    owner: Any,
    parent: Mapping[str, Any],
    source: Mapping[str, Any],
    proofs: Mapping[str, Any],
    sample: Mapping[str, Any],
    quiet_certificate: Mapping[str, Any] | None,
    terminal_certificate: Mapping[str, Any] | None,
) -> dict[str, Any]:
    raw = sample["raw"]
    previous = sample["previous"]
    old_external = sample["old_external"]
    old_compact = sample["old_compact"]
    witness = _parent_long_witness(sample["parent_sample"])
    long_atoms = list(witness["radius_weighted_C0_C1_log_bounds"])
    if len(long_atoms) != 2:
        raise ValueError("Exactly two parent long native0 C0/Z1 atoms are required")
    terminal_zero = terminal_certificate is not None
    quiet_transport = quiet_certificate is not None
    positive_inputs = [owner.decode(atom) for atom in long_atoms]
    if not terminal_zero:
        positive_inputs.append(owner.decode(sample["compact_atom"]))
    bound = owner.positive_union(positive_inputs)
    replacement = copy.deepcopy(old_external)
    replacement["name"] = NEW_EXTERNAL_TERM
    replacement["log_scale_terms"] = [owner.record(bound)]
    replacement["actual_function_binding"] = {
        "error_enclosure_only": True,
        "function": _trace_chain(terminal_zero=terminal_zero, quiet_transport=quiet_transport)["function"],
        "delta_definition": "R*(Qy_modified-Qy_leading)",
        "ordinary_Z1": True,
        "physical_R_once": True,
        "external_endpoints_not_zeroed": not terminal_zero,
        "surviving_long0_inlet_not_zeroed": True,
        "endpoint_zero_only_by_source_identity": True,
        "signed_centers_evaluated": False,
        "no_error_cap_subtraction": True,
        "no_interval_midpoint_or_overlap_closure": True,
        "directed_positive_union_long_native0_then_compact210": not terminal_zero,
        "parent_long_native0_budget_retained": True,
        "compact210_aggregate_retained_as_provenance": True,
        "compact210_aggregate_used_as_error_only": not terminal_zero,
        "old_outer_full_202_aggregate_used": False,
        "terminal_D_zero_certificate_used": terminal_zero,
        "quiet_constant_D_transport_used": quiet_transport,
        "source_functions_values_or_caps_not_selected": True,
        "power_compact_source_join_proved": True,
        "compact_quiet_source_join_proved": True,
        "source_seam_proofs": copy.deepcopy(dict(proofs)),
        "terminal_D_zero_certificate": copy.deepcopy(dict(terminal_certificate)) if terminal_certificate is not None else None,
        "quiet_transport_certificate": copy.deepcopy(dict(quiet_certificate)) if quiet_certificate is not None else None,
    }

    removed_names: list[str] = []
    reduced: list[dict[str, Any]] = []
    for term in previous:
        name = term.get("name")
        if name in {OLD_EXTERNAL_TERM, COMPACT_VISCOUS_TERM}:
            removed_names.append(name)
        else:
            reduced.append(copy.deepcopy(term))
    # Place the replacement at the old external term's position while keeping
    # every other 188 entry byte-for-byte unchanged and in order.
    insertion_index = next(i for i, term in enumerate(previous) if term.get("name") == OLD_EXTERNAL_TERM)
    reduced.insert(insertion_index, replacement)
    if sorted(removed_names) != sorted((OLD_EXTERNAL_TERM, COMPACT_VISCOUS_TERM)) or len(reduced) != NEW_REDUCED_COUNT:
        raise ValueError("Exact new 189-term view required")
    unchanged_old = [term for term in previous if term.get("name") not in {OLD_EXTERNAL_TERM, COMPACT_VISCOUS_TERM}]
    unchanged_new = [term for term in reduced if term.get("name") != NEW_EXTERNAL_TERM]
    if len(unchanged_old) != UNCHANGED_COUNT or unchanged_new != unchanged_old:
        raise ValueError("Other 188 correlated sectors must remain verbatim")
    if reduced[:6] != previous[:6]:
        raise ValueError("Six signed inner/lower terms changed")
    if sum(term.get("name") == NEW_EXTERNAL_TERM for term in reduced) != 1:
        raise ValueError("Exactly one new compact2Rc external replacement required")
    if sum(term.get("name") == QUIET_TERM for term in reduced) != 1:
        raise ValueError("Quiet mean-square term must remain once")
    if sum(term.get("name") == RESERVED_TERM for term in reduced) != 1:
        raise ValueError("Reserved own-swirl term must remain once")

    target = copy.deepcopy(sample.get("target"))
    if target is None:
        target = copy.deepcopy(parent.get("independent_target"))
    if target is None:
        raise ValueError("Independent target must remain available")
    trace = _trace_chain(terminal_zero=terminal_zero, quiet_transport=quiet_transport)
    packet = sample["compact_packet_witness"]
    # Keep the compact upper packet in the same canonical shape as the
    # earlier numerical adapters.  The runner uses ``aggregate_log_atom``
    # and ``aggregate`` to recompute the directed union; provenance-only
    # copies below deliberately use ``raw_term_index`` so a terminal
    # D-zero certificate cannot accidentally authorize the packet as a
    # positive union component.
    compact_aggregate_packet = {
        "join": "compact_quiet",
        "chart": "compact",
        "endpoint_coordinate": COMPACT_LOWER,
        "aggregate_coordinate": COMPACT_UPPER,
        "prefix_count": COMPACT_PREFIX_COUNT,
        "term_count": COMPACT_TERM_COUNT,
        "raw_index": COMPACT_PREFIX_COUNT + 3,
        "term_name": COMPACT_VISCOUS_TERM,
        "exact_prefix_witness": copy.deepcopy(packet),
        "aggregate": copy.deepcopy(sample["compact_aggregate"]),
        "aggregate_log_atom": copy.deepcopy(sample["compact_atom"]),
        "selected_integer": copy.deepcopy(sample["selected_N"]),
        "source_report": COMPACT_QUIET_REPORT,
    }
    compact_provenance = {
        "join": "compact_quiet",
        "chart": "compact",
        "endpoint_coordinate": COMPACT_LOWER,
        "aggregate_coordinate": COMPACT_UPPER,
        "raw_term_index": COMPACT_PREFIX_COUNT + 3,
        "term_name": COMPACT_VISCOUS_TERM,
        "exact_prefix_witness": copy.deepcopy(packet),
        "term": copy.deepcopy(sample["compact_aggregate"]),
        "log_atom": copy.deepcopy(sample["compact_atom"]),
        "selected_integer": copy.deepcopy(sample["selected_N"]),
        "source_report": COMPACT_QUIET_REPORT,
    }
    budget: dict[str, Any] = {
        "own_longreshape0": copy.deepcopy(dict(witness)),
        "positive_union": owner.record(bound),
        "positive_union_inputs": copy.deepcopy(long_atoms) + ([] if terminal_zero else [copy.deepcopy(sample["compact_atom"])]),
        "source_binding": copy.deepcopy(dict(source)),
        "compact210_aggregate_used": not terminal_zero,
        "compact210_aggregate_retained_as_provenance": True,
        "safe_no_quiet_transport_fallback_at_compact2Rc": not quiet_transport,
        "old_outer_full_202_aggregate_used": False,
        "no_cap_subtraction_value_midpoint_or_overlap": True,
    }
    if terminal_zero:
        budget["compact_upper_aggregate_provenance"] = compact_provenance
    else:
        budget["compact_upper_aggregate"] = compact_aggregate_packet
    return {
        "Z": sample["Z"],
        "source_family": getattr(owner, "family", parent.get("source_family")),
        "raw_m15_C0_C1_terms": copy.deepcopy(raw),
        "correlated_reduced_m15_C0_C1_terms": reduced,
        "original_raw_terms": RAW_COUNT,
        "previous_reduced_terms": PREVIOUS_REDUCED_COUNT,
        "correlated_reduced_terms": NEW_REDUCED_COUNT,
        "unchanged_other_sectors": UNCHANGED_COUNT,
        "signed_trace_chain": trace,
        "source_seam_proofs": copy.deepcopy(dict(proofs)),
        "compact_quiet_packet_witness": copy.deepcopy(packet),
        "external_endpoint_budget_witness": {
            **budget,
            "terminal_D_zero_certificate": copy.deepcopy(dict(terminal_certificate)) if terminal_certificate is not None else None,
            "quiet_transport_certificate": copy.deepcopy(dict(quiet_certificate)) if quiet_certificate is not None else None,
        },
        "old_removed_tail_terms": [copy.deepcopy(old_external), copy.deepcopy(old_compact)],
        "previous_external_replacement_retained_as_provenance": copy.deepcopy(old_external),
        "previous_compact_viscous_retained_as_provenance": copy.deepcopy(old_compact),
        "replacement_tail_endpoints": copy.deepcopy(replacement),
        "independent_target": target,
        "selected_N_evidence": copy.deepcopy(sample["selected_N"]),
        "largest_remaining_individual_log_caps": owner.largest_caps(reduced, target),
        "contract": {
            "raw_213_verbatim": True,
            "previous_190_to_new_189": True,
            "other_188_verbatim": True,
            "six_signed_inner_lower_retained": True,
            "independent_target_unchanged": True,
            "reserved_retained_once": True,
            "quiet_mean_square_retained_once": True,
            "ordinary_Z1_and_physical_R_once": True,
            "leading_subtraction_retained": True,
            "surviving_long0_inlet_not_zeroed": True,
            "endpoint_zero_only_by_source_identity": True,
            "power_compact_source_join_proved": True,
            "compact_quiet_source_join_proved": True,
            "compact210_aggregate_retained": True,
            "compact210_aggregate_used_as_error_only": not terminal_zero,
            "no_error_cap_as_value": True,
            "cap_is_not_signed_function_value": True,
            "all_other_volume_width_numerical_errors_retained": True,
            "all_other_errors_retained": True,
            "terminal_D_zero_certificate_used": terminal_zero,
            "quiet_constant_D_transport_used": quiet_transport,
            "old_outer_full_202_aggregate_used": False,
            "source_functions_not_cap_values": True,
            "ancestors_not_reexecuted": True,
        },
        "limits": {
            "compact_quiet_internal_trace_resolved": True,
            "external_signed_values_resolved": False,
            "independent_numeric_Rp_seam": False,
            "global_five_moment_common_N_closure": False,
            "exact_heat_and_energy": False,
            "admissible_stress_flat_remainder": False,
            "temporal_n_recursion_and_summation": False,
            "independent_corrected_Cartesian_NS": False,
        },
    }


def build_tail_numeric(owner: Any, seam_proofs: Mapping[str, Any]) -> dict[str, Any]:
    """Return the saved compact/quiet 189-term contraction packet."""
    proofs, quiet_certificate, terminal_certificate = _proofs(seam_proofs)
    parent = _parent_report(owner)
    if parent.get("parent_report_sha256") != "3f936cbb5f91fd05677f0c292ed783ad72788b0971ed869298f43be9dfc0850d":
        raise ValueError("Outer-chain parent provenance changed")
    family = getattr(owner, "family", parent.get("source_family"))
    if parent.get("source_family") != family:
        raise ValueError("Tail owner and outer-chain parent source family differ")
    compact_report = _load_json(owner, COMPACT_QUIET_REPORT, COMPACT_QUIET_SHA)
    if compact_report.get("checks", {}).get("passed") is not True:
        raise ValueError("Accepted compact/quiet numerical report must pass")
    geometry_binding = _geometry_source_binding(owner)
    samples: list[dict[str, Any]] = []
    validated: list[dict[str, Any]] = []
    for Z in SAVED_Z:
        parent_sample = _find(parent.get("samples"), label="outer-chain parent sample", Z=Z)
        validated.append(_validate_sample(owner, parent_sample, compact_report, family))
    selected = [sample["selected_N"] for sample in validated]
    if not selected or any(value != selected[0] for value in selected[1:]):
        raise ValueError("All tail samples must use one selected N")
    source = {
        "passed": True,
        "function": _trace_chain(terminal_zero=terminal_certificate is not None,
                                  quiet_transport=quiet_certificate is not None)["function"],
        "fallback_source_function": _trace_chain(terminal_zero=terminal_certificate is not None,
                                                  quiet_transport=quiet_certificate is not None)["function"],
        "delta_definition": "D=R*(Qy_modified-Qy_leading)",
        "physical_radius_power": 1,
        "ordinary_Z_orders": [0, 1],
        "physical_R_once": True,
        "external_endpoints_not_zeroed": terminal_certificate is None,
        "surviving_long0_inlet_not_zeroed": True,
        "endpoint_zero_only_by_source_identity": True,
        "leading_subtraction_retained": True,
        "cap_is_not_signed_function_value": True,
        "all_other_errors_retained": True,
        "source_seam_proofs": copy.deepcopy(proofs),
        "source_proof_sha256": {name: _digest(proof) for name, proof in proofs.items()},
        "compact_quiet_geometry_source_binding": geometry_binding,
        "compact_quiet_report": {"path": COMPACT_QUIET_REPORT, "sha256": COMPACT_QUIET_SHA, "passed": True},
        "parent_report": {"path": PARENT, "sha256": PARENT_SHA, "passed": True},
        "compact_packet_route": {
            "power_compact": "O3_power native2Rc -> compact native1Rc",
            "compact_quiet": "compact native2Rc -> quiet native0 2Rc",
            "compact_lower_coordinate": COMPACT_LOWER,
            "compact_upper_coordinate": COMPACT_UPPER,
            "quiet_lower_fraction": QUIET_LOWER,
            "quiet_Rp_fraction": QUIET_UPPER,
            "quiet_source_coordinate_node": QUIET_SOURCE_COORDINATE_NODE,
            "quiet_source_coordinate": "s=L*q",
            "quiet_lower_source_s": "0",
            "quiet_Rp_source_s": "L",
            "quiet_Rp_fraction_is_not_native_s_one": True,
            "compact_prefix_count": COMPACT_PREFIX_COUNT,
            "compact_complete_count": COMPACT_TERM_COUNT,
            "quiet_prefix_count": QUIET_PREFIX_COUNT,
            "quiet_complete_count": QUIET_TERM_COUNT,
            "compact_viscous_raw_index": 210,
            "quiet_mean_square_raw_index": 211,
            "reserved_raw_index": 212,
        },
        "terminal_D_zero_certificate": copy.deepcopy(dict(terminal_certificate)) if terminal_certificate is not None else None,
        "quiet_transport_certificate": copy.deepcopy(dict(quiet_certificate)) if quiet_certificate is not None else None,
        "terminal_D_zero_used": terminal_certificate is not None,
        "quiet_transport_used": quiet_certificate is not None,
        "safe_no_quiet_transport_fallback_at_compact2Rc": quiet_certificate is None,
        "source_functions_values_or_caps_not_selected": True,
        "old_outer_full_202_aggregate_used": False,
        "accepted_ancestor_constructors_not_run": True,
    }
    for sample in validated:
        samples.append(_evaluate_sample(owner, parent, source, proofs, sample, quiet_certificate, terminal_certificate))
    return {
        "source_family": family,
        "parent_report": PARENT,
        "parent_report_sha256": PARENT_SHA,
        "compact_quiet_report": COMPACT_QUIET_REPORT,
        "compact_quiet_report_sha256": COMPACT_QUIET_SHA,
        "source_binding": source,
        "source_seam_proofs": {
            "passed": True,
            "joins": copy.deepcopy(proofs),
            "quiet_transport_certificate": copy.deepcopy(dict(quiet_certificate)) if quiet_certificate is not None else None,
            "terminal_D_zero_certificate": copy.deepcopy(dict(terminal_certificate)) if terminal_certificate is not None else None,
        },
        "samples": samples,
        "checks": {
            "passed": True,
            "raw_213_verbatim": True,
            "previous_reduced_190": True,
            "new_reduced_189": True,
            "unchanged_other_188": True,
            "all_six_signed_inner_lower_retained": True,
            "reserved_once": True,
            "independent_target_unchanged": True,
            "compact_lower_upper_quiet_full_dictionary_prefixes": len(samples),
            "selected_N_exact_across_parent_compact_quiet": True,
            "compact210_retained": True,
            "compact210_used_as_error_only": terminal_certificate is None,
            "old_outer_full_202_aggregate_used": False,
            "safe_no_quiet_transport_fallback_at_compact2Rc": quiet_certificate is None,
            "terminal_D_zero_certificate_used": terminal_certificate is not None,
            "quiet_transport_certificate_used": quiet_certificate is not None,
            "cap_as_value": False,
        },
        "contract": {
            "raw_213_retained": True,
            "previous_190_to_new_189": True,
            "unchanged_other_188": True,
            "raw_213_and_target_retained": True,
            "parent_190_to_new_189": True,
            "other_188_verbatim": True,
            "six_signed_inner_lower_unchanged": True,
            "reserved_n1_own_swirl_retained_once": True,
            "ordinary_Z1_and_physical_R_once": True,
            "compact210_positive_aggregate_retained": True,
            "surviving_long0_inlet_not_zeroed": True,
            "endpoint_zero_only_by_source_identity": True,
            "quiet211_mean_square_not_viscous": True,
            "power_compact_and_compact_quiet_source_joins_required": True,
            "terminal_D_zero_only_with_explicit_certificate": True,
            "quiet_constant_D_only_with_explicit_certificate": True,
            "no_old_full_202_aggregate": True,
            "no_error_cap_as_value": True,
            "accepted_ancestor_constructors_not_run": True,
        },
        "limits": {
            "external_signed_values_resolved": False,
            "independent_numeric_Rp_seam": False,
            "global_five_moment_common_N_closure": False,
            "exact_heat_and_energy": False,
            "admissible_stress_flat_remainder": False,
            "temporal_n_recursion_and_summation": False,
            "independent_corrected_Cartesian_NS": False,
        },
    }


__all__ = [
    "PARENT",
    "PARENT_SHA",
    "COMPACT_QUIET_REPORT",
    "COMPACT_QUIET_SHA",
    "build_tail_numeric",
]
