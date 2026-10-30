"""Bounded source-only audit of the saved F ``Q^2`` Z0/Z1 pair.

The entry point is intentionally one-shot and self-contained.  It authenticates
all manifest bytes before parsing the two prior packets or decompressing the F
archive, then reads only the two saved Q-squared roots and at most the two
non-protected product children named by the Z1 sum.  Q0, Q1, protected roots,
and the radius handle are opaque endpoints.  No constructor import, source
value selection, quadrature, BFS, or recursive graph reconstruction occurs.
"""

from __future__ import annotations

import copy
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


SCHEMA = "f_gap_qsquare_transport_v1"
F_GRAPH_NAMESPACE = "F_exact_graph"

ROLE_RESULT = Path("tmp/f_pressure_path_integral_roles_v1/result.json")
ROLE_RESULT_SHA256 = "ecbaeda54bdfa82952e09ed1639498ea0bb057716ffc087ab078aa31673c780b"
ROLE_COMPONENT = "f_pressure_path_integral_roles_v1"
ROLE_WRAPPER_SCHEMA = ROLE_COMPONENT + "_result_v1"
ROLE_INNER_SCHEMA = ROLE_COMPONENT
LEADING_RESULT = Path("tmp/pulse_gap_actual_source_program_v1/result.json")
LEADING_RESULT_SHA256 = "e611cdc7e049ce40e6667275b4e9706c3478044880b92c93baebf8ab8d06ebf2"
LEADING_COMPONENT = "pulse_gap_actual_source_program_v1"
F_ARCHIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz"
)
F_ARCHIVE_SHA256 = "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe"

Q_SQUARED_ROOTS = {0: 279381, 1: 279477}
Q_SQUARED_INTEGRANDS = {0: 279382, 1: 279478}
Q_SQUARED_INTEGRALS = {0: 279383, 1: 279479}
SHARED_MEASURE = 34721
PROTECTED_NODE_IDS = frozenset(
    {0, 1, 17, 33572, 25220, 263344, 263348, 263349, 262918, 205147, 263346, 33982}
)
REQUIRED_INPUTS = {
    ROLE_RESULT.as_posix(): ROLE_RESULT_SHA256,
    LEADING_RESULT.as_posix(): LEADING_RESULT_SHA256,
    F_ARCHIVE.as_posix(): F_ARCHIVE_SHA256,
}


class FGapQSquareTransportInputError(ValueError):
    """A frozen manifest or saved source packet is absent or malformed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise FGapQSquareTransportInputError(message)


def _sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _canon(value: Any) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True,
                      separators=(",", ":"), allow_nan=False).encode("utf-8")


def _path(value: Any) -> str:
    _require(isinstance(value, str) and value, "manifest path missing")
    return value.replace("\\", "/")


def _digest(value: Any, label: str) -> str:
    _require(isinstance(value, str) and len(value) == 64, label + " SHA malformed")
    value = value.lower()
    _require(all(char in "0123456789abcdef" for char in value), label + " SHA malformed")
    return value


def _manifest_hashes(source: Mapping[str, Any], hashes: Mapping[str, Any] | None) -> dict[str, str]:
    entries: dict[str, str] = {}

    def add(name: Any, digest: Any, label: str) -> None:
        key, value = _path(name), _digest(digest, label)
        _require(key not in entries or entries[key] == value, "manifest hash conflict: " + key)
        entries[key] = value

    files = source.get("files")
    if isinstance(files, list):
        for index, item in enumerate(files):
            _require(isinstance(item, Mapping), "manifest files entry malformed")
            add(item.get("path"), item.get("sha256"), f"manifest.files[{index}]")
    for field in ("sha256", "input_hashes", "hashes"):
        table = source.get(field)
        if isinstance(table, Mapping):
            for name, digest in table.items():
                add(name, digest, f"manifest.{field}.{name}")
    if isinstance(hashes, Mapping):
        for name, digest in hashes.items():
            add(name, digest, f"caller hashes.{name}")
    _require(entries, "caller frozen manifest has no hashes")
    return entries


def _required_digest(entries: Mapping[str, str], required: str) -> str:
    if required in entries:
        return entries[required]
    key = required.replace("\\", "/")
    matches = [digest for name, digest in entries.items()
               if name.replace("\\", "/").endswith("/" + key)]
    _require(len(matches) == 1, "required manifest hash missing: " + required)
    return matches[0]


def _safe(root: Path, name: str) -> Path:
    relative = Path(name)
    _require(not relative.is_absolute() and ".." not in relative.parts,
             "manifest path escapes root")
    result = (root / relative).resolve()
    _require(result.is_relative_to(root.resolve()), "manifest path escapes root")
    return result


def _canonical_manifest_name(name: str) -> str:
    name = name.replace("\\", "/")
    if name.startswith("/") or ":/" in name:
        aliases = [required for required in REQUIRED_INPUTS
                   if name.endswith("/" + required)]
        _require(len(aliases) == 1, "absolute manifest alias is not required input")
        return aliases[0]
    return name


def _authenticate(root: Path, source: Mapping[str, Any], hashes: Mapping[str, Any] | None,
                  ) -> tuple[dict[str, bytes], dict[str, str], int]:
    """Hash every supplied byte before any JSON parse or gzip decompression."""
    entries = _manifest_hashes(source, hashes)
    for required, expected in REQUIRED_INPUTS.items():
        _require(_required_digest(entries, required) == expected,
                 "required SHA changed: " + required)
    raw: dict[str, bytes] = {}
    actual: dict[str, str] = {}
    for raw_name, expected in entries.items():
        name = _canonical_manifest_name(raw_name)
        data = _safe(root, name).read_bytes()
        observed = _sha(data)
        _require(observed == expected, "raw SHA changed: " + name)
        _require(name not in actual or actual[name] == observed, "manifest alias conflict: " + name)
        raw[name], actual[name] = data, observed
    for required, expected in REQUIRED_INPUTS.items():
        if required in raw:
            continue
        data = _safe(root, required).read_bytes()
        observed = _sha(data)
        _require(observed == expected, "required raw SHA changed: " + required)
        raw[required], actual[required] = data, observed
    return raw, actual, len(entries)


def _json(raw: bytes, label: str) -> Mapping[str, Any]:
    try:
        value = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise FGapQSquareTransportInputError("invalid JSON: " + label) from exc
    _require(isinstance(value, Mapping), label + " is not an object")
    return value


def _gzip_json(raw: bytes, label: str) -> Mapping[str, Any]:
    try:
        decoded = gzip.decompress(raw)
    except (OSError, EOFError) as exc:
        raise FGapQSquareTransportInputError("invalid gzip: " + label) from exc
    return _json(decoded, label)


def _unwrap(packet: Mapping[str, Any], schema: str, label: str) -> Mapping[str, Any]:
    value: Any = packet
    if packet.get("schema") != schema:
        value = packet.get("result")
        _require(isinstance(value, Mapping), label + " result missing")
    _require(value.get("schema") == schema, label + " schema changed")
    return value


def _unwrap_role_result(packet: Mapping[str, Any], label: str) -> Mapping[str, Any]:
    """Unwrap the saved pressure runner envelope, not its outer result object.

    The pressure runner stores both ``synthetic_fixture`` and the actual leaf
    under ``result``.  Only ``result.integral_roles`` is the role packet used
    below.  Every outer discriminator is checked before returning that inner
    packet; a leading-role-shaped wrapper is therefore rejected here.
    """

    _require(packet.get("status") == "returned", label + " wrapper status is not returned")
    _require(packet.get("component") == ROLE_COMPONENT, label + " wrapper component changed")
    _require(packet.get("schema") == ROLE_WRAPPER_SCHEMA, label + " wrapper schema changed")
    if "error" in packet:
        _require(packet.get("error") is None, label + " wrapper contains an error")
    result = packet.get("result")
    _require(isinstance(result, Mapping), label + " wrapper result missing")
    inner = result.get("integral_roles")
    _require(isinstance(inner, Mapping), label + " result.integral_roles missing")
    _require(inner.get("schema") == ROLE_INNER_SCHEMA, label + " inner schema changed")
    _require(inner.get("component") == ROLE_COMPONENT, label + " inner component changed")
    return inner


def _ref(value: Any, role: str, namespace: str = F_GRAPH_NAMESPACE) -> dict[str, Any]:
    valid = type(value) is int and value >= 0
    return {
        "status": "opaque_source_ref" if valid else "unresolved_source_ref",
        "namespace": namespace,
        "node": value if valid else None,
        "source_node_id": value if valid else None,
        "source_role": role,
        "node_id_is_not_a_numeric_value": True,
        "source_value_selected": False,
    }


def _node_from_mapping(value: Any) -> int | None:
    if not isinstance(value, Mapping):
        return None
    for key in ("node", "source_node_id"):
        item = value.get(key)
        if type(item) is int and item >= 0:
            return item
    nested = value.get("source_ref")
    return _node_from_mapping(nested)


def _edge(edge: Any) -> tuple[str | None, int | None, bool, str | None]:
    if not isinstance(edge, Mapping):
        return None, None, False, None
    field, ref = edge.get("field"), edge.get("ref")
    if not isinstance(ref, Mapping):
        return field if isinstance(field, str) else None, None, False, None
    return (
        field if isinstance(field, str) else None,
        _node_from_mapping(ref),
        edge.get("typed_edge") is True,
        ref.get("namespace") if isinstance(ref.get("namespace"), str) else None,
    )


def _has_edge(edges: Any, field: str, node: int) -> bool:
    if not isinstance(edges, list):
        return False
    return any(
        ef == field and en == node and typed and namespace == F_GRAPH_NAMESPACE
        for ef, en, typed, namespace in (_edge(item) for item in edges)
    )


def _prior_roles(prior: Mapping[str, Any]) -> dict[str, Any]:
    """Check only saved role metadata and typed edges; no F graph reads."""
    issues: list[str] = []
    bindings = prior.get("integral_bindings")
    rows_i, rows_f = prior.get("integral_rows"), prior.get("integrand_rows")
    if not isinstance(bindings, list):
        bindings, issues = [], ["integral_bindings_missing"]
    if not isinstance(rows_i, Mapping):
        rows_i, issues = {}, issues + ["integral_rows_missing"]
    if not isinstance(rows_f, Mapping):
        rows_f, issues = {}, issues + ["integrand_rows_missing"]
    records: list[dict[str, Any]] = []
    for order in (0, 1):
        integral, integrand, qroot = Q_SQUARED_INTEGRALS[order], Q_SQUARED_INTEGRANDS[order], Q_SQUARED_ROOTS[order]
        matches = [item for item in bindings if isinstance(item, Mapping) and item.get("root_order") == order]
        local: list[str] = []
        if len(matches) != 1:
            local.append("missing_or_ambiguous_integral_binding")
            binding: Mapping[str, Any] = {}
        else:
            binding = matches[0]
            for key, expected in (("chart", "gap"), ("channel", "I1_q_squared"),
                                  ("formal_integrand", "Q^2"), ("formal_m15_coefficient", "-1/4")):
                if binding.get(key) != expected:
                    local.append(key + "_mismatch")
        metadata = binding.get("integral_metadata")
        metadata_ok = isinstance(metadata, Mapping) and metadata.get("measure") == "(R^2/mu)*dxi" and \
            isinstance(metadata.get("variable"), str) and metadata.get("variable", "").endswith("gap_xi") and \
            metadata.get("ordinary_Z_order") == order and metadata.get("endpoint_independent_of_Z_declared") is True and \
            metadata.get("physical_Pstar_not_reapplied_declared") is True
        if not metadata_ok:
            local.append("measure_or_order_metadata_mismatch")
        ir, fr = rows_i.get(str(order)), rows_f.get(str(order))
        typed_ok = isinstance(ir, Mapping) and ir.get("operation") == "definite_integral" and \
            isinstance(fr, Mapping) and fr.get("operation") == "product"
        if typed_ok:
            typed_ok = _has_edge(ir.get("typed_edges"), "integrand", integrand) and \
                _has_edge(fr.get("typed_edges"), "arguments[0]", SHARED_MEASURE) and \
                _has_edge(fr.get("typed_edges"), "arguments[1]", qroot)
        if not typed_ok:
            local.append("typed_edges_to_measure_or_Q_squared_root_missing")
        records.append({
            "ordinary_Z_order": order,
            "integral_ref": _ref(integral, f"prior.integral_rows[{order}]"),
            "integrand_ref": _ref(integrand, f"prior.integrand_rows[{order}]"),
            "Q_squared_root_ref": _ref(qroot, f"prior.integrand_rows[{order}].arguments[1]"),
            "unique_role_record": len(matches) == 1,
            "metadata_matches": metadata_ok and not local,
            "typed_edges_match": typed_ok,
            "issues": sorted(set(local)),
        })
        issues.extend(f"Z{order}:{item}" for item in local)
    return {
        "status": "matched" if not issues else "unresolved",
        "actual_two_gap_I1_q_squared_roles": not issues,
        "records": records,
        "issues": sorted(set(issues)),
        "source_values_selected": False,
        "integrals_evaluated": False,
    }


def _leading_prior(packet: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    issues: list[str] = []
    sides = packet.get("sides")
    frame = sides.get("F") if isinstance(sides, Mapping) else None
    gap = frame.get("leading_fields_gap") if isinstance(frame, Mapping) else None
    q = gap.get("Q") if isinstance(gap, Mapping) else None
    rows = q.get("ordinary_Z0_to_Z2") if isinstance(q, Mapping) else None
    if not isinstance(q, Mapping) or not isinstance(rows, list) or len(rows) < 2:
        return {
            "status": "unresolved",
            "Q0": _ref(None, "prior.leading.sides.F.leading_fields_gap.Q[0]"),
            "Q1": _ref(None, "prior.leading.sides.F.leading_fields_gap.Q[1]"),
            "source_values_selected": False,
        }, ["F_gap_Q_roles_missing"]
    if q.get("status") != "direct_source_rows":
        issues.append("F_gap_Q_roles_not_direct_source_rows")
    refs = []
    for order in (0, 1):
        item, node = rows[order], _node_from_mapping(rows[order])
        if not isinstance(item, Mapping) or item.get("node_id_is_not_a_numeric_value") is not True:
            issues.append(f"Q{order}:opaque_role_marker_missing")
        if node is None:
            issues.append(f"Q{order}:source_node_id_missing")
        refs.append(_ref(node, f"prior.leading.sides.F.leading_fields_gap.Q[{order}]"))
    if refs[0]["node"] == refs[1]["node"]:
        issues.append("Q0_Q1_same_opaque_handle")
    return {
        "status": "matched" if not issues else "unresolved",
        "Q0": refs[0],
        "Q1": refs[1],
        "source_role": "prior.leading.sides.F.leading_fields_gap.Q",
        "source_values_selected": False,
    }, sorted(set(issues))


def _archive_q_roles(archive: Mapping[str, Any], prior: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    issues: list[str] = []
    bindings = archive.get("Rp_native_C5_pulse_source_bindings")
    leading = bindings.get("leading_fields") if isinstance(bindings, Mapping) else None
    gap = leading.get("gap") if isinstance(leading, Mapping) else None
    qrows = gap.get("Q") if isinstance(gap, Mapping) else None
    if not isinstance(qrows, list) or len(qrows) < 2:
        return {
            "status": "unresolved",
            "Q0": _ref(None, "F_archive.leading_fields.gap.Q[0]"),
            "Q1": _ref(None, "F_archive.leading_fields.gap.Q[1]"),
            "source_values_selected": False,
        }, ["F_archive_gap_Q_roles_missing"]
    refs = [_ref(qrows[order], f"F_archive.leading_fields.gap.Q[{order}]") for order in (0, 1)]
    for order, item in enumerate(refs):
        if item["status"] != "opaque_source_ref":
            issues.append(f"archive_Q{order}:opaque_ref_malformed")
        if item["node"] in PROTECTED_NODE_IDS:
            issues.append(f"archive_Q{order}:protected_Q_ref")
    p0, p1 = prior.get("Q0"), prior.get("Q1")
    correspondence = isinstance(p0, Mapping) and isinstance(p1, Mapping) and \
        p0.get("status") == "opaque_source_ref" and p1.get("status") == "opaque_source_ref" and \
        p0.get("namespace") == F_GRAPH_NAMESPACE and p1.get("namespace") == F_GRAPH_NAMESPACE and \
        p0.get("node") == refs[0]["node"] and p1.get("node") == refs[1]["node"]
    if not correspondence:
        issues.append("prior_to_archive_Q_role_ref_mismatch")
    return {
        "status": "matched" if not issues else "unresolved",
        "Q0": refs[0],
        "Q1": refs[1],
        "same_F_snapshot_role_ref_correspondence": correspondence,
        "source_values_selected": False,
    }, sorted(set(issues))


def _direct_row(nodes: Any, node: Any, role: str, read_ids: list[int]) -> dict[str, Any] | None:
    if type(node) is not int or node < 0 or node in PROTECTED_NODE_IDS:
        return None
    _require(isinstance(nodes, list), "F exact_graph_nodes missing")
    _require(node < len(nodes), "direct F row outside archive")
    row = nodes[node]
    _require(isinstance(row, Mapping), "direct F row malformed")
    read_ids.append(node)
    copied = copy.deepcopy(dict(row))
    return {
        "node_ref": _ref(node, role),
        "operation": copied.get("operation"),
        "row": copied,
        "row_sha256": _sha(_canon(copied)),
        "children_not_indexed": True,
        "protected_roots_not_indexed": True,
    }


def _args(record: Mapping[str, Any] | None) -> list[Any] | None:
    row = record.get("row") if isinstance(record, Mapping) else None
    args = row.get("arguments") if isinstance(row, Mapping) else None
    return args if isinstance(args, list) else None


def _orientation(args: Any, q0: int, q1: int) -> str | None:
    if not isinstance(args, list) or len(args) != 2:
        return None
    if args == [q1, q0]:
        return "Q1_Q0"
    if args == [q0, q1]:
        return "Q0_Q1"
    return None


def _shape_audit(nodes: Any, q0: int | None, q1: int | None) -> dict[str, Any]:
    """Read roots 279381/279477 and at most their two direct Z1 product rows."""
    read_ids: list[int] = []
    reasons: list[str] = []
    roots = {
        str(order): _direct_row(nodes, Q_SQUARED_ROOTS[order], f"F.Q_squared_root[{order}]", read_ids)
        for order in (0, 1)
    }
    if q0 is None or q1 is None or q0 in PROTECTED_NODE_IDS or q1 in PROTECTED_NODE_IDS:
        return {
            "status": "unresolved",
            "actual_row_shape_matched": False,
            "shape_kind": "unresolved",
            "root_rows": roots,
            "Z1_product_children": [],
            "direct_graph_row_ids": sorted(set(read_ids)),
            "direct_graph_row_count": len(set(read_ids)),
            "maximum_indexed_rows": 4,
            "reasons": ["Q0_or_Q1_role_unresolved_or_protected"],
            "Q_refs_are_endpoints_only": True,
            "protected_roots_not_indexed": True,
            "radius_target_not_indexed": True,
        }
    r0, r1 = roots.get("0"), roots.get("1")
    r0row, r1row = r0.get("row") if isinstance(r0, Mapping) else None, r1.get("row") if isinstance(r1, Mapping) else None
    r0ok = isinstance(r0row, Mapping) and r0row.get("operation") == "product" and _args(r0) == [q0, q0]
    if not r0ok:
        reasons.append("Z0:expected_product_Q0_Q0_mismatch")
    r1args = _args(r1)
    r1ok = isinstance(r1row, Mapping) and r1row.get("operation") == "sum"
    if not r1ok:
        reasons.append("Z1:expected_sum_root_mismatch")
    child_ids: list[int] = []
    if r1ok and not isinstance(r1args, list):
        reasons.append("Z1:missing_or_ambiguous_product_children")
    elif r1ok and len(r1args) != 2:
        reasons.append("Z1:missing_or_ambiguous_product_children")
    elif r1ok and any(type(item) is not int for item in r1args):
        reasons.append("Z1:typed_child_reference_malformed")
    elif r1ok:
        child_ids = list(r1args)
        if len(set(child_ids)) != 2:
            reasons.append("Z1:missing_or_ambiguous_product_children")
        if any(item in PROTECTED_NODE_IDS or item in (q0, q1) for item in child_ids):
            reasons.append("Z1:child_is_protected_or_Q_endpoint")
    can_read = r1ok and len(child_ids) == 2 and len(set(child_ids)) == 2 and \
        all(type(item) is int and item not in PROTECTED_NODE_IDS and item not in (q0, q1) for item in child_ids)
    children: list[dict[str, Any]] = []
    if can_read:
        for index, child_id in enumerate(child_ids):
            child = _direct_row(nodes, child_id, f"F.Q_squared_Z1_product_child[{index}]", read_ids)
            if child is None:
                reasons.append(f"Z1:child_{index}_unreadable")
            else:
                children.append(child)
        if len(children) == 2:
            orientations = [_orientation(_args(child), q0, q1) if child.get("operation") == "product" else None for child in children]
            if any(child.get("operation") != "product" for child in children):
                reasons.append("Z1:nonproduct_child")
            if (
                orientations.count("Q0_Q1") != 1
                or orientations.count("Q1_Q0") != 1
                or any(item not in ("Q0_Q1", "Q1_Q0") for item in orientations)
            ):
                reasons.append("Z1:product_children_do_not_match_Q1Q0_plus_Q0Q1")
        else:
            reasons.append("Z1:missing_or_ambiguous_product_children")
    mismatch = any(item in reasons for item in ("Z0:expected_product_Q0_Q0_mismatch", "Z1:expected_sum_root_mismatch", "Z1:nonproduct_child", "Z1:product_children_do_not_match_Q1Q0_plus_Q0Q1"))
    unresolved = any("missing_or_ambiguous" in item or "unreadable" in item or "malformed" in item or "protected" in item for item in reasons)
    matched = r0ok and len(children) == 2 and not reasons
    return {
        "status": "matched" if matched else "unresolved" if unresolved else "mismatch",
        "actual_row_shape_matched": matched,
        "shape_kind": "matched" if matched else "unresolved" if unresolved else "mismatch",
        "root_rows": roots,
        "Z1_product_children": children,
        "direct_graph_row_ids": sorted(set(read_ids)),
        "direct_graph_row_count": len(set(read_ids)),
        "maximum_indexed_rows": 4,
        "reasons": sorted(set(reasons)),
        "structural_mismatch": mismatch,
        "Q_refs_are_endpoints_only": True,
        "protected_roots_not_indexed": True,
        "radius_target_not_indexed": True,
    }


def _symbol(name: str) -> dict[str, Any]:
    return {"kind": "formal_symbol", "name": name}


def _literal(value: str) -> dict[str, Any]:
    return {"kind": "exact_scalar_literal", "value": value}


def _expr(operator: str, *arguments: Any, **metadata: Any) -> dict[str, Any]:
    result = {"operator": operator, "arguments": list(arguments)}
    result.update(metadata)
    return result


def _formal_candidate(q0_ref: Mapping[str, Any], q1_ref: Mapping[str, Any]) -> dict[str, Any]:
    q0: Any = dict(q0_ref) if q0_ref.get("status") == "opaque_source_ref" else _symbol("Q0")
    q1: Any = dict(q1_ref) if q1_ref.get("status") == "opaque_source_ref" else _symbol("Q1")
    xi, mu, mu_z = _symbol("xi"), _symbol("mu_F"), _symbol("mu_F_Z")
    log_rp, log_rp_z = _symbol("logRp_F"), _symbol("logRp_F_Z")
    xi_l, xi_u = _symbol("xi_l(Z)"), _symbol("xi_u(Z)")
    xi_l_z, xi_u_z = _symbol("xi_l_Z"), _symbol("xi_u_Z")
    q_partial = _symbol("partial_Z_Q0_fixed_xi")
    d_q = _expr("subtract", q_partial, q1, candidate="D_Q=(partial_Z Q0)_xi-Q1", source_bound=False)
    e_q = _expr("subtract", q1, q_partial, candidate="E_Q=Q1-(partial_Z Q0)_xi=-D_Q", source_bound=False)
    radius = _expr("exp", _expr("add", log_rp, _expr("divide", xi, mu)), candidate="R=exp(logRp_F+xi/mu_F)")
    measure = _expr("divide", _expr("power", radius, _literal("2")), mu, candidate="M=R^2/mu_F")
    q2 = _expr("power", q0, _literal("2"))
    qcross = _expr("add", _expr("multiply", q1, q0), _expr("multiply", q0, q1), candidate="Q1Q0+Q0Q1")
    qpartialcross = _expr("add", _expr("multiply", q_partial, q0), _expr("multiply", q0, q_partial))
    i0_integrand, i1_integrand = _expr("multiply", measure, q2), _expr("multiply", measure, qcross)
    i0 = _expr("definite_integral", i0_integrand, xi_l, xi_u, variable="xi", candidate="I0=integral M Q0^2 dxi")
    i1 = _expr("definite_integral", i1_integrand, xi_l, xi_u, variable="xi", candidate="I1=integral M(Q1Q0+Q0Q1) dxi")
    mz = _expr(
        "multiply", measure,
        _expr("add", _expr("multiply", _literal("2"), log_rp_z),
               _expr("negate", _expr("divide", _expr("multiply", _literal("2"), xi, mu_z), _expr("power", mu, _literal("2")))),
               _expr("negate", _expr("divide", mu_z, mu))),
        candidate="M_Z=M*(2logRp_Z-2xi*mu_Z/mu^2-mu_Z/mu)", source_bound=False,
    )
    f_u, f_l = _expr("evaluate_at", i0_integrand, xi_u), _expr("evaluate_at", i0_integrand, xi_l)
    bulk = _expr("definite_integral", _expr("add", _expr("multiply", mz, q2), _expr("multiply", measure, _expr("add", _expr("multiply", d_q, q0), _expr("multiply", q0, d_q)))), xi_l, xi_u, variable="xi", candidate="integral[M_Z Q0^2+M(D_Q Q0+Q0 D_Q)] dxi")
    discrepancy = _expr("multiply", _literal("-1/4"), _expr("add", bulk, _expr("multiply", f_u, xi_u_z), _expr("negate", _expr("multiply", f_l, xi_l_z))), candidate="d_Z(-I0/4)-(-I1/4); D_Q retained", formal_conditional=True)
    d_form = _expr(
        "add",
        _expr("multiply", _literal("-1/4"), _expr("definite_integral", _expr("multiply", mz, q2), xi_l, xi_u, variable="xi")),
        _expr("multiply", _literal("-1/2"), _expr("definite_integral", _expr("multiply", measure, _expr("multiply", q0, d_q)), xi_l, xi_u, variable="xi")),
        _expr("multiply", _literal("-1/4"), _expr("multiply", f_u, xi_u_z)),
        _expr("multiply", _literal("1/4"), _expr("multiply", f_l, xi_l_z)),
        candidate="commutative D_Q form; D_Q not set to zero", formal_conditional=True,
    )
    e_form = _expr(
        "add",
        _expr("multiply", _literal("-1/4"), _expr("definite_integral", _expr("multiply", mz, q2), xi_l, xi_u, variable="xi")),
        _expr("multiply", _literal("1/2"), _expr("definite_integral", _expr("multiply", measure, _expr("multiply", q0, e_q)), xi_l, xi_u, variable="xi")),
        _expr("multiply", _literal("-1/4"), _expr("multiply", f_u, xi_u_z)),
        _expr("multiply", _literal("1/4"), _expr("multiply", f_l, xi_l_z)),
        candidate="E_Q form; formal conditional only", formal_conditional=True,
    )
    return {
        "status": "conditional_formal_expression_not_evaluated",
        "Q0": q0, "Q1": q1, "D_Q": d_q, "E_Q": e_q,
        "radius_candidate": radius, "measure_candidate": measure, "M_Z_candidate": mz,
        "I0_candidate": i0, "I1_candidate": i1,
        "m15_Z0_candidate": _expr("multiply", _literal("-1/4"), i0),
        "m15_Z1_saved_candidate": _expr("multiply", _literal("-1/4"), i1),
        "missing_correction_candidate": discrepancy,
        "missing_correction_commutative_D_Q_form": d_form,
        "missing_correction_E_Q_form": e_form,
        "endpoint_candidates": {
            "xi_l": _literal("11"), "xi_u": _expr("subtract", _literal("13"), _expr("multiply", _literal("4"), mu)),
            "xi_l_Z_candidate": _literal("0"), "xi_u_Z_candidate": _expr("negate", _expr("multiply", _literal("4"), mu_z)),
            "xi_l_Z": xi_l_z, "xi_u_Z": xi_u_z,
            "lower_Z_zero_source_bound": False, "upper_Z_source_bound": False,
            "moving_upper_bound_retained": True,
        },
        "conditions": [
            "Q1=partial_Z Q0 at fixed xi only after independently proving D_Q=0",
            "lower and upper bounds are differentiable source-owned functions",
            "I1 is exactly the saved ordinary-Z source row for I0",
            "same source-owned M, coordinate, and capture-free derivative tree",
            "endpoint owners and fixed-xi partials are source-authenticated",
        ],
        "Q_derivative_defect_zero_proved": False,
        "mu_Z_source_bound": False, "logRp_Z_source_bound": False,
        "M_Z_zero_assumed": False, "upper_Z_zero_assumed": False, "lower_Z_zero_assumed": False,
        "source_values_selected": False, "integrals_evaluated": False,
        "numeric_use_admissible": False, "actual_force_or_residual": False,
    }


def _synthetic_rows(q0: int, q1: int) -> dict[int, dict[str, Any]]:
    return {301: {"operation": "product", "arguments": [q0, q0]},
            302: {"operation": "sum", "arguments": [401, 402]},
            401: {"operation": "product", "arguments": [q1, q0]},
            402: {"operation": "product", "arguments": [q0, q1]}}


def _synthetic_shape(rows: Mapping[int, Mapping[str, Any]], q0: int | None, q1: int | None) -> dict[str, Any]:
    if q0 is None or q1 is None:
        return {"status": "unresolved", "actual_row_shape_matched": False}
    r0, r1 = rows.get(301), rows.get(302)
    if not isinstance(r0, Mapping) or not isinstance(r1, Mapping):
        return {"status": "unresolved", "actual_row_shape_matched": False}
    if r0.get("operation") != "product" or r0.get("arguments") != [q0, q0]:
        return {"status": "mismatch", "actual_row_shape_matched": False}
    args = r1.get("arguments")
    if r1.get("operation") != "sum" or not isinstance(args, list) or len(args) != 2 or len(set(args)) != 2:
        return {"status": "unresolved", "actual_row_shape_matched": False}
    children = [rows.get(item) for item in args]
    if any(not isinstance(item, Mapping) for item in children):
        return {"status": "unresolved", "actual_row_shape_matched": False}
    orientations = [_orientation(item.get("arguments"), q0, q1) if item.get("operation") == "product" else None for item in children]
    matched = (
        orientations.count("Q0_Q1") == 1
        and orientations.count("Q1_Q0") == 1
        and all(item in ("Q0_Q1", "Q1_Q0") for item in orientations)
    )
    return {"status": "matched" if matched else "mismatch", "actual_row_shape_matched": matched,
            "orientations": orientations, "source_values_selected": False}


def _synthetic_fixtures() -> list[dict[str, Any]]:
    q0, q1 = 501, 502
    correct, wrong, nonproduct = _synthetic_rows(q0, q1), _synthetic_rows(q0, q1), _synthetic_rows(q0, q1)
    wrong[401] = {"operation": "product", "arguments": [q0, q0]}
    nonproduct[401] = {"operation": "sum", "arguments": [q1, q0]}
    missing, ambiguous = _synthetic_rows(q0, q1), _synthetic_rows(q0, q1)
    missing[302], ambiguous[302] = {"operation": "sum", "arguments": [401]}, {"operation": "sum", "arguments": [401, 401]}
    return [
        {"name": "correct_pair", "kind": "structural_fixture", "result": _synthetic_shape(correct, q0, q1), "mathematical_identity_proved": False},
        {"name": "Q1_wrong_replacement", "kind": "structural_fixture", "result": _synthetic_shape(wrong, q0, q1), "mathematical_identity_proved": False},
        {"name": "nonproduct_child", "kind": "structural_fixture", "result": _synthetic_shape(nonproduct, q0, q1), "mathematical_identity_proved": False},
        {"name": "missing_or_ambiguous_child", "kind": "structural_fixture", "variants": {"missing": _synthetic_shape(missing, q0, q1), "ambiguous": _synthetic_shape(ambiguous, q0, q1)}, "mathematical_identity_proved": False},
        {"name": "unbound_role", "kind": "role_fixture", "result": _synthetic_shape(correct, None, q1), "mathematical_identity_proved": False},
        {"name": "moving_endpoint_nonzero", "kind": "formal_boundary_fixture", "result": {"status": "conditional_only", "xi_u_Z": "-4mu_F_Z", "mu_F_Z_nonzero": True, "boundary_term_retained": True, "ordinary_Z_transport_proved": False}, "mathematical_identity_proved": False},
        {"name": "M_Z_nonzero", "kind": "formal_measure_fixture", "result": {"status": "conditional_only", "M_Z_zero_assumed": False, "M_Z_candidate_retained": True, "ordinary_Z_transport_proved": False}, "mathematical_identity_proved": False},
    ]


def _synthetic_wrapper_checks() -> dict[str, bool]:
    """Exercise the two saved-result envelope shapes without source reads."""

    role_inner = {
        "schema": ROLE_INNER_SCHEMA,
        "component": ROLE_COMPONENT,
        "status": "bounded_F_integral_roles_partial",
    }
    role_wrapper = {
        "schema": ROLE_WRAPPER_SCHEMA,
        "component": ROLE_COMPONENT,
        "status": "returned",
        "error": None,
        "result": {
            "synthetic_fixture": {"synthetic_only": True},
            "integral_roles": role_inner,
        },
    }
    leading_inner = {
        "schema": LEADING_COMPONENT,
        "component": LEADING_COMPONENT,
        "status": "source_program_returned",
    }
    leading_wrapper = {
        "schema": LEADING_COMPONENT + "_result_v1",
        "component": LEADING_COMPONENT,
        "status": "returned",
        "result": leading_inner,
    }

    def role_accepts(packet: Mapping[str, Any]) -> bool:
        try:
            value = _unwrap_role_result(packet, "synthetic role wrapper")
        except FGapQSquareTransportInputError:
            return False
        return value.get("schema") == ROLE_INNER_SCHEMA

    def leading_accepts(packet: Mapping[str, Any]) -> bool:
        try:
            value = _unwrap(packet, LEADING_COMPONENT, "synthetic leading wrapper")
        except FGapQSquareTransportInputError:
            return False
        return value.get("schema") == LEADING_COMPONENT

    bad_role = copy.deepcopy(role_wrapper)
    bad_role["result"]["integral_roles"]["schema"] = "wrong_role_schema"
    bad_leading = copy.deepcopy(leading_wrapper)
    bad_leading["result"]["schema"] = "wrong_leading_schema"
    return {
        "role_result_integral_roles_shape": role_accepts(role_wrapper),
        "leading_result_direct_shape": leading_accepts(leading_wrapper),
        "role_wrapper_mismatch_fail_closed": not role_accepts(bad_role),
        "leading_wrapper_mismatch_fail_closed": not leading_accepts(bad_leading),
    }


def verify_synthetic_fixture() -> dict[str, Any]:
    """Check the seven bounded synthetic guards without reading any source."""

    fixtures = {item["name"]: item for item in _synthetic_fixtures()}
    checks: dict[str, Any] = {}

    def check(name: str, observed: Any, expected: Any, description: str) -> None:
        checks[name] = {
            "description": description,
            "observed": observed,
            "expected": expected,
            "passed": observed == expected,
        }

    correct = fixtures["correct_pair"]["result"]
    check(
        "correct_pair",
        correct.get("actual_row_shape_matched"),
        True,
        "Q0*Q0 and the two oppositely ordered Q1/Q0 products match",
    )

    wrong = fixtures["Q1_wrong_replacement"]["result"]
    check(
        "Q1_wrong_replacement",
        wrong.get("actual_row_shape_matched"),
        False,
        "replacing Q1 with Q0 must not match the derivative shape",
    )

    nonproduct = fixtures["nonproduct_child"]["result"]
    check(
        "nonproduct_child",
        nonproduct.get("actual_row_shape_matched"),
        False,
        "a non-product Z1 child must not match",
    )

    child_variants = fixtures["missing_or_ambiguous_child"]["variants"]
    missing = child_variants["missing"]
    ambiguous = child_variants["ambiguous"]
    check(
        "missing_child",
        missing.get("actual_row_shape_matched"),
        False,
        "a missing Z1 product child remains unresolved",
    )
    check(
        "ambiguous_child",
        ambiguous.get("actual_row_shape_matched"),
        False,
        "duplicate Z1 child references remain unresolved",
    )

    unbound = fixtures["unbound_role"]["result"]
    check(
        "unbound_role",
        unbound.get("actual_row_shape_matched"),
        False,
        "unbound Q roles cannot establish a structural match",
    )

    moving = fixtures["moving_endpoint_nonzero"]["result"]
    check(
        "moving_endpoint_nonzero",
        (
            moving.get("boundary_term_retained") is True
            and moving.get("mu_F_Z_nonzero") is True
            and moving.get("ordinary_Z_transport_proved") is False
        ),
        True,
        "a nonzero moving-endpoint derivative remains in the formal candidate",
    )

    measure = fixtures["M_Z_nonzero"]["result"]
    check(
        "M_Z_nonzero",
        (
            measure.get("M_Z_candidate_retained") is True
            and measure.get("M_Z_zero_assumed") is False
            and measure.get("ordinary_Z_transport_proved") is False
        ),
        True,
        "the M_Z correction remains formal and is not set to zero",
    )

    wrapper_checks = _synthetic_wrapper_checks()
    check(
        "role_result_integral_roles_shape",
        wrapper_checks["role_result_integral_roles_shape"],
        True,
        "ROLE_RESULT must unwrap result.integral_roles rather than outer result",
    )
    check(
        "leading_result_direct_shape",
        wrapper_checks["leading_result_direct_shape"],
        True,
        "LEADING_RESULT keeps its direct result wrapper shape",
    )
    check(
        "role_wrapper_mismatch_fail_closed",
        wrapper_checks["role_wrapper_mismatch_fail_closed"],
        True,
        "a mismatched nested role schema must remain rejected",
    )
    check(
        "leading_wrapper_mismatch_fail_closed",
        wrapper_checks["leading_wrapper_mismatch_fail_closed"],
        True,
        "a mismatched leading result schema must remain rejected",
    )

    passed = all(item["passed"] for item in checks.values())
    return {
        "schema": SCHEMA + "_synthetic_fixture_v1",
        "component": SCHEMA,
        "status": "synthetic_checks_passed" if passed else "synthetic_checks_failed",
        "passed": passed,
        "checks": checks,
        "check_count": len(checks),
        "fixture_count": len(fixtures),
        "wrapper_fixture_checks": wrapper_checks,
        "source_values_computed": False,
        "source_values_selected": False,
        "source_graph_read": False,
        "archive_decoded": False,
        "integrals_evaluated": False,
        "numeric_use_admissible": False,
        "actual_force_or_residual": False,
        "scope": {
            "synthetic_only": True,
            "source_only": True,
            "no_source_numeric_values": True,
            "no_prior_result_read": True,
            "no_F_archive_read": True,
            "no_graph_traversal": True,
            "ordinary_Z_transport_proved": False,
            "Q_derivative_defect_zero_proved": False,
            "capture_free_substitution_proved": False,
            "Leibniz_moving_boundary_proved": False,
            "M_Z_zero_proved": False,
        },
    }


def build_qsquare_transport(root: str | Path = ".", source: Mapping[str, Any] | None = None,
                            hashes: Mapping[str, Any] | None = None) -> dict[str, Any]:
    """Build one bounded source-role packet for the F Q-squared pair."""
    _require(isinstance(source, Mapping), "caller frozen manifest is required")
    raw, actual, manifest_count = _authenticate(Path(root).resolve(), source, hashes)
    # No packet is parsed until every manifest byte has passed its raw SHA check.
    role = _unwrap_role_result(_json(raw[ROLE_RESULT.as_posix()], ROLE_RESULT.as_posix()), "pressure role prior")
    leading = _unwrap(_json(raw[LEADING_RESULT.as_posix()], LEADING_RESULT.as_posix()), LEADING_COMPONENT, "F leading role prior")
    archive = _gzip_json(raw[F_ARCHIVE.as_posix()], F_ARCHIVE.as_posix())
    prior = _prior_roles(role)
    q_prior, q_prior_issues = _leading_prior(leading)
    q_archive, q_archive_issues = _archive_q_roles(archive, q_prior)
    q0 = q_archive.get("Q0", {}).get("node") if isinstance(q_archive.get("Q0"), Mapping) else None
    q1 = q_archive.get("Q1", {}).get("node") if isinstance(q_archive.get("Q1"), Mapping) else None
    nodes = archive.get("exact_graph_nodes")
    _require(isinstance(nodes, list), "F exact_graph_nodes missing")
    shape = _shape_audit(nodes, q0, q1)
    formal = _formal_candidate(q_archive.get("Q0", _ref(None, "F.Q0")), q_archive.get("Q1", _ref(None, "F.Q1")))
    correspondence = bool(q_archive.get("same_F_snapshot_role_ref_correspondence") is True and not q_prior_issues and not q_archive_issues)
    status = "bounded_F_qsquare_structural_match_conditional" if shape["actual_row_shape_matched"] and correspondence and not prior["issues"] else "bounded_F_qsquare_structural_unresolved"
    return {
        "schema": SCHEMA, "component": SCHEMA, "status": status, "source_only": True,
        "runner_interface": {"entry": "build_qsquare_transport", "arguments": ["root", "source", "hashes"], "one_shot": True, "coordination": "/root/bound_mu"},
        "raw_hash_authentication": {"all_manifest_raw_hashes_checked_before_decode": True, "manifest_entry_count": manifest_count, "required_inputs": {name: {"path": name, "expected_sha256": digest, "actual_sha256": actual[name]} for name, digest in REQUIRED_INPUTS.items()}},
        "prior_scope": {"pressure_roles": prior, "leading_Q_prior": q_prior, "leading_Q_prior_issues": q_prior_issues, "role_result_path": ROLE_RESULT.as_posix(), "leading_result_path": LEADING_RESULT.as_posix()},
        "F_snapshot": {"archive_path": F_ARCHIVE.as_posix(), "leading_Q_roles": q_archive, "leading_Q_issues": q_archive_issues, "namespace": F_GRAPH_NAMESPACE, "whole_json_decoded": True, "BFS_or_recursive_graph_reconstruction": False},
        "direct_graph_scope": {"requested_Q_squared_root_ids": [279381, 279477], "indexed_rows": shape["direct_graph_row_ids"], "indexed_row_count": shape["direct_graph_row_count"], "maximum_indexed_rows": 4, "Q_refs_only_endpoints": True, "protected_root_ids": sorted(PROTECTED_NODE_IDS), "protected_roots_indexed": False, "protected_nodes_indexed": False, "radius_target_33982_indexed": False, "target_node_dereferenced": False, "children_not_recursed": True, "graph_edges_followed": False, "source_graph_traversed": False},
        "actual_row_shape_matched": shape["actual_row_shape_matched"], "same_F_snapshot_role_ref_correspondence": correspondence, "shape_audit": shape,
        "formal_symbolic_candidate": formal,
        "source_provenance": {
            "saved_role_consumer": {"path": "src/openai_ns_reconstruction/f_pressure_path_integral_roles.py", "role": "saved gap I1_q_squared Z0/Z1 metadata and typed-edge contract"},
            "Q_squared_source_formula": {"path": "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_physical_nonradial_pressure_fields.py", "lines": "259-279", "role": "_z1_mul and I1_q_squared source construction"},
            "F_rewrite_entry": {"path": "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.py", "lines": "98-102", "role": "saved leading_fields/pressure rewrite boundary"},
            "identity_scope": "operation-tree shape only; IDs, row hashes, and metadata are not function semantics",
        },
        "synthetic_fixtures": _synthetic_fixtures(),
        "source_gates": {"actual_row_shape_matched": bool(shape["actual_row_shape_matched"]), "same_F_snapshot_role_ref_correspondence": correspondence, "ordinary_Z_derivative_proved": False, "ordinary_Z_transport_proved": False, "actual_ordinary_Z_transport_proved": False, "Q_derivative_defect_zero_proved": False, "capture_free_substitution_proved": False, "F_radius_owner_bound": False, "F_field_radius_owner_bound": False, "F_index_semantics_proved": False, "index_semantics_proved": False, "common_measure_proved": False, "measure_owner_proved": False, "endpoint_transport_proved": False, "Leibniz_moving_boundary_proved": False, "source_values_selected": False, "source_values_computed": False, "integrals_evaluated": False, "numeric_use_admissible": False, "m15_closure_proved": False, "error_transport_proved": False},
        "open_gates": ["ordinary_Z_derivative_proved", "Q_derivative_defect_zero_proved", "capture_free_substitution_proved", "F_radius_owner_bound", "F_index_semantics_proved", "common_measure_proved", "endpoint_transport_proved", "Leibniz_moving_boundary_proved", "source_values_selected", "integrals_evaluated", "numeric_use_admissible", "m15_closure_proved", "error_transport_proved"],
        "retained_errors": {"old_errors_removed": False, "prior_ledgers_and_corrections_reinterpreted": False, "formal_missing_correction_deleted": False, "numeric_values_computed": False, "source_values_selected": False, "actual_force_or_residual": False},
    }


__all__ = ["FGapQSquareTransportInputError", "build_qsquare_transport", "verify_synthetic_fixture"]
