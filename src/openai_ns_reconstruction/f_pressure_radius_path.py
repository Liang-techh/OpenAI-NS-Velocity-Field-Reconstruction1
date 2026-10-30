"""Bounded structural F-pressure-to-radius path query.

This leaf consumes one authenticated F archive and only its saved manifest.
It reads the two saved ``pressure.m15`` roots at ordinary Z orders zero and
one, then traverses the archive's exact graph with a frozen operand schema.
The saved ``rewritten_nodes['33982']`` value is an opaque radius boundary:
the boundary row is never dereferenced.

The result is structural existence evidence only.  A path to the boundary is
not a numerical radius value, a proof of mathematical dependence, a proof of
capture-free substitution, or a proof that the boundary is the radius used by
the current pulse-coordinate metadata.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from collections import deque
from pathlib import Path
from typing import Any, Iterable, Mapping


SCHEMA = "f_pressure_radius_path_v1"
COMPONENT = SCHEMA
F_ARCHIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz"
)
F_ARCHIVE_SHA256 = "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe"

PRESSURE_BINDING_PATH = "$.Rp_native_C5_pulse_source_bindings.pressure.m15"
REWRITTEN_RADIUS_KEY = "33982"
F_GRAPH_NAMESPACE = "F_exact_graph"
EXPECTED_PRESSURE_M15_ROWS = 4
TRAVERSED_PRESSURE_ROWS = (0, 1)
MAX_NODES_PER_ROOT = 4096
MAX_DEPTH = 128

# These are saved handles and are never indexed.  A protected handle is a
# terminal frontier, not an operator row in the bounded query.
PROTECTED_NODE_IDS = frozenset(
    {
        0,
        1,
        17,
        33572,
        25220,
        263344,
        263348,
        263349,
        262918,
        205147,
        263346,
    }
)

# Frozen from the saved source producer's reference schema.  Only these
# fields are interpreted as graph operands; all other row metadata is opaque.
FROZEN_OPERAND_SCHEMA: dict[str, tuple[str, ...]] = {
    "mathematical_positive_infinity": (),
    "exact_original_source_heat_phi": ("argument",),
    "exact_rational": (),
    "bound_variable": (),
    "mathematical_pi": (),
    "current_original_source_parameter": (),
    "original_source_parameter": (),
    "exact_positive_integer_expression": (),
    "shared_positive_integer": (),
    "sum": ("arguments",),
    "product": ("arguments",),
    "negative": ("argument",),
    "analytic_unary": ("argument",),
    "positive_quotient": ("numerator", "denominator"),
    "exact_real_comparison": ("left", "right"),
    "logical_or": ("arguments",),
    "compact_raw_beta": ("argument",),
    "exact_original_source_flat_sigma": ("argument",),
    "current_original_leading_function_recipe": ("coordinate", "Z_variable"),
    "current_original_physical_source_scale": (),
    "current_original_physical_source_quantity": ("coordinate", "Z_variable"),
    "current_original_physical_radius": ("log_radius", "free_query"),
    "substitute_original_native_coordinate": ("body", "query_coordinate"),
    "substitute_free_query_at_native_coordinate": ("body", "free_query"),
    "substitute_original_inverse_angle": ("body", "inverse_angle"),
    "function_substitution": ("expression", "variable", "value"),
    "definite_integral": ("integrand", "lower", "upper"),
    "original_lazy_flat_branch": (
        "active_body",
        "Delta",
        "eta",
        "flat_predicate",
        "original_collar_coordinate",
        "original_collar_upper",
    ),
    "original_monotone_phase_inverse": (
        "phase_function",
        "target_phase",
        "angle_lower",
        "angle_upper",
        "Delta",
        "eta",
    ),
}


class FPressureRadiusPathInputError(ValueError):
    """A manifest, archive, or direct F root is malformed or changed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise FPressureRadiusPathInputError(message)


def _sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _path(value: Any) -> str:
    _require(isinstance(value, str) and value, "manifest path missing")
    return value.replace("\\", "/")


def _digest(value: Any, label: str) -> str:
    _require(isinstance(value, str) and len(value) == 64, label + " SHA malformed")
    result = value.lower()
    _require(
        all(character in "0123456789abcdef" for character in result),
        label + " SHA malformed",
    )
    return result


def _safe(root: Path, name: str) -> Path:
    relative = Path(name)
    _require(
        not relative.is_absolute() and ".." not in relative.parts,
        "manifest path escapes root",
    )
    result = (root / relative).resolve()
    _require(result.is_relative_to(root.resolve()), "manifest path escapes root")
    return result


def _manifest_hashes(
    source: Mapping[str, Any], hashes: Mapping[str, str] | None
) -> dict[str, str]:
    entries: dict[str, str] = {}

    def add(name_value: Any, digest_value: Any, label: str) -> None:
        name = _path(name_value)
        digest = _digest(digest_value, label)
        _require(
            name not in entries or entries[name] == digest,
            "manifest hash conflict: " + name,
        )
        entries[name] = digest

    files = source.get("files")
    if isinstance(files, list):
        for index, item in enumerate(files):
            _require(isinstance(item, Mapping), "manifest files entry malformed")
            add(item.get("path"), item.get("sha256"), f"manifest.files[{index}]")

    for field in ("sha256", "input_hashes", "hashes"):
        table = source.get(field)
        if isinstance(table, Mapping):
            for name, digest in table.items():
                add(name, digest, "manifest." + field)

    if isinstance(hashes, Mapping):
        for name, digest in hashes.items():
            add(name, digest, "caller hashes")

    _require(entries, "caller frozen manifest has no hashes")
    return entries


def _find_hash(entries: Mapping[str, str], name: str) -> str:
    if name in entries:
        return entries[name]
    normalized = name.replace("\\", "/")
    candidates = [
        value
        for key, value in entries.items()
        if key.replace("\\", "/").endswith("/" + normalized)
    ]
    _require(len(candidates) == 1, "required manifest hash missing: " + name)
    return candidates[0]


def _canonical_manifest_name(raw_name: str) -> str:
    """Map an absolute alias only when it names the required F archive."""

    name = raw_name.replace("\\", "/")
    if name.startswith("/") or ":/" in name:
        aliases = [
            item
            for item in (F_ARCHIVE.as_posix(),)
            if name.endswith("/" + item)
        ]
        _require(len(aliases) == 1, "absolute manifest path is not the F archive")
        return aliases[0]
    return name


def _authenticate(
    root: Path, source: Mapping[str, Any], hashes: Mapping[str, str] | None
) -> tuple[bytes, dict[str, str], int]:
    """Hash every manifest file before the F archive is decompressed."""

    entries = _manifest_hashes(source, hashes)
    _require(
        _find_hash(entries, F_ARCHIVE.as_posix()) == F_ARCHIVE_SHA256,
        "required F archive SHA changed",
    )

    actual: dict[str, str] = {}
    archive_raw: bytes | None = None
    for raw_name, expected in entries.items():
        name = _canonical_manifest_name(raw_name)
        data = _safe(root, name).read_bytes()
        observed = _sha(data)
        _require(observed == expected, "raw SHA changed: " + name)
        actual[name] = observed
        if name == F_ARCHIVE.as_posix():
            archive_raw = data

    # A manifest may identify the required input by a unique suffix rather
    # than its canonical path.  Authenticate/read it once more only if the
    # first pass did not retain the required bytes.
    if archive_raw is None:
        canonical_expected = _find_hash(entries, F_ARCHIVE.as_posix())
        data = _safe(root, F_ARCHIVE.as_posix()).read_bytes()
        observed = _sha(data)
        _require(observed == canonical_expected == F_ARCHIVE_SHA256, "required F archive SHA changed")
        archive_raw = data
        actual[F_ARCHIVE.as_posix()] = observed

    return archive_raw, actual, len(entries)


def _decode_archive(raw: bytes) -> Mapping[str, Any]:
    try:
        decoded = gzip.decompress(raw)
        packet = json.loads(decoded.decode("utf-8"))
    except (OSError, EOFError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise FPressureRadiusPathInputError("F archive is not valid gzip JSON") from exc
    _require(isinstance(packet, Mapping), "F archive is not a JSON object")
    return packet


def _node_id(value: Any) -> bool:
    return type(value) is int and value >= 0


def _opaque_ref(node: Any, role: str | None = None) -> dict[str, Any]:
    result: dict[str, Any] = {
        "namespace": F_GRAPH_NAMESPACE,
        "node": node,
        "node_id_is_not_a_numeric_value": True,
    }
    if role is not None:
        result["source_role"] = role
    return result


def _extract_direct_inputs(
    packet: Mapping[str, Any],
) -> tuple[list[int], int, list[Mapping[str, Any]]]:
    bindings = packet.get("Rp_native_C5_pulse_source_bindings")
    _require(isinstance(bindings, Mapping), "F source bindings missing")

    pressure = bindings.get("pressure")
    _require(isinstance(pressure, Mapping), "F pressure bindings missing")
    m15 = pressure.get("m15")
    _require(
        isinstance(m15, list) and len(m15) == EXPECTED_PRESSURE_M15_ROWS,
        "F pressure.m15 must contain exactly four C3 rows",
    )
    roots: list[int] = []
    for order, value in enumerate(m15):
        _require(
            _node_id(value),
            f"F pressure.m15[{order}] is not an opaque integer reference",
        )
        roots.append(value)

    rewritten = bindings.get("rewritten_nodes")
    _require(isinstance(rewritten, Mapping), "F rewritten_nodes mapping missing")
    _require(REWRITTEN_RADIUS_KEY in rewritten, "F rewritten_nodes radius key missing")
    target = rewritten.get(REWRITTEN_RADIUS_KEY)
    _require(_node_id(target), "F rewritten radius target is not an opaque integer reference")

    nodes = packet.get("exact_graph_nodes")
    _require(isinstance(nodes, list), "F exact_graph_nodes is not a list")
    return roots, target, nodes


def _step(
    node: int,
    operator: Any,
    operand_label: str | None,
    depth: int,
    dereferenced: bool,
) -> dict[str, Any]:
    return {
        "node_ref": _opaque_ref(node),
        "operator": operator if isinstance(operator, str) else None,
        "operand_label": operand_label,
        "depth": depth,
        "dereferenced": dereferenced,
    }


def _operand_items(value: Any, field: str) -> Iterable[tuple[str, Any]]:
    if isinstance(value, list):
        for index, item in enumerate(value):
            yield f"{field}[{index}]", item
    else:
        yield field, value


def _walk_root(
    root: int,
    target: int,
    nodes: list[Any],
    *,
    max_nodes: int = MAX_NODES_PER_ROOT,
    max_depth: int = MAX_DEPTH,
    role: str = "F.pressure.m15",
    row_order: int | None = None,
) -> dict[str, Any]:
    record: dict[str, Any] = {
        "role": role,
        "row_order": row_order,
        "root_ref": _opaque_ref(root, f"{PRESSURE_BINDING_PATH}[{row_order}]" if row_order is not None else None),
        "target_ref": _opaque_ref(target, "saved rewritten_nodes['33982'] target"),
        "existence": False,
        "status": "unresolved",
        "stopped_on_target": False,
        "complete_traversal": False,
        "path": [],
        "visited_node_count": 0,
        "scheduled_node_count": 0,
        "visited_node_ids": [],
        "edges_examined": 0,
        "edge_counts": {
            "operand_edges": 0,
            "target_edges": 0,
            "protected_terminations": 0,
            "shared_revisits": 0,
            "queued_edges": 0,
            "unresolved_edges": 0,
        },
        "unresolved": [],
        "protected_node_ids": sorted(PROTECTED_NODE_IDS),
        "protected_boundary_seen": False,
        "node_budget": max_nodes,
        "depth_budget": max_depth,
        "target_row_dereferenced": False,
    }

    def unresolved(code: str, **detail: Any) -> None:
        item = {"code": code}
        item.update(detail)
        record["unresolved"].append(item)
        record["edge_counts"]["unresolved_edges"] += 1

    if not _node_id(root):
        unresolved("root_ref_malformed", root=root)
        record["status"] = "unresolved"
        return record
    if root == target:
        record["existence"] = True
        record["status"] = "hit"
        record["stopped_on_target"] = True
        record["path"] = [_step(target, None, "<root>", 0, False)]
        record["target_row_dereferenced"] = False
        return record
    if root in PROTECTED_NODE_IDS:
        record["status"] = "unresolved_protected_boundary"
        record["protected_boundary_seen"] = True
        record["edge_counts"]["protected_terminations"] += 1
        return record
    if root >= len(nodes):
        unresolved("root_ref_missing", node=root, node_count=len(nodes))
        return record

    queue: deque[tuple[int, list[dict[str, Any]], frozenset[int], int, str | None]] = deque()
    queue.append((root, [], frozenset({root}), 0, None))
    scheduled: set[int] = {root}
    visited: set[int] = set()

    while queue:
        node, prefix, active, depth, via = queue.popleft()
        if node == target:
            # This should normally be handled at edge creation, but keeping
            # this guard makes the target terminal even for a root queue item.
            record["existence"] = True
            record["status"] = "hit"
            record["stopped_on_target"] = True
            record["path"] = prefix + [_step(target, None, via or "<root>", depth, False)]
            break
        if node in PROTECTED_NODE_IDS:
            record["protected_boundary_seen"] = True
            record["edge_counts"]["protected_terminations"] += 1
            continue
        if node in visited:
            record["edge_counts"]["shared_revisits"] += 1
            continue
        if node >= len(nodes):
            unresolved("missing_node", node=node, node_count=len(nodes))
            continue
        row = nodes[node]
        visited.add(node)
        record["visited_node_count"] = len(visited)
        record["visited_node_ids"] = sorted(visited)
        if not isinstance(row, Mapping):
            unresolved("node_row_malformed", node=node)
            continue
        operation = row.get("operation")
        current_path = prefix + [_step(node, operation, via, depth, True)]
        if not isinstance(operation, str) or operation not in FROZEN_OPERAND_SCHEMA:
            unresolved("unknown_operation", node=node, operation=operation)
            continue

        for field in FROZEN_OPERAND_SCHEMA[operation]:
            if field not in row or row.get(field) is None:
                continue
            for operand_label, operand in _operand_items(row.get(field), field):
                record["edges_examined"] += 1
                record["edge_counts"]["operand_edges"] += 1
                if not _node_id(operand):
                    unresolved(
                        "undeclared_or_malformed_ref",
                        node=node,
                        operand=operand_label,
                        value_kind=type(operand).__name__,
                    )
                    continue
                next_depth = depth + 1
                if next_depth > max_depth:
                    unresolved(
                        "depth_budget_exhausted",
                        node=node,
                        operand=operand_label,
                        target=operand,
                        depth=next_depth,
                    )
                    continue
                if operand == target:
                    record["edge_counts"]["target_edges"] += 1
                    record["existence"] = True
                    record["status"] = "hit"
                    record["stopped_on_target"] = True
                    record["path"] = current_path + [
                        _step(operand, None, operand_label, next_depth, False)
                    ]
                    record["target_row_dereferenced"] = False
                    break
                if operand in PROTECTED_NODE_IDS:
                    record["protected_boundary_seen"] = True
                    record["edge_counts"]["protected_terminations"] += 1
                    continue
                if operand in active:
                    unresolved(
                        "cycle",
                        node=node,
                        operand=operand_label,
                        target=operand,
                    )
                    continue
                if operand in scheduled:
                    record["edge_counts"]["shared_revisits"] += 1
                    continue
                if len(scheduled) >= max_nodes:
                    unresolved(
                        "node_budget_exhausted",
                        node=node,
                        operand=operand_label,
                        target=operand,
                        scheduled_node_count=len(scheduled),
                    )
                    continue
                if operand >= len(nodes):
                    unresolved(
                        "missing_node",
                        node=node,
                        operand=operand_label,
                        target=operand,
                        node_count=len(nodes),
                    )
                    continue
                scheduled.add(operand)
                record["edge_counts"]["queued_edges"] += 1
                queue.append(
                    (
                        operand,
                        current_path,
                        active | {node},
                        next_depth,
                        operand_label,
                    )
                )
            if record["stopped_on_target"]:
                break
        if record["stopped_on_target"]:
            break

    record["scheduled_node_count"] = len(scheduled)
    if record["stopped_on_target"]:
        # A witness remains valid even if another branch was intentionally
        # blocked at a protected boundary.  Keep existence=True but expose
        # that the whole root was not structurally closed.
        record["status"] = (
            "hit_with_protected_boundary"
            if record["protected_boundary_seen"]
            else "hit"
        )
        record["complete_traversal"] = False
    elif record["protected_boundary_seen"]:
        # Protected rows are never indexed.  Their presence prevents a
        # closed no-path conclusion, even when no other unresolved edge was
        # recorded.
        record["status"] = "unresolved_protected_boundary"
        record["complete_traversal"] = False
    elif record["unresolved"]:
        record["status"] = "unresolved"
        record["complete_traversal"] = False
    else:
        record["status"] = "no_path"
        record["complete_traversal"] = True
    return record


def _archive_result(
    packet: Mapping[str, Any],
    *,
    archive_sha256: str,
    manifest_entry_count: int,
    actual_hashes: Mapping[str, str],
) -> dict[str, Any]:
    roots, target, nodes = _extract_direct_inputs(packet)
    root_results = [
        _walk_root(
            roots[order],
            target,
            nodes,
            role="F.pressure.m15",
            row_order=order,
        )
        for order in TRAVERSED_PRESSURE_ROWS
    ]
    any_hit = any(item["existence"] for item in root_results)
    all_hit = all(item["existence"] for item in root_results)
    unresolved_count = sum(len(item["unresolved"]) for item in root_results)
    blocked_roots = [
        item["row_order"]
        for item in root_results
        if item.get("protected_boundary_seen") is True
        or item.get("status") == "unresolved_protected_boundary"
        or item.get("status") == "hit_with_protected_boundary"
    ]
    protected_boundary_seen = bool(blocked_roots)
    all_roots_closed_no_path = bool(root_results) and all(
        item.get("status") == "no_path"
        and item.get("complete_traversal") is True
        for item in root_results
    )
    all_roots_clean_hit = bool(root_results) and all(
        item.get("existence") is True
        and item.get("status") == "hit"
        and item.get("protected_boundary_seen") is False
        and item.get("complete_traversal") is False
        for item in root_results
    ) and unresolved_count == 0
    status = (
        "bounded_F_pressure_radius_path_found"
        if all_hit and all_roots_clean_hit
        else "bounded_F_pressure_radius_path_partial"
        if any_hit
        else "bounded_F_pressure_radius_path_no_path"
        if all_roots_closed_no_path
        else "bounded_F_pressure_radius_path_unresolved"
    )
    return {
        "schema": SCHEMA,
        "component": COMPONENT,
        "status": status,
        "source_only": True,
        "existence": any_hit,
        "all_traversed_roots_exist": all_hit,
        "root_count": len(root_results),
        "pressure_m15_row_count": len(roots),
        "traversed_pressure_rows": list(TRAVERSED_PRESSURE_ROWS),
        "roots": root_results,
        "blocked_root_count": len(blocked_roots),
        "blocked_roots": blocked_roots,
        "protected_boundary_seen": protected_boundary_seen,
        "all_roots_closed_no_path": all_roots_closed_no_path,
        "target": {
            "source_role": "$.Rp_native_C5_pulse_source_bindings.rewritten_nodes['33982']",
            "old_mapping_key": REWRITTEN_RADIUS_KEY,
            "target_ref": _opaque_ref(target),
            "from_saved_mapping": True,
            "operator_row_dereferenced": False,
        },
        "archive_scope": {
            "archive_path": F_ARCHIVE.as_posix(),
            "archive_sha256": archive_sha256,
            "exact_graph_nodes_count": len(nodes),
            "F_archive_exact_graph_traversed": True,
            "bounded_BFS": True,
            "max_nodes_per_root": MAX_NODES_PER_ROOT,
            "max_depth": MAX_DEPTH,
            "protected_node_ids": sorted(PROTECTED_NODE_IDS),
            "protected_nodes_indexed": False,
            "target_node_dereferenced": False,
            "only_pressure_m15_Z0_Z1_roots_traversed": True,
            "prior_20_root_query_reentered": False,
        },
        "scope": {
            "source_graph_traversed": True,
            "F_archive_exact_graph_traversed": True,
            "source_values_computed": False,
            "source_values_selected": False,
            "numeric_values_computed": False,
            "numeric_use_admissible": False,
            "integrals_evaluated": False,
            "old_producer_replayed": False,
            "ancestor_replayed": False,
            "fresh_graph_generated": False,
            "prior_20_root_query_reentered": False,
            "P_F_G_compared": False,
            "target_row_dereferenced": False,
            "protected_root_indexed": False,
        },
        "source_gates": {
            "F_pressure_radius_structural_path_exists": any_hit,
            "F_field_radius_owner_bound": False,
            "actual_field_radius_identity_proved": False,
            "mathematical_nonzero_dependence_proved": False,
            "capture_free_substitution_proved": False,
            "ordinary_Z_transport_proved": False,
            "five_moment_closure_proved": False,
            "m15_integral_evaluated": False,
            "P_equals_F_proved": False,
            "G_equals_F_proved": False,
            "numeric_use_admissible": False,
        },
        "raw_hash_authentication": {
            "all_manifest_raw_hashes_checked_before_decompression": True,
            "manifest_entry_count": manifest_entry_count,
            "required_archive": {
                "path": F_ARCHIVE.as_posix(),
                "expected_sha256": F_ARCHIVE_SHA256,
                "actual_sha256": actual_hashes.get(F_ARCHIVE.as_posix(), archive_sha256),
                "decompressed_after_raw_hash": True,
            },
        },
    }


def build_f_pressure_radius_path(
    root: str | Path = ".",
    source: Mapping[str, Any] | None = None,
    hashes: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Authenticate a one-shot manifest, then traverse only F pressure roots."""

    _require(isinstance(source, Mapping), "caller frozen manifest is required")
    repo = Path(root).resolve()
    archive_raw, actual, manifest_entry_count = _authenticate(repo, source, hashes)
    packet = _decode_archive(archive_raw)
    return _archive_result(
        packet,
        archive_sha256=actual.get(F_ARCHIVE.as_posix(), F_ARCHIVE_SHA256),
        manifest_entry_count=manifest_entry_count,
        actual_hashes=actual,
    )


def _synthetic_nodes_hit() -> list[dict[str, Any]]:
    return [
        {"operation": "exact_rational"},
        {"operation": "exact_rational"},
        {"operation": "analytic_unary", "argument": 9},
    ]


def verify_synthetic_fixture() -> dict[str, Any]:
    """Run in-memory path/negative controls without reading repository input."""

    # IDs zero and one are protected in the production archive.  Keep them as
    # inert list padding only; every synthetic root and traversed node starts
    # at two or above.
    def padded(*rows: dict[str, Any]) -> list[dict[str, Any]]:
        return [
            {"operation": "exact_rational"},
            {"operation": "exact_rational"},
            *rows,
        ]

    hit = _walk_root(
        2,
        9,
        padded({"operation": "analytic_unary", "argument": 9}),
        role="synthetic.hit",
        row_order=0,
    )

    unknown_boundary = _walk_root(
        2,
        9,
        padded({"operation": "analytic_unary", "argument": "external_ref"}),
        role="synthetic.unknown_boundary",
        row_order=0,
    )

    protected = _walk_root(
        2,
        9,
        padded({"operation": "analytic_unary", "argument": 17}),
        role="synthetic.protected",
        row_order=0,
    )

    unknown_operation = _walk_root(
        2,
        9,
        padded({"operation": "unknown_operation", "argument": 3}),
        role="synthetic.unknown_operation",
        row_order=0,
    )

    cycle = _walk_root(
        2,
        9,
        padded(
            {"operation": "analytic_unary", "argument": 3},
            {"operation": "analytic_unary", "argument": 2},
        ),
        role="synthetic.cycle",
        row_order=0,
    )

    budget = _walk_root(
        2,
        9,
        padded(
            {"operation": "analytic_unary", "argument": 3},
            {"operation": "analytic_unary", "argument": 4},
            {"operation": "exact_rational"},
        ),
        max_nodes=2,
        role="synthetic.budget",
        row_order=0,
    )

    protected_packet = {
        "Rp_native_C5_pulse_source_bindings": {
            "pressure": {"m15": [2, 2, 2, 2]},
            "rewritten_nodes": {REWRITTEN_RADIUS_KEY: 9},
        },
        "exact_graph_nodes": padded({"operation": "analytic_unary", "argument": 17}),
    }
    protected_top_level = _archive_result(
        protected_packet,
        archive_sha256=F_ARCHIVE_SHA256,
        manifest_entry_count=1,
        actual_hashes={F_ARCHIVE.as_posix(): F_ARCHIVE_SHA256},
    )

    checks = {
        "hit_records_complete_path_without_target_row": (
            hit["existence"] is True
            and hit["status"] == "hit"
            and hit["target_row_dereferenced"] is False
            and hit["path"][-1]["dereferenced"] is False
        ),
        "unknown_boundary_is_unresolved_not_no_path": (
            unknown_boundary["existence"] is False
            and unknown_boundary["status"] == "unresolved"
            and any(item["code"] == "undeclared_or_malformed_ref" for item in unknown_boundary["unresolved"])
        ),
        "protected_boundary_is_terminal_without_indexing": (
            protected["existence"] is False
            and protected["status"] == "unresolved_protected_boundary"
            and protected["complete_traversal"] is False
            and protected["protected_boundary_seen"] is True
            and protected["edge_counts"]["protected_terminations"] == 1
            and protected["unresolved"] == []
        ),
        "protected_boundary_prevents_no_path": (
            protected["status"] != "no_path"
            and protected["complete_traversal"] is False
        ),
        "unknown_operation_is_unresolved": (
            unknown_operation["existence"] is False
            and unknown_operation["status"] == "unresolved"
            and any(item["code"] == "unknown_operation" for item in unknown_operation["unresolved"])
        ),
        "cycle_is_unresolved": (
            cycle["existence"] is False
            and cycle["status"] == "unresolved"
            and any(item["code"] == "cycle" for item in cycle["unresolved"])
        ),
        "node_budget_is_unresolved": (
            budget["existence"] is False
            and budget["status"] == "unresolved"
            and any(item["code"] == "node_budget_exhausted" for item in budget["unresolved"])
        ),
        "top_level_protected_boundary_prevents_no_path": (
            protected_top_level["existence"] is False
            and protected_top_level["status"] == "bounded_F_pressure_radius_path_unresolved"
            and protected_top_level["protected_boundary_seen"] is True
            and protected_top_level["blocked_root_count"] == 2
            and protected_top_level["all_roots_closed_no_path"] is False
        ),
    }
    return {
        "schema": SCHEMA + "_synthetic_fixture_v1",
        "component": COMPONENT,
        "status": "synthetic_negative_controls_passed"
        if all(checks.values())
        else "synthetic_negative_controls_failed",
        "synthetic_only": True,
        "actual_source_read": False,
        "F_archive_exact_graph_traversed": False,
        "source_graph_traversed": False,
        "source_values_computed": False,
        "numeric_use_admissible": False,
        "checks": checks,
        "cases": {
            "hit": hit,
            "unknown_boundary": unknown_boundary,
            "protected": protected,
            "unknown_operation": unknown_operation,
            "cycle": cycle,
            "budget": budget,
            "protected_top_level": protected_top_level,
        },
    }


__all__ = [
    "COMPONENT",
    "F_ARCHIVE",
    "F_ARCHIVE_SHA256",
    "F_GRAPH_NAMESPACE",
    "FROZEN_OPERAND_SCHEMA",
    "FPressureRadiusPathInputError",
    "MAX_DEPTH",
    "MAX_NODES_PER_ROOT",
    "PROTECTED_NODE_IDS",
    "SCHEMA",
    "build_f_pressure_radius_path",
    "verify_synthetic_fixture",
]
