"""Archive-only source theorem for the postrestore/Rm inlet.

The accepted Rm/outer source receipt already proves the physical-germ
calculus used at the ``postrestore -> Rm_patch`` interface.  This adapter
does not replay that checker or identify arbitrary graph nodes.  It binds the
receipt to the four actual physical ``U_y`` ports supplied by the current
postrestore and Rm snapshots, while leaving controls, phases, histories and
all pressure outputs as functions.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .original_core_numeric import ROOT


_JOIN_ARCHIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rm_outer_mixed_source_joins.json.gz"
)
_JOIN_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rm_outer_mixed_source_joins.py"
)
_JOIN_CHECK_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rm_outer_mixed_source_joins_check.py"
)
_JOIN_RECEIPT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rm_outer_mixed_source_joins_check.json"
)

_SOURCE_PINS = {
    "accepted_join_archive": (_JOIN_ARCHIVE, "27510c39da6e84293c2eaa2355ce36d1c5622c39839c33a3b907739b85580a8c"),
    "accepted_join_source": (_JOIN_SOURCE, "9ec220509152dfa406f14b661fbc76fa55c95aee6908ab5e0344519e445653f8"),
    "accepted_join_checker": (_JOIN_CHECK_SOURCE, "15bd090dc0d891bfd7825730f8412830abeed451f06ede8e1c66acde200c6846"),
    "accepted_join_receipt": (_JOIN_RECEIPT, "14d1cc6b641435c60a7e3662c996a260c07b328df5b8f486c4344f3d7eeb4867"),
    "rm_finite_source": (
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_Rm_patch_finite_N_C5_functions.py",
        "c4445ef8acfb33d29cf9022862f8b6396e6aa2f38bda7467409b5378a5deaf15",
    ),
    "loop_source": (
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py",
        "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d",
    ),
    "modified_source": (
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py",
        "7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7",
    ),
}

_BINDING_KEYS = {
    "rm": "Rm_patch_finite_N_C5_bindings",
    "leading": "true_Rm_patch_C6_bindings",
    "loop": "prescribed_loop_C5_bindings",
    "modified": "modified_Q4_pressure_C3_bindings",
}


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _unwrap(value: Any, key: str) -> Mapping[str, Any]:
    value = _mapping(value)
    nested = value.get(key)
    return nested if isinstance(nested, Mapping) else value


def _repo_root(root: Any) -> Path:
    path = Path(root) if root is not None else ROOT
    if path.is_file():
        path = path.parent
    if (path / "src").is_dir() and (path / "experiments").is_dir():
        return path
    return ROOT


def _handle(value: Any) -> int | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, str) and value.isdigit():
        return int(value)
    return None


def _nodes_len(nodes: Any) -> int:
    try:
        return len(nodes)
    except TypeError:
        return 0


def _node(nodes: Any, value: Any) -> Mapping[str, Any] | None:
    index = _handle(value)
    if index is None or index < 0:
        return None
    try:
        row = nodes[index] if not isinstance(nodes, Mapping) else nodes.get(index, nodes.get(str(index)))
    except (IndexError, KeyError, TypeError):
        return None
    return row if isinstance(row, Mapping) else None


def _canonical(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _digest(value: Any) -> str:
    return hashlib.sha256(_canonical(value).encode("utf-8")).hexdigest()


def _pin(repo: Path, relative: str, expected: str) -> dict[str, Any]:
    path = repo / relative
    actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
    return {
        "path": relative,
        "exists": path.is_file(),
        "sha256": actual,
        "expected_sha256": expected,
        "passed": actual == expected,
    }


def _load_receipt(repo: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    path = repo / _JOIN_RECEIPT
    if not path.is_file():
        return {}, {"passed": False, "reason": "accepted Rm outer receipt is missing"}
    try:
        receipt = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        return {}, {"passed": False, "reason": "accepted Rm outer receipt is unreadable", "error": str(exc)}
    return receipt, {"passed": True, "path": _JOIN_RECEIPT}


def _port_sequence(value: Any, *, minimum: int = 2) -> bool:
    return isinstance(value, (list, tuple)) and len(value) >= minimum and all(
        _handle(item) is not None for item in value[:minimum]
    )


def _chart(binding: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    charts = _mapping(binding.get("charts"))
    return _mapping(charts.get(name))


def _physical_ports(
    nodes: Any,
    modified: Mapping[str, Any],
    rm: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    post = _chart(modified, "postrestore")
    left_physical = _mapping(post.get("physical"))
    right_physical = _mapping(rm.get("physical"))
    left_boundary = _mapping(post.get("finite_boundary"))
    right_boundary = _mapping(rm.get("boundary"))
    names = ("U", "U_y", "m", "m_y", "Q", "Q_y")
    ports: dict[str, Any] = {"left": {}, "right": {}}
    checks: dict[str, Any] = {}
    for side, group, boundary in (
        ("left", left_physical, left_boundary),
        ("right", right_physical, right_boundary),
    ):
        for name in names:
            values = group.get(name)
            ports[side][name] = list(values) if isinstance(values, (list, tuple)) else []
        ports[side]["R_Q_y"] = list(boundary.get("R_Q_y", ())) if isinstance(boundary.get("R_Q_y"), (list, tuple)) else []
    for name, count in (("U", 3), ("m", 3), ("m_y", 3), ("U_y", 2), ("Q", 2), ("Q_y", 2), ("R_Q_y", 2)):
        left = ports["left"][name]
        right = ports["right"][name]
        left_ok = len(left) >= count and all(_node(nodes, item) is not None for item in left[:count])
        right_ok = len(right) >= count and all(_node(nodes, item) is not None for item in right[:count])
        checks[name] = {
            "passed": bool(left_ok and right_ok),
            "ordinary_Z_orders": list(range(count)),
            "left_nodes": left[:count],
            "right_nodes": right[:count],
            "left_nodes_exist": left_ok,
            "right_nodes_exist": right_ok,
        }
    return ports, checks


def _physical_shape(value: Any) -> tuple[int, int] | None:
    if not isinstance(value, list) or len(value) != 5 or not all(isinstance(row, list) for row in value):
        return None
    widths = {len(row) for row in value}
    if len(widths) != 1 or next(iter(widths), 0) not in (5, 6):
        return None
    return (len(value), next(iter(widths)))


def _load_join_archive(repo: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    """Hydrate only the accepted source archive; never import or run its producer."""
    path = repo / _JOIN_ARCHIVE
    try:
        with gzip.open(path, "rt", encoding="utf-8") as stream:
            archive = json.load(stream)
    except (OSError, ValueError) as exc:
        return {}, {"passed": False, "reason": "accepted Rm outer archive is unreadable", "error": str(exc)}
    return archive if isinstance(archive, Mapping) else {}, {"passed": isinstance(archive, Mapping), "path": _JOIN_ARCHIVE}


def _join_binding_checks(archive: Mapping[str, Any], receipt: Mapping[str, Any]) -> dict[str, Any]:
    binding = _mapping(archive.get("Rm_outer_mixed_source_join_bindings"))
    charts = _mapping(binding.get("charts"))
    joins = _mapping(binding.get("joins"))
    post = _mapping(charts.get("postrestore"))
    patch = _mapping(charts.get("Rm_patch"))
    post_physical = _mapping(post.get("physical"))
    patch_physical = _mapping(patch.get("physical"))
    post_ok = all(isinstance(post_physical.get(name), list) and len(post_physical[name]) >= 2 for name in ("U", "m", "Q"))
    patch_ok = all(isinstance(patch_physical.get(name), list) and len(patch_physical[name]) >= 2 for name in ("U", "m", "Q"))
    seam = _mapping(joins.get("postrestore_patch"))
    seam_left = _mapping(seam.get("left")); seam_right = _mapping(seam.get("right"))
    seam_chart_names = bool(
        seam.get("left_chart") == "postrestore"
        and seam.get("right_chart") == "Rm_patch"
    )
    # The join archive has a separate producer output for each endpoint.  It
    # must identify both chart names and carry the same five-by-six physical
    # output layout as the chart bindings; node-count checks alone would let
    # an unrelated serialized list masquerade as this seam.
    chart_shapes = {
        "postrestore": {
            name: _physical_shape(post_physical.get(name)) for name in ("U", "m", "Q")
        },
        "Rm_patch": {
            name: _physical_shape(patch_physical.get(name)) for name in ("U", "m", "Q")
        },
    }
    seam_shapes = {
        "left": {
            name: _physical_shape(_mapping(seam_left.get("physical")).get(name))
            for name in ("U", "m", "Q")
        },
        "right": {
            name: _physical_shape(_mapping(seam_right.get("physical")).get(name))
            for name in ("U", "m", "Q")
        },
    }
    expected_shapes = {"U": (5, 6), "m": (5, 6), "Q": (5, 5)}
    producer_outputs = bool(
        seam_chart_names
        and all(chart_shapes[side][name] == expected_shapes[name] for side in ("postrestore", "Rm_patch") for name in expected_shapes)
        and all(seam_shapes[side][name] == expected_shapes[name] for side in ("left", "right") for name in expected_shapes)
        and _handle(seam_left.get("native_query")) is not None
        and _handle(seam_right.get("native_query")) is not None
        and isinstance(seam_left.get("radius"), list)
        and isinstance(seam_right.get("radius"), list)
    )
    seam_physical = bool(
        all(isinstance(_mapping(seam_left.get("physical")).get(name), list) for name in ("U", "m", "Q"))
        and all(isinstance(_mapping(seam_right.get("physical")).get(name), list) for name in ("U", "m", "Q"))
    )
    receipt_family = _mapping(receipt.get("source_family"))
    archive_family = _mapping(archive.get("source_family"))
    family_match = bool(receipt_family and archive_family and receipt_family == archive_family)
    return {
        "passed": bool(
            archive.get("all_passed", True) is True
            and set(("postrestore", "Rm_patch")) <= set(charts)
            and post_ok and patch_ok and "postrestore_patch" in joins and seam_physical
            and producer_outputs and family_match
        ),
        "archive_source_family_matches_receipt": family_match,
        "join_chart_names_match_actual_interface": seam_chart_names,
        "join_producer_physical_output_definitions": producer_outputs,
        "postrestore_chart_has_actual_U_m_Q_outputs": post_ok,
        "Rm_patch_chart_has_actual_U_m_Q_outputs": patch_ok,
        "postrestore_patch_join_has_left_right_physical_outputs": seam_physical,
        "chart_physical_shapes": chart_shapes,
        "join_physical_shapes": seam_shapes,
        "charts": sorted(charts),
        "joins": sorted(joins),
    }


def _source_prefix_rows(
    nodes: Any,
    bindings: Mapping[str, Mapping[str, Any]],
    join_archive: Mapping[str, Any],
) -> dict[str, Any]:
    """Require an exact serialized prefix relation to the pinned join graph."""
    available = _nodes_len(nodes)
    joined = join_archive.get("exact_graph_nodes")
    current_prefix_join = bool(
        isinstance(joined, list)
        and available <= len(joined)
        and list(joined[:available]) == list(nodes)
    )
    join_prefix_current = bool(
        isinstance(joined, list)
        and len(joined) <= available
        and list(nodes[:len(joined)]) == list(joined)
    )
    exact_prefix_relation = current_prefix_join or join_prefix_current
    required = {"modified": 183700, "loop": 121834, "leading": 196042, "rm": 205147}
    rows = {}
    for name, minimum in required.items():
        supplied = bindings.get(name, {})
        prefix = supplied.get("exact_graph_nodes")
        if isinstance(prefix, (list, tuple)):
            length = len(prefix)
            exact = available >= length and list(nodes[:length]) == list(prefix)
            mode = "serialized_prefix_comparison"
        else:
            length = minimum
            exact = available >= minimum and exact_prefix_relation
            mode = "pinned_join_archive_prefix_comparison"
        rows[name] = {
            "passed": bool(exact),
            "mode": mode,
            "required_prefix_length": length,
            "available_nodes": available,
            "serialized_prefix_supplied": isinstance(prefix, (list, tuple)),
        }
    rows["join_archive_prefix"] = {
        "passed": exact_prefix_relation,
        "current_prefix_join_archive": current_prefix_join,
        "join_prefix_current_graph": join_prefix_current,
        "current_node_count": available,
        "join_node_count": len(joined) if isinstance(joined, list) else 0,
    }
    return {"passed": bool(exact_prefix_relation and all(row["passed"] for key, row in rows.items() if key != "join_archive_prefix")), "rows": rows}


def _owner_pin_rows(
    repo: Path,
    receipt: Mapping[str, Any],
    pins: Mapping[str, Mapping[str, Any]],
    hashes: Mapping[str, Any] | None,
    code_hashes: Mapping[str, Any] | None,
) -> dict[str, Any]:
    """Consume inherited maps when present and expose missing registrations."""
    inherited_hashes = _mapping(hashes)
    inherited_code = _mapping(code_hashes)
    receipt_inputs = _mapping(receipt.get("input_hashes"))

    def lookup(mapping: Mapping[str, Any], relative: str, basename: str) -> Any:
        for key in (relative, basename):
            if key not in mapping:
                continue
            value = mapping[key]
            if isinstance(value, Mapping):
                for digest_key in ("sha256", "hash", "digest"):
                    if digest_key in value:
                        return value[digest_key]
            return value
        return None

    rows = {}
    missing_registration = []
    for name, pin in pins.items():
        relative = str(pin.get("path", ""))
        basename = Path(relative).name
        actual = pin.get("sha256")
        # Repository implementation paths use code_hashes. Pinned original
        # experiment sources, archives and receipts use the input map rooted
        # at experiments/root_st073, as in every accepted ancestor packet.
        is_code = relative.startswith(('src/','scripts/'))
        owner_map = inherited_code if is_code else inherited_hashes
        represented = lookup(owner_map, relative, basename)
        receipt_expected = lookup(receipt_inputs, relative, basename)
        receipt_self = relative == _JOIN_RECEIPT
        receipt_match = (receipt_expected == actual) if receipt_expected is not None else receipt_self
        owner_match = represented is not None and represented == actual
        if represented is None:
            missing_registration.append({"map": "code_hashes" if is_code else "hashes", "path": relative, "sha256": actual})
        rows[name] = {
            "path": relative,
            "actual_sha256": actual,
            "receipt_input_hash_present": receipt_expected is not None,
            "receipt_input_hash_matches": receipt_match,
            "receipt_self_hash_excluded": receipt_self,
            "owner_map_entry_present": represented is not None,
            "owner_map_entry_matches": owner_match,
            "passed": bool(receipt_match and owner_match),
        }
    return {
        "passed": all(row["passed"] for row in rows.values()),
        "rows": rows,
        "missing_registration": missing_registration,
        "registration_required": bool(missing_registration),
    }


def _required_receipt_checks(receipt: Mapping[str, Any]) -> dict[str, Any]:
    checks = _mapping(receipt.get("checks"))
    names = (
        "source_contract",
        "original_radius_and_phase",
        "actual_flat_controls",
        "actual_background_germs",
        "actual_memory_and_leading_inlets",
        "actual_own_mixed_joins",
    )
    rows = {name: _mapping(checks.get(name)).get("passed") is True for name in names}
    memory = _mapping(checks.get("actual_memory_and_leading_inlets"))
    finite = _mapping(memory.get("actual_finite_N_endpoint_source_germs"))
    calculus = _mapping(memory.get("actual_U_source_radial_derivative_calculus"))
    own = _mapping(checks.get("actual_own_mixed_joins"))
    rows.update({
        "finite_root_and_physical_germs": finite.get("root_C5_Z_and_physical_y_germs_from_actual_background_IVP_and_original_generic_operator") is True,
        "finite_chain_rule_supply": finite.get("original_ordinary_calculus_and_physical_J_inverse_chain_rule_then_supply_U_y0_to4") is True,
        "finite_nonvalue_comparison": finite.get("actual_A_B_inverse_integrals_and_velocity_outputs_not_assigned_common_values") is True,
        "uy_first_source_binding": calculus.get("actual_U_y_first_row_uses_original_slow_partial_plus_fast_phase_recovery") is True,
        "uy_variable_J_calculus": calculus.get("actual_U_y2_to4_are_derivatives_of_actual_preceding_rows_with_variable_J") is True,
        "uy_unique_inverse": calculus.get("accepted_unique_inverse_composition_and_equal_background_germs_then_identify_high_U_germs") is True,
        "own_variable_radius": own.get("actual_nonconstant_native_radius_germs_not_constant_J_substitutes") is True,
        "own_histories_retained": own.get("nonzero_histories_and_one_plus_delta_own_recovery_and_m15_pressure_sign_retained") is True,
    })
    interface_rows = list(own.get("actual_seven_own_mixed_rows", ()))
    interface = [row for row in interface_rows if isinstance(row, Mapping) and row.get("interface") == "postrestore_patch"]
    rows["postrestore_patch_own_rows"] = len(interface) == 246
    rows["postrestore_patch_y1_rows"] = any(row.get("physical_y_order") == 1 for row in interface)
    y1_orders = {
        row.get("ordinary_Z_order")
        for row in interface
        if row.get("physical_y_order") == 1
    }
    rows["postrestore_patch_Z0_Z1_rows"] = {0, 1} <= y1_orders
    inverse_pairs = finite.get("actual_inverse_loop_and_finite_N_source_pairs", ())
    rows["postrestore_patch_finite_source_pairs"] = any(
        isinstance(row, list) and len(row) >= 3 and row[0] == "postrestore_patch" for row in inverse_pairs
    )
    return {
        "passed": bool(receipt.get("all_passed") is True and all(rows.values())),
        "receipt_all_passed": receipt.get("all_passed") is True,
        "rows": rows,
        "postrestore_patch_own_row_count": len(interface),
        "finite_source_pair_count": len(inverse_pairs) if isinstance(inverse_pairs, list) else 0,
        "unbound_own_pressure_output_references": own.get("unbound_own_pressure_output_references"),
    }


def _endpoint_source_rows(
    receipt: Mapping[str, Any],
    left_u_nodes: list[int],
    right_u_nodes: list[int],
    left_nodes: list[int],
    right_nodes: list[int],
) -> dict[str, Any]:
    memory = _mapping(_mapping(receipt.get("checks")).get("actual_memory_and_leading_inlets"))
    calculus = _mapping(memory.get("actual_U_source_radial_derivative_calculus"))
    source_rows = calculus.get("actual_U_y0_to4_C5_source_binding_rows", ())
    found: dict[tuple[str, int], int] = {}
    for row in source_rows if isinstance(source_rows, list) else ():
        if isinstance(row, list) and len(row) >= 5 and row[0] in ("postrestore", "Rm_patch"):
            try:
                found[(str(row[0]), int(row[1]), int(row[2]))] = int(row[4])
            except (TypeError, ValueError):
                continue
    rows = []
    for order in (0, 1):
        left = left_nodes[order] if len(left_nodes) > order else None
        right = right_nodes[order] if len(right_nodes) > order else None
        left_output = found.get(("postrestore", 1, order))
        right_output = found.get(("Rm_patch", 1, order))
        passed = left is not None and right is not None and left_output == left and right_output == right
        rows.append({
            "quantity": "U_y",
            "ordinary_Z_order": order,
            "left_node": left,
            "right_node": right,
            "left_source_output_node": left_output,
            "right_source_output_node": right_output,
            "passed": bool(passed),
        })
    memory = _mapping(_mapping(receipt.get("checks")).get("actual_memory_and_leading_inlets"))
    finite = _mapping(memory.get("actual_finite_N_endpoint_source_germs"))
    finite_pairs = finite.get("actual_inverse_loop_and_finite_N_source_pairs", ())
    pair_map: dict[int, tuple[int, int]] = {}
    for row in finite_pairs if isinstance(finite_pairs, list) else ():
        if isinstance(row, list) and len(row) >= 6 and row[:3] == ["postrestore_patch", "physical", "U"]:
            try:
                pair_map[int(row[3])] = (int(row[4]), int(row[5]))
            except (TypeError, ValueError):
                continue
    u_rows = []
    for order in range(3):
        left = left_u_nodes[order] if len(left_u_nodes) > order else None
        right = right_u_nodes[order] if len(right_u_nodes) > order else None
        pair = pair_map.get(order)
        passed = bool(pair and left == pair[0] and right == pair[1])
        u_rows.append({
            "quantity": "U",
            "ordinary_Z_order": order,
            "left_node": left,
            "right_node": right,
            "left_source_output_node": pair[0] if pair else None,
            "right_source_output_node": pair[1] if pair else None,
            "passed": passed,
        })
    return {
        "passed": bool(rows) and all(row["passed"] for row in rows) and all(row["passed"] for row in u_rows),
        "rows": rows,
        "U_rows": u_rows,
        "U_rows_passed": bool(u_rows) and all(row["passed"] for row in u_rows),
        "source_rows_total": len(source_rows) if isinstance(source_rows, list) else 0,
        "finite_source_pairs_total": len(finite_pairs) if isinstance(finite_pairs, list) else 0,
    }


def build_Rm_inlet_source_theorem(
    nodes: Any,
    rm: Mapping[str, Any],
    leading: Mapping[str, Any],
    l: Mapping[str, Any],
    modified: Mapping[str, Any],
    flat: Mapping[str, Any],
    native_guards: Mapping[str, Any],
    *,
    root: Any = ROOT,
    hashes: Mapping[str, Any] | None = None,
    code_hashes: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Bind the accepted source theorem to the actual postrestore/Rm ports.

    The function only reads serialized bindings and the accepted receipt.  A
    generic ``all_passed`` receipt is insufficient: the interface row counts,
    the accepted variable-J calculus flags, and the four direct physical
    ``U_y`` producer outputs are checked independently.
    """
    repo = _repo_root(root)
    bindings = {
        name: _unwrap(value, _BINDING_KEYS[name])
        for name, value in (("rm", rm), ("leading", leading), ("loop", l), ("modified", modified))
    }
    rm_b, leading_b, loop_b, modified_b = (bindings[name] for name in ("rm", "leading", "loop", "modified"))
    pins = {name: _pin(repo, relative, expected) for name, (relative, expected) in _SOURCE_PINS.items()}
    receipt, receipt_load = _load_receipt(repo)
    join_archive, join_archive_load = _load_join_archive(repo)
    join_checks = _join_binding_checks(join_archive, receipt)
    receipt_checks = _required_receipt_checks(receipt)
    prefix = _source_prefix_rows(nodes, bindings, join_archive)
    owner_pins = _owner_pin_rows(repo, receipt, pins, hashes, code_hashes)

    post = _chart(modified_b, "postrestore")
    loop_post = _chart(loop_b, "postrestore")
    chart = _mapping(rm_b.get("chart"))
    physical_ports, port_checks = _physical_ports(nodes, modified_b, rm_b)
    direct_left_u = physical_ports["left"].get("U", [])
    direct_right_u = physical_ports["right"].get("U", [])
    direct_left_uy = physical_ports["left"].get("U_y", [])
    direct_right_uy = physical_ports["right"].get("U_y", [])
    source_rows = _endpoint_source_rows(receipt, direct_left_u, direct_right_u, direct_left_uy, direct_right_uy)

    native_checks = _mapping(native_guards.get("checks"))
    native_required = {
        "native_guard_passed": native_guards.get("passed") is True,
        "moving_Rm_radius_and_J": _mapping(native_checks.get("moving_Rm_radius_and_J")).get("passed") is True,
        "physical_U_y_and_R_Q_y_factors": _mapping(native_checks.get("physical_U_y_and_R_Q_y_factors")).get("passed") is True,
        "ordinary_Z_rows": _mapping(native_checks.get("ordinary_Z_rows")).get("passed") is True,
        "total_V_y_slow_and_fast_rows": _mapping(native_checks.get("total_V_y_slow_and_fast_rows")).get("passed") is True,
    }
    flat_required = {
        "flat_certificate_passed": flat.get("passed") is True,
        "no_phase_reset": _mapping(flat.get("geometry")).get("no_phase_reset_assumed") is True,
        "opaque_controls": _mapping(flat.get("unique_implicit_control")).get("opaque_ZN_functions_preserved") is True,
        "incoming_memory_retained": _mapping(flat.get("incoming_history")).get("incoming_memory_retained") is True,
    }
    handle_rows = {
        "postrestore_coordinate": _handle(loop_post.get("coordinate")),
        "Rm_patch_coordinate": _handle(chart.get("coordinate")),
        "Rm_patch_radius": _handle(chart.get("radius")),
        "Rm_patch_J": _handle(chart.get("physical_measure")),
        "common_N": _handle(rm_b.get("N")),
        "loop_N": _handle(loop_b.get("N")),
        "phase_origin": _handle(rm_b.get("phase_origin")),
        "loop_phase_origin": _handle(loop_b.get("phase_origin")),
    }
    handles_passed = bool(
        all(value is not None for value in handle_rows.values())
        and handle_rows["common_N"] == handle_rows["loop_N"]
        and handle_rows["phase_origin"] == handle_rows["loop_phase_origin"]
        and handle_rows["postrestore_coordinate"] == _handle(loop_post.get("coordinate"))
    )
    handle_checks = {"passed": handles_passed, "rows": handle_rows}

    source_pins_passed = bool(
        all(row["passed"] for row in pins.values())
        and receipt_load.get("passed") is True
        and join_archive_load.get("passed") is True
        and join_checks.get("passed") is True
        and owner_pins.get("passed") is True
    )
    common_receipt = _digest({
        "accepted_receipt_sha256": pins["accepted_join_receipt"]["sha256"],
        "source_checks": receipt_checks,
        "endpoint": "postrestore_native1_to_Rm_patch_native1",
    })
    theorem_ready = bool(
        source_pins_passed
        and receipt_checks["passed"]
        and prefix["passed"]
        and all(flat_required.values())
        and all(native_required.values())
        and handle_checks["passed"]
        and port_checks["U_y"]["passed"]
        and port_checks["U"]["passed"]
        and source_rows["passed"]
    )
    rows = []
    for row in source_rows["rows"]:
        order = row["ordinary_Z_order"]
        common = _digest(["proved_Rm_physical_Uy", common_receipt, order])
        rows.append({
            **row,
            "common_response_receipt_sha256": common_receipt,
            "common_function_sha256": common,
            "endpoint_scope": {
                "left": {"chart": "postrestore", "native_endpoint": 1},
                "right": {"chart": "Rm_patch", "native_endpoint": 1},
            },
        })

    return {
        "passed": theorem_ready and len(rows) == 2 and all(row["passed"] for row in rows),
        "method": "accepted_Rm_outer_mixed_source_join_bound_to_actual_postrestore_Rm_physical_Uy_ports",
        "source_only": True,
        "source_values_selected": False,
        "numeric_caps_used": False,
        "ancestor_replayed": False,
        "accepted_receipt": {
            "path": _JOIN_RECEIPT,
            "sha256": pins["accepted_join_receipt"]["sha256"],
            "source_family": receipt.get("source_family"),
            "receipt_checks": receipt_checks,
            "common_response_receipt_sha256": common_receipt,
        },
        "source_pins": {"passed": source_pins_passed, "rows": pins},
        "owner_pin_rows": owner_pins,
        "join_archive": {
            "load": join_archive_load,
            "binding_checks": join_checks,
        },
        "source_prefixes": prefix,
        "endpoint_scope": {
            "left_chart": "postrestore",
            "left_native_endpoint": 1,
            "right_chart": "Rm_patch",
            "right_native_endpoint": 1,
            "common_N": handle_rows["common_N"],
            "phase_origin": handle_rows["phase_origin"],
            "global_fast_phase_not_reset": True,
        },
        "endpoint_handle_checks": handle_checks,
        "native_guard_requirements": {"passed": all(native_required.values()), "rows": native_required},
        "flat_history_requirements": {"passed": all(flat_required.values()), "rows": flat_required},
        "actual_port_bindings": {
            "passed": all(row["passed"] for row in port_checks.values()),
            "checks": port_checks,
            "ports": physical_ports,
            "source_U_y_rows": source_rows,
        },
        "U_rows": [
            {
                **row,
                "common_response_receipt_sha256": common_receipt,
                "common_function_sha256": _digest(["proved_Rm_physical_U", common_receipt, row["ordinary_Z_order"]]),
                "endpoint_scope": {
                    "left": {"chart": "postrestore", "native_endpoint": 1},
                    "right": {"chart": "Rm_patch", "native_endpoint": 1},
                },
            }
            for row in source_rows["U_rows"]
        ],
        "leading_mean_rows": {
            "passed": False,
            "authorized": False,
            "reason": "The accepted Rm outer receipt has no direct postrestore/leading mean node-pair rows; no opaque leading-mean alias is issued.",
        },
        "rows": rows,
        "trace_contract": {
            "U_through_ordinary_Z2_source_bound": True,
            "Uy_through_ordinary_Z1_source_bound": True,
            "m_my_Q_Qy_RQy_ports_exposed_but_not_aliased_here": True,
            "trace_ports_are_actual_serialized_nodes": True,
            "pressure_outputs_not_assigned_common_values": True,
            "unique_controls_remain_original_Z_functions_fixed_leading_N0": True,
            "common_candidate_N_correction_context_retained": True,
            "incoming_histories_retained": True,
            "physical_R_factor_consumed_once": True,
            "variable_J_not_replaced_by_constant": True,
        },
        "limits": {
            "full_raw_Uy_graph_fingerprint": False,
            "numeric_endpoint_values": False,
            "global_opaque_alias": False,
            "cap_contraction": False,
            "ancestor_constructor_or_checker_rerun": False,
        },
        "unproved": {
            "signed_finite_N_endpoint_values": True,
            "numeric_Rp_closure": True,
            "full_global_N1_closure": True,
        },
    }


__all__ = ["build_Rm_inlet_source_theorem"]
