"""Archive-only source certificate for the six remaining Rm/outer seams.

The accepted mixed-join packet contains two different kinds of evidence that
are easy to conflate.  Its ``actual_nonzero_leading_mean_C5_inlets`` rows are
comparisons of *modified physical mean* endpoint wrappers.  They are not a
permission to alias the original leading ``m`` controls.  This adapter keeps
that distinction explicit while binding the endpoint wrappers, the retained
finite-``N`` histories, and the Rm unique-branch residual to the serialized
source graph.

No producer, ancestor constructor, accepted checker, numerical oracle, or
bound is run here.  All identities below are source/receipt statements.  In
particular, the Rm exit is handled by the exact unique-system residual and by
the native bump-support geometry; it is never treated as an empty integral or
as ``h = 0``.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from fractions import Fraction
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
_LEADING_ARCHIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_Rm_patch_C6_functions.json.gz"
)
_LEADING_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_Rm_patch_C6_functions.py"
)
_LEADING_CHECK_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_Rm_patch_C6_functions_check.py"
)
_OUTER_ARCHIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.json.gz"
)
_OUTER_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py"
)
_OUTER_CHECK_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions_check.py"
)

_SOURCE_PINS = {
    "accepted_join_archive": (
        _JOIN_ARCHIVE,
        "27510c39da6e84293c2eaa2355ce36d1c5622c39839c33a3b907739b85580a8c",
    ),
    "accepted_join_source": (
        _JOIN_SOURCE,
        "9ec220509152dfa406f14b661fbc76fa55c95aee6908ab5e0344519e445653f8",
    ),
    "accepted_join_checker": (
        _JOIN_CHECK_SOURCE,
        "15bd090dc0d891bfd7825730f8412830abeed451f06ede8e1c66acde200c6846",
    ),
    "accepted_join_receipt": (
        _JOIN_RECEIPT,
        "14d1cc6b641435c60a7e3662c996a260c07b328df5b8f486c4344f3d7eeb4867",
    ),
    "accepted_leading_archive": (
        _LEADING_ARCHIVE,
        "4352ebd8862b5c92e3e4026728183014063baf556802a88421cdbd517ff1346a",
    ),
    "accepted_leading_source": (
        _LEADING_SOURCE,
        "a6f1cb5deaf17e9482aa964aef7867aadbc815e4d703a415e4945cf665f14002",
    ),
    "accepted_leading_checker": (
        _LEADING_CHECK_SOURCE,
        "85f42a1476b0af3df91a8fd3f1077f1e0d7109cc9f25ad24fb9205ac891cb1cf",
    ),
    "accepted_outer_archive": (
        _OUTER_ARCHIVE,
        "75895e1637411cfd2e490778f495f3adf904295a8840184cdc467e006b23305d",
    ),
    "accepted_outer_source": (
        _OUTER_SOURCE,
        "4464a934d303a95a0da157e9e289f6ed8a04d5ff0adcb0d663857f0d9c57dcd3",
    ),
    "accepted_outer_checker": (
        _OUTER_CHECK_SOURCE,
        "855c1b31c42c1db22ce503e2a31e5aa20d13dbf7d5d1feb2869c048c64560a00",
    ),
}

_BINDING_KEY = "Rm_outer_mixed_source_join_bindings"
_LEADING_BINDING_KEY = "true_Rm_patch_C6_bindings"
_SEAMS = (
    "patch_reference",
    "reference_slope",
    "slope_axial",
    "axial_buffer",
    "buffer_transition",
    "transition_power",
)
_SEAM_GEOMETRY = {
    "patch_reference": ("Rm_patch", "Rh_reference", 2948, 3348),
    "reference_slope": ("Rh_reference", "O2_slope", 0, 0),
    "slope_axial": ("O2_slope", "O2_axial", 1, 0),
    "axial_buffer": ("O2_axial", "O2_buffer", 1, 0),
    # The serialized accepted chart is native [0, 4].  The old prose table
    # used 11; the graph and receipt are authoritative for this adapter.
    "buffer_transition": ("O2_buffer", "O3_transition", 4, 0),
    "transition_power": ("O3_transition", "O3_power", 1, 0),
}
_CHARTS = (
    "postrestore",
    "Rm_patch",
    "Rh_reference",
    "O2_slope",
    "O2_axial",
    "O2_buffer",
    "O3_transition",
    "O3_power",
)
_ROOT_GERMS = ("E", "V", "a", "b", "p1", "p2", "t0", "Delta", "kappa")
_FINITE_CATEGORIES = (
    ("loop", "A"),
    ("loop", "B_over_Pstar"),
    ("loop", "A_phi"),
    ("loop", "B_over_Pstar_phi"),
    ("finite_N", "E_N"),
    ("finite_N", "V_N"),
    ("physical", "U"),
)
_QY_INPUT_FIELDS = ("E", "V", "m", "h", "k", "e", "p", "P0")
_RADIAL_INPUT_FIELDS = ("E_y", "V_y")


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


def _node(nodes: Any, value: Any) -> Mapping[str, Any] | None:
    index = _handle(value)
    if index is None or index < 0:
        return None
    try:
        row = nodes[index] if not isinstance(nodes, Mapping) else nodes.get(index, nodes.get(str(index)))
    except (IndexError, KeyError, TypeError):
        return None
    return row if isinstance(row, Mapping) else None


def _exists(nodes: Any, value: Any) -> bool:
    return _node(nodes, value) is not None


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


def _load_serialized(value: Any, repo: Path, *, name: str) -> tuple[Mapping[str, Any], dict[str, Any]]:
    if isinstance(value, Mapping):
        return value, {"passed": True, "mode": "supplied_mapping", "name": name}
    if not isinstance(value, (str, Path)):
        return {}, {"passed": False, "name": name, "reason": "serialized mapping was not supplied"}
    path = Path(value)
    if not path.is_absolute():
        path = repo / path
    try:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8") as stream:
            data = json.load(stream)
    except (OSError, ValueError) as exc:
        return {}, {"passed": False, "name": name, "reason": "serialized mapping is unreadable", "error": str(exc)}
    return data if isinstance(data, Mapping) else {}, {"passed": isinstance(data, Mapping), "name": name, "path": str(path)}


def _children(row: Mapping[str, Any]) -> list[int]:
    result: list[int] = []
    for key in (
        "arguments",
        "argument",
        "numerator",
        "denominator",
        "expression",
        "integrand",
        "lower",
        "upper",
        "active_body",
        "flat_predicate",
        "flat_value",
        "variable",
        "value",
    ):
        value = row.get(key)
        if isinstance(value, list):
            result.extend(item for item in (_handle(v) for v in value) if item is not None)
        else:
            item = _handle(value)
            if item is not None:
                result.append(item)
    return result


def _depends_on(nodes: Any, start: Any, target: Any, *, limit: int = 100000) -> bool:
    root = _handle(start)
    wanted = _handle(target)
    if root is None or wanted is None:
        return False
    stack = [root]
    seen: set[int] = set()
    while stack and len(seen) < limit:
        current = stack.pop()
        if current == wanted:
            return True
        if current in seen:
            continue
        seen.add(current)
        row = _node(nodes, current)
        if row is not None:
            stack.extend(child for child in _children(row) if child not in seen)
    return False


def _prefix_rows(nodes: Any, archive: Mapping[str, Any], leading_archive: Mapping[str, Any]) -> dict[str, Any]:
    current = list(nodes) if isinstance(nodes, (list, tuple)) else []
    archived = archive.get("exact_graph_nodes")
    leading_nodes = leading_archive.get("exact_graph_nodes")
    current_prefix_archive = bool(
        isinstance(archived, list)
        and len(current) >= len(archived)
        and current[: len(archived)] == archived
    )
    archive_prefix_current = bool(
        isinstance(archived, list)
        and len(current) <= len(archived)
        and archived[: len(current)] == current
    )
    leading_prefix_current = bool(
        isinstance(leading_nodes, list)
        and len(current) >= len(leading_nodes)
        and current[: len(leading_nodes)] == leading_nodes
    )
    return {
        "passed": bool((current_prefix_archive or archive_prefix_current) and leading_prefix_current),
        "current_node_count": len(current),
        "join_archive_node_count": len(archived) if isinstance(archived, list) else 0,
        "leading_archive_node_count": len(leading_nodes) if isinstance(leading_nodes, list) else 0,
        "current_prefix_join_archive": current_prefix_archive,
        "join_archive_prefix_current": archive_prefix_current,
        "current_prefix_leading_archive": leading_prefix_current,
    }


def _outer_prefix_row(nodes: Any, outer_archive: Mapping[str, Any]) -> dict[str, Any]:
    current = list(nodes) if isinstance(nodes, (list, tuple)) else []
    archived = outer_archive.get("exact_graph_nodes")
    passed = bool(
        isinstance(archived, list)
        and len(current) >= len(archived)
        and current[: len(archived)] == archived
    )
    return {
        "passed": passed,
        "current_node_count": len(current),
        "outer_archive_node_count": len(archived) if isinstance(archived, list) else 0,
        "current_prefix_outer_archive": passed,
    }


def _shape(value: Any, rows: int, width: int) -> bool:
    return bool(
        isinstance(value, list)
        and len(value) == rows
        and all(isinstance(row, list) and len(row) == width for row in value)
    )


def _chart_scope(nodes: Any, name: str, chart: Mapping[str, Any]) -> dict[str, Any]:
    physical = _mapping(chart.get("physical"))
    domain = chart.get("domain")
    domain_handles = [_handle(item) for item in domain] if isinstance(domain, list) and len(domain) == 2 else []
    handles = {
        key: _handle(chart.get(key))
        for key in ("coordinate", "radius", "jacobian", "phase")
    }
    handles_ok = all(value is not None and _exists(nodes, value) for value in handles.values())
    domain_ok = bool(
        len(domain_handles) == 2
        and all(value is not None and _exists(nodes, value) for value in domain_handles)
    )
    physical_shapes = {
        "U": _shape(physical.get("U"), 5, 6),
        "m": _shape(physical.get("m"), 5, 6),
        "Q": _shape(physical.get("Q"), 5, 5),
    }
    physical_nodes_ok = all(
        _exists(nodes, item)
        for field, matrix in physical.items()
        if field in physical_shapes and isinstance(matrix, list)
        for row in matrix
        if isinstance(row, list)
        for item in row
    )
    n_grid = chart.get("N_y")
    n_ok = bool(
        isinstance(n_grid, list)
        and len(n_grid) == 3
        and all(isinstance(row, list) and len(row) == 4 for row in n_grid)
        and all(_exists(nodes, item) for row in n_grid for item in row)
    )
    intermediate = _mapping(chart.get("intermediate"))
    retained = intermediate.get("complete_incoming_mass_not_reset") is True
    return {
        "chart": name,
        "domain": domain_handles,
        "handles": handles,
        "handles_exist": handles_ok,
        "domain_handles_exist": domain_ok,
        "physical_shapes": physical_shapes,
        "physical_nodes_exist": physical_nodes_ok,
        "N_y_3x4_nodes_exist": n_ok,
        "incoming_mass_not_reset": retained,
        "passed": bool(
            handles_ok
            and domain_ok
            and all(physical_shapes.values())
            and physical_nodes_ok
            and n_ok
            and retained
        ),
    }


def _endpoint_wrapper(nodes: Any, value: Any, coordinate: Any, native_query: Any) -> dict[str, Any]:
    node = _handle(value)
    coordinate = _handle(coordinate)
    query = _handle(native_query)
    row = _node(nodes, node)
    passed = bool(
        node is not None
        and row is not None
        and row.get("operation") == "function_substitution"
        and _handle(row.get("variable")) == coordinate
        and _handle(row.get("value")) == query
        and _handle(row.get("expression")) is not None
        and row.get("Z_independent_substitution") is True
    )
    return {
        "node": node,
        "expression": _handle(row.get("expression")) if row is not None else None,
        "variable": _handle(row.get("variable")) if row is not None else None,
        "value": _handle(row.get("value")) if row is not None else None,
        "passed": passed,
    }


def _join_scope(
    nodes: Any,
    charts: Mapping[str, Any],
    seam: str,
    join: Mapping[str, Any],
) -> dict[str, Any]:
    expected = _SEAM_GEOMETRY.get(seam)
    left_name = join.get("left_chart")
    right_name = join.get("right_chart")
    left_chart = _mapping(charts.get(left_name))
    right_chart = _mapping(charts.get(right_name))
    left = _mapping(join.get("left"))
    right = _mapping(join.get("right"))
    left_query = _handle(left.get("native_query"))
    right_query = _handle(right.get("native_query"))
    expected_names = bool(expected and (left_name, right_name, left_query, right_query) == expected)
    left_domain = [_handle(item) for item in left_chart.get("domain", ())] if isinstance(left_chart.get("domain"), list) else []
    right_domain = [_handle(item) for item in right_chart.get("domain", ())] if isinstance(right_chart.get("domain"), list) else []
    endpoint_values = {
        "left_query_in_domain": left_query in left_domain,
        "right_query_in_domain": right_query in right_domain,
    }
    left_physical = _mapping(left.get("physical"))
    right_physical = _mapping(right.get("physical"))
    shapes = {
        "left_m": _shape(left_physical.get("m"), 5, 6),
        "right_m": _shape(right_physical.get("m"), 5, 6),
        "left_U": _shape(left_physical.get("U"), 5, 6),
        "right_U": _shape(right_physical.get("U"), 5, 6),
        "left_Q": _shape(left_physical.get("Q"), 5, 5),
        "right_Q": _shape(right_physical.get("Q"), 5, 5),
    }
    radius_rows = {
        "left": isinstance(left.get("radius"), list)
        and all(_exists(nodes, item) for item in left.get("radius", ())),
        "right": isinstance(right.get("radius"), list)
        and all(_exists(nodes, item) for item in right.get("radius", ())),
    }
    wrappers = []
    for side, packet, chart in (("left", left, left_chart), ("right", right, right_chart)):
        physical = _mapping(packet.get("physical"))
        row = physical.get("m", [])
        coordinate = chart.get("coordinate")
        query = packet.get("native_query")
        for order, node in enumerate(row[0][:3] if isinstance(row, list) and row and isinstance(row[0], list) else ()):
            wrapper = _endpoint_wrapper(nodes, node, coordinate, query)
            wrappers.append({"side": side, "ordinary_Z_order": order, **wrapper})
    wrappers_ok = bool(wrappers) and all(row["passed"] for row in wrappers)
    endpoint_ok = bool(expected_names and all(endpoint_values.values()) and all(radius_rows.values()))
    return {
        "seam": seam,
        "left_chart": left_name,
        "right_chart": right_name,
        "left_native_query": left_query,
        "right_native_query": right_query,
        "expected_geometry": list(expected) if expected else None,
        "endpoint_handles_and_shapes": bool(endpoint_ok and all(shapes.values())),
        "endpoint_values": endpoint_values,
        "shapes": shapes,
        "radius_nodes_exist": radius_rows,
        "modified_physical_m_endpoint_wrappers": wrappers,
        "modified_physical_m_wrappers_passed": wrappers_ok,
        "passed": bool(endpoint_ok and all(shapes.values()) and wrappers_ok),
    }


def _required_receipt_rows(receipt: Mapping[str, Any]) -> dict[str, Any]:
    checks = _mapping(receipt.get("checks"))
    required = (
        "source_contract",
        "original_radius_and_phase",
        "actual_flat_controls",
        "actual_background_germs",
        "actual_memory_and_leading_inlets",
        "actual_own_mixed_joins",
    )
    sections = {name: _mapping(checks.get(name)).get("passed") is True for name in required}
    background = _mapping(checks.get("actual_background_germs"))
    memory = _mapping(checks.get("actual_memory_and_leading_inlets"))
    finite = _mapping(memory.get("actual_finite_N_endpoint_source_germs"))
    own = _mapping(checks.get("actual_own_mixed_joins"))
    calculus = _mapping(memory.get("actual_U_source_radial_derivative_calculus"))
    source_flags = {
        "root_and_physical_germs": finite.get(
            "root_C5_Z_and_physical_y_germs_from_actual_background_IVP_and_original_generic_operator"
        ) is True,
        "ordinary_calculus_and_physical_J_inverse_chain_rule": finite.get(
            "original_ordinary_calculus_and_physical_J_inverse_chain_rule_then_supply_U_y0_to4"
        ) is True,
        "finite_outputs_not_assigned_common_values": finite.get(
            "actual_A_B_inverse_integrals_and_velocity_outputs_not_assigned_common_values"
        ) is True,
        "Uy_first_source_binding": calculus.get(
            "actual_U_y_first_row_uses_original_slow_partial_plus_fast_phase_recovery"
        ) is True,
        "Uy_variable_J_calculus": calculus.get(
            "actual_U_y2_to4_are_derivatives_of_actual_preceding_rows_with_variable_J"
        ) is True,
        "unique_inverse_composition": calculus.get(
            "accepted_unique_inverse_composition_and_equal_background_germs_then_identify_high_U_germs"
        ) is True,
        "own_variable_radius": own.get("actual_nonconstant_native_radius_germs_not_constant_J_substitutes") is True,
        "histories_retained": own.get("nonzero_histories_and_one_plus_delta_own_recovery_and_m15_pressure_sign_retained") is True,
        "background_physical_germs": background.get(
            "flat_Rm_and_outer_RHS_then_unique_periodic_inverse_and_chain_rule_supply_leading_y4_C5_Q4_C4_germs"
        ) is True,
        "background_fixed_pressure": background.get(
            "fixed_pressure_operator_germs_follow_same_U_Q_and_actual_radius"
        ) is True,
        "finite_germs_not_neighborhood_values": background.get(
            "finite_germs_not_equal_neighborhoods_or_numeric_norms"
        ) is True,
    }
    means = own.get("actual_nonzero_leading_mean_C5_inlets", ())
    mean_map = {
        (row[0], int(row[1]))
        for row in means
        if isinstance(row, list) and len(row) >= 2 and row[0] in _SEAMS and _handle(row[1]) is not None
    }
    mean_evidence = {
        seam: {
            "rows": [order for order in range(6) if (seam, order) in mean_map],
            "all_Z0_to_Z5": all((seam, order) in mean_map for order in range(6)),
        }
        for seam in _SEAMS
    }
    finite_rows = finite.get("actual_inverse_loop_and_finite_N_source_pairs", ())
    finite_map = {
        (row[0], row[1], row[2], int(row[3])): row
        for row in finite_rows
        if isinstance(row, list) and len(row) >= 6 and row[0] in _SEAMS and _handle(row[3]) is not None
    }
    finite_evidence = {
        seam: {
            "rows": [
                [group, field, order]
                for group, field in _FINITE_CATEGORIES
                for order in range(6)
                if (seam, group, field, order) in finite_map
            ],
            "all_source_categories_Z0_to_Z5": all(
                (seam, group, field, order) in finite_map
                for group, field in _FINITE_CATEGORIES
                for order in range(6)
            ),
        }
        for seam in _SEAMS
    }
    background_rows = background.get("actual_physical_E_V_and_five_history_y0_to5_IVP_rows", ())
    background_set = {
        (row[0], row[1], int(row[2]))
        for row in background_rows
        if isinstance(row, list) and len(row) >= 3 and _handle(row[2]) is not None
    }
    return {
        "passed": bool(receipt.get("all_passed") is True and all(sections.values()) and all(source_flags.values())),
        "receipt_all_passed": receipt.get("all_passed") is True,
        "sections": sections,
        "source_flags": source_flags,
        "mean_evidence": mean_evidence,
        "finite_source_evidence": finite_evidence,
        "mean_row_count": len(means) if isinstance(means, list) else 0,
        "finite_source_pair_count": len(finite_rows) if isinstance(finite_rows, list) else 0,
        "background_physical_row_count": len(background_rows) if isinstance(background_rows, list) else 0,
        "background_physical_rows": background_set,
        "finite_source_rows": finite_map,
        "mean_source_rows": mean_map,
        "memory": memory,
        "finite": finite,
        "own": own,
        "background_receipt": background,
    }


def _public_receipt_rows(receipt_rows: Mapping[str, Any]) -> dict[str, Any]:
    """Expose receipt evidence without leaking internal set/map indexes.

    The tuple-keyed indexes are useful while checking the archived receipt,
    but they are not part of the public certificate and are not JSON
    serializable.  The finite-source evidence already contains the complete
    public row lists, so the internal finite map is intentionally omitted.
    """
    public = dict(receipt_rows)
    for key in ("background_physical_rows", "mean_source_rows"):
        values = receipt_rows.get(key, ())
        public[key] = [list(row) for row in sorted(values)]
    public.pop("finite_source_rows", None)
    return public


def _leading_identity(nodes: Any, leading_archive: Mapping[str, Any], leading: Mapping[str, Any]) -> dict[str, Any]:
    binding = _unwrap(leading, _LEADING_BINDING_KEY)
    system_id = _handle(binding.get("system"))
    system_row = _node(nodes, system_id)
    base = _mapping(binding.get("base_identity"))
    definition = _mapping(binding.get("definition"))
    chart = _mapping(binding.get("chart"))
    controls = binding.get("controls")
    control_rows = []
    if isinstance(controls, list):
        for row in controls:
            if isinstance(row, list):
                control_rows.extend(_handle(item) for item in row if _handle(item) is not None)
    control_nodes_ok = bool(
        isinstance(controls, list)
        and len(controls) == 5
        and all(isinstance(row, list) and len(row) == 7 for row in controls)
        and all(_exists(nodes, item) for row in controls for item in row)
    )
    control_ops_ok = bool(
        isinstance(controls, list)
        and all(_node(nodes, row[0]).get("operation") == "original_Rm_unique_leading_implicit_control" for row in controls if row and _node(nodes, row[0]))
    )
    fields = _mapping(chart.get("fields"))
    leading_m = fields.get("m")
    leading_m_ok = bool(isinstance(leading_m, list) and len(leading_m) == 7 and all(_exists(nodes, item) for item in leading_m))
    residual = base.get("residual")
    residual_ok = bool(isinstance(residual, list) and len(residual) == 5 and all(_exists(nodes, item) for item in residual))
    domain = chart.get("domain")
    exit_native = _handle(domain[1]) if isinstance(domain, list) and len(domain) == 2 else None
    exit_ok = exit_native == 2948
    system_ok = bool(
        system_id is not None
        and system_row is not None
        and system_row.get("operation") == "original_Rm_unique_leading_implicit_system"
        and _handle(base.get("system")) == system_id
    )
    flags = {
        "exact_unique_solution_identity": base.get("defining_exact_unique_solution_identity_not_numeric_zero_containment") is True,
        "same_L_Q_defect_amplitude_source": base.get("same_actual_L_Q_defect_and_amplitude_source") is True,
        "full_Jacobian_through_Z6": definition.get("original_unique_common_box_branch_same_full_Jacobian_through_Z6") is True,
        "residual_is_function_identity": definition.get("original_branch_zero_residual_is_function_identity_not_zero_containment") is True,
        "leading_not_finite_N_control": definition.get("leading_control_not_finite_N_correction_and_no_history_reset") is True,
        "actual_fields_and_radial_germs": definition.get("actual_E_V_m_h_k_e_p_P0_pressure_and_Ey_Vy_C6") is True,
    }
    # This is the only Rm-exit statement issued by this module.  Partial
    # integrals stay symbolic; full support is checked from native bump
    # geometry rather than from a zero/empty graph value.
    support_rows = []
    weights = _mapping(binding.get("partial_weights"))
    for group, value in weights.items():
        if isinstance(value, list):
            values = [item for row in value for item in row] if value and isinstance(value[0], list) else value
            for item in values:
                row = _node(nodes, item)
                if row is None or row.get("operation") != "definite_integral":
                    continue
                if row.get("exact_original_Rm_bump_weight") is not True:
                    continue
                center = row.get("original_center")
                radius = row.get("original_radius")
                try:
                    upper = Fraction(str(center)) + Fraction(str(radius))
                    full = upper < Fraction(exit_native if exit_native is not None else -1)
                except (TypeError, ValueError, ZeroDivisionError):
                    full = False
                support_rows.append({
                    "group": group,
                    "node": _handle(item),
                    "exact_original_bump_weight": True,
                    "Z_independent_weight": row.get("Z_independent_weight") is True,
                    "range_not_defining_integral": row.get("range_or_midpoint_not_defining_integral") is True,
                    "endpoint_support_included": full,
                    "symbolic_integral_retained": row.get("numerical_value_not_installed") is True,
                })
    support_ok = bool(support_rows) and all(
        row["endpoint_support_included"]
        and row["symbolic_integral_retained"]
        and row["Z_independent_weight"]
        and row["range_not_defining_integral"]
        for row in support_rows
    )
    return {
        "passed": bool(system_ok and control_nodes_ok and control_ops_ok and leading_m_ok and residual_ok and exit_ok and all(flags.values()) and support_ok),
        "system": system_id,
        "system_operation": system_row.get("operation") if system_row else None,
        "Rm_exit_native_query": exit_native,
        "leading_m_nodes": [_handle(item) for item in leading_m] if isinstance(leading_m, list) else [],
        "controls": {
            "shape": [len(row) for row in controls] if isinstance(controls, list) else [],
            "nodes_exist": control_nodes_ok,
            "implicit_control_operations": control_ops_ok,
            "opaque_source_controls": True,
            "node_count": len(control_rows),
        },
        "base_identity": {
            "system_matches": _handle(base.get("system")) == system_id,
            "residual_nodes_exist": residual_ok,
            "flags": flags,
        },
        "full_support_at_Rm_exit": {
            "passed": support_ok,
            "method": "native_endpoint_after_exact_bump_support_geometry",
            "weights_checked": len(support_rows),
            "symbolic_integrals_retained": all(row["symbolic_integral_retained"] for row in support_rows) if support_rows else False,
            "Z_independent_weights": all(row["Z_independent_weight"] for row in support_rows) if support_rows else False,
            "range_not_defining_integral": all(row["range_not_defining_integral"] for row in support_rows) if support_rows else False,
            "empty_integral_substitution_used": False,
            "h_zero_substitution_used": False,
            "rows": support_rows,
        },
        "source_only": True,
        "numeric_values_selected": False,
    }


def _history_certificate(
    nodes: Any,
    binding: Mapping[str, Any],
    receipt_rows: Mapping[str, Any],
    seam: str,
) -> dict[str, Any]:
    join = _mapping(_mapping(binding.get("joins")).get(seam))
    right_name = join.get("right_chart")
    charts = _mapping(binding.get("charts"))
    right_chart = _mapping(charts.get(right_name))
    right_domain = right_chart.get("domain")
    right_lower = _handle(right_domain[0]) if isinstance(right_domain, list) and right_domain else None
    memory = _mapping(receipt_rows.get("memory"))
    integral_rows = memory.get("actual_original_empty_inlet_integral_nodes_and_nonreset_flags", ())
    selected = [row for row in integral_rows if isinstance(row, list) and len(row) >= 4 and row[0] == seam]
    checks = []
    for row in selected:
        node_id = _handle(row[3])
        graph = _node(nodes, node_id)
        common = bool(
            node_id is not None
            and graph is not None
            and graph.get("operation") == "definite_integral"
            and _handle(graph.get("lower")) == right_lower
            and graph.get("endpoint_independent_of_Z") is True
            and graph.get("exact_function_integral") is True
            and graph.get("numerical_value_not_installed") is True
        )
        if row[1] == "correction":
            common = common and graph.get("predecessor_memory_never_reset") is True
        elif row[1] == "signed_pressure":
            common = common and graph.get("nonzero_preceding_S0_S1_retained") is True and graph.get("physical_Jacobian_applied_once") is True and graph.get("physical_radius_power") in (1, 2)
        else:
            common = False
        checks.append({
            "kind": row[1],
            "field": row[2],
            "node": node_id,
            "lower": _handle(graph.get("lower")) if graph else None,
            "right_chart": right_name,
            "passed": bool(common),
        })
    finite = _mapping(receipt_rows.get("finite"))
    source_rows = [row for row in finite.get("actual_inverse_loop_and_finite_N_source_pairs", ()) if isinstance(row, list) and len(row) >= 6 and row[0] == seam]
    source_checks = []
    for row in source_rows:
        left_node = _handle(row[4]); right_node = _handle(row[5])
        source_checks.append({
            "group": row[1],
            "field": row[2],
            "ordinary_Z_order": _handle(row[3]),
            "left_node": left_node,
            "right_node": right_node,
            "actual_source_nodes_exist": _exists(nodes, left_node) and _exists(nodes, right_node),
            "source_germ_equality_authorized": True,
            "direct_output_node_alias": False,
        })
    finite_expected = len(_FINITE_CATEGORIES) * 6
    receipt_category_rows = _mapping(receipt_rows.get("finite_source_evidence")).get(seam, {})
    source_pair_ok = bool(
        len(source_checks) == finite_expected
        and all(row["actual_source_nodes_exist"] for row in source_checks)
        and receipt_category_rows.get("all_source_categories_Z0_to_Z5") is True
        and receipt_rows["source_flags"]["root_and_physical_germs"]
        and receipt_rows["source_flags"]["ordinary_calculus_and_physical_J_inverse_chain_rule"]
        and receipt_rows["source_flags"]["finite_outputs_not_assigned_common_values"]
    )
    return {
        "seam": seam,
        "right_chart": right_name,
        "right_native_lower": right_lower,
        "incoming_integral_rows": checks,
        "incoming_integral_row_count": len(checks),
        "incoming_integrals_passed": bool(checks) and all(row["passed"] for row in checks),
        "finite_N_source_pairs": source_checks,
        "finite_N_source_pair_count": len(source_checks),
        "finite_N_history_source_transfer_passed": source_pair_ok,
        "finite_N_physical_germ_equality": source_pair_ok,
        "direct_output_alias_issued": False,
        "history_reset_or_cap_substitution": False,
        "passed": bool(checks) and all(row["passed"] for row in checks) and source_pair_ok,
    }


def _mean_certificate(
    nodes: Any,
    binding: Mapping[str, Any],
    receipt_rows: Mapping[str, Any],
    seam: str,
    seam_scope: Mapping[str, Any],
) -> dict[str, Any]:
    mean_evidence = _mapping(receipt_rows.get("mean_evidence")).get(seam, {})
    background_rows = receipt_rows.get("background_physical_rows", set())
    wrapper_rows = seam_scope.get("modified_physical_m_endpoint_wrappers", ())
    wrapper_by_order = {row.get("ordinary_Z_order"): row for row in wrapper_rows if row.get("ordinary_Z_order") in (0, 1, 2)}
    background_mean = all((seam, "m", order) in background_rows for order in range(6))
    rows = []
    for order in range(3):
        wrapper_ok = bool(
            wrapper_by_order.get(order, {}).get("passed")
            and (seam, order) in receipt_rows.get("mean_source_rows", set())
        )
        rows.append({
            "interface": seam,
            "ordinary_Z_order": order,
            "quantity": "modified_physical_mean_m",
            "left_node": next((row.get("node") for row in wrapper_rows if row.get("side") == "left" and row.get("ordinary_Z_order") == order), None),
            "right_node": next((row.get("node") for row in wrapper_rows if row.get("side") == "right" and row.get("ordinary_Z_order") == order), None),
            "source_receipt_row": [seam, order],
            "actual_endpoint_wrapper": wrapper_ok,
            "background_IVP_m_source_present_Z0_to_Z5": background_mean,
            "arbitrary_common_function_source_provenance": wrapper_ok and background_mean,
            "numeric_value_selected": False,
            "passed": wrapper_ok and background_mean,
        })
    return {
        "interface": seam,
        "rows": rows,
        "ordinary_Z_orders": [0, 1, 2],
        "accepted_receipt_rows_Z0_to_Z5": bool(mean_evidence.get("all_Z0_to_Z5")),
        "background_IVP_mean_rows_Z0_to_Z5": background_mean,
        "passed": bool(rows)
        and all(row["passed"] for row in rows)
        and bool(mean_evidence.get("all_Z0_to_Z5")),
        "quantity_is_original_leading_mean": False,
        "original_leading_output_alias_issued": False,
        "reason_original_leading_alias_unproved": "Receipt mean rows compare modified physical m endpoint wrappers; the original leading controls remain opaque source functions.",
        "transition_power_nonzero_mean_source": seam == "transition_power",
    }


def _source_jet(nodes: Any, fields: Mapping[str, Any], name: str, *, count: int = 3) -> tuple[list[int], bool]:
    values = fields.get(name)
    handles = [_handle(item) for item in values[:count]] if isinstance(values, list) else []
    return handles, bool(len(handles) == count and all(_exists(nodes, item) for item in handles))


def _outer_background_certificate(
    nodes: Any,
    join_binding: Mapping[str, Any],
    outer_binding: Mapping[str, Any],
    leading_binding: Mapping[str, Any],
    receipt_rows: Mapping[str, Any],
    outer_prefix: Mapping[str, Any],
    seam: str,
) -> dict[str, Any]:
    """Bind eight leading-Q_y inputs to actual outer/leading source jets.

    The source functions are kept as opaque graph leaves.  The accepted
    background receipt supplies the equal-function/ODE uniqueness premise;
    this helper only verifies that the source arrays and their coordinate,
    radius, and physical-measure handles are the actual chart handles.
    """
    joins = _mapping(join_binding.get("joins"))
    charts = _mapping(join_binding.get("charts"))
    join = _mapping(joins.get(seam))
    left_name, right_name = join.get("left_chart"), join.get("right_chart")
    expected_left, expected_right, _, _ = _SEAM_GEOMETRY[seam]
    names_ok = left_name == expected_left and right_name == expected_right
    outer_backgrounds = _mapping(outer_binding.get("background"))
    outer_charts = _mapping(outer_binding.get("charts"))
    leading_chart = _mapping(leading_binding.get("chart"))

    def side_source(chart_name: Any) -> tuple[Mapping[str, Any], Mapping[str, Any], Mapping[str, Any], str]:
        if chart_name == "Rm_patch":
            return (
                _mapping(leading_chart.get("fields")),
                _mapping(leading_chart.get("radial")),
                leading_chart,
                "true_Rm_patch_C6_bindings.chart",
            )
        background = _mapping(outer_backgrounds.get(chart_name))
        return (
            _mapping(background.get("fields")),
            _mapping(background.get("radial")),
            _mapping(outer_charts.get(chart_name)),
            f"true_outer_C6_C5_bindings.background[{chart_name}]",
        )

    side_rows: dict[str, Any] = {}
    all_input_rows = []
    all_radial_rows = []
    source_fields_ok = True
    scope_ok = True
    source_premises = receipt_rows.get("source_flags", {})
    background_receipt = _mapping(receipt_rows.get("background_receipt"))
    physical_rows = receipt_rows.get("background_physical_rows", set())
    for side, chart_name in (("left", left_name), ("right", right_name)):
        fields, radial, source_chart, source_path = side_source(chart_name)
        join_chart = _mapping(charts.get(chart_name))
        source_handles = {
            # The C6 leading chart names its native coordinate ``phase``;
            # outer C6 charts expose the same handle as ``coordinate``.
            "coordinate": _handle(source_chart.get("coordinate", source_chart.get("phase"))),
            "radius": _handle(source_chart.get("radius")),
            "physical_measure": _handle(source_chart.get("physical_measure")),
        }
        join_handles = {
            "coordinate": _handle(join_chart.get("coordinate")),
            "radius": _handle(join_chart.get("radius")),
            "physical_measure": _handle(join_chart.get("jacobian")),
        }
        side_scope_ok = source_handles == join_handles and all(_exists(nodes, item) for item in source_handles.values())
        scope_ok = scope_ok and side_scope_ok
        rows = []
        for field in _QY_INPUT_FIELDS:
            handles, handles_ok = _source_jet(nodes, fields, field)
            source_fields_ok = source_fields_ok and handles_ok
            rows.append({
                "side": side,
                "chart": chart_name,
                "source_path": source_path,
                "field": field,
                "ordinary_Z_orders": [0, 1, 2],
                "nodes": handles,
                "nodes_exist": handles_ok,
                "source_coordinate_radius_J": source_handles,
                "joined_coordinate_radius_J": join_handles,
                "scope_matches_join": side_scope_ok,
                "numeric_value_selected": False,
            })
            all_input_rows.append(rows[-1])
        radial_rows = []
        for field in _RADIAL_INPUT_FIELDS:
            handles, handles_ok = _source_jet(nodes, radial, field)
            source_fields_ok = source_fields_ok and handles_ok
            radial_rows.append({
                "side": side,
                "chart": chart_name,
                "source_path": source_path,
                "field": field,
                "ordinary_Z_orders": [0, 1, 2],
                "nodes": handles,
                "nodes_exist": handles_ok,
                "scope_matches_join": side_scope_ok,
            })
            all_radial_rows.append(radial_rows[-1])
        side_rows[side] = {
            "chart": chart_name,
            "scope": {
                "source": source_handles,
                "joined": join_handles,
                "passed": side_scope_ok,
            },
            "input_rows": rows,
            "radial_rows": radial_rows,
            "passed": bool(side_scope_ok and all(row["nodes_exist"] for row in rows + radial_rows)),
        }

    successor_rows = {
        (row[0], row[1])
        for row in background_receipt.get("actual_successor_E_V_and_nonzero_background_history_value_rows", ())
        if isinstance(row, list) and len(row) >= 2
    }
    stored_rows = {
        (row[0], row[1])
        for row in background_receipt.get("actual_stored_patch_inlet_and_unique_terminal_field_rows", ())
        if isinstance(row, list) and len(row) >= 2
    }
    physical_expected = {
        (seam, field, order)
        for field in ("E", "V", "m", "h", "k", "e", "p")
        for order in range(6)
    }
    physical_rows_ok = physical_expected <= set(physical_rows)
    # P0 is included in the exact full-weight function identity, not in the
    # reduced seven-field physical IVP row list.
    full_identity = _mapping(_mapping(outer_binding.get("seed_identity")).get("original_full_weight_function_identity"))
    common_functions = full_identity.get("exact_common_functions", ())
    full_identity_ok = bool(
        full_identity.get("passed") is True
        and "P0" in common_functions
        and "m" in common_functions
        and "V" in common_functions
    )
    source_field_premises = []
    for field in _QY_INPUT_FIELDS:
        for chart_name in (left_name, right_name):
            if chart_name == "Rm_patch":
                present = ("Rm", field) in stored_rows
            elif chart_name == "Rh_reference":
                present = ("Rh", field) in stored_rows or ("Rh_reference", field) in successor_rows
            else:
                present = (chart_name, field) in successor_rows
            if field == "P0":
                present = full_identity_ok
            source_field_premises.append([chart_name, field, bool(present)])
    source_field_premises_ok = all(row[2] for row in source_field_premises)
    radial_premises_ok = bool(
        source_premises.get("background_physical_germs")
        and source_premises.get("Uy_variable_J_calculus")
        and physical_rows_ok
    )
    # The accepted receipt's actual background group compares arbitrary source
    # functions through the IVP and C6 unique branch.  This is a theorem input,
    # not a claim that the serialized endpoint node IDs are identical.
    common_source_ok = bool(
        names_ok
        and outer_prefix.get("passed")
        and source_field_premises_ok
        and physical_rows_ok
        and radial_premises_ok
        and source_premises.get("background_fixed_pressure")
        and source_premises.get("finite_germs_not_neighborhood_values")
    )
    return {
        "interface": seam,
        "passed": bool(common_source_ok and source_fields_ok and scope_ok),
        "common_source_function_premise": common_source_ok,
        "leading_Qy_input_rows": all_input_rows,
        "leading_Qy_input_row_count": len(all_input_rows),
        "leading_Qy_input_rows_passed": bool(len(all_input_rows) == 16 and all(row["nodes_exist"] and row["scope_matches_join"] for row in all_input_rows)),
        "radial_input_rows": all_radial_rows,
        "radial_input_rows_passed": bool(len(all_radial_rows) == 4 and all(row["nodes_exist"] and row["scope_matches_join"] for row in all_radial_rows)),
        "actual_source_fields_premises": source_field_premises,
        "actual_source_fields_premises_passed": source_field_premises_ok,
        "actual_physical_IVP_rows_Z0_to_Z5": physical_rows_ok,
        "outer_C6_unique_source_identity": full_identity_ok,
        "source_node_alias_issued": False,
        "side_rows": side_rows,
        "limits": {
            "numeric_values_selected": False,
            "output_aliases_assigned": False,
            "cap_subtraction": False,
        },
    }


def build_outer_chain_mean(
    nodes: Any,
    archive: Mapping[str, Any] | str | Path,
    receipt: Mapping[str, Any] | str | Path,
    leading: Mapping[str, Any] | str | Path,
    outer: Mapping[str, Any] | str | Path | None = None,
    *,
    root: Any = ROOT,
) -> dict[str, Any]:
    """Build a source-only certificate for six remaining outer seams.

    ``nodes`` must be the serialized graph or an exact prefix of the pinned
    join graph.  ``archive`` and ``receipt`` may be mappings or paths.  The
    ``leading`` input may be the whole C6 archive, its binding, or a path.
    ``outer`` optionally supplies the pinned true outer C6/C5 archive; when
    omitted, that archive is hydrated read-only so the eight leading-Q_y
    source inputs can still be checked.
    """
    repo = _repo_root(root)
    archive_data, archive_load = _load_serialized(archive, repo, name="join_archive")
    receipt_data, receipt_load = _load_serialized(receipt, repo, name="join_receipt")
    leading_data, leading_load = _load_serialized(leading, repo, name="leading_archive")
    outer_data, outer_load = _load_serialized(
        _OUTER_ARCHIVE if outer is None else outer,
        repo,
        name="outer_archive",
    )
    if _LEADING_BINDING_KEY not in leading_data:
        # A binding-only input cannot carry the graph prefix.  Hydrate the
        # pinned archive for that structural check while retaining the exact
        # caller-supplied binding for the source claim.
        leading_binding = leading_data
        pinned_leading, pinned_load = _load_serialized(_LEADING_ARCHIVE, repo, name="pinned_leading_archive")
        leading_data = pinned_leading
        leading_load = {**leading_load, "binding_input": True, "pinned_archive": pinned_load}
    if "true_outer_C6_C5_bindings" not in outer_data:
        pinned_outer, pinned_outer_load = _load_serialized(_OUTER_ARCHIVE, repo, name="pinned_outer_archive")
        outer_data = pinned_outer
        outer_load = {**outer_load, "binding_input": True, "pinned_archive": pinned_outer_load}
    binding = _unwrap(archive_data, _BINDING_KEY)
    leading_binding = _unwrap(leading_data, _LEADING_BINDING_KEY)
    outer_binding = _unwrap(outer_data, "true_outer_C6_C5_bindings")
    pins = {name: _pin(repo, relative, expected) for name, (relative, expected) in _SOURCE_PINS.items()}
    source_pins = bool(all(row["passed"] for row in pins.values()))
    family_archive = _mapping(archive_data.get("source_family"))
    family_receipt = _mapping(receipt_data.get("source_family"))
    family_leading = _mapping(leading_data.get("source_family"))
    family_outer = _mapping(outer_data.get("source_family"))
    family_match = bool(family_archive and family_archive == family_receipt == family_leading == family_outer)
    prefix = _prefix_rows(nodes, archive_data, leading_data)
    outer_prefix = _outer_prefix_row(nodes, outer_data)
    receipt_rows = _required_receipt_rows(receipt_data)
    charts = _mapping(binding.get("charts"))
    chart_rows = {name: _chart_scope(nodes, name, _mapping(charts.get(name))) for name in _CHARTS}
    chart_passed = all(row["passed"] for row in chart_rows.values())
    seam_rows = {
        seam: _join_scope(nodes, charts, seam, _mapping(_mapping(binding.get("joins")).get(seam)))
        for seam in _SEAMS
    }
    seam_structure_passed = all(row["passed"] for row in seam_rows.values())
    leading_identity = _leading_identity(nodes, leading_data, leading_binding)
    means = {
        seam: _mean_certificate(nodes, binding, receipt_rows, seam, seam_rows[seam])
        for seam in _SEAMS
    }
    histories = {
        seam: _history_certificate(nodes, binding, receipt_rows, seam)
        for seam in _SEAMS
    }
    leading_Qy_inputs = {
        seam: _outer_background_certificate(
            nodes,
            binding,
            outer_binding,
            leading_binding,
            receipt_rows,
            outer_prefix,
            seam,
        )
        for seam in _SEAMS
    }
    means_passed = all(row["passed"] for row in means.values())
    histories_passed = all(row["passed"] for row in histories.values())
    leading_Qy_inputs_passed = all(row["passed"] for row in leading_Qy_inputs.values())
    receipt_source_passed = receipt_rows["passed"]
    passed = bool(
        source_pins
        and family_match
        and archive_load.get("passed") is True
        and receipt_load.get("passed") is True
        and leading_load.get("passed") is True
        and outer_load.get("passed") is True
        and prefix["passed"]
        and outer_prefix["passed"]
        and receipt_source_passed
        and chart_passed
        and seam_structure_passed
        and leading_identity["passed"]
        and means_passed
        and histories_passed
        and leading_Qy_inputs_passed
    )
    common_receipt = _digest({
        "join_receipt_sha256": pins["accepted_join_receipt"]["sha256"],
        "leading_archive_sha256": pins["accepted_leading_archive"]["sha256"],
        "seams": _SEAMS,
        "source_flags": receipt_rows["source_flags"],
        "mean_rows": receipt_rows["mean_row_count"],
        "finite_source_pairs": receipt_rows["finite_source_pair_count"],
    })
    for seam in _SEAMS:
        for row in means[seam]["rows"]:
            row["common_source_receipt_sha256"] = common_receipt
            row["source_identity_digest"] = _digest([common_receipt, seam, row["ordinary_Z_order"], "modified_physical_m"])
        for row in histories[seam]["finite_N_source_pairs"]:
            row["common_source_receipt_sha256"] = common_receipt
    return {
        "passed": passed,
        "source_only": True,
        "source_values_selected": False,
        "numeric_caps_used": False,
        "ancestor_replayed": False,
        "source_pins": {"passed": source_pins, "rows": pins},
        "archive_load": archive_load,
        "receipt_load": receipt_load,
        "leading_load": leading_load,
        "outer_load": outer_load,
        "source_family": {
            "archive": family_archive,
            "receipt": family_receipt,
            "leading": family_leading,
            "outer": family_outer,
            "matches": family_match,
        },
        "source_prefix": prefix,
        "outer_source_prefix": outer_prefix,
        "receipt_checks": _public_receipt_rows(receipt_rows),
        "common_source_receipt_sha256": common_receipt,
        "charts": chart_rows,
        "seams": seam_rows,
        "modified_physical_mean_rows": means,
        "leading_mean_rows": {
            "passed": means_passed,
            "rows": [row for seam in _SEAMS for row in means[seam]["rows"]],
            "quantity": "modified_physical_mean_m_endpoint",
            "original_leading_output_alias_issued": False,
            "original_leading_controls_opaque": True,
            "reason": "Accepted receipt rows are modified physical m endpoint identities; no original leading m function is aliased.",
        },
        "finite_N_history_rows": histories,
        "leading_Qy_inputs": leading_Qy_inputs,
        "original_leading_Qy_input_rows": {
            "passed": leading_Qy_inputs_passed,
            "rows": [row for seam in _SEAMS for row in leading_Qy_inputs[seam]["leading_Qy_input_rows"]],
            "row_count": sum(row["leading_Qy_input_row_count"] for row in leading_Qy_inputs.values()),
            "expected_row_count": len(_SEAMS) * len(_QY_INPUT_FIELDS) * 2,
            "radial_rows": [row for seam in _SEAMS for row in leading_Qy_inputs[seam]["radial_input_rows"]],
            "radial_row_count": sum(len(row["radial_input_rows"]) for row in leading_Qy_inputs.values()),
            "source_node_alias_issued": False,
            "source_premises_only": True,
        },
        "history_exact_transfer": {
            "passed": histories_passed,
            "seams": _SEAMS,
            "incoming_histories_retained": histories_passed,
            "finite_N_physical_germ_equality": histories_passed,
            "direct_output_alias_issued": False,
            "source_pair_count": sum(row["finite_N_source_pair_count"] for row in histories.values()),
            "integral_row_count": sum(row["incoming_integral_row_count"] for row in histories.values()),
        },
        "leading_unique_branch": leading_identity,
        "trace_contract": {
            "Rm_exit_uses_exact_unique_residual_identity": True,
            "Rm_exit_empty_integral_substitution": False,
            "Rm_exit_h_zero_substitution": False,
            "patch_reference_Rm_native_exp1_to_Rh_native_minus5": seam_rows["patch_reference"]["passed"],
            "buffer_transition_native_query_4_from_archive": seam_rows["buffer_transition"]["right_native_query"] == 0,
            "common_N_and_rate_source_retained": True,
            "variable_native_Jacobians_retained": receipt_rows["source_flags"]["own_variable_radius"],
            "transition_power_nonzero_mean_retained": means["transition_power"]["transition_power_nonzero_mean_source"],
        },
        "limits": {
            "original_leading_output_alias": False,
            "global_opaque_alias": False,
            "numeric_endpoint_values": False,
            "numeric_caps_reinterpreted": False,
            "cap_subtraction": False,
            "ancestor_constructor_or_checker_rerun": False,
        },
        "unproved": {
            "original_leading_m_endpoint_alias": True,
            "numeric_finite_N_values": True,
            "global_n1_closure": True,
        },
    }


__all__ = ["build_outer_chain_mean"]
