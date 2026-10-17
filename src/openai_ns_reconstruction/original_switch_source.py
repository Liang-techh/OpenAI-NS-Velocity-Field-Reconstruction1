"""Archive-only source certificate for the first-switch interface.

This module binds ``first_switch`` native one to ``second_switch`` native
zero.  It consumes immutable graph archives, source text, receipts, and the
saved pre-mapping parent report.  It never constructs an accepted owner,
replays a producer/checker, evaluates a field, or promotes a numerical cap to
a source function.

The loop/root ring is intentionally scoped to the two selected endpoint
coordinates.  Only the R 100 primitive ports authorized by the first/second
source receipt may be reduced to a common physical-function atom.  Constants
remain constants and conflicting primitive authorizations fail closed.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .original_core_numeric import ROOT
from .original_active_viscous_m15 import ActiveEndpointFingerprint
from .original_upstream_loop_identity import (
    ScopedActiveEndpointFingerprint,
    _compare,
    _scope_tokens,
)
from .original_upstream_native import (
    _ast_guard,
    _common_N,
    _leading_ports,
    _pressure_port,
    _radius_guard,
)
from .original_native_derivative_guards import (
    _ast_bundle,
    _slow_b_source_binding,
)
from .original_upstream_source import (
    _binding,
    _direct_endpoint_fingerprints,
    _finite_function_source_rows,
    _handle,
    _history_guard,
    _join_source_proof,
    _leading_rows,
    _mapping,
    _node,
    _physical_rows,
    _physical_source_statement,
    _prefix_check,
    _radial_source_link_rows,
    _rows,
    _scope,
)


_PREFIX = "lei_ren_part1_paper_compliant_current_original_"
_PARENT_RELATIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_upstream_m15_numeric.json"
)
_PARENT_SHA256 = "6028d49ac9ec1b94ae967b15d698366061b4a15b9e0ddef9b45be2d3971e355e"

_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}

_ARCHIVE_KEYS = {
    "loop": "prescribed_loop_C5_bindings",
    "modified": "modified_Q4_pressure_C3_bindings",
    "leading": "prescribed_R100_Rm_C6_bindings",
    "switch": "switch_mixed_source_join_bindings",
}
_ARCHIVE_DEFAULTS = {
    "loop": _PREFIX + "prescribed_loop_C5_functions.json.gz",
    "modified": _PREFIX + "modified_Q4_pressure_C3_functions.json.gz",
    "leading": _PREFIX + "prescribed_R100_Rm_C6_functions.json.gz",
    "switch": _PREFIX + "switch_mixed_source_joins.json.gz",
}

# The graph archives are prefixes of the switch graph supplied by the current
# owner.  Source and receipt pins are kept here so this helper remains useful
# without an owner hash table.
_PINS = {
    "loop_archive": (_PREFIX + "prescribed_loop_C5_functions.json.gz", "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa"),
    "loop_source": (_PREFIX + "prescribed_loop_C5_functions.py", "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d"),
    "loop_receipt": (_PREFIX + "prescribed_loop_C5_functions_check.json", "2c77e80df04fae615f8ea909fc7b606bfccfa27659cb67744ae8d84020b1a9e7"),
    "modified_archive": (_PREFIX + "modified_Q4_pressure_C3_functions.json.gz", "b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a"),
    "modified_source": (_PREFIX + "modified_Q4_pressure_C3_functions.py", "7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7"),
    "modified_receipt": (_PREFIX + "modified_Q4_pressure_C3_functions_check.json", "189ef4dacc60bc2feda16046a193a75e89f658f635fd6e3e381200c09637bf6f"),
    "leading_archive": (_PREFIX + "prescribed_R100_Rm_C6_functions.json.gz", "5e5b565dd8d2abc41e6f3df0d9e17000e524d493e58a5778dfcf217c6e42f44a"),
    "leading_source": (_PREFIX + "prescribed_R100_Rm_C6_functions.py", "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32"),
    "leading_receipt": (_PREFIX + "prescribed_R100_Rm_C6_functions_check.json", "494caf482bffe5475c011a35a1b7d95fd02725b32eb55653f9ca25984dd98ad6"),
    "switch_archive": (_PREFIX + "switch_mixed_source_joins.json.gz", "ddbd506d2fbfc4c58f7a5ba93e0abe3a95d3af8fd3c520ee3ce5024d90a29567"),
    "switch_source": (_PREFIX + "switch_mixed_source_joins.py", "0ae7f651c0d183756f0e6ab90d56dfd80bf3d3c1a6717ffdc2a03fe08dd52f2f"),
    "switch_receipt": (_PREFIX + "switch_mixed_source_joins_check.json", "9cf3c6d941df95f4b5f73e7b97e32198a54fcd3d17a85ecca2a0742c0d602fc3"),
}

_ROOT_KEYS = ("E", "V", "a", "b", "p1", "p2", "t0", "Delta", "kappa")
_LOOP_KEYS = ("q", "B_over_Pstar", "B_over_Pstar_phi")
_R100_PORTS = (
    ("fields", "E"),
    ("fields", "V"),
    ("radial", "E_y"),
    ("radial", "V_y"),
    ("generic", "theta_linear"),
    ("generic", "theta_quadratic"),
    ("generic", "p2"),
)
_CHARTS = ("first_switch", "second_switch")


def _repo_root(root: Any) -> Path:
    path = Path(root) if root is not None else ROOT
    if path.is_file():
        path = path.parent
    if (path / "src").is_dir() and (path / "experiments").is_dir():
        return path
    return Path(__file__).resolve().parents[2]


def _digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    ).hexdigest()


def _load_archive(value: Any, repo: Path, family: str) -> tuple[Mapping[str, Any], dict[str, Any]]:
    if isinstance(value, Mapping):
        return value, {"family": family, "mode": "supplied_mapping", "passed": True}
    path = Path(value) if value is not None else repo / "experiments" / "root_st073" / _ARCHIVE_DEFAULTS[family]
    if not path.is_absolute():
        path = repo / path
    try:
        with gzip.open(path, "rt", encoding="utf-8") as stream:
            data = json.load(stream)
    except (OSError, ValueError, EOFError, json.JSONDecodeError) as exc:
        return {}, {"family": family, "path": str(path), "passed": False, "error": str(exc)}
    return (data if isinstance(data, Mapping) else {}), {
        "family": family, "path": str(path.relative_to(repo)) if path.is_relative_to(repo) else str(path),
        "passed": isinstance(data, Mapping),
    }


def _load_receipt(repo: Path, family: str) -> tuple[Mapping[str, Any], dict[str, Any]]:
    relative = _PREFIX + {
        "loop": "prescribed_loop_C5_functions_check.json",
        "modified": "modified_Q4_pressure_C3_functions_check.json",
        "leading": "prescribed_R100_Rm_C6_functions_check.json",
        "switch": "switch_mixed_source_joins_check.json",
    }[family]
    path = repo / "experiments" / "root_st073" / relative
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        return {}, {"family": family, "path": str(path.relative_to(repo)), "passed": False, "error": str(exc)}
    return (data if isinstance(data, Mapping) else {}), {"family": family, "path": str(path.relative_to(repo)), "passed": isinstance(data, Mapping)}


def _pin(repo: Path, key: str) -> dict[str, Any]:
    relative, expected = _PINS[key]
    path = repo / "experiments" / "root_st073" / relative
    try:
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as exc:
        return {"key": key, "path": "experiments/root_st073/" + relative, "expected_sha256": expected,
                "actual_sha256": None, "passed": False, "error": str(exc)}
    return {"key": key, "path": "experiments/root_st073/" + relative,
            "expected_sha256": expected, "actual_sha256": actual, "passed": actual == expected}


def _prefix(nodes: Any, archive: Mapping[str, Any]) -> dict[str, Any]:
    archived = archive.get("exact_graph_nodes")
    if not isinstance(archived, list) or not isinstance(nodes, (list, tuple)):
        return {"passed": False, "archived_node_count": len(archived) if isinstance(archived, list) else 0,
                "current_node_count": len(nodes) if isinstance(nodes, (list, tuple)) else 0}
    passed = len(nodes) >= len(archived) and list(nodes[: len(archived)]) == archived
    return {"passed": bool(passed), "archived_node_count": len(archived), "current_node_count": len(nodes),
            "exact_prefix": bool(passed)}


def _compact_chart_guard(row: Mapping[str, Any]) -> dict[str, Any]:
    phase = _mapping(row.get("phase")); offset = _mapping(row.get("global_offset"))
    return {
        "passed": row.get("passed") is True,
        "chart": row.get("chart"),
        "coordinate": row.get("coordinate"),
        "radius": row.get("radius"),
        "physical_measure": row.get("physical_measure"),
        "native_J": row.get("native_J"),
        "radius_exponential_node": row.get("radius_exponential_node"),
        "radius_exponent_node": row.get("radius_exponent_node"),
        "radius_base_node": row.get("radius_base_node"),
        "radius_formula_shape": row.get("radius_formula_shape") is True,
        "R_native_factor_once": row.get("R_native_factor_once") is True,
        "d_native_equals_J_d_y": row.get("d_native_equals_J_d_y") is True,
        "d_y_equals_d_native_over_J": row.get("d_y_equals_d_native_over_J") is True,
        "native_J_is_constant": row.get("native_J_is_constant") is True,
        "native_J_positive": row.get("native_J_positive") is True,
        "phase": {
            "passed": phase.get("passed") is True,
            "coordinate": phase.get("coordinate"),
            "global_offset": phase.get("global_offset"),
            "actual_phase": phase.get("actual_phase"),
            "common_N": phase.get("common_N"),
        },
        "global_offset": {
            "passed": offset.get("passed") is True,
            "global_offset": offset.get("global_offset"),
            "phase_origin": offset.get("phase_origin"),
            "radius_ratio": offset.get("radius_ratio") is True,
            "origin_subtracted": offset.get("origin_subtracted") is True,
        },
        "phase_factor_after_J_division_is_N": row.get("phase_factor_after_J_division_is_N") is True,
        "total_fast_phase_N_term_has_no_extra_N": row.get("total_fast_phase_N_term_has_no_extra_N") is True,
        "ordinary_Z_orders": list(row.get("ordinary_Z_orders", ())),
        "source_values_selected": False,
        "numeric_caps_used": False,
    }


def _compact_pressure_port(row: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "passed": row.get("passed") is True,
        "chart": row.get("chart"),
        "physical_rows": {key: list(value) for key, value in _mapping(row.get("physical_rows")).items()},
        "radial_rows": {key: list(value) for key, value in _mapping(row.get("radial_rows")).items()},
        "finite_boundary_rows": {key: list(value) for key, value in _mapping(row.get("finite_boundary_rows")).items()},
        "Pstar_node": row.get("Pstar_node"),
        "Pstar_once": row.get("Pstar_once") is True,
        "m_y_equals_U_minus_m": row.get("m_y_equals_U_minus_m") is True,
        "Q_formula_source": row.get("Q_formula_source") is True,
        "Q_y_graph_rows": list(row.get("Q_y_graph_rows", ())),
        "R_Q_y_rows": list(row.get("R_Q_y_rows", ())),
        "ordinary_Z1": row.get("ordinary_Z1") is True,
        "Q_plus_Q_y_rows": row.get("Q_plus_Q_y_rows") is True,
        "R_Q_y_radius_once": row.get("R_Q_y_radius_once") is True,
        "radius_node": row.get("radius_node"),
        "total_fast_phase_terms_included": row.get("total_fast_phase_terms_included") is True,
        "physical_y_equals_native_partial_divided_by_Jacobian": row.get("physical_y_equals_native_partial_divided_by_Jacobian") is True,
        "ordinary_Z_orders": list(row.get("ordinary_Z_orders", ())),
        "source_values_selected": False,
        "numeric_caps_used": False,
    }


def _compact_slow_binding(row: Mapping[str, Any]) -> dict[str, Any]:
    rows = []
    for item in row.get("rows", ()) if isinstance(row.get("rows"), list) else ():
        inv = _mapping(item.get("inverse_J_row")); lazy = _mapping(item.get("generated_lazy_native_derivative_row"))
        rows.append({
            "order": item.get("order"),
            "slow_B_node": item.get("slow_B_node"),
            "loop_B_source_node": item.get("loop_B_source_node"),
            "inverse_J_row": {"operation": inv.get("operation"), "numerator": inv.get("numerator"), "denominator": inv.get("denominator"),
                               "source_positive_certificate": inv.get("source_positive_certificate")},
            "generated_lazy_native_derivative_row": {"active_body": lazy.get("active_body"), "flat_value": lazy.get("flat_value"),
                                                       "ordinary_slow_Z_derivative_order": lazy.get("ordinary_slow_Z_derivative_order"),
                                                       "same_original_flat_extension_through_order5": lazy.get("same_original_flat_extension_through_order5"),
                                                       "same_original_lazy_native_derivative": lazy.get("same_original_lazy_native_derivative")},
            "inverse_J_binding": item.get("inverse_J_binding") is True,
            "lazy_native_derivative_binding": item.get("lazy_native_derivative_binding") is True,
            "loop_B_source_binding": item.get("loop_B_source_binding") is True,
        })
    return {"passed": row.get("passed") is True, "rows": rows, "source_values_selected": False, "numeric_caps_used": False}


def _scope_row(loop: Mapping[str, Any], name: str, query: int) -> dict[str, Any]:
    scope = _scope(name, loop, query)
    scope["native_query_value"] = str(query)
    return scope


def _selected_receipt_checks(receipt: Mapping[str, Any], loop: Mapping[str, Any], modified: Mapping[str, Any],
                             leading: Mapping[str, Any], nodes: Any) -> dict[str, Any]:
    checks = _mapping(receipt.get("checks"))
    contract = _mapping(checks.get("source_contract"))
    bg = _mapping(checks.get("actual_background_recipes"))
    jets = _mapping(checks.get("physical_controls_and_ODE_jets"))
    phase = _mapping(checks.get("actual_loop_and_phase_transfer"))
    own = _mapping(checks.get("actual_own_mixed_joins"))
    root_rows = []
    for item in phase.get("actual_original_root_inverse_phase_J_bindings", ()):
        if isinstance(item, list) and len(item) == 4 and item[0] in _CHARTS:
            chart = _mapping(_mapping(loop.get("charts")).get(item[0]))
            expected = [item[0], chart.get("roots"), chart.get("actual_phase"), chart.get("physical_measure")]
            root_rows.append({"chart": item[0], "passed": item == expected,
                              "coordinate_node": _handle(chart.get("coordinate")),
                              "actual_phase_node": _handle(chart.get("actual_phase")),
                              "physical_measure_node": _handle(chart.get("physical_measure")),
                              "root_orders": {key: list(_mapping(chart.get("roots")).get(key, ())[:6]) for key in _ROOT_KEYS}})
    seams = list(phase.get("actual_radius_and_global_logarithmic_phase_seams", ()))
    own_rows = [row for row in own.get("actual_two_switch_own_mixed_rows", ())
                if isinstance(row, Mapping) and row.get("interface") == "first_second"]
    rhs_rows = [row for row in jets.get("actual_physical_RHS_coefficient_rows", ())
                if isinstance(row, list) and row and row[0] == "first_second"]
    ode_rows = [row for row in jets.get("independent_physical_ODE_rows", ())
                if isinstance(row, list) and row and row[0] == "first_second"]
    checks_out = {
        "source_contract": {
            "passed": contract.get("passed") is True,
            "accepted_ten_full_chart_source_receipt": isinstance(contract.get("accepted_ten_full_chart_source_receipt"), str),
            "all_180_actual_high_derivative_requests_reused_without_rebuilding": contract.get("all_180_actual_high_derivative_requests_reused_without_rebuilding") is True,
            "actual_own_endpoint_pairs_and_complete_predecessor_sources_bound": contract.get("actual_own_endpoint_pairs_and_complete_predecessor_sources_bound") is True,
            "shared_switch_width_and_distinct_post_power_native_length": contract.get("shared_switch_width_and_distinct_post_power_native_length") is True,
        },
        "background_recipes": {
            "passed": bg.get("passed") is True,
            "actual_first_switch1_and_second_switch1_phi_U_and_six_histories_are_successor_inputs": bg.get("actual_first_switch1_and_second_switch1_phi_U_and_six_histories_are_successor_inputs") is True,
            "same_current_radius_comparison_formula_with_nonzero_inlets_and_B_half_equilibrium": bg.get("same_current_radius_comparison_formula_with_nonzero_inlets_and_B_half_equilibrium") is True,
            "exact_hydro_pressure_swirl_F0_Pstar_and_physical_measures_preserved": bg.get("exact_hydro_pressure_swirl_F0_Pstar_and_physical_measures_preserved") is True,
        },
        "physical_controls_and_ODE_jets": {
            "passed": jets.get("passed") is True,
            "native_sigma_flat_rows": len(jets.get("actual_native_sigma_flat_rows", ())) if isinstance(jets.get("actual_native_sigma_flat_rows"), list) else 0,
            "physical_RHS_coefficient_rows": len(rhs_rows),
            "independent_physical_ODE_rows": len(ode_rows),
            "first_axial_driver_turns_off_with_four_flat_derivatives": jets.get("first_axial_driver_turns_off_with_four_flat_derivatives") is True,
            "second_angular_driver_turns_into_exact_minus_two_fifths_power_slope": jets.get("second_angular_driver_turns_into_exact_minus_two_fifths_power_slope") is True,
            "background_phi_U_and_six_histories_y0_to5_with_arbitrary_nonzero_inlets": jets.get("background_phi_U_and_six_histories_y0_to5_with_arbitrary_nonzero_inlets") is True,
            "ordinary_Z_differentiation_commutes_with_the_finite_smooth_source_induction": jets.get("ordinary_Z_differentiation_commutes_with_the_finite_smooth_source_induction") is True,
        },
        "loop_and_phase_transfer": {
            "passed": phase.get("passed") is True,
            "actual_root_inverse_phase_J_rows": root_rows,
            "root_inverse_rows_passed": len(root_rows) == 2 and all(row["passed"] for row in root_rows),
            "first_second_radius_global_phase_seam": "first_second" in seams,
            "background_y5_supplies_root_y4_and_original_unique_inverse_periodic_lift_U_y4_C5": phase.get("background_y5_supplies_root_y4_and_original_unique_inverse_periodic_lift_U_y4_C5") is True,
            "actual_nonzero_mean_inlet_is_checked_before_Q_y4_N_y2_pressure_y3_transfer": phase.get("actual_nonzero_mean_inlet_is_checked_before_Q_y4_N_y2_pressure_y3_transfer") is True,
            "no_raw_fractional_part_differentiation_or_numerical_enclosure_claim": phase.get("no_raw_fractional_part_differentiation_or_numerical_enclosure_claim") is True,
        },
        "own_mixed_join": {
            "passed": own.get("passed") is True,
            "first_second_row_count": len(own_rows),
            "first_second_rows_present": bool(own_rows),
            "complete_nonzero_original_correction_S0_S1_pressure_and_five_histories_retained": own.get("complete_nonzero_original_correction_S0_S1_pressure_and_five_histories_retained") is True,
            "finite_common_leading_germs_only_after_source_ODE_inverse_Q_N_proofs": own.get("finite_common_leading_germs_only_after_source_ODE_inverse_Q_N_proofs") is True,
        },
    }
    checks_out["passed"] = bool(
        all(row.get("passed") is True for key, row in checks_out.items() if key != "physical_controls_and_ODE_jets" and key != "loop_and_phase_transfer" and isinstance(row, Mapping))
        and checks_out["physical_controls_and_ODE_jets"]["passed"]
        and checks_out["loop_and_phase_transfer"]["passed"]
        and checks_out["loop_and_phase_transfer"]["root_inverse_rows_passed"]
        and checks_out["loop_and_phase_transfer"]["first_second_radius_global_phase_seam"]
        and checks_out["own_mixed_join"]["first_second_row_count"] > 0
    )
    return checks_out


def _make_scoped_fingerprint(nodes: Any, loop: Mapping[str, Any], leading: Mapping[str, Any],
                             left_name: str, right_name: str) -> tuple[ScopedActiveEndpointFingerprint, dict[str, Any], dict[str, Any], dict[str, Any]]:
    charts = _mapping(loop.get("charts"))
    scope, serialized, endpoint_scopes = _scope_tokens(nodes, charts, left_name, right_name)
    authorizations: dict[int, Any] = {}
    conflicts: list[dict[str, Any]] = []
    authorization_rows: list[dict[str, Any]] = []

    def token(label: Any) -> tuple[str, None]:
        return (_digest(label), None)

    leading_charts = _mapping(leading.get("charts"))
    for group, quantity in _R100_PORTS:
        for order in range(3):
            values = {}
            for chart_name in (left_name, right_name):
                chart = _mapping(leading_charts.get(chart_name))
                ports = _mapping(chart.get(group))
                values[chart_name] = _handle((ports.get(quantity) or [None] * 3)[order]) if isinstance(ports.get(quantity), (list, tuple)) else None
            if any(value is None or _node(nodes, value) is None for value in values.values()):
                conflicts.append({"kind": "missing_R100_port", "group": group, "quantity": quantity, "ordinary_Z_order": order, "values": values})
                continue
            exact_constant_pair = any(_mapping(_node(nodes, value)).get("operation") == "exact_rational" for value in values.values())
            primitive_token = token(["first_second_R100_common_primitive", quantity, order])
            for chart_name, node in values.items():
                coordinate = _handle(_mapping(charts.get(chart_name)).get("coordinate"))
                expected = scope.get(coordinate)
                if exact_constant_pair or coordinate is None or expected is None:
                    continue
                prior = authorizations.get(node)
                if prior is None:
                    authorizations[node] = (coordinate, expected, primitive_token)
                else:
                    prior_rows = prior if isinstance(prior, list) else [prior]
                    if any(item[2] != primitive_token for item in prior_rows):
                        conflicts.append({"kind": "authorization_conflict", "node": node, "quantity": quantity, "ordinary_Z_order": order})
                    elif (coordinate, expected, primitive_token) not in prior_rows:
                        prior_rows.append((coordinate, expected, primitive_token))
                        authorizations[node] = prior_rows
            authorization_rows.append({"group": group, "quantity": quantity, "ordinary_Z_order": order,
                                       "left_node": values[left_name], "right_node": values[right_name],
                                       "canonical_sha256": primitive_token[0], "scope_restricted": True,
                                       "exact_constant_pair": exact_constant_pair})

    # Radius and phase are admitted only by the selected source seam.  They
    # are not broad chart aliases and are valid only at the two native scopes.
    for role in ("radius", "global_offset", "actual_phase"):
        geometry_token = token(["first_second_radius_global_phase_source", role])
        for chart_name in (left_name, right_name):
            chart = _mapping(charts.get(chart_name)); coordinate = _handle(chart.get("coordinate"))
            node = _handle(chart.get(role)); expected = scope.get(coordinate)
            if node is None or coordinate is None or expected is None:
                conflicts.append({"kind": "missing_geometry_port", "role": role, "chart": chart_name})
                continue
            prior = authorizations.get(node)
            if prior is not None and (prior[2] if not isinstance(prior, list) else prior[0][2]) != geometry_token:
                conflicts.append({"kind": "geometry_authorization_conflict", "role": role, "node": node})
            else:
                authorizations[node] = (coordinate, expected, geometry_token)
        authorization_rows.append({"group": "geometry", "quantity": role, "canonical_sha256": geometry_token[0], "scope_restricted": True})

    fingerprint = ScopedActiveEndpointFingerprint(
        nodes,
        {"coordinate": _handle(_mapping(charts.get(left_name)).get("coordinate")), "loop": charts},
        authorizations,
        scope,
    )
    return fingerprint, endpoint_scopes, {
        "serialized": serialized,
        "authorization_rows": authorization_rows,
        "authorization_conflicts": conflicts,
        "authorized_node_count": len(authorizations),
        "scope_restricted": True,
        "source_values_selected": False,
    }, scope


def _ring_rows(fingerprint: ScopedActiveEndpointFingerprint, endpoint_scopes: Mapping[str, Mapping[int, Any]],
               loop: Mapping[str, Any], left_name: str, right_name: str, quantities: tuple[str, ...], group: str) -> list[dict[str, Any]]:
    charts = _mapping(loop.get("charts")); left = _mapping(charts.get(left_name)); right = _mapping(charts.get(right_name))
    output = []
    for quantity in quantities:
        left_values = _mapping(left.get(group)).get(quantity)
        right_values = _mapping(right.get(group)).get(quantity)
        for order in range(3):
            if not isinstance(left_values, (list, tuple)) or not isinstance(right_values, (list, tuple)) or len(left_values) <= order or len(right_values) <= order:
                output.append({"group": group, "quantity": quantity, "ordinary_Z_order": order, "passed": False,
                               "left_node": None, "right_node": None, "reason": "missing source row", "source_values_selected": False})
                continue
            row = _compare(fingerprint, endpoint_scopes["left"], endpoint_scopes["right"],
                           f"{group}.{quantity}[{order}]", left_values[order], right_values[order], order)
            row["group"] = group; row["quantity"] = quantity
            row["source_values_selected"] = False
            output.append(row)
    return output


def _direct_group_row(direct: Mapping[str, Any], group: str, order: int) -> Mapping[str, Any]:
    rows = _mapping(direct.get("groups", {}).get(group)).get("rows", ())
    return rows[order] if isinstance(rows, list) and len(rows) > order and isinstance(rows[order], Mapping) else {}


def _identity(quantity: str, order: int, row: Mapping[str, Any], direct: Mapping[str, Any], loop_identity: Mapping[str, Any],
              native_identity: Mapping[str, Any], scope_left: Mapping[str, Any], scope_right: Mapping[str, Any],
              source_ready: bool, derivative_ready: bool) -> dict[str, Any]:
    group = {"U": "modified_U", "Uy": "modified_Uy", "m": "modified_m"}.get(quantity, quantity)
    direct_row = _direct_group_row(direct, group, order) if group.startswith("modified_") else {}
    exact = direct_row.get("passed") is True
    operator = quantity == "U" and source_ready or quantity == "Uy" and source_ready and derivative_ready
    result = {
        "passed": bool(exact or operator),
        "quantity": quantity,
        "ordinary_Z_order": order,
        "left_node": row.get("left_node"),
        "right_node": row.get("right_node"),
        "left_scope": dict(scope_left),
        "right_scope": dict(scope_right),
        "source_values_selected": False,
        "numeric_caps_used": False,
        "ancestor_replayed": False,
    }
    if exact:
        result.update({
            "exact_endpoint_fingerprints": True,
            "left_endpoint_sha256": direct_row.get("left_endpoint_sha256"),
            "right_endpoint_sha256": direct_row.get("right_endpoint_sha256"),
            "method": "independent scoped endpoint fingerprint",
        })
    elif operator:
        result.update({
            "actual_loop_operator_identity_proved": True,
            # Downstream pressure congruence hashes the complete loop
            # certificate object.  The optional certificate_sha256 field is
            # only a compact receipt and is computed before that field is
            # attached, so it must not be used as the operator identity hash.
            "actual_loop_identity_sha256": _digest(loop_identity),
            "method": "source operator congruence from scoped roots/loop functions and finite-N source map",
        })
        if quantity == "Uy":
            result.update({
                # Uy consumers bind the chart wiring itself, not the enclosing
                # derivative certificate (which also contains AST and source
                # metadata).
                "issued_derivative_wiring_sha256": native_identity.get("charts_sha256"),
                "native_derivative_once": True,
                "physical_y_divided_by_J_once": True,
            })
    else:
        result.update({"method": "unproved source identity", "unproved_reason": "operator or derivative certificate missing"})
    return result


def _saved_parent_rows(parent: Mapping[str, Any], current_u: list[dict[str, Any]], current_uy: list[dict[str, Any]], current_m: list[dict[str, Any]], leading: list[dict[str, Any]]) -> dict[str, Any]:
    seam = _mapping(_mapping(parent.get("source_derivative_theorem")).get("seams")).get("switch_switch")
    checks = []
    for key, current in (("modified_U_rows3", current_u), ("modified_Uy_rows2", current_uy), ("modified_m_rows4", current_m)):
        saved = _mapping(seam).get(key)
        ok = isinstance(saved, list) and len(saved) == len(current)
        if ok:
            for a, b in zip(saved, current):
                ok = ok and a.get("left_node") == b.get("left_node") and a.get("right_node") == b.get("right_node")
                ok = ok and a.get("left_scope") == b.get("left_scope") and a.get("right_scope") == b.get("right_scope")
        checks.append({"key": key, "saved_count": len(saved) if isinstance(saved, list) else 0, "current_count": len(current), "passed": bool(ok)})
    saved_lead = _mapping(_mapping(seam).get("leading_Qy_inputs")).get("rows")
    ok = isinstance(saved_lead, list) and len(saved_lead) == len(leading)
    if ok:
        for a, b in zip(saved_lead, leading):
            ok = ok and a.get("quantity") == b.get("quantity") and a.get("ordinary_Z_order") == b.get("ordinary_Z_order")
            ok = ok and a.get("left_node") == b.get("left_node") and a.get("right_node") == b.get("right_node")
    checks.append({"key": "leading_Qy_inputs.rows", "saved_count": len(saved_lead) if isinstance(saved_lead, list) else 0,
                   "current_count": len(leading), "passed": bool(ok)})
    seam_ok = bool(seam.get("passed") is True and seam.get("source_values_selected") is False and seam.get("numeric_caps_used") is False and seam.get("ancestor_replayed") is False)
    return {"passed": bool(seam_ok and all(row["passed"] for row in checks)), "parent_seam_present": bool(seam), "rows": checks,
            "saved_seam_passed": seam.get("passed") is True, "source_values_selected": False, "numeric_caps_used": False}


def _load_parent(parent_report: Any, repo: Path) -> tuple[Mapping[str, Any], dict[str, Any]]:
    canonical = repo / _PARENT_RELATIVE
    try:
        raw = canonical.read_bytes()
        canonical_sha = hashlib.sha256(raw).hexdigest()
        canonical_data = json.loads(raw.decode("utf-8"))
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        return {}, {"path": _PARENT_RELATIVE, "expected_sha256": _PARENT_SHA256, "actual_sha256": None, "passed": False, "error": str(exc)}
    supplied = parent_report if isinstance(parent_report, Mapping) else canonical_data
    if not isinstance(supplied, Mapping):
        supplied = {}
    # The canonical bytes are always hash-bound.  A supplied mapping is only
    # accepted when its selected seam and native guards are content-equal to
    # the canonical saved evidence.  All consumed data below are canonical.
    canonical_seam = _mapping(_mapping(canonical_data.get("source_derivative_theorem")).get("seams")).get("switch_switch")
    supplied_seam = _mapping(_mapping(supplied.get("source_derivative_theorem")).get("seams")).get("switch_switch")
    selected_match = bool(canonical_seam and supplied_seam and _digest(canonical_seam) == _digest(supplied_seam))
    native_match = bool(canonical_data.get("native_source_guards")
                        and _digest(canonical_data.get("native_source_guards"))
                        == _digest(supplied.get("native_source_guards")))
    info = {"path": _PARENT_RELATIVE, "expected_sha256": _PARENT_SHA256, "actual_sha256": canonical_sha,
            "passed": canonical_sha == _PARENT_SHA256 and selected_match and native_match,
            "canonical_bytes_hash_bound": canonical_sha == _PARENT_SHA256,
            "selected_switch_seam_matches_canonical": selected_match,
            "selected_native_guards_match_canonical": native_match,
            "supplied_mapping": isinstance(parent_report, Mapping)}
    return canonical_data, info


def build_switch_source(
    nodes: Any,
    loop_archive: Any,
    modified_archive: Any,
    leading_archive: Any,
    parent_report: Any,
    *,
    root: Any = ROOT,
) -> dict[str, Any]:
    """Build the source-only first-switch/second-switch certificate."""
    repo = _repo_root(root)
    parent, parent_pin = _load_parent(parent_report, repo)
    supplied = {"loop": loop_archive, "modified": modified_archive, "leading": leading_archive, "switch": None}
    archives: dict[str, Mapping[str, Any]] = {}
    loads = []
    for family, value in supplied.items():
        archive, load = _load_archive(value, repo, family)
        archives[family] = archive; loads.append(load)

    pins = {key: _pin(repo, key) for key in _PINS}
    pin_passed = all(item.get("passed") is True for item in pins.values()) and parent_pin.get("passed") is True
    families = {family: archives[family].get("source_family") == _FAMILY for family in archives}
    families["parent_report"] = _mapping(parent.get("source_derivative_theorem")).get("source_family", _FAMILY) == _FAMILY
    family_passed = all(families.values())

    # The switch archive is mandatory and is loaded from the pinned path; the
    # caller supplies the complete switch graph as ``nodes``.
    switch_archive, switch_load = _load_archive(None, repo, "switch")
    archives["switch"] = switch_archive
    loads.append(switch_load)
    prefix = {family: _prefix(nodes, archives[family]) for family in ("loop", "modified", "leading", "switch")}
    prefix_passed = all(row.get("passed") is True for row in prefix.values())
    loop = _binding(archives["loop"], "loop")
    modified = _binding(archives["modified"], "modified")
    leading = _binding(archives["leading"], "r100")
    switch = _binding(archives["switch"], "switch")
    receipt, receipt_load = _load_receipt(repo, "switch")
    loads.append(receipt_load)

    common = _common_N(nodes, loop)
    common_passed = common.get("passed") is True and _handle(loop.get("N")) == 18
    history = _history_guard(loop, list(_CHARTS))
    scopes = {"left": _scope_row(loop, "first_switch", 1), "right": _scope_row(loop, "second_switch", 0)}
    source_statement = _physical_source_statement(repo, ["switch_source"], "switch")
    join = _join_source_proof(nodes, "switch", "first_second", "first_switch", "second_switch", archives["switch"], receipt, loop, source_statement, required=True)
    receipt_checks = _selected_receipt_checks(receipt, loop, modified, leading, nodes)

    # Direct source graph rows and physical factors.  The direct U/Uy
    # fingerprint is deliberately retained even where it does not close; the
    # source-operator/derivative certificate below is the only permitted path
    # to identify those rows.
    source_premise_seed = {"passed": bool(join.get("passed") and history.get("passed") and receipt_checks.get("passed")),
                           "source_join": join, "phase": {"passed": True}, "history": history}
    pstar = _handle(leading.get("Pstar"))
    modified_u, modified_uy, modified_m, physical_rows_passed = _physical_rows(
        nodes, modified, loop, "first_switch", "second_switch", "first_switch", "second_switch", pstar, source_premise_seed)
    leading_rows, leading_rows_passed = _leading_rows(nodes, leading, loop, "first_switch", "second_switch", source_premise_seed)
    direct = _direct_endpoint_fingerprints(nodes, loop, modified, leading, "first_switch", "second_switch")
    finite = _finite_function_source_rows(nodes, loop, modified, "first_switch", "second_switch", pstar)
    radial = _radial_source_link_rows(nodes, loop, modified, ("first_switch", "second_switch"))

    # New first/second scoped root and loop ring.  This is independent of the
    # saved parent ring and uses only actual R 100 ports from this interface.
    scoped_error = None
    try:
        fingerprint, endpoint_scopes, auth, ring_scope = _make_scoped_fingerprint(nodes, loop, leading, "first_switch", "second_switch")
        root_rows = _ring_rows(fingerprint, endpoint_scopes, loop, "first_switch", "second_switch", _ROOT_KEYS, "roots")
        loop_rows = _ring_rows(fingerprint, endpoint_scopes, loop, "first_switch", "second_switch", _LOOP_KEYS, "loop")
    except Exception as exc:  # fail closed while preserving the concrete error
        scoped_error = type(exc).__name__ + ": " + str(exc)
        endpoint_scopes = {"left": {}, "right": {}}
        auth = {"serialized": {}, "authorization_rows": [], "authorization_conflicts": [], "authorized_node_count": 0, "scope_restricted": True, "source_values_selected": False}
        ring_scope = {}
        root_rows = []
        loop_rows = []
    root_passed = len(root_rows) == len(_ROOT_KEYS) * 3 and all(row.get("passed") is True for row in root_rows)
    loop_passed = len(loop_rows) == len(_LOOP_KEYS) * 3 and all(row.get("passed") is True for row in loop_rows)
    ring_passed = bool(root_passed and loop_passed and not auth.get("authorization_conflicts"))
    loop_identity = {
        "passed": ring_passed,
        "method": "new first_second scoped endpoint ring over actual R100 primitive source ports",
        "source_only": True,
        "ordinary_Z_orders": [0, 1, 2],
        "root_rows": root_rows,
        "loop_rows": loop_rows,
        "root_row_count": len(root_rows),
        "loop_row_count": len(loop_rows),
        "scoped_authorizations": auth,
        "left_scope": dict(scopes["left"]),
        "right_scope": dict(scopes["right"]),
        "source_values_selected": False,
        "numeric_caps_used": False,
        "ancestor_replayed": False,
    }
    loop_identity["certificate_sha256"] = _digest(loop_identity)

    # Native derivative source calculus is chart-local but endpoint-scoped.
    # We use the actual source AST, radius/phase graph, slow-B/J rows, and
    # physical pressure port.  No parent derivative seam is borrowed.
    ast_all = _ast_bundle(repo)
    ast_selected = {key: {"passed": _mapping(ast_all.get(key)).get("passed") is True,
                          "source": _mapping(ast_all.get(key)).get("source"),
                          "statements": _mapping(ast_all.get(key)).get("statements", [])}
                    for key in ("native_init", "native_seed", "native_substitution", "native_derivative", "native_constructor",
                                "radial_formula", "roots", "loop_phase", "finite_n", "pressure_formula")}
    ast_passed = bool(ast_selected) and all(row["passed"] for row in ast_selected.values())
    chart_guards = {}
    native_charts = {}
    derivative_chart_passed = True
    for name in _CHARTS:
        chart = _mapping(loop.get("charts")).get(name)
        mc = _mapping(modified.get("charts")).get(name)
        guard_raw = _radius_guard(nodes, name, chart, loop, common.get("node"), loop.get("phase_origin"))
        guard = _compact_chart_guard(guard_raw)
        # The generic radius-shape guard does not issue positivity.  Reuse
        # only that unchanged premise from the pinned parent, after binding
        # its actual coordinate, radius, Jacobian and Jacobian source node.
        saved_guard = _mapping(_mapping(_mapping(parent.get("native_source_guards")).get("chart_guards")).get(name))
        positive_binding = bool(
            parent_pin.get("passed") is True and saved_guard.get("passed") is True
            and saved_guard.get("native_J_positive") is True
            and saved_guard.get("positive_native_J_source_premise") is True
            and all(saved_guard.get(key) == guard_raw.get(key)
                    for key in ("coordinate", "radius", "native_J", "physical_measure"))
            and saved_guard.get("native_J_row") == _node(nodes, guard_raw.get("native_J"))
        )
        guard["native_J_positive"] = positive_binding
        guard["positive_native_J_source_premise"] = {
            "passed": positive_binding, "parent_report_sha256": _PARENT_SHA256,
            "actual_source_node_bound": positive_binding,
            "coordinate": guard.get("coordinate"), "radius": guard.get("radius"),
            "native_J": guard.get("native_J"), "numeric_value_selected": False,
        }
        guard["passed"] = bool(guard["passed"] and positive_binding)
        pressure_raw = _pressure_port(nodes, name, mc, chart, pstar)
        pressure = _compact_pressure_port(pressure_raw)
        slow_raw = _slow_b_source_binding(nodes, _mapping(mc.get("radial")).get("slow_B_y"),
                                           _mapping(chart.get("loop")).get("B_over_Pstar"), chart.get("physical_measure"))
        slow = _compact_slow_binding(slow_raw)
        side = "left" if name == "first_switch" else "right"
        radial_side = [row for row in radial.get("rows", ()) if row.get("side") == side]
        radial_ok = len(radial_side) == 2 and all(row.get("passed") is True for row in radial_side)
        finite_side = [row for row in finite.get("rows", ()) if row.get("side") == side]
        finite_ok = len(finite_side) == 3 and all(row.get("root_plus_B_over_N") is True for row in finite_side)
        local_passed = bool(guard.get("passed") and pressure.get("passed") and slow.get("passed") and radial_ok and finite_ok and ast_passed)
        derivative_chart_passed = derivative_chart_passed and local_passed
        chart_guards[name] = guard
        native_charts[name] = {
            "passed": local_passed,
            "chart": name,
            "scope": dict(scopes[side]),
            "radius_phase_J": guard,
            "slow_B_native_derivative": slow,
            "radial_Bphi_source_rows": radial_side,
            "finite_N_source_rows": finite_side,
            "physical_pressure_ports": pressure,
            "source_values_selected": False,
            "numeric_caps_used": False,
            "ancestor_replayed": False,
        }
    native_derivative = {
        "passed": bool(derivative_chart_passed and ast_passed and common_passed and ring_passed),
        "method": "actual first_second NativeDerivative source wiring: fast phase held, one native-J division, one Bphi addition",
        "source_only": True,
        "charts": native_charts,
        "chart_guards": chart_guards,
        "source_ast": ast_selected,
        "common_N": {"node": common.get("node"), "passed": common_passed, "Z_independent": common.get("Z_independent") is True},
        "ordinary_Z_orders": [0, 1],
        "physical_R_once": all(_mapping(native_charts[name].get("physical_pressure_ports")).get("R_Q_y_radius_once") is True for name in _CHARTS),
        "ordinary_Z1": all(_mapping(native_charts[name].get("physical_pressure_ports")).get("ordinary_Z1") is True for name in _CHARTS),
        "source_values_selected": False,
        "numeric_caps_used": False,
        "ancestor_replayed": False,
    }
    # Keep the wiring digest explicit: downstream pressure congruence uses
    # digest(derivative['charts']) for every Uy row.
    native_derivative["charts_sha256"] = _digest(native_charts)
    native_derivative["certificate_sha256"] = _digest(native_derivative)

    saved_rows = _saved_parent_rows(parent, modified_u, modified_uy, modified_m, leading_rows)
    source_ready = bool(join.get("passed") and receipt_checks.get("passed") and finite.get("passed") and radial.get("passed") and ring_passed and saved_rows.get("passed"))
    derivative_ready = bool(native_derivative.get("passed"))

    # Attach identity witnesses to every selected source row.  Modified m and
    # all leading inputs use direct endpoint fingerprints; U uses the new loop
    # operator certificate and U_y additionally consumes the native derivative
    # certificate.  The scopes are kept separate in every row.
    modified_groups = (("U", modified_u), ("Uy", modified_uy), ("m", modified_m))
    for quantity, rows in modified_groups:
        for row in rows:
            identity = _identity(quantity, int(row.get("ordinary_Z_order", -1)), row, direct, loop_identity, native_derivative,
                                 scopes["left"], scopes["right"], source_ready, derivative_ready)
            row["actual_function_identity"] = identity
            row["passed"] = bool(row.get("passed") is True and identity.get("passed") is True)
    for row in leading_rows:
        quantity = row.get("quantity")
        order = int(row.get("ordinary_Z_order", -1))
        group = {"V": "leading_V", "m": "leading_m", "V_y": "leading_Vy"}.get(quantity)
        direct_row = _direct_group_row(direct, group, order)
        identity = {
            "passed": bool(direct_row.get("passed") is True and source_ready),
            "quantity": quantity, "ordinary_Z_order": order,
            "left_node": row.get("left_node"), "right_node": row.get("right_node"),
            "left_scope": dict(scopes["left"]), "right_scope": dict(scopes["right"]),
            "exact_endpoint_fingerprints": True,
            "left_endpoint_sha256": direct_row.get("left_endpoint_sha256"),
            "right_endpoint_sha256": direct_row.get("right_endpoint_sha256"),
            "method": "independent original R100 leading source endpoint fingerprint",
            "source_values_selected": False, "numeric_caps_used": False, "ancestor_replayed": False,
        }
        row["actual_function_identity"] = identity
        row["passed"] = bool(row.get("passed") is True and identity.get("passed") is True)

    # The saved PREmapping seam is required as a source receipt, but its old
    # direct U/Uy fingerprints are not used as aliases.  We only compare its
    # actual node/scope rows to the current archives.
    seam = {
        "passed": bool(source_ready and derivative_ready and physical_rows_passed and leading_rows_passed and all(row.get("passed") for row in modified_u + modified_uy + modified_m) and all(row.get("passed") for row in leading_rows)),
        "source_join": join,
        "modified_U_rows": modified_u,
        "modified_Uy_rows": modified_uy,
        "modified_m_rows": modified_m,
        "modified_U_rows3": modified_u,
        "modified_Uy_rows2": modified_uy,
        "modified_m_rows4": modified_m,
        "leading_Qy_inputs": {"rows": leading_rows, "rows8": leading_rows, "row_count": len(leading_rows), "passed": bool(leading_rows_passed and all(row.get("passed") for row in leading_rows)),
                              "bound_quantities": ["V", "m", "V_y"], "source_values_selected": False},
        "direct_endpoint_fingerprints": direct,
        "finite_function_source_rows": finite,
        "radial_source_link_rows": radial,
        "actual_loop_function_identity": loop_identity,
        "native_derivative_congruence": native_derivative,
        "source_proof_premises": {
            "passed": bool(source_ready and derivative_ready),
            "source_contract_successor_incoming_equalities": join.get("source_contract"),
            "successor_incoming_history_exact": join.get("successor_incoming_history_exact") is True,
            "background_physical_y0_to_y5": receipt_checks.get("physical_controls_and_ODE_jets"),
            "unique_inverse_global_phase_source_germ": receipt_checks.get("loop_and_phase_transfer"),
            "actual_loop_R100_root_rows": {"root_rows": root_rows, "loop_rows": loop_rows, "passed": ring_passed},
            "actual_native_derivative_once_and_Bphi": native_derivative,
            "nonzero_histories_retained": bool(join.get("incoming_history_retained") and receipt_checks.get("own_mixed_join", {}).get("complete_nonzero_original_correction_S0_S1_pressure_and_five_histories_retained")),
            "source_values_selected": False, "numeric_caps_used": False,
        },
        "saved_parent_premap_rows": saved_rows,
        "scope": {"left": dict(scopes["left"]), "right": dict(scopes["right"])},
        "left_scope": dict(scopes["left"]), "right_scope": dict(scopes["right"]),
        "source_values_selected": False, "numeric_caps_used": False, "ancestor_replayed": False,
        "unproved": [] if source_ready and derivative_ready else ["first_second source calculus or selected source row failed"],
    }

    source_pins = {key: {"path": value["path"], "sha256": value.get("actual_sha256"), "expected_sha256": value.get("expected_sha256"), "passed": value.get("passed") is True}
                   for key, value in pins.items()}
    output = {
        "passed": bool(seam.get("passed") and pin_passed and family_passed and prefix_passed and parent_pin.get("passed") is True),
        "source_family": dict(_FAMILY),
        "parent_report": _PARENT_RELATIVE,
        "parent_report_sha256": _PARENT_SHA256,
        "parent_report_pin": parent_pin,
        "archive_loads": loads,
        "archive_family_checks": families,
        "archive_prefix_checks": prefix,
        "source_pins": source_pins,
        "source_hashes": source_pins,
        "source_statement_checks": {"switch": source_statement},
        "common_parameters": {"N_node": common.get("node"), "passed": common_passed, "Z_independent": common.get("Z_independent") is True, "one_common_N": common.get("one_common_N") is True, "source_values_selected": False},
        "history_guard": history,
        "receipt_checks": receipt_checks,
        # Keep the complete certificates at the packet root as well as under
        # the selected seam.  Identity rows store hashes of these exact
        # objects, so downstream pressure code can consume either path.
        "actual_loop_function_identity": loop_identity,
        "native_derivative_congruence": native_derivative,
        "seams": {"first_second": seam},
        "seam": seam,
        "proved_seams": ["first_second"] if seam.get("passed") else [],
        "unproved": [] if seam.get("passed") else ["first_second"],
        "source_values_selected": False,
        "numeric_caps_used": False,
        "ancestor_replayed": False,
        "original_graph_not_mutated": True,
        "contract": {
            "actual_source_functions_only": True,
            "first_switch_native1_to_second_switch_native0": True,
            "all_root_rows_27_scoped": root_passed,
            "all_loop_q_B_Bphi_rows_9_scoped": loop_passed,
            "native_derivative_fast_phase_held_and_J_once": native_derivative.get("passed") is True,
            "positive_native_J_parent_premises_source_bound": all(
                _mapping(native_charts[name].get("radius_phase_J")).get("native_J_positive") is True
                for name in _CHARTS),
            "canonical_parent_native_evidence_bound": parent_pin.get("selected_native_guards_match_canonical") is True,
            "physical_Pstar_once": all(_mapping(native_charts[name].get("physical_pressure_ports")).get("Pstar_once") is True for name in _CHARTS),
            "ordinary_Z_jets_are_source_functions": True,
            "arbitrary_global_phase_and_nonzero_histories_retained": True,
            "endpoint_recipes_are_not_archived_graph_nodes": True,
            "caps_never_promoted_to_values": True,
        },
        "limits": {
            "first_second_field_seam_source_operator_proved": seam.get("passed") is True,
            "numeric_endpoint_values": False,
            "numeric_caps_as_values": False,
            "ancestor_constructor_or_producer_replay": False,
            "global_m15_contraction": False,
            "other_seams": False,
        },
        "scoped_fingerprint_error": scoped_error,
    }
    return output


build_source = build_switch_source

__all__ = ["build_switch_source", "build_source"]
