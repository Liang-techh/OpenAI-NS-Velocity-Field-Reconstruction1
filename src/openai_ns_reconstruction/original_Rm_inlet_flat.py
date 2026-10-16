"""Archive-only x=1 guards for the prescribed Rm finite-N inlet.

The Rm patch contains two operations which the older endpoint fingerprint did
not admit: the exact compact bump and the original unique implicit leading
control.  This extension admits only their source-bound endpoint semantics.
Compact bumps are zero when their exact endpoint argument is outside the
support, and empty ``dt/t`` histories reduce to their incoming memories.  The
implicit controls remain opaque functions of Z and their defining source system
with fixed leading N0. The common correction N is carried as context; it is not
substituted for that leading input or used as a control value.

This module does not construct a producer, replay an accepted calculation, or
claim the postrestore modified ``Q_y``/``R Q_y`` seam.  It only returns compact
source evidence for the two x=1 scopes.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path
from typing import Any, Mapping

from .original_active_viscous_m15 import ActiveEndpointFingerprint, ROOT


_RM_ARCHIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_Rm_patch_finite_N_C5_functions.json.gz"
)
_RM_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_Rm_patch_finite_N_C5_functions.py"
)
_RM_C6_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_Rm_patch_C6_functions.py"
)
_BETA_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_five_moment_repair.py"
)
_RM_ARCHIVE_SHA256 = "8f3bc343e9dbffbdda465ce94d6183735184645a70e291cb900d3cde726769ad"
_RM_SOURCE_SHA256 = "c4445ef8acfb33d29cf9022862f8b6396e6aa2f38bda7467409b5378a5deaf15"
_RM_C6_SOURCE_SHA256 = "a6f1cb5deaf17e9482aa964aef7867aadbc815e4d703a415e4945cf665f14002"
_BETA_SOURCE_SHA256 = "978de37342cc2a63aa7c1bfda573e512f344e0157aed28b6017c7dae2ff8eba5"


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _unwrap(value: Any, key: str) -> Mapping[str, Any]:
    value = _mapping(value)
    nested = value.get(key)
    return nested if isinstance(nested, Mapping) else value


def _node(nodes: Any, index: Any) -> Mapping[str, Any] | None:
    if isinstance(index, bool) or not isinstance(index, int):
        return None
    try:
        row = nodes[index] if not isinstance(nodes, Mapping) else nodes.get(index, nodes.get(str(index)))
    except (IndexError, KeyError, TypeError):
        return None
    return row if isinstance(row, Mapping) else None


def _repo_root(root: Any) -> Path:
    path = Path(root) if root is not None else ROOT
    if path.is_file():
        path = path.parent
    if (path / "src").is_dir() and (path / "experiments").is_dir():
        return path
    return ROOT


def _pin(root: Path, relative: str, expected: str) -> dict[str, Any]:
    path = root / relative
    actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
    return {
        "path": relative,
        "exists": path.is_file(),
        "sha256": actual,
        "expected_sha256": expected,
        "passed": actual == expected,
    }


def _children(row: Mapping[str, Any]) -> list[int]:
    """Return graph references from a serialized operation row only."""
    refs: list[int] = []
    for key in (
        "arguments", "argument", "numerator", "denominator", "expression", "value",
        "phase_function", "target_phase", "angle_lower", "angle_upper", "flat_predicate",
        "active_body", "flat_value", "system", "body", "lower", "upper", "integrand",
        "outer", "dx", "direct", "chain",
    ):
        value = row.get(key)
        if isinstance(value, list):
            refs.extend(x for x in value if isinstance(x, int) and not isinstance(x, bool))
        elif isinstance(value, int) and not isinstance(value, bool):
            refs.append(value)
    return refs


def _reachable(nodes: Any, roots: list[int]) -> set[int]:
    seen: set[int] = set()
    stack = list(roots)
    while stack:
        index = stack.pop()
        if not isinstance(index, int) or isinstance(index, bool) or index in seen:
            continue
        row = _node(nodes, index)
        if row is None:
            continue
        seen.add(index)
        stack.extend(_children(row))
    return seen


def _depends_on(
    nodes: Any,
    root: Any,
    target: int,
    memo: dict[tuple[int, int], bool] | None = None,
    active: set[tuple[int, int]] | None = None,
) -> bool:
    """Structural dependency on the selected native coordinate only."""
    if not isinstance(root, int) or isinstance(root, bool):
        return False
    if root == target:
        return True
    memo = memo if memo is not None else {}
    active = active if active is not None else set()
    key = (root, target)
    if key in memo:
        return memo[key]
    # Serialized source graphs may retain a non-semantic self-reference in a
    # metadata-like child slot.  It cannot establish coordinate dependence.
    if key in active:
        return False
    active.add(key)
    row = _node(nodes, root)
    if row is None:
        active.remove(key)
        memo[key] = False
        return False
    answer = any(_depends_on(nodes, child, target, memo, active) for child in _children(row))
    active.remove(key)
    memo[key] = answer
    return answer


def _ints(value: Any) -> list[int]:
    if isinstance(value, bool):
        return []
    if isinstance(value, int):
        return [value]
    if isinstance(value, list):
        out: list[int] = []
        for item in value:
            out.extend(_ints(item))
        return out
    if isinstance(value, Mapping):
        out = []
        for item in value.values():
            out.extend(_ints(item))
        return out
    return []


def _rational_json(value: Any) -> list[int] | None:
    """Serialize EndpointFingerprint's exact Fraction without losing type."""
    if value is None:
        return None
    numerator = getattr(value, "numerator", None)
    denominator = getattr(value, "denominator", None)
    if numerator is None or denominator is None:
        raise TypeError(f"Expected an exact rational, got {type(value)!r}")
    return [int(numerator), int(denominator)]


class RmInletFingerprint(ActiveEndpointFingerprint):
    """Active endpoint fingerprint with scoped Rm compact-beta semantics."""

    def __init__(self, nodes: Any, bindings: Mapping[str, Any]):
        super().__init__(nodes, bindings)
        self.flat_beta_nodes: set[int] = set()
        self.flat_beta_arguments: dict[int, Any] = {}
        self.flat_derivative_nodes: set[int] = set()
        self.unique_control_rows: list[dict[str, Any]] = []
        self.opaque_unique_control_digests: dict[int, str] = {}

    def _opaque(self, tag: str, node: int, row: Mapping[str, Any], scope: Mapping[int, Any]) -> tuple[str, None]:
        # Do not include a numeric field value.  The source-system identity,
        # component, Z handle, common N and endpoint scope are the identity of
        # this opaque function; different systems cannot alias accidentally.
        scope_keys = sorted(int(key) for key in scope)
        payload = {
            "tag": tag,
            "source_system": row.get("system"),
            "component": row.get("component"),
            "Z_variable": row.get("Z_variable"),
            "N_handle": self.b.get("N"),
            "fixed_leading_input_N0": (_node(self.nodes,row.get("system")) or {}).get("fixed_leading_input_N0"),
            "source_family": row.get("source_family"),
            "scope_coordinates": scope_keys,
        }
        return self.digest([tag, payload]), None

    def at(self, node: int, scope: Mapping[int, Any] | None = None):
        scope = scope if scope is not None else {self.b["coordinate"]: self.one}
        row = _node(self.nodes, node)
        if row is None:
            raise ValueError(f"Unknown graph node {node}")
        operation = row.get("operation")

        if operation in {
            "compact_raw_beta",
            "compact_raw_beta_derivative",
            "compact_raw_beta_second_derivative",
        }:
            key = (node, tuple(sorted((k, v[0]) for k, v in scope.items())))
            if key in self.cache:
                return self.cache[key]
            argument = row.get("argument")
            value = self.at(argument, scope)
            rational = value[1]
            if rational is not None and (rational <= -1 or rational >= 1):
                out = self.zero
                self.flat_beta_nodes.add(node)
                self.flat_beta_arguments[node] = rational
                if operation != "compact_raw_beta":
                    self.flat_derivative_nodes.add(node)
            else:
                # The active beta itself is left opaque.  In particular, this
                # never turns a nonzero interior beta or its Z derivative into
                # zero merely because it is inconvenient to fingerprint.
                out = self.digest([
                    operation,
                    row.get("exact_definition"),
                    row.get("original_source"),
                    value[0],
                ]), None
            self.cache[key] = out
            return out

        if operation == "original_Rm_unique_leading_implicit_control":
            key = (node, tuple(sorted((k, v[0]) for k, v in scope.items())))
            if key in self.cache:
                return self.cache[key]
            out = self._opaque("original_Rm_unique_leading_implicit_control_ZN", node, row, scope)
            self.cache[key] = out
            self.unique_control_rows.append({
                "node": node,
                "system": row.get("system"),
                "component": row.get("component"),
                "Z_variable": row.get("Z_variable"),
                "N_handle": self.b.get("N"),
                "source_family": row.get("source_family"),
                "exact_same_unique_original_common_box_solution": row.get(
                    "exact_value_of_same_unique_original_common_box_solution"
                ) is True,
                "rational_value": out[1],
                "opaque_digest": out[0],
                "nonzero_Z_function_preserved": out[1] is None,
            })
            self.opaque_unique_control_digests[node] = out[0]
            return out

        if operation == "original_Rm_unique_leading_implicit_system":
            key = (node, tuple(sorted((k, v[0]) for k, v in scope.items())))
            if key not in self.cache:
                self.cache[key] = self._opaque("original_Rm_unique_leading_implicit_system", node, row, scope)
            return self.cache[key]

        return super().at(node, scope)


def _binding_rows(rm_bindings: Any, loop_bindings: Any) -> tuple[Mapping[str, Any], Mapping[str, Any]]:
    rm = _unwrap(rm_bindings, "Rm_patch_finite_N_C5_bindings")
    loop = _unwrap(loop_bindings, "prescribed_loop_C5_bindings")
    if not isinstance(rm.get("chart"), Mapping):
        raise ValueError("Rm patch chart binding required")
    if not isinstance(loop.get("charts"), Mapping):
        raise ValueError("prescribed loop chart bindings required")
    return rm, loop


def _endpoint_scope(fingerprint: RmInletFingerprint, rm: Mapping[str, Any], loop: Mapping[str, Any]):
    rm_coordinate = rm["chart"]["coordinate"]
    post_coordinate = loop["charts"]["postrestore"]["coordinate"]
    one = fingerprint.rational(1)
    return {
        "rm_patch_x1": rm_coordinate,
        "postrestore_native1": post_coordinate,
        "scope": {rm_coordinate: one, post_coordinate: one},
    }


def _compare_nodes(fingerprint: RmInletFingerprint, left: int, right: int, scope: Mapping[int, Any], quantity: str):
    l_digest, l_rational = fingerprint.at(left, scope)
    r_digest, r_rational = fingerprint.at(right, scope)
    return {
        "quantity": quantity,
        "left_node": left,
        "right_node": right,
        "passed": l_digest == r_digest and l_rational == r_rational,
        "canonical_endpoint_sha256": l_digest if l_digest == r_digest else None,
        "left_rational": _rational_json(l_rational),
        "right_rational": _rational_json(r_rational),
        "left_endpoint_sha256": l_digest,
        "right_endpoint_sha256": r_digest,
    }


def _compare_modified_fields(
    fingerprint: RmInletFingerprint,
    post: Mapping[str, Any],
    chart: Mapping[str, Any],
    scope: Mapping[int, Any],
) -> dict[str, Any]:
    """Compare the actual finite-N modified-field ports in both charts.

    The Rm archive exposes no postrestore radial ``E_y``/``V_y`` ports, and
    this helper deliberately reports that gap instead of identifying a field
    with a derivative or manufacturing one from a numeric value.
    """
    names = ("E", "V", "m", "h", "k", "e", "p", "P0", "pressure")
    post_fields = _mapping(post.get("modified_fields"))
    rm_fields = _mapping(chart.get("modified_fields"))
    rows: list[dict[str, Any]] = []
    missing: list[dict[str, Any]] = []
    for name in names:
        left = post_fields.get(name)
        right = rm_fields.get(name)
        if not isinstance(left, list) or not isinstance(right, list) or len(left) < 3 or len(right) < 3:
            missing.append({
                "quantity": name,
                "postrestore_orders": len(left) if isinstance(left, list) else None,
                "Rm_patch_orders": len(right) if isinstance(right, list) else None,
            })
            continue
        for order in range(3):
            rows.append(_compare_nodes(fingerprint, left[order], right[order], scope, f"modified_field_{name}_Z{order}"))
    radial_gap = {
        "quantity": "leading_radial_E_y_V_y",
        "postrestore_ports": sorted(k for k in post if k in ("radial", "radial_background")),
        "Rm_patch_ports": sorted(k for k in chart if k in ("radial", "radial_background")),
        "passed": False,
        "status": "unproved_missing_postrestore_radial_ports",
    }
    return {
        "passed": bool(rows) and not missing and all(row["passed"] for row in rows),
        "ordinary_Z_orders": [0, 1, 2],
        "rows": rows,
        "missing": missing,
        "radial_derivative_gap": radial_gap,
        "source_only": True,
        "numeric_values_used": False,
        "actual_source_group": "modified_fields",
        "not_a_leading_only_comparison": True,
    }


def build_Rm_inlet_flat_certificate(
    nodes: Any,
    rm_bindings: Mapping[str, Any],
    loop_bindings: Mapping[str, Any],
    root: Any = ROOT,
) -> dict[str, Any]:
    """Build strict archive-only evidence at postrestore native 1/Rm x=1."""
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 20000))
    rm, loop = _binding_rows(rm_bindings, loop_bindings)
    repo = _repo_root(root)
    pins = {
        "rm_archive": _pin(repo, _RM_ARCHIVE, _RM_ARCHIVE_SHA256),
        "rm_source": _pin(repo, _RM_SOURCE, _RM_SOURCE_SHA256),
        "rm_c6_source": _pin(repo, _RM_C6_SOURCE, _RM_C6_SOURCE_SHA256),
        "beta_source": _pin(repo, _BETA_SOURCE, _BETA_SOURCE_SHA256),
    }

    chart = rm["chart"]
    post = loop["charts"]["postrestore"]
    bindings = {"coordinate": chart["coordinate"], "loop": loop["charts"], "N": rm["N"]}
    fingerprint = RmInletFingerprint(nodes, bindings)
    endpoints = _endpoint_scope(fingerprint, rm, loop)
    scope = endpoints["scope"]

    geometry_rows = []
    for key in ("radius", "physical_measure", "global_offset", "actual_phase"):
        geometry_rows.append(_compare_nodes(fingerprint, post[key], chart[key], scope, key))
    geometry_rows.append({
        "quantity": "postrestore_coordinate_and_Rm_patch_coordinate_are_both_native_one",
        "postrestore_coordinate": post["coordinate"],
        "Rm_patch_coordinate": chart["coordinate"],
        "passed": scope[post["coordinate"]][1] == 1 and scope[chart["coordinate"]][1] == 1,
    })

    roots = _ints(rm.get("physical")) + _ints(rm.get("boundary")) + _ints(rm.get("radial"))
    reachable = _reachable(nodes, roots)
    compact_nodes = sorted(i for i in reachable if _node(nodes, i).get("operation") == "compact_raw_beta")
    all_lazy_nodes = sorted(i for i in reachable if _node(nodes, i).get("operation") == "original_lazy_flat_branch")
    # The Rm graph also contains unrelated lazy sigma/collar branches.  Only
    # beta derivative branches carrying the source's exact-beta flag belong to
    # this endpoint certificate.
    lazy_nodes = sorted(
        i for i in all_lazy_nodes
        if _node(nodes, i).get("original_beta_exact_flat_derivative") is True
    )
    coordinate = chart["coordinate"]
    dependency_memo: dict[tuple[int, int], bool] = {}
    control_nodes = sorted(i for i in reachable if _node(nodes, i).get("operation") == "original_Rm_unique_leading_implicit_control")
    compact_rows = []
    for index in compact_nodes:
        digest, rational = fingerprint.at(index, scope)
        row = _node(nodes, index)
        compact_rows.append({
            "node": index,
            "argument": row.get("argument"),
            "endpoint_dependent": _depends_on(nodes, row.get("argument"), coordinate, dependency_memo),
            "passed_source_definition": row.get("exact_definition") == "exp(-1/(1-x^2)) for |x|<1; exact0 elsewhere",
            "passed_source_module": _mapping(row.get("original_source")).get("module") == "lei_ren_part1_paper_compliant_five_moment_repair.py",
            "endpoint_rational": _rational_json(rational),
            "endpoint_zero": rational == 0,
            "canonical_endpoint_sha256": digest,
        })
    lazy_rows = []
    for index in lazy_nodes:
        digest, rational = fingerprint.at(index, scope)
        row = _node(nodes, index)
        lazy_rows.append({
            "node": index,
            "Delta": row.get("Delta"),
            "endpoint_dependent": _depends_on(nodes, index, coordinate, dependency_memo),
            "flat_predicate": row.get("flat_predicate"),
            "passed_original_beta_exact_flat_derivative": row.get("original_beta_exact_flat_derivative") is True,
            "endpoint_zero": rational == 0,
            "canonical_endpoint_sha256": digest,
        })

    control_rows = []
    for index in control_nodes:
        digest, rational = fingerprint.at(index, scope)
        row = _node(nodes, index)
        control_rows.append({
            "node": index,
            "system": row.get("system"),
            "component": row.get("component"),
            "Z_variable": row.get("Z_variable"),
            "N_handle": rm.get("N"),
            "source_family": row.get("source_family"),
            "exact_same_unique_original_common_box_solution": row.get(
                "exact_value_of_same_unique_original_common_box_solution"
            ) is True,
            "opaque_digest": digest,
            "rational_value": _rational_json(rational),
            "nonzero_Z_function_preserved": rational is None,
        })

    history_rows = []
    history = chart.get("history", {})
    incoming = rm.get("incoming", {})
    for name in ("m", "h", "k", "e", "p"):
        for order, (history_node, incoming_node) in enumerate(zip(history.get(name, ()), incoming.get(name, ()) )):
            row = _compare_nodes(fingerprint, history_node, incoming_node, scope, f"incoming_history_{name}_Z{order}")
            integral_id = rm.get("chart", {}).get("contributions", {}).get(name, [None] * 6)[order]
            integral = _node(nodes, integral_id)
            # The serialized endpoints are the lower bound ``1`` and the
            # native-coordinate node.  Endpoint emptiness is proved by the
            # scoped fingerprint after that coordinate is set to one, not by
            # comparing their raw graph IDs.
            row["empty_integral"] = integral is not None and integral_id in fingerprint.empty
            history_rows.append(row)

    integral_nodes = sorted(fingerprint.empty)
    definition = _mapping(rm.get("definition"))
    required_definition = (
        "actual_same_N_Section11_loop_inverse_density_and_history_C5",
        "nonzero_true_Rm_correction_and_signed_pressure_inlets_retained",
        "exact_log_coordinate_transport_and_rate_zero_pressure_memory",
        "original_global_phase_Ra_hbB_sc_half_not_new_patch_phase_origin",
        "original_collar_not_applied_to_Rm_patch",
        "full_fast_and_slow_total_radial_C5_and_modified_generic_C4",
        "physical_U_Uy_mean_C5_Q_Qy_C4_and_signed_rminus_to_Rh_pressure_C3",
        "exact_physical_Rm_and_Rm_squared_x_pressure_measures_once",
    )
    definition_rows = {key: definition.get(key) is True for key in required_definition}

    geometry_passed = all(row.get("passed") for row in geometry_rows)
    endpoint_compact_rows = [row for row in compact_rows if row["endpoint_dependent"]]
    endpoint_lazy_rows = [row for row in lazy_rows if row["endpoint_dependent"]]
    compact_passed = bool(endpoint_compact_rows) and all(
        row["passed_source_definition"] and row["passed_source_module"] and row["endpoint_zero"]
        for row in endpoint_compact_rows
    )
    lazy_passed = bool(endpoint_lazy_rows) and all(
        row["passed_original_beta_exact_flat_derivative"] and row["endpoint_zero"]
        for row in endpoint_lazy_rows
    )
    controls_passed = bool(control_nodes) and all(
        row["exact_same_unique_original_common_box_solution"] and row["nonzero_Z_function_preserved"]
        for row in control_rows
    )
    history_passed = bool(history_rows) and all(row["passed"] and row["empty_integral"] for row in history_rows)
    modified_probe = _compare_modified_fields(fingerprint, post, chart, scope)
    source_passed = all(row["passed"] for row in pins.values()) and all(definition_rows.values())

    return {
        "passed": bool(source_passed and geometry_passed and compact_passed and lazy_passed and controls_passed and history_passed),
        "method": "pinned_Rm_patch_x1_compact_beta_and_memory_endpoint_fingerprint",
        "source_only": True,
        "numeric_caps_used": False,
        "global_opaque_alias_used": False,
        "source_pins": pins,
        "source_definition_flags": definition_rows,
        "scope": {
            "postrestore_native_endpoint": 1,
            "Rm_patch_native_x": 1,
            "postrestore_coordinate": post["coordinate"],
            "Rm_patch_coordinate": chart["coordinate"],
            "common_N_handle": rm.get("N"),
            "phase_origin_handle": rm.get("phase_origin"),
        },
        "geometry": {
            "passed": geometry_passed,
            "rows": geometry_rows,
            "no_phase_reset_assumed": True,
            "physical_measure_compared_once": True,
        },
        "raw_modified_field_probe": modified_probe,
        "actual_source_ports": {
            "postrestore_loop": {
                key: post.get(key)
                for key in ("coordinate", "radius", "physical_measure", "global_offset", "actual_phase", "roots", "loop", "history", "incoming")
            },
            "Rm_patch": {
                "chart": {key: chart.get(key) for key in ("coordinate", "radius", "physical_measure", "global_offset", "actual_phase")},
                "physical": {key: rm.get("physical", {}).get(key) for key in ("U", "U_y", "m", "m_y", "Q", "Q_y")},
                "boundary_R_Q_y": rm.get("boundary", {}).get("R_Q_y"),
                "radial": rm.get("radial"),
            },
        },
        "compact_raw_beta": {
            "passed": compact_passed and lazy_passed,
            "reachable_compact_beta_nodes": len(compact_nodes),
            "endpoint_zero_nodes": sum(row["endpoint_zero"] for row in endpoint_compact_rows),
            "endpoint_dependent_nodes": len(endpoint_compact_rows),
            "non_endpoint_normalization_nodes": len(compact_rows) - len(endpoint_compact_rows),
            "rows": compact_rows,
            "flat_derivative_rows": lazy_rows,
            "all_lazy_flat_branch_nodes": len(all_lazy_nodes),
            "endpoint_dependent_flat_derivative_nodes": len(endpoint_lazy_rows),
            "source_theorem": "exact compact support |argument|>=1 gives beta=0; beta first/second derivative lazy branches carry original_beta_exact_flat_derivative=True",
            "products_with_actual_beta_vanish_at_x1": compact_passed and lazy_passed,
            "nonzero_interior_beta_not_replaced": True,
        },
        "unique_implicit_control": {
            "passed": controls_passed,
            "rows": control_rows,
            "reachable_control_nodes": len(control_nodes),
            "opaque_ZN_functions_preserved": controls_passed,
            "control_values_not_set_to_zero": True,
            "leading_control_definition": "h(Z; fixed leading N0), in common finite-N correction context",
            "candidate_N_not_substituted_for_leading_N0": True,
            "source_systems_are_scoped": sorted({row["system"] for row in control_rows}),
        },
        "incoming_history": {
            "passed": history_passed,
            "rows": history_rows,
            "empty_integral_nodes": integral_nodes,
            "empty_integral_count": len(integral_nodes),
            "incoming_memory_retained": history_passed,
            "rate_zero_pressure_memory_retained": True,
            "physical_measure": "dt/t",
            "no_history_reset": True,
        },
        "limits": {
            "full_modified_C0_Z1_postrestore_Rm_Qy_identity": False,
            "postrestore_modified_physical_ports_supplied_to_this_api": False,
            "postrestore_radial_Ey_Vy_ports_supplied_to_this_api": False,
            "raw_modified_fields_E_V_m_h_k_e_p_P0_pressure_C0_C1_C2": modified_probe["passed"],
            "native_derivative_theorem": False,
            "numeric_Rm_or_Rp_value": False,
            "cap_as_function_value": False,
            "global_opaque_control_alias": False,
        },
    }


__all__ = ["RmInletFingerprint", "build_Rm_inlet_flat_certificate"]
