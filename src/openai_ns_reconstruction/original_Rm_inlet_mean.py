"""Source-bound x=1 leading-mean identity for the Rm inlet.

The finite-N Rm patch stores the leading mean in an exact ordinary-
``Z`` graph.  At its x=1 endpoint every partial bump weight is an empty
integral.  After those terms are removed, the remaining numerator is the
source leading mean multiplied by the same positive ``Pstar`` as the
normalizing denominator.  This module proves that identity in the exact
source polynomial ring while treating the post-restore mean ports as
independent opaque functions of ``Z``.

No source value, interval cap, or ancestor producer is selected here.  The
certificate is endpoint scoped and does not alias the modified mean or any
other field.
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path
from typing import Any, Mapping

from .original_core_numeric import ROOT
from .original_Rm_inlet_flat import RmInletFingerprint, _node, _reachable, _unwrap
from .original_reference_viscous_m15 import PhysicalEndpointFingerprint


_LEADING_ARCHIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_Rm_patch_C6_functions.json.gz"
)
_OUTER_ARCHIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.json.gz"
)
_LOOP_ARCHIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz"
)
_MODIFIED_ARCHIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz"
)
_RM_ARCHIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_Rm_patch_finite_N_C5_functions.json.gz"
)
_ARCHIVE_SHA256 = {
    _LEADING_ARCHIVE: "4352ebd8862b5c92e3e4026728183014063baf556802a88421cdbd517ff1346a",
    _OUTER_ARCHIVE: "5e5b565dd8d2abc41e6f3df0d9e17000e524d493e58a5778dfcf217c6e42f44a",
    _LOOP_ARCHIVE: "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa",
    _MODIFIED_ARCHIVE: "b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a",
    _RM_ARCHIVE: "8f3bc343e9dbffbdda465ce94d6183735184645a70e291cb900d3cde726769ad",
}


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


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


def _flatten_ints(value: Any) -> list[int]:
    if isinstance(value, bool):
        return []
    if isinstance(value, int):
        return [value]
    if isinstance(value, list):
        out: list[int] = []
        for item in value:
            out.extend(_flatten_ints(item))
        return out
    return []


def _row_op(nodes: Any, node: Any) -> str | None:
    row = _node(nodes, node)
    return row.get("operation") if row is not None else None


def _scoped_key(scope: Mapping[int, Any]) -> tuple[tuple[int, Any], ...]:
    return tuple(sorted((int(k), v[0]) for k, v in scope.items()))


class _MeanFingerprint(PhysicalEndpointFingerprint, RmInletFingerprint):
    """Exact ring with source post-restore means as opaque atoms.

    ``PhysicalEndpointFingerprint`` evaluates products eagerly.  Partial
    weights are exact zero factors, so checking their direct integral factors
    before the other product factors prevents expansion of unrelated implicit
    controls.  The short circuit is valid only for serialized definite
    integrals, whose endpoints are separately evaluated at native one.
    """

    def __init__(self, nodes: Any, bindings: Mapping[str, Any], opaque_means: list[int]):
        super().__init__(nodes, bindings)
        self.opaque_means = set(opaque_means)
        self.opaque_mean_digests: dict[int, str] = {}

    def _opaque_mean(self, node: int, scope: Mapping[int, Any]):
        key = (node, _scoped_key(scope))
        if key in self.cache:
            return self.cache[key]
        row = _node(self.nodes, node) or {}
        out = self.digest(
            [
                "original_Rm_postrestore_mean_opaque_Z_function",
                node,
                row.get("operation"),
                row.get("original_source"),
                _scoped_key(scope),
                self.b.get("N"),
            ]
        ), None
        self.cache[key] = out
        self.opaque_mean_digests[node] = out[0]
        return out

    def at(self, node: int, scope: Mapping[int, Any] | None = None):
        scope = scope if scope is not None else {self.b["coordinate"]: self.one}
        if node in self.opaque_means:
            if scope.get(self.b['postrestore_coordinate'])!=self.one:
                raise ValueError('Opaque inlet mean leaf is restricted to actual postrestore native1')
            return self._opaque_mean(node, scope)

        row = _node(self.nodes, node)
        if row is not None and row.get("operation") == "product":
            key = (node, _scoped_key(scope))
            if key in self.cache:
                return self.cache[key]
            # Check only cheap direct integral/negative-integral factors.  Do
            # not inspect or assign any value to the other source factors.
            for child in row.get("arguments", ()):
                child_row = _node(self.nodes, child)
                if child_row is None:
                    continue
                candidate = child
                if child_row.get("operation") == "negative":
                    candidate = child_row.get("argument")
                    child_row = _node(self.nodes, candidate)
                if child_row is None or child_row.get("operation") != "definite_integral":
                    continue
                value = self.at(candidate, scope)
                if value[1] == 0:
                    self.cache[key] = self.zero
                    return self.zero
        return super().at(node, scope)


def _binding(value: Any, key: str) -> Mapping[str, Any]:
    return _unwrap(value, key)


def _endpoint_scope(fp: _MeanFingerprint, leading: Mapping[str, Any], loop: Mapping[str, Any], rm: Mapping[str, Any]):
    chart = _mapping(leading.get("chart"))
    post = _mapping(loop.get("charts", {}).get("postrestore"))
    chart_coordinate = chart.get("phase")
    post_coordinate = post.get("coordinate")
    if not isinstance(chart_coordinate, int) or not isinstance(post_coordinate, int):
        raise ValueError("Leading/Rm endpoint coordinates are required")
    one = fp.one
    return {
        chart_coordinate: one,
        post_coordinate: one,
    }, chart_coordinate, post_coordinate


def _wrapper_rows(nodes: Any, leading: Mapping[str, Any], outer: Mapping[str, Any], fp: _MeanFingerprint, scope):
    inlet = _mapping(leading.get("inlet"))
    post = _mapping(outer.get("charts", {}).get("post_restore"))
    inlet_fields = _mapping(inlet.get("fields"))
    outer_fields = _mapping(post.get("fields"))
    inlet_means = inlet_fields.get("m")
    outer_means = outer_fields.get("m")
    post_phase = post.get("phase")
    rows: list[dict[str, Any]] = []
    if not isinstance(inlet_means, list) or not isinstance(outer_means, list) or len(inlet_means) < 3 or len(outer_means) < 3:
        return rows, False
    for order in range(3):
        wrapper_node = inlet_means[order]
        source_node = outer_means[order]
        row = _node(nodes, wrapper_node) or {}
        structure_passed = (
            row.get("operation") == "function_substitution"
            and row.get("expression") == source_node
            and row.get("variable") == post_phase
            and row.get("value") == 1
            and row.get("Z_independent_substitution") is True
        )
        wrapper_value = fp.at(wrapper_node, scope)
        source_value = fp.at(source_node, scope)
        rows.append(
            {
                "ordinary_Z_order": order,
                "wrapper_node": wrapper_node,
                "source_expression_node": source_node,
                "variable_node": row.get("variable"),
                "substitution_value": row.get("value"),
                "passed_structure": structure_passed,
                "passed_endpoint_identity": wrapper_value[0] == source_value[0] and wrapper_value[1] == source_value[1],
                "source_values_selected": False,
                "numeric_value": None,
            }
        )
    return rows, bool(rows) and all(row["passed_structure"] and row["passed_endpoint_identity"] for row in rows)


def _partial_rows(nodes: Any, leading: Mapping[str, Any], fp: _MeanFingerprint, scope):
    partial = _mapping(leading.get("partial_weights"))
    roots: list[tuple[str, str, int]] = []

    def collect(value: Any, path: str):
        if isinstance(value, int) and not isinstance(value, bool):
            if value:
                roots.append((path, path.rsplit("/", 1)[-1], value))
        elif isinstance(value, list):
            for index, item in enumerate(value):
                collect(item, f"{path}/{index}")

    for name, value in partial.items():
        collect(value, name)
    rows: list[dict[str, Any]] = []
    integral_ids: set[int] = set()

    def weight_integrals(root: int) -> set[int]:
        """Find serialized weight integrals without entering an integrand.

        The integrand itself contains the source loop's many unrelated
        history integrals.  They are not partial bump weights and must not be
        mistaken for endpoint tails in this certificate.
        """
        found: set[int] = set()
        seen: set[int] = set()
        stack = [root]
        while stack:
            node = stack.pop()
            if not isinstance(node, int) or node in seen:
                continue
            seen.add(node)
            row = _node(nodes, node)
            if row is None:
                continue
            if row.get("operation") == "definite_integral":
                found.add(node)
                continue
            for key in ("arguments", "argument", "numerator", "denominator"):
                value = row.get(key)
                if isinstance(value, int):
                    stack.append(value)
                elif isinstance(value, list):
                    stack.extend(item for item in value if isinstance(item, int) and not isinstance(item, bool))
        return found

    for path, _, node in roots:
        integral_ids.update(weight_integrals(node))
        value = fp.at(node, scope)
        rows.append(
            {
                "weight_path": path,
                "weight_node": node,
                "operation": _row_op(nodes, node),
                "passed_exact_zero": value[1] == 0,
                "endpoint_rational": [int(value[1].numerator), int(value[1].denominator)] if value[1] is not None else None,
                "source_values_selected": False,
            }
        )
    integral_rows: list[dict[str, Any]] = []
    for node in sorted(integral_ids):
        row = _node(nodes, node) or {}
        flags = {
            "exact_function_integral": row.get("exact_function_integral") is True,
            "numerical_value_not_installed": row.get("numerical_value_not_installed") is True,
            "exact_original_Rm_bump_weight": row.get("exact_original_Rm_bump_weight") is True,
            "full_original_support": row.get("full_original_support") is False,
            "ordinary_measure_dx": row.get("ordinary_measure") == "dx",
            "Z_independent_weight": row.get("Z_independent_weight") is True,
            "range_or_midpoint_not_defining_integral": row.get("range_or_midpoint_not_defining_integral") is True,
        }
        integral_rows.append(
            {
                "integral_node": node,
                "lower_node": row.get("lower"),
                "upper_node": row.get("upper"),
                "empty_at_endpoint": node in fp.empty,
                "source_definition_flags": flags,
                "passed": node in fp.empty and all(flags.values()),
            }
        )
    passed = bool(rows) and all(row["passed_exact_zero"] for row in rows) and bool(integral_rows) and all(
        row["passed"] for row in integral_rows
    )
    return rows, integral_rows, passed


def _control_rows(nodes: Any, leading: Mapping[str, Any], fp: _MeanFingerprint, scope):
    chart = _mapping(leading.get("chart"))
    roots = _flatten_ints(chart.get("fields", {}).get("m", [])) + _flatten_ints(leading.get("partial_weights"))
    reachable = _reachable(nodes, roots)
    # The system node is the defining implicit map, not a selected control
    # value.  Only the five actual control leaves are reported here.
    controls = sorted(i for i in reachable if _row_op(nodes, i) == "original_Rm_unique_leading_implicit_control")
    rows = []
    for node in controls:
        row = _node(nodes, node) or {}
        digest, rational = fp.at(node, scope)
        rows.append(
            {
                "node": node,
                "operation": row.get("operation"),
                "system": row.get("system"),
                "component": row.get("component"),
                "Z_variable": row.get("Z_variable"),
                "N_handle": fp.b.get("N"),
                "same_unique_original_common_box_solution": row.get(
                    "exact_value_of_same_unique_original_common_box_solution"
                ) is True,
                "source_family": row.get("source_family"),
                "opaque_function_sha256": digest,
                "numeric_value": None,
                "rational_value": None,
                "passed_opaque": rational is None
                and row.get("exact_value_of_same_unique_original_common_box_solution") is True,
                "source_values_selected": False,
            }
        )
    return rows, bool(rows) and all(row["passed_opaque"] for row in rows) and len({row["opaque_function_sha256"] for row in rows}) == len(rows)


def build_Rm_inlet_mean_certificate(
    nodes: Any,
    leading: Mapping[str, Any],
    rm: Mapping[str, Any],
    l: Mapping[str, Any],
    modified: Mapping[str, Any],
    outer: Mapping[str, Any],
    *,
    root: Any = ROOT,
) -> dict[str, Any]:
    """Return the source-bound ordinary-Z leading mean identity at native one."""
    del modified  # The modified mean is intentionally outside this certificate.
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 20000))
    leading = _binding(leading, "true_Rm_patch_C6_bindings")
    rm = _binding(rm, "Rm_patch_finite_N_C5_bindings")
    l = _binding(l, "prescribed_loop_C5_bindings")
    outer = _binding(outer, "prescribed_R100_Rm_C6_bindings")
    repo = _repo_root(root)
    pins = {path.rsplit("/", 1)[-1]: _pin(repo, path, expected) for path, expected in _ARCHIVE_SHA256.items()}

    family_row = _node(nodes, leading.get("system"))
    family = family_row.get("source_family") if family_row is not None else leading.get("source_family")
    families = [
        family,
        rm.get("source_family"),
        l.get("source_family"),
        outer.get("source_family"),
    ]
    source_family_passed = (
        isinstance(family,Mapping) and bool(family)
        and all(value == family for value in families if value is not None)
    )
    chart = _mapping(leading.get("chart"))
    post = _mapping(outer.get("charts", {}).get("post_restore"))
    loop_post = _mapping(l.get("charts", {}).get("postrestore"))
    chart_coordinate = chart.get("phase")
    post_coordinate = post.get("phase")
    rm_coordinate = _mapping(rm.get("chart")).get("coordinate")
    N_handle = rm.get("N")
    if not all(isinstance(value, int) and not isinstance(value, bool) for value in (chart_coordinate, post_coordinate, rm_coordinate, N_handle)):
        raise ValueError("Actual leading, loop, and Rm endpoint handles required")
    if post_coordinate != loop_post.get("coordinate") or chart_coordinate != _mapping(rm.get("chart")).get("coordinate"):
        raise ValueError("Leading and Rm endpoint coordinate crosswalk changed")

    opaque_means = list(_mapping(post.get("fields")).get("m", ()))[:3]
    if len(opaque_means) != 3 or any(not isinstance(node, int) for node in opaque_means):
        raise ValueError("Three actual postrestore leading mean ports required")
    bindings = {
        "coordinate": chart_coordinate,
        "loop": l.get("charts", {}),
        "N": N_handle,
        "postrestore_coordinate": post_coordinate,
    }
    fp = _MeanFingerprint(nodes, bindings, opaque_means)
    scope, _, _ = _endpoint_scope(fp, leading, l, rm)
    wrapper_rows, wrapper_passed = _wrapper_rows(nodes, leading, outer, fp, scope)
    partial_rows, integral_rows, partial_passed = _partial_rows(nodes, leading, fp, scope)
    controls, controls_passed = _control_rows(nodes, leading, fp, scope)

    chart_means = list(_mapping(chart.get("fields")).get("m", ()))[:3]
    rows: list[dict[str, Any]] = []
    if len(chart_means) == 3:
        for order, (left, right) in enumerate(zip(opaque_means, chart_means)):
            left_value = fp.at(left, scope)
            right_value = fp.at(right, scope)
            equal = left_value[0] == right_value[0] and left_value[1] == right_value[1]
            rows.append(
                {
                    "ordinary_Z_order": order,
                    "left_node": left,
                    "right_node": right,
                    "left_function_sha256": left_value[0],
                    "right_function_sha256": right_value[0],
                    "common_function_sha256": left_value[0] if equal else None,
                    "passed": equal,
                    "endpoint_scope": {
                        "left": {"chart": "postrestore", "native_endpoint": 1},
                        "right": {"chart": "Rm_patch", "native_endpoint": 1},
                    },
                    "source_values_selected": False,
                    "numeric_value": None,
                }
            )
    pstar = outer.get("Pstar")
    pstar_row = _node(nodes, pstar)
    pstar_passed = (
        isinstance(pstar, int)
        and isinstance(pstar_row, Mapping)
        and pstar_row.get("operation") == "analytic_unary"
        and pstar_row.get("name") == "exp"
        and isinstance(pstar_row.get("argument"), int)
    )
    root_pstar_rows = []
    for node in chart_means:
        row = _node(nodes, node) or {}
        factors = row.get("arguments", ()) if row.get("operation") == "product" else ()
        pstar_factor = next(
            (
                factor
                for factor in factors
                if isinstance(factor, int)
                and (_node(nodes, factor) or {}).get("operation") == "positive_quotient"
                and (_node(nodes, factor) or {}).get("numerator") == 1
                and (_node(nodes, factor) or {}).get("denominator") == pstar
            ),
            None,
        )
        root_pstar_rows.append(
            {
                "node": node,
                "operation": row.get("operation"),
                "contains_Pstar_inverse": pstar_factor is not None,
                "Pstar_inverse_node": pstar_factor,
                "pstar_handle": pstar,
            }
        )
    formula_passed = bool(rows) and all(row["passed"] for row in rows) and pstar_passed and all(
        row["contains_Pstar_inverse"] for row in root_pstar_rows
    )

    return {
        "passed": bool(
            source_family_passed
            and all(pin["passed"] for pin in pins.values())
            and wrapper_passed
            and partial_passed
            and controls_passed
            and formula_passed
        ),
        "source_only": True,
        "source_values_selected": False,
        "numeric_values_selected": False,
        "numeric_caps_used": False,
        "assigned_numeric_values": False,
        "global_opaque_alias_used": False,
        "source_family": family,
        "source_family_consistent": source_family_passed,
        "source_family_from_original_unique_system": True,
        "source_pins": pins,
        "N_handle": N_handle,
        "endpoint_scope": {
            "left": {"chart": "postrestore", "native_endpoint": 1, "coordinate": post_coordinate},
            "right": {"chart": "Rm_patch", "native_endpoint": 1, "coordinate": rm_coordinate},
        },
        "ordinary_Z_orders": [0, 1, 2],
        "wrapper_source_identity": {
            "passed": wrapper_passed,
            "rows": wrapper_rows,
            "source_expression_is_actual_outer_post_restore_m": True,
            "substitution_variable_is_actual_post_restore_phase": post_coordinate,
            "substitution_value_is_one": True,
            "Z_independent_substitution": True,
        },
        "partial_weights": {
            "passed": partial_passed,
            "rows": partial_rows,
            "empty_integral_rows": integral_rows,
            "nonzero_weight_count": len(partial_rows),
            "empty_integral_count": len(integral_rows),
            "all_nonzero_partial_weights_exact_zero": partial_passed,
            "no_midpoint_or_numeric_weight_selected": True,
        },
        "controls": {
            "passed": controls_passed,
            "rows": controls,
            "opaque_source_functions": True,
            "distinct_control_leaves": len({row["node"] for row in controls}) == len(controls),
            "numeric_values_selected": False,
            "control_values_not_set_to_zero": True,
        },
        "rows": rows,
        "leading_mean_identity": {
            "passed": formula_passed,
            "ordinary_Z_orders": [0, 1, 2],
            "left_nodes_are_outer_post_restore_m": opaque_means,
            "right_nodes_are_Rm_patch_chart_m": chart_means,
            "Pstar_handle": pstar,
            "Pstar_log_handle": pstar_row.get("argument") if isinstance(pstar_row, Mapping) else None,
            "Pstar_source_formula": "exp(logP)",
            "Pstar_source_positive": pstar_passed,
            "root_rows": root_pstar_rows,
            "identity": "(4Z jet + Pstar*inlet_m - 4Z jet + empty partial weights)/Pstar = postrestore M_k",
            "exact_polynomial_ring_denominator_cancellation": True,
            "opaque_mean_leaves": True,
        },
        "limits": {
            "modified_mean_or_mean_y_alias": False,
            "modified_Qy_or_RQy_identity": False,
            "full_ancestry_expansion": False,
            "numeric_source_value": False,
            "numeric_error_cap": False,
            "global_opaque_control_alias": False,
        },
    }


__all__ = ["build_Rm_inlet_mean_certificate"]
