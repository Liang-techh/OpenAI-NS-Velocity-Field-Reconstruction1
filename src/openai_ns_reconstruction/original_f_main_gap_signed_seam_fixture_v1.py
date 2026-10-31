"""Independent in-memory fixture for :mod:`original_f_main_gap_signed_seam_v1`.

The fixture creates two fresh indexed graphs and supplies them directly to a
caller-provided ``build`` function.  It never imports an older seam, source
constructor, graph evaluator, or C3 endpoint helper.  The checks inspect only
the symbolic records returned by the caller and deliberately do not assign a
meaning to an integer graph handle.
"""

from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
from typing import Any, Callable


SCHEMA = "original_f_main_gap_signed_seam_fixture_v1"
F_SOURCE_KEY = "Rp_native_C5_pulse_source_bindings"
OWN_SOURCE_KEY = "own_main_gap_end_fields"
FAMILY_KEYS = (
    "actual_five_defect_family_sha256",
    "implicit_source_sha256",
    "datum_enclosure_sha256",
)


class _ProtectedGraph(Mapping[int, Any]):
    """Mapping that fails if a caller dereferences the protected handle."""

    def __init__(self, rows: Mapping[int, Any]):
        self._rows = dict(rows)
        self.dereferenced_protected = False

    def __getitem__(self, key: int) -> Any:
        if key == 33982:
            self.dereferenced_protected = True
            raise AssertionError("fixture: protected node 33982 was dereferenced")
        return self._rows[key]

    def __iter__(self):
        return iter(self._rows)

    def __len__(self) -> int:
        return len(self._rows)


def _family(tag: str = "synthetic-family") -> dict[str, str]:
    return {
        "actual_five_defect_family_sha256": tag + ":five",
        "implicit_source_sha256": tag + ":implicit",
        "datum_enclosure_sha256": tag + ":datum",
    }


def _graph_rows() -> dict[int, dict[str, Any]]:
    """Build a small fresh graph with direct radius-shape rows only."""

    rows: dict[int, dict[str, Any]] = {
        2: {"operation": "bound_variable", "name": "xi"},
        3: {"operation": "bound_variable", "name": "xi"},
        # Geometry rows are intentionally opaque to the seam.  Their handles
        # are used only as names in formal source atoms.
        4: {"operation": "source_parameter", "name": "delta"},
        5: {"operation": "source_parameter", "name": "mu"},
        6: {"operation": "source_parameter", "name": "Rp"},
        # R = Rp * exp(query / mu), for each independent chart.
        10: {"operation": "product", "arguments": [6, 11]},
        11: {"operation": "analytic_unary", "name": "exp", "argument": 12},
        12: {"operation": "positive_quotient", "numerator": 2, "denominator": 5},
        15: {"operation": "product", "arguments": [6, 16]},
        16: {"operation": "analytic_unary", "name": "exp", "argument": 17},
        17: {"operation": "positive_quotient", "numerator": 3, "denominator": 5},
        # Field ports are source expressions, not scalar values.
        20: {"operation": "source_field", "name": "main.U.Z0"},
        21: {"operation": "source_field", "name": "main.U.Z1"},
        22: {"operation": "source_field", "name": "main.Q.Z0"},
        23: {"operation": "source_field", "name": "main.Q.Z1"},
        24: {"operation": "source_field", "name": "main.Q_y.Z0"},
        25: {"operation": "source_field", "name": "main.Q_y.Z1"},
        26: {"operation": "source_field", "name": "gap.U.Z0"},
        27: {"operation": "source_field", "name": "gap.U.Z1"},
        28: {"operation": "source_field", "name": "gap.Q.Z0"},
        29: {"operation": "source_field", "name": "gap.Q.Z1"},
        30: {"operation": "source_field", "name": "gap.Q_y.Z0"},
        31: {"operation": "source_field", "name": "gap.Q_y.Z1"},
    }
    return rows


def _archives() -> tuple[dict[str, Any], dict[str, Any]]:
    family = _family()
    f_nodes = _graph_rows()
    own_nodes = _graph_rows()
    f_source = {
        "leading_fields": {
            "main": {"U": [20, 21], "Q": [22, 23], "Q_y": [24, 25]},
            "gap": {"U": [26, 27], "Q": [28, 29], "Q_y": [30, 31]},
        },
        "pulse_coordinates": {
            "main": {"query": 2, "radius": 10},
            "gap": {"query": 3, "radius": 15},
        },
        "geometry": {"delta": 4, "mu": 5, "Rp": 6},
        # This is deliberately retained by the consumer only as a candidate
        # reference and can never replace the F leading rows above.
        "inputs": {
            "lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C3_functions.json.gz": {
                "endpoints": {
                    "main": [{}, {name: [40 + 2*i, 41 + 2*i] for i, name in enumerate(('R_squared_q', 'R_squared_q_squared', 'R_squared_U_squared', 'R_q_y'))}],
                    "gap": [{name: [50 + 2*i, 51 + 2*i] for i, name in enumerate(('R_squared_q', 'R_squared_q_squared', 'R_squared_U_squared', 'R_q_y'))}, {}],
                }
            }
        },
    }
    own_source = {
        "main": {"query": 2, "domain": {"lower": 0, "upper": 18}},
        "gap": {"query": 3, "domain": {"lower": 19, "upper": 0}},
    }
    f_archive = {
        "source_family": family,
        F_SOURCE_KEY: f_source,
        "exact_graph_nodes": f_nodes,
    }
    own_archive = {
        "source_family": dict(family),
        OWN_SOURCE_KEY: own_source,
        "exact_graph_nodes": own_nodes,
    }
    # The selected own edges are the only exact-rational rows admitted by the
    # main leaf.  They are not used as arithmetic values in this fixture.
    own_nodes[18] = {"operation": "exact_rational", "numerator": 11, "denominator": 1}
    own_nodes[19] = {"operation": "exact_rational", "numerator": 11, "denominator": 1}
    return f_archive, own_archive


def _endpoint(result: Mapping[str, Any], key: str) -> Mapping[str, Any]:
    endpoints = result.get("endpoints")
    if not isinstance(endpoints, Mapping):
        return {}
    value = endpoints.get(key)
    return value if isinstance(value, Mapping) else {}


def _source_at_exact11(endpoint: Mapping[str, Any], field: str) -> bool:
    value = endpoint.get(field)
    if not isinstance(value, Mapping) or value.get("operator") != "source_at":
        return False
    args = value.get("arguments")
    if not isinstance(args, list) or len(args) != 3:
        return False
    return isinstance(args[2], Mapping) and args[2].get("typed_exact_rational_11") is True


def _contains_text(value: Any, *needles: str) -> bool:
    text = repr(value)
    return all(needle in text for needle in needles)


def _call(build: Callable[[Mapping[str, Any], Mapping[str, Any]], Mapping[str, Any]],
          F: Mapping[str, Any], own: Mapping[str, Any]) -> Mapping[str, Any] | None:
    try:
        value = build(F, own)
    except Exception:
        return None
    return value if isinstance(value, Mapping) else None


def run_fixture(build: Callable[[Mapping[str, Any], Mapping[str, Any]], Mapping[str, Any]]) -> dict[str, Any]:
    """Run independent positive and fail-closed synthetic checks.

    The returned ``checks`` values are plain booleans so a runner can combine
    them without importing a test framework.  No check asserts a numerical
    source value or a main/gap identity.
    """

    checks: dict[str, bool] = {}
    F, own = _archives()
    good = _call(build, F, own)
    checks["source_endpoints"] = bool(
        good is not None
        and good.get("status") == "constructed"
        and all(
            _source_at_exact11(_endpoint(good, key), field)
            for key in ("main.upper", "gap.lower")
            for field in ("R", "U", "Q", "Q_y")
        )
    )
    main = _endpoint(good or {}, "main.upper")
    gap = _endpoint(good or {}, "gap.lower")
    main_radius = main.get("R", {})
    gap_radius = gap.get("R", {})
    checks["independent_sides"] = bool(
        isinstance(main_radius, Mapping)
        and isinstance(gap_radius, Mapping)
        and main_radius.get("expression_ref", {}).get("source_node_id")
        != gap_radius.get("expression_ref", {}).get("source_node_id")
        and isinstance((good or {}).get("signed_seam"), Mapping)
        and (good or {}).get("signed_seam", {}).get("K", {}).get("operator") == "subtract"
        and (good or {}).get("signed_seam", {}).get("K_Z", {}).get("operator") == "subtract"
    )
    checks["full_R_Z_delta_Z_terms"] = bool(
        isinstance(main.get("R_Z"), Mapping)
        and main["R_Z"].get("unresolved") is True
        and main["R_Z"].get("shape_candidate") is not None
        and _contains_text(main["R_Z"], "Rp_Z", "mu_Z")
        and _contains_text(main.get("geometry", {}).get("L_Z"), "delta_Z")
        and main.get("geometry", {}).get("delta_Z", {}).get("operator") == "partial_Z_at_fixed_query"
    )
    checks["old_C3_reference_only"] = bool(
        isinstance(good, Mapping)
        and good.get("candidate_old_C3_endpoints", {}).get("present") is True
        and all(len(good.get("candidate_old_C3_endpoints", {}).get("selected_endpoint_refs", {}).get(chart, {})) == 4 for chart in ('main', 'gap'))
        and good.get("candidate_old_C3_endpoints", {}).get("used_as_endpoint_source") is False
        and good.get("candidate_old_C3_endpoints", {}).get("equal_to_new_F_endpoint_source_proved") is False
    )
    checks["product_quotient_sources_retained"] = bool(
        isinstance(main.get("K"), Mapping)
        and _contains_text(main["K"], "R^2*Q/(2*L)", "R^2*Q^2/4", "Z^2*R^2*U^2/L^2")
        and isinstance(main.get("K_Z"), Mapping)
        and _contains_text(main["K_Z"], "L_Z", "C_Z", "D_Z")
    )
    checks["no_numeric_or_interval_result"] = bool(
        isinstance(good, Mapping)
        and good.get("numeric_values_computed") is False
        and good.get("interval_result") is False
        and good.get("graph_handles_used_as_values") is False
    )

    missing = deepcopy(F)
    missing[F_SOURCE_KEY]["leading_fields"]["main"]["Q"] = [22]
    missing_result = _call(build, missing, deepcopy(own))
    checks["missing_Z1_unresolved"] = bool(
        missing_result is not None
        and missing_result.get("status") == "unresolved"
        and missing_result.get("passed") is False
    )

    wrong_family = deepcopy(own)
    wrong_family["source_family"][FAMILY_KEYS[0]] = "synthetic-family:wrong"
    wrong_family_result = _call(build, deepcopy(F), wrong_family)
    checks["wrong_family_rejected"] = bool(
        wrong_family_result is not None
        and wrong_family_result.get("status") == "unresolved"
        and wrong_family_result.get("source_family", {}).get("comparison", {}).get("passed") is False
    )

    wrong_edge = deepcopy(own)
    wrong_edge["exact_graph_nodes"][18] = {
        "operation": "exact_rational",
        "numerator": 10,
        "denominator": 1,
    }
    wrong_edge_result = _call(build, deepcopy(F), wrong_edge)
    checks["wrong_edge_rejected"] = bool(
        wrong_edge_result is not None
        and wrong_edge_result.get("status") == "unresolved"
        and wrong_edge_result.get("audit", {}).get("endpoint_checks", {}).get("main.upper", {}).get("edge_check", {}).get("passed") is False
    )

    protected_own = deepcopy(own)
    guarded = _ProtectedGraph(protected_own["exact_graph_nodes"])
    protected_own["exact_graph_nodes"] = guarded
    protected_own[OWN_SOURCE_KEY]["main"]["domain"]["upper"] = 33982
    protected_result = _call(build, deepcopy(F), protected_own)
    protected_check = (
        protected_result is not None
        and protected_result.get("status") == "unresolved"
        and protected_result.get("audit", {}).get("endpoint_checks", {}).get("main.upper", {}).get("edge_check", {}).get("protected_node_not_dereferenced") is True
        and guarded.dereferenced_protected is False
    )
    checks["protected_33982_not_dereferenced"] = bool(protected_check)

    checks = {name: bool(value) for name, value in checks.items()}
    return {
        "schema": SCHEMA,
        "synthetic_only": True,
        "passed": bool(all(checks.values())),
        "checks": checks,
        "source_values_computed": False,
        "numeric_values_computed": False,
        "old_helpers_called": False,
    }


__all__ = ["SCHEMA", "run_fixture"]
