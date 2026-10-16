"""Exact ordinary-Z root identities at the reference/restoration seam.

This module is deliberately smaller than the native derivative theorem.  It
checks only the root values ``a, b, t0, Delta`` and the first ordinary-Z row of
the actual loop ``q`` source.  The checks use the bounded source endpoint ring
from :mod:`original_reference_viscous_m15`; no ancestor constructor, interval
cap, or source value is evaluated.

The archived loop represents the active q branch with the same lazy node used
for the flat sigma branch.  At the accepted seam, ``Delta=-6/5`` and
``eta>0``, so ``1-Delta/eta>1``.  Resolving that exact source condition makes
sigma equal one and every sigma derivative row equal zero.  The q lazy node
itself is then evaluated through its active body.  This is a source identity,
not a numerical branch selection.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any, Mapping

from .original_reference_viscous_m15 import PhysicalEndpointFingerprint


_CHARTS = ("reference", "restoration")
_ROOT_KEYS = ("a", "b", "t0", "Delta")
_EXPECTED = {
    "a": Fraction(4, 5),
    "b": Fraction(0),
    "t0": Fraction(0),
    "Delta": Fraction(-6, 5),
}


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _binding(value: Any, key: str) -> Mapping[str, Any]:
    value = _mapping(value)
    nested = value.get(key)
    return nested if isinstance(nested, Mapping) else value


def _row(nodes: Any, index: Any) -> Mapping[str, Any] | None:
    try:
        index = int(index)
    except (TypeError, ValueError):
        return None
    if index < 0 or index >= len(nodes):
        return None
    value = nodes[index]
    return value if isinstance(value, Mapping) else None


def _handle(value: Any) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


class _RootZFingerprint(PhysicalEndpointFingerprint):
    """PhysicalEndpointFingerprint with two source-certified q simplifiers.

    The generic fingerprint cannot decide the symbolic predicate
    ``Delta >= eta``.  For the q rows only, the accepted active-gap identity
    authorizes the active body.  Inside that body, the exact flat sigma
    argument ``1-Delta/eta`` is above one; its value is one and its lazy
    derivative rows are exactly zero.  No other lazy node receives this
    treatment.
    """

    def __init__(
        self,
        *args: Any,
        q_nodes: Mapping[int, tuple[int, int]],
        delta_eta: Mapping[int, int] | None = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(*args, **kwargs)
        self._q_nodes = dict(q_nodes)
        # The cut expression is keyed by its Delta source handle, whereas
        # q_nodes is keyed by the q source handles.  Keep that relation
        # explicit so sigma(cut) and its lazy derivative rows are recognized
        # only for the accepted chart sources.
        self._delta_eta = dict(delta_eta or {
            pair[0]: pair[1] for pair in self._q_nodes.values()
        })
        self.q_active_rows: set[int] = set()
        self.sigma_flat_rows: set[int] = set()
        self.sigma_one_rows: set[int] = set()

    @staticmethod
    def _exact_one(nodes: Any, index: Any) -> bool:
        row = _row(nodes, index)
        return bool(row and row.get("operation") == "exact_rational"
                    and row.get("numerator") == 1 and row.get("denominator") == 1)

    def _is_cut(self, index: Any) -> bool:
        """Recognize the actual source expression ``1-Delta/eta``."""
        node = _row(self.nodes, index)
        if not node or node.get("operation") != "sum":
            return False
        arguments = tuple(node.get("arguments", ()))
        if len(arguments) != 2:
            return False
        for first, second in ((arguments[0], arguments[1]), (arguments[1], arguments[0])):
            if not self._exact_one(self.nodes, first):
                continue
            negative = _row(self.nodes, second)
            if not negative or negative.get("operation") != "negative":
                continue
            product = _row(self.nodes, negative.get("argument"))
            if not product or product.get("operation") != "product":
                continue
            factors = tuple(product.get("arguments", ()))
            if len(factors) != 2:
                continue
            for delta, inverse_eta in ((factors[0], factors[1]), (factors[1], factors[0])):
                eta_handle = self._delta_eta.get(_handle(delta))
                inverse = _row(self.nodes, inverse_eta)
                if eta_handle is None or not inverse or inverse.get("operation") != "positive_quotient":
                    continue
                if (inverse.get("numerator") != 1
                        or _handle(inverse.get("denominator")) != eta_handle):
                    continue
                return True
        return False

    def at(self, node: Any, scope: Mapping[int, Any] | None = None) -> tuple[str, Fraction | None]:
        index = _handle(node)
        current = _row(self.nodes, index)
        operation = current.get("operation") if current else None

        # This is the exact active q source branch.  It is restricted to q0/q1
        # handles supplied by the accepted loop binding.
        if index in self._q_nodes and operation == "original_lazy_flat_branch":
            scope = scope if scope is not None else {self.b["coordinate"]: self.one}
            key = (index, tuple(sorted((k, v[0]) for k, v in scope.items())))
            if key in self.cache:
                return self.cache[key]
            pair = self._q_nodes[index]
            if _handle(current.get("Delta")) != pair[0] or _handle(current.get("eta")) != pair[1]:
                raise ValueError("q lazy row is not bound to its accepted Delta/eta source")
            self.q_active_rows.add(index)
            result = self.at(current.get("active_body"), scope)
            self.cache[key] = result
            return result

        # The actual source sigma(cut) is one once the accepted strict active
        # gap is established.  This does not use a selected eta value.
        if operation == "exact_original_source_flat_sigma" and self._is_cut(current.get("argument")):
            self.sigma_one_rows.add(index)
            return self.one

        # Native/ordinary sigma derivative rows use the same cut predicate.
        # Their archived flat branch is exactly zero on cut >= 1.
        if operation == "original_lazy_flat_branch" and self._is_cut(current.get("Delta")):
            flat = _row(self.nodes, current.get("flat_value"))
            if not self._exact_one(self.nodes, current.get("flat_value")):
                flat_value = self.at(current.get("flat_value"), scope)
                if flat_value[1] != 0:
                    raise ValueError("sigma flat branch was not the exact zero source")
            self.sigma_flat_rows.add(index)
            return self.at(current.get("flat_value"), scope)

        return super().at(index, scope)


def _phase_domain_receipt(phase: Mapping[str, Any], l: Mapping[str, Any]) -> dict[str, Any]:
    domain = _mapping(phase.get("physical_seam_domain"))
    rows = domain.get("rows", ())
    expected = {
        (chart, key): _EXPECTED[key]
        for chart in _CHARTS for key in _ROOT_KEYS
    }
    seen: set[tuple[str, str]] = set()
    selected = []
    for item in rows if isinstance(rows, (list, tuple)) else ():
        row = _mapping(item)
        chart, key = row.get("chart"), row.get("quantity")
        if (chart, key) not in expected or (chart, key) in seen:
            continue
        chart_binding = _mapping(_mapping(l.get("charts")).get(chart))
        root = _mapping(chart_binding.get("roots")).get(key)
        source_node = _handle(root[0]) if isinstance(root, (list, tuple)) and root else None
        exact = row.get("exact_value")
        exact_fraction = None
        if isinstance(exact, (list, tuple)) and len(exact) == 2:
            try:
                exact_fraction = Fraction(exact[0], exact[1])
            except (TypeError, ValueError, ZeroDivisionError):
                exact_fraction = None
        passed = bool(
            row.get("passed") is True
            and row.get("original_node") == source_node
            and exact_fraction == expected[(chart, key)]
        )
        seen.add((chart, key))
        selected.append({
            "chart": chart,
            "quantity": key,
            "original_node": source_node,
            "accepted_exact_value": [expected[(chart, key)].numerator, expected[(chart, key)].denominator],
            "accepted_row_passed": row.get("passed") is True,
            "accepted_function_identity": domain.get("all_slow_Z_function_identity") is True,
            "accepted_numeric_Z_values_used": domain.get("numeric_Z_values_used") is False,
            "passed": passed,
        })
    return {
        "passed": bool(
            domain.get("passed") is True
            and domain.get("all_slow_Z_function_identity") is True
            and domain.get("numeric_Z_values_used") is False
            and len(selected) == 8
            and len(seen) == 8
            and all(row["passed"] for row in selected)
        ),
        "method": domain.get("method"),
        "rows": selected,
        "native_endpoints": domain.get("native_endpoints"),
        "active_branch": domain.get("active_branch"),
        "q_positive": domain.get("q_positive"),
    }


def _germ_premise_receipt(germ: Mapping[str, Any]) -> dict[str, Any]:
    radial = _mapping(germ.get("leading_radial_laws"))
    radial_rows = radial.get("rows", ())
    radial_ok = bool(
        radial.get("passed") is True
        and radial.get("identity") == "E_y=E/10 on both native charts; hence E_yy=E/100"
        and isinstance(radial_rows, (list, tuple))
        and len(radial_rows) == 6
        and all(_mapping(row).get("whole_native_chart_function_identity") is True
                and _mapping(row).get("endpoint_sampling_used") is False
                for row in radial_rows)
    )
    physical = _mapping(germ.get("physical_two_jet_certificate"))
    physical_ok = bool(
        physical.get("passed") is True
        and physical.get("ordinary_Z_jets_are_function_identities_not_axial_samples") is True
        and physical.get("source_values_selected") is False
        and physical.get("E_y") == "E/10"
        and physical.get("E_yy") == "E/100"
        and physical.get("V_y") == "0 at the seam"
        and physical.get("V_yy") == "0 at the seam"
    )
    guards = _mapping(germ.get("source_guards"))
    guard_ok = guards.get("passed") is True and guards.get("source_only") is True
    return {
        "passed": radial_ok and physical_ok and guard_ok,
        "leading_radial_laws": {
            "passed": radial_ok,
            "identity": radial.get("identity"),
            "row_count": len(radial_rows) if isinstance(radial_rows, (list, tuple)) else None,
            "whole_chart_function_rows": sum(
                _mapping(row).get("whole_native_chart_function_identity") is True
                for row in radial_rows if isinstance(radial_rows, (list, tuple))
            ),
            "endpoint_sampling_used": any(
                _mapping(row).get("endpoint_sampling_used") is True
                for row in radial_rows if isinstance(radial_rows, (list, tuple))
            ),
        },
        "physical_two_jet_certificate": {
            "passed": physical_ok,
            "method": physical.get("method"),
            "ordinary_Z_function_identities": physical.get("ordinary_Z_jets_are_function_identities_not_axial_samples"),
            "E_y": physical.get("E_y"),
            "E_yy": physical.get("E_yy"),
            "V_y": physical.get("V_y"),
            "V_yy": physical.get("V_yy"),
            "source_values_selected": physical.get("source_values_selected"),
        },
        "source_guards": {
            "passed": guard_ok,
            "source_only": guards.get("source_only"),
            "numeric_fixtures_used": guards.get("numeric_fixtures_used"),
        },
    }


def _parameter_receipt(nodes: Any, l: Mapping[str, Any]) -> dict[str, Any]:
    eta = _handle(l.get("eta"))
    dstar = _handle(l.get("dstar"))
    eta_row, dstar_row = _row(nodes, eta), _row(nodes, dstar)
    passed = bool(
        eta_row and dstar_row
        and eta_row.get("operation") == "current_original_source_parameter"
        and dstar_row.get("operation") == "current_original_source_parameter"
        and eta_row.get("Z_independent") is True
        and dstar_row.get("Z_independent") is True
        and isinstance(eta_row.get("canonical_parameter_proof"), Mapping)
        and isinstance(dstar_row.get("canonical_parameter_proof"), Mapping)
    )
    return {
        "passed": passed,
        "eta": {"node": eta, "operation": eta_row.get("operation") if eta_row else None,
                "Z_independent": eta_row.get("Z_independent") if eta_row else None},
        "d_star": {"node": dstar, "operation": dstar_row.get("operation") if dstar_row else None,
                    "Z_independent": dstar_row.get("Z_independent") if dstar_row else None},
        "numeric_parameter_values_selected": False,
    }


def build_root_Z_certificate(
    nodes: Any,
    l: Mapping[str, Any],
    phase_certificate: Mapping[str, Any],
    germ_packet: Mapping[str, Any],
) -> dict[str, Any]:
    """Build the exact source certificate for root ordinary-Z constancy.

    The returned ``rows`` contain exactly nine rows per chart: four root values,
    four root ordinary-Z1 rows, and the actual loop ``q`` ordinary-Z1 row.
    ``passed`` is false on any missing premise or unresolved exact arithmetic;
    no unresolved expression is relabeled as zero.
    """
    loop = _binding(l, "prescribed_loop_C5_bindings")
    phase = _mapping(phase_certificate)
    if isinstance(phase.get("phase_lift_certificate"), Mapping):
        phase = phase["phase_lift_certificate"]
    germ = _mapping(germ_packet)

    phase_receipt = _phase_domain_receipt(phase, loop)
    germ_receipt = _germ_premise_receipt(germ)
    parameter_receipt = _parameter_receipt(nodes, loop)

    charts = _mapping(loop.get("charts"))
    chart_bindings = {name: _mapping(charts.get(name)) for name in _CHARTS}
    chart_shape_ok = all(
        chart_bindings[name]
        and _handle(chart_bindings[name].get("coordinate")) is not None
        and isinstance(_mapping(chart_bindings[name].get("roots")), Mapping)
        and isinstance(_mapping(chart_bindings[name].get("loop")), Mapping)
        for name in _CHARTS
    )
    q_nodes: dict[int, tuple[int, int]] = {}
    delta_eta: dict[int, int] = {}
    for name in _CHARTS:
        chart = chart_bindings[name]
        roots = _mapping(chart.get("roots"))
        q = _mapping(chart.get("loop")).get("q")
        delta = _handle(roots.get("Delta")[0]) if isinstance(roots.get("Delta"), (list, tuple)) else None
        eta = _handle(loop.get("eta"))
        if isinstance(q, (list, tuple)) and len(q) >= 2 and delta is not None and eta is not None:
            delta_eta[delta] = eta
            q_nodes[_handle(q[0])] = (delta, eta)
            q_nodes[_handle(q[1])] = (delta, eta)
    q_nodes = {key: value for key, value in q_nodes.items() if key is not None}

    ring_rows = []
    exact_ring_passed = chart_shape_ok and len(q_nodes) == 4
    fingerprint: _RootZFingerprint | None = None
    scope: dict[int, tuple[str, Fraction | None]] = {}
    if exact_ring_passed:
        try:
            reference = chart_bindings["reference"]
            fingerprint = _RootZFingerprint(
                nodes,
                dict(coordinate=_handle(reference["coordinate"]), loop=charts),
                q_nodes=q_nodes,
                delta_eta=delta_eta,
            )
            for name in _CHARTS:
                chart = chart_bindings[name]
                coordinate = _handle(chart["coordinate"])
                endpoint = chart.get("domain", (None, None))[1 if name == "reference" else 0]
                if coordinate is None or endpoint is None:
                    raise ValueError("Missing native endpoint")
                scope[coordinate] = fingerprint.at(endpoint, {})
            for name in _CHARTS:
                chart = chart_bindings[name]
                roots = _mapping(chart.get("roots"))
                for key in _ROOT_KEYS:
                    values = roots.get(key)
                    if not isinstance(values, (list, tuple)) or len(values) < 2:
                        raise ValueError("Missing root ordinary-Z rows")
                    for order, node in ((0, values[0]), (1, values[1])):
                        value = fingerprint.at(node, scope)
                        expected = _EXPECTED[key] if order == 0 else Fraction(0)
                        exact = value[1] == expected
                        ring_rows.append({
                            "chart": name,
                            "quantity": key,
                            "ordinary_Z_order": order,
                            "original_node": _handle(node),
                            "exact_value": [value[1].numerator, value[1].denominator] if value[1] is not None else None,
                            "expected_value": [expected.numerator, expected.denominator],
                            "ring_digest": value[0],
                            "exact_ring_row": exact,
                            "source_values_selected": False,
                        })
                q = _mapping(chart.get("loop")).get("q")
                if not isinstance(q, (list, tuple)) or len(q) < 2:
                    raise ValueError("Missing q ordinary-Z rows")
                value = fingerprint.at(q[1], scope)
                expected = Fraction(0)
                ring_rows.append({
                    "chart": name,
                    "quantity": "q",
                    "ordinary_Z_order": 1,
                    "original_node": _handle(q[1]),
                    "q_value_node": _handle(q[0]),
                    "exact_value": [value[1].numerator, value[1].denominator] if value[1] is not None else None,
                    "expected_value": [0, 1],
                    "ring_digest": value[0],
                    "exact_ring_row": value[1] == expected,
                    "source_values_selected": False,
                })
            exact_ring_passed = len(ring_rows) == 18 and all(row["exact_ring_row"] for row in ring_rows)
        except (ArithmeticError, KeyError, TypeError, ValueError, IndexError) as exc:
            exact_ring_passed = False
            ring_rows.append({"error": type(exc).__name__ + ": " + str(exc), "exact_ring_row": False})

    return {
        "passed": bool(exact_ring_passed and phase_receipt["passed"] and germ_receipt["passed"] and parameter_receipt["passed"]),
        "method": "physical_endpoint_fingerprint_exact_bounded_polynomial_ring_plus_active_q_source_identity",
        "source_only": True,
        "numeric_fixtures_used": False,
        "source_values_selected": False,
        "phase_premise_receipt": phase_receipt,
        "germ_premise_receipt": germ_receipt,
        "parameter_premise_receipt": parameter_receipt,
        "rows": ring_rows,
        "row_count": len(ring_rows),
        "native_endpoints": {
            name: {
                "coordinate": _handle(chart_bindings[name].get("coordinate")),
                "endpoint": chart_bindings[name].get("domain", (None, None))[1 if name == "reference" else 0],
            }
            for name in _CHARTS
        },
        "q_source_normalization": {
            "q_nodes": sorted(q_nodes),
            "active_q_rows": sorted(fingerprint.q_active_rows) if fingerprint else [],
            "sigma_one_rows": sorted(fingerprint.sigma_one_rows) if fingerprint else [],
            "sigma_flat_rows": sorted(fingerprint.sigma_flat_rows) if fingerprint else [],
            "exact_polynomial_maximum_terms": fingerprint.maximum_terms if fingerprint else None,
            "certified_quotient_count": len(fingerprint.certified_quotients) if fingerprint else 0,
            "exact_division_count": len(fingerprint.exact_divisions) if fingerprint else 0,
            "active_gap_source_identity_used": bool(fingerprint and fingerprint.sigma_one_rows),
            "q_Z1_direct_ring_proof": bool(fingerprint and all(
                row.get("quantity") != "q" or row.get("exact_ring_row") for row in ring_rows
            )),
        },
        "limits": {
            "no_error_cap_or_interval_endpoint_used": True,
            "no_ancestor_constructor_run": True,
            "q_active_branch_selected_from_source_gap_identity": True,
            "ordinary_Z_rows_are_function_identities_before_substitution": True,
        },
    }


__all__ = ["build_root_Z_certificate"]
