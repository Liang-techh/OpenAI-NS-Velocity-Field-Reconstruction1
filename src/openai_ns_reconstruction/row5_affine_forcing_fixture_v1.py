"""Independent synthetic checks for :mod:`row5_affine_forcing_v1`.

The fixture owns a small exact-rational source packet and caller-owned state
jets.  It receives ``build`` and ``evaluate`` explicitly, so it does not
import a producer, read a receipt, or execute an ancestor leaf.  Expected
values are accumulated with :class:`fractions.Fraction`; the implementation
under test is checked only for directed containment and provenance structure.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any, Mapping


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
STATE_COMPONENTS = (1, 2, 3, 5)
SENTINEL_SOURCE = "fixture:row5:source"
SENTINEL_AXIS = "fixture:row5:axis:wholeZ"


__all__ = [
    "SOURCE_POWER_BASIS",
    "STATE_COMPONENTS",
    "SENTINEL_SOURCE",
    "SENTINEL_AXIS",
    "run_fixture",
    "run",
]


def _clone(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _clone(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_clone(item) for item in value]
    return value


def _q(c: Any, value: Fraction | int) -> Any:
    value = Fraction(value)
    return c.mpf(value.numerator) / value.denominator


def _term(
    c: Any,
    value: Fraction | int,
    source_id: str,
    powers: tuple[int, int, int, int],
    path: Any,
    error_class: str,
    *,
    uncertainty: bool = False,
    source_ids: list[str] | None = None,
) -> dict[str, Any]:
    """Make one exact point source term while retaining fixture bookkeeping."""

    rational = Fraction(value)
    identifiers = list(source_ids if source_ids is not None else [source_id])
    return {
        "coefficient": _q(c, rational),
        "powers": list(powers),
        "source_ids": identifiers,
        "source_path": _clone(path),
        "error_classes": [error_class],
        "source_error_classes": [error_class],
        "uncertainty": bool(uncertainty),
        "error_only": False,
        "signed_field": not uncertainty,
        "fixture_value": rational,
        "fixture_source_id": source_id,
    }


def _axis(c: Any, value: Fraction | int, path: Any, error_class: str) -> dict[str, Any]:
    return {
        "coefficient": _q(c, value),
        "powers": [0, 0, 0, 0],
        "source_ids": [SENTINEL_AXIS],
        "source_path": _clone(path),
        "error_classes": [error_class],
        "source_error_classes": [error_class],
        "uncertainty": True,
        "whole_Z": True,
        "scope": "same-owner whole-Z axis condition",
        "fixture_value": Fraction(value),
    }


def _state(c: Any) -> tuple[dict[str, Any], dict[tuple[int, int], Fraction]]:
    atoms = {
        "hb": Fraction(2),
        "F0_anchor": Fraction(3),
        "Pstar_squared": Fraction(5),
        "Lambda": Fraction(7),
    }
    values: dict[tuple[int, int], Fraction] = {}
    rows: dict[int, dict[str, Any]] = {}
    for component in STATE_COMPONENTS:
        rows[component] = {}
        for order in (0, 1, 2):
            value = Fraction(10 * component + order + 1, component + order + 2)
            values[(component, order)] = value
            rows[component][f"Z{order}"] = _q(c, value)
    return {
        "formal_atoms": {name: _q(c, value) for name, value in atoms.items()},
        "W": rows,
    }, values


def _formal_value(value: Fraction, powers: list[int], atoms: Mapping[str, Fraction]) -> Fraction:
    result = value
    for name, power in zip(SOURCE_POWER_BASIS, powers):
        result *= atoms[name] ** int(power)
    return result


def _source_value(term: Mapping[str, Any]) -> Fraction:
    value = term.get("fixture_value")
    if not isinstance(value, Fraction):
        raise TypeError("fixture source term lost exact rational value")
    return value


def _expected(packet: Mapping[str, Any], chi: Mapping[str, Any], chi_z: Mapping[str, Any], values: Mapping[tuple[int, int], Fraction]) -> tuple[Fraction, Fraction]:
    atoms = {
        "hb": Fraction(2),
        "F0_anchor": Fraction(3),
        "Pstar_squared": Fraction(5),
        "Lambda": Fraction(7),
    }

    def contribution(term: Mapping[str, Any], component: int | None, order: int | None, scalar: int = 1) -> Fraction:
        value = _formal_value(_source_value(term), list(term["powers"]), atoms) * scalar
        if component is not None and order is not None:
            value *= values[(component, order)]
        return value

    h5 = sum((contribution(term, None, None) for term in packet["g5"]["Z0"]), Fraction(0))
    for name, keys, state_order in (("C0", ("5,1", "5,2", "5,3", "5,5"), 0), ("B1", ("5,1", "5,2", "5,3"), 1)):
        for key in keys:
            component = int(key.split(",")[1])
            h5 += sum((contribution(term, component, state_order) for term in packet[name]["Z0"][key]), Fraction(0))
    h5 += sum((contribution(term, 1, 0, 2) for term in [chi]), Fraction(0))

    h5_z1 = sum((contribution(term, None, None) for term in packet["g5"]["Z1"]), Fraction(0))
    for key in ("5,1", "5,2", "5,3", "5,5"):
        component = int(key.split(",")[1])
        h5_z1 += sum((contribution(term, component, 0) for term in packet["C0"]["Z1"][key]), Fraction(0))
        h5_z1 += sum((contribution(term, component, 1) for term in packet["C0"]["Z0"][key]), Fraction(0))
    for key in ("5,1", "5,2", "5,3"):
        component = int(key.split(",")[1])
        h5_z1 += sum((contribution(term, component, 1) for term in packet["B1"]["Z1"][key]), Fraction(0))
        h5_z1 += sum((contribution(term, component, 2) for term in packet["B1"]["Z0"][key]), Fraction(0))
    h5_z1 += sum((contribution(chi_z, 1, 0, 2), contribution(chi, 1, 1, 2)), Fraction(0))
    return h5, h5_z1


def _inside(c: Any, interval: Any, value: Fraction) -> bool:
    try:
        point = _q(c, value)
        # ``point`` is itself a directed enclosure after c.mpf/division.
        # Check containment by endpoints rather than asking an interval
        # object to perform a chained comparison against another interval;
        # this remains stable under the context's outward rounding.
        return bool(interval.a <= point.a and point.b <= interval.b)
    except (AttributeError, TypeError, ValueError, ArithmeticError):
        return False


def _callable(value: Any, name: str) -> Any:
    candidate = getattr(value, name, value)
    if not callable(candidate):
        raise TypeError(f"{name} must be callable or expose {name}")
    return candidate


def _rejected(call: Any) -> bool:
    try:
        call()
    except (TypeError, ValueError, KeyError, ZeroDivisionError):
        return True
    return False


def run_fixture(build: Any, evaluate: Any, c: Any) -> dict[str, Any]:
    """Run exact-rational product-rule, provenance, scope and guard checks."""

    build_fn = _callable(build, "build")
    evaluate_fn = _callable(evaluate, "evaluate")
    packet: dict[str, Any] = {
        "g5": {
            "Z0": [
                _term(c, Fraction(2, 3), "g5:Z0:a", (1, 0, 0, 0), ["g5", "Z0", 0], "g5-Z0-a"),
                _term(c, Fraction(-1, 5), "g5:Z0:b", (0, 1, 0, 0), ["g5", "Z0", 1], "g5-Z0-b", uncertainty=True,
                      source_ids=["g5:Z0:b", "g5:Z0:b"]),
            ],
            "Z1": [
                _term(c, Fraction(3, 7), "g5:Z1:a", (0, 0, 1, 0), ["g5", "Z1", 0], "g5-Z1-a"),
            ],
        },
        "C0": {"Z0": {}, "Z1": {}},
        "B1": {"Z0": {}, "Z1": {}},
    }
    for index, key in enumerate(("5,1", "5,2", "5,3", "5,5"), start=1):
        packet["C0"]["Z0"][key] = [_term(c, Fraction(index, 11), f"C0:Z0:{key}", (0, 0, 0, index - 2), ["C0", "Z0", key], f"C0-Z0-{key}")]
        packet["C0"]["Z1"][key] = [_term(c, Fraction(-index, 13), f"C0:Z1:{key}", (0, 0, 0, 2 - index), ["C0", "Z1", key], f"C0-Z1-{key}")]
    for index, key in enumerate(("5,1", "5,2", "5,3"), start=1):
        packet["B1"]["Z0"][key] = [_term(c, Fraction(index + 1, 17), f"B1:Z0:{key}", (0, 1, 0, 0), ["B1", "Z0", key], f"B1-Z0-{key}")]
        packet["B1"]["Z1"][key] = [_term(c, Fraction(-index - 1, 19), f"B1:Z1:{key}", (1, 0, 0, 0), ["B1", "Z1", key], f"B1-Z1-{key}")]

    chi = _axis(c, Fraction(1, 4), ["axis", "chi"], "axis-chi")
    chi_z = _axis(c, Fraction(-2, 9), ["axis", "chi_Z"], "axis-chi-Z1")
    state, state_values = _state(c)
    expected_h5, expected_h5_z1 = _expected(packet, chi, chi_z, state_values)

    checks: dict[str, bool] = {}
    try:
        expression = build_fn(c, packet, chi, chi_z)
        evaluation = evaluate_fn(c, expression, state)
        checks["build_completed"] = isinstance(expression, Mapping)
        checks["evaluate_completed"] = isinstance(evaluation, Mapping)
    except Exception as exc:
        return {
            "schema": "row5_affine_forcing_fixture_v1",
            "passed": False,
            "checks": {"build_completed": False},
            "check_count": 1,
            "error": f"{type(exc).__name__}: {exc}",
            "synthetic_only": True,
            "implementation_imported": False,
        }

    checks["schema"] = expression.get("schema") == "row5_affine_forcing_v1"
    checks["status_is_conditional"] = expression.get("status") == "conditional_source_bound_row5_affine_forcing"
    checks["formal_basis"] = tuple(expression.get("source_power_basis", ())) == SOURCE_POWER_BASIS
    checks["ordinary_fixed_x"] = expression.get("fixed_x_ordinary_Z") is True and expression.get("physical_dx_applied") is False
    checks["not_physical_h5"] = expression.get("physical_h5_materialized") is False and expression.get("ns_validation") is False
    checks["not_solved_state"] = expression.get("complete_W1") is False and expression.get("state_solved") is False
    checks["axis_owner_is_condition"] = expression.get("axis_same_owner_required") is True and expression.get("axis_same_owner_proved") is False
    checks["build_did_not_sum"] = expression.get("interval_sum_in_build") is False
    checks["h5_contains_nonzero_chi"] = any(term.get("role") == "two_chi_times_W1" and term.get("scalar") == 2 for term in expression["h5"])
    checks["h5_Z1_contains_chi_Z"] = any(term.get("role") == "two_chi_Z_times_W1" for term in expression["h5_Z1"])
    checks["h5_Z1_contains_W_ZZ"] = any(term.get("state_reference") == {"component": 1, "ordinary_Z_order": 2} for term in expression["h5_Z1"])
    checks["derivative_counts"] = len(expression["h5_Z1"]) == 1 + 2 * 4 + 2 * 3 + 2
    checks["source_duplicate_retained"] = any(term.get("source_ids") == ["g5:Z0:b", "g5:Z0:b"] for term in expression["h5"])
    checks["source_paths_retained"] = all(term.get("source_paths") for term in expression["h5"] + expression["h5_Z1"])
    checks["original_source_retained"] = all("original_source_term" in term and "original_term" in term for term in expression["h5"] + expression["h5_Z1"])
    checks["error_classes_not_mixed"] = all(
        len(term.get("error_classes", ())) == 1
        and term.get("error_classes") == term.get("source_error_classes")
        for term in expression["h5"] + expression["h5_Z1"]
    )
    checks["powers_copied"] = all(
        term.get("powers") == term.get("formal_powers") and len(term.get("powers", ())) == 4
        for term in expression["h5"] + expression["h5_Z1"]
    )
    checks["source_only_g5"] = all(term.get("source_only") is True and term.get("state_reference") is None for term in expression["h5"][:2])
    checks["evaluation_contains_h5"] = _inside(c, evaluation["h5"], expected_h5)
    checks["evaluation_contains_h5_Z1"] = _inside(c, evaluation["h5_Z1"], expected_h5_z1)
    checks["term_values_preserved"] = len(evaluation.get("term_values", {}).get("h5", ())) == len(expression["h5"])
    checks["state_atom_values_are_caller_owned"] = evaluation.get("formal_atoms_caller_owned") is True and evaluation.get("state_caller_owned") is True
    checks["evaluation_not_physical"] = evaluation.get("physical_h5_materialized") is False and evaluation.get("complete_W1") is False

    missing_state = _clone(state)
    del missing_state["W"][3]["Z2"]
    checks["missing_W_ZZ_rejected"] = _rejected(lambda: evaluate_fn(c, expression, missing_state))
    missing_source = _clone(packet)
    del missing_source["C0"]["Z0"]["5,5"]
    checks["missing_matrix_entry_rejected"] = _rejected(lambda: build_fn(c, missing_source, chi, chi_z))
    bad_axis = _clone(chi)
    bad_axis["powers"] = [1, 0, 0, 0]
    checks["nonzero_axis_powers_rejected"] = _rejected(lambda: build_fn(c, packet, bad_axis, chi_z))

    return {
        "schema": "row5_affine_forcing_fixture_v1",
        "passed": all(checks.values()),
        "check_count": len(checks),
        "checks": checks,
        "synthetic_only": True,
        "exact_fraction_expected": True,
        "ordinary_Z1_product_rule_checked": True,
        "cross_terms_checked": True,
        "W_ZZ_checked": True,
        "missing_state_rejected": True,
        "source_provenance_checked": True,
        "error_classes_kept_separate": True,
        "same_owner_not_proved": True,
        "physical_h5_not_claimed": True,
        "ns_validation_not_claimed": True,
        "implementation_imported": False,
    }


def run(build: Any, evaluate: Any, c: Any) -> dict[str, Any]:
    return run_fixture(build, evaluate, c)

