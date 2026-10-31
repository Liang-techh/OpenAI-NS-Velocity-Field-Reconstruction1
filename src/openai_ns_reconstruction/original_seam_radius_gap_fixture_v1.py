"""Independent synthetic dual-number fixture for the radius-gap seam.

The fixture never imports the seam constructor and never uses its factorized
formula to form an expected answer.  ``run_fixture(build)`` supplies twelve
opaque synthetic tokens to a caller-provided builder, parses the returned
restricted JSON DAG, and evaluates it with local floating-point dual
arithmetic.  The independent expectation is formed directly from
``Ro=Rf*exp(D)``, ``K_old=Ro**2*Ho-Ro*qyo``, and
``K_new=Rf**2*Hf-Rf*qyf``.  A second scenario sets ``D=0`` while keeping
``Ho`` and ``Hf`` different, so the endpoint field difference cannot silently
disappear.

This module contains no source paths, graph handles, archive readers, gzip
logic, or module-level execution.  Its values are synthetic doubles only;
they are not source values and provide no certification of owner identity,
coordinate or unit semantics, error transport, or seam cancellation.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence
from typing import Any


SCHEMA = "original_seam_radius_gap_fixture_v1"
PORT_NAMES = (
    "Rf",
    "Rf_Z",
    "D",
    "D_Z",
    "Ho",
    "Ho_Z",
    "Hf",
    "Hf_Z",
    "qyo",
    "qyo_Z",
    "qyf",
    "qyf_Z",
)
_BASE_NAMES = ("Rf", "D", "Ho", "Hf", "qyo", "qyf")


class _Dual:
    """A value/ordinary-Z1 pair used only for synthetic arithmetic."""

    __slots__ = ("value", "derivative")

    def __init__(self, value: float, derivative: float = 0.0) -> None:
        self.value = float(value)
        self.derivative = float(derivative)

    def __add__(self, other: Any) -> "_Dual":
        right = _as_dual(other)
        return _Dual(self.value + right.value, self.derivative + right.derivative)

    def __radd__(self, other: Any) -> "_Dual":
        return self + other

    def __neg__(self) -> "_Dual":
        return _Dual(-self.value, -self.derivative)

    def __sub__(self, other: Any) -> "_Dual":
        return self + (-_as_dual(other))

    def __rsub__(self, other: Any) -> "_Dual":
        return _as_dual(other) - self

    def __mul__(self, other: Any) -> "_Dual":
        right = _as_dual(other)
        return _Dual(
            self.value * right.value,
            self.derivative * right.value + self.value * right.derivative,
        )

    def __rmul__(self, other: Any) -> "_Dual":
        return self * other

    def power(self, exponent: int) -> "_Dual":
        if exponent < 0:
            raise ValueError("fixture only admits nonnegative integer powers")
        if exponent == 0:
            return _Dual(1.0)
        value = self.value**exponent
        derivative = exponent * (self.value ** (exponent - 1)) * self.derivative
        return _Dual(value, derivative)

    def exp(self) -> "_Dual":
        value = math.exp(self.value)
        return _Dual(value, value * self.derivative)

    def expm1(self) -> "_Dual":
        value = math.expm1(self.value)
        return _Dual(value, math.exp(self.value) * self.derivative)


def _as_dual(value: Any) -> _Dual:
    if isinstance(value, _Dual):
        return value
    return _Dual(float(value))


# Every derivative entry is nonzero and independent in the primary scenario.
# These numbers are fixture-local doubles, never values copied from a source
# graph or from a source token.
_PRIMARY: dict[str, float] = {
    "Rf": 1.73,
    "Rf_Z": 0.29,
    "D": 0.43,
    "D_Z": -0.17,
    "Ho": 2.41,
    "Ho_Z": 0.38,
    "Hf": -0.76,
    "Hf_Z": -0.23,
    "qyo": 0.67,
    "qyo_Z": -0.31,
    "qyf": -0.44,
    "qyf_Z": 0.19,
}

# D is exactly zero here; equal q fields isolate the retained Ho-Hf value
# difference while D_Z remains nonzero for the derivative path.
_ZERO_D: dict[str, float] = {
    "Rf": 1.21,
    "Rf_Z": -0.36,
    "D": 0.0,
    "D_Z": 0.27,
    "Ho": 2.8,
    "Ho_Z": -0.41,
    "Hf": -0.6,
    "Hf_Z": 0.22,
    "qyo": 0.35,
    "qyo_Z": 0.16,
    "qyf": 0.35,
    "qyf_Z": 0.16,
}


def _synthetic_ports() -> dict[str, Any]:
    """Return opaque token records; numeric fixture data stay in this module."""

    return {
        name: {
            "fixture_token": name,
            "token_is_opaque": True,
            "source_values_selected": False,
        }
        for name in PORT_NAMES
    }


def _port_dual(port: str, values: Mapping[str, float]) -> _Dual:
    if port in _BASE_NAMES:
        return _Dual(values[port], values[port + "_Z"])
    if port in PORT_NAMES and port.endswith("_Z"):
        return _Dual(values[port])
    raise ValueError("fixture encountered an unknown source port")


def _direct_gap(values: Mapping[str, float]) -> _Dual:
    """Compute the independent old-minus-new expectation using exp directly."""

    Rf = _port_dual("Rf", values)
    D = _port_dual("D", values)
    Ho = _port_dual("Ho", values)
    Hf = _port_dual("Hf", values)
    qyo = _port_dual("qyo", values)
    qyf = _port_dual("qyf", values)

    Ro = Rf * D.exp()
    old_K = Ro * Ro * Ho - Ro * qyo
    new_K = Rf * Rf * Hf - Rf * qyf
    return old_K - new_K


def _literal(node: Mapping[str, Any]) -> _Dual:
    text = node.get("text")
    if not isinstance(text, str):
        raise ValueError("exact_literal text is missing")
    return _Dual(float(text))


def _arguments(node: Mapping[str, Any]) -> list[Any]:
    arguments = node.get("arguments")
    if not isinstance(arguments, Sequence) or isinstance(arguments, (str, bytes, bytearray)):
        raise ValueError("DAG arguments are not a sequence")
    return list(arguments)


def _eval(node: Any, values: Mapping[str, float]) -> _Dual:
    """Evaluate only the small operator vocabulary emitted by the constructor."""

    if not isinstance(node, Mapping):
        raise ValueError("DAG node is not a mapping")
    if node.get("source_values_selected") is not False:
        raise ValueError("fixture refuses a node with selected source values")
    operator = node.get("operator")
    args = _arguments(node)

    if operator == "opaque_input_token":
        port = node.get("port")
        if not isinstance(port, str):
            raise ValueError("opaque input port is missing")
        return _port_dual(port, values)
    if operator == "exact_literal":
        return _literal(node)
    if operator in ("add", "sum"):
        total = _Dual(0.0)
        for argument in args:
            total = total + _eval(argument, values)
        return total
    if operator in ("multiply", "product"):
        total = _Dual(1.0)
        for argument in args:
            total = total * _eval(argument, values)
        return total
    if operator in ("subtract",):
        if len(args) < 2:
            raise ValueError("subtract expects at least two arguments")
        result = _eval(args[0], values)
        for argument in args[1:]:
            result = result - _eval(argument, values)
        return result
    if operator in ("negate", "negative"):
        if len(args) != 1:
            raise ValueError("negate expects one argument")
        return -_eval(args[0], values)
    if operator == "power":
        if len(args) != 2:
            raise ValueError("power expects base and exponent")
        base = _eval(args[0], values)
        exponent = _eval(args[1], values)
        if exponent.derivative != 0.0 or not exponent.value.is_integer():
            raise ValueError("fixture power exponent must be constant integer")
        return base.power(int(exponent.value))
    if operator == "exp":
        if len(args) != 1:
            raise ValueError("exp expects one argument")
        return _eval(args[0], values).exp()
    if operator == "expm1":
        if len(args) != 1:
            raise ValueError("expm1 expects one argument")
        return _eval(args[0], values).expm1()
    raise ValueError("unsupported operator in restricted fixture parser")


def _close(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=4e-12, abs_tol=4e-12)


def _output_node(built: Mapping[str, Any], name: str) -> Any:
    direct = built.get(name)
    if isinstance(direct, Mapping):
        return direct
    outputs = built.get("outputs")
    if isinstance(outputs, Mapping) and isinstance(outputs.get(name), Mapping):
        return outputs[name]
    ordinary = built.get("ordinary_Z1")
    if isinstance(ordinary, Mapping) and isinstance(ordinary.get(name), Mapping):
        return ordinary[name]
    raise ValueError("builder result is missing " + name)


def _scenario_checks(
    delta_node: Any,
    delta_z_node: Any,
    values: Mapping[str, float],
    *,
    require_nonzero_D: bool,
) -> dict[str, bool]:
    actual = _eval(delta_node, values)
    actual_z1 = _eval(delta_z_node, values)
    expected = _direct_gap(values)
    checks = {
        "delta_K_value_matches_direct_exp_dual": _close(actual.value, expected.value),
        "delta_K_dual_derivative_matches_direct_exp_dual": _close(
            actual.derivative, expected.derivative
        ),
        "delta_K_Z_matches_direct_exp_dual": _close(actual_z1.value, expected.derivative),
        "delta_K_Z_matches_delta_K_dual_derivative": _close(
            actual_z1.value, actual.derivative
        ),
    }
    if require_nonzero_D:
        checks.update(
            {
                "D_nonzero": values["D"] != 0.0,
                "D_Z_nonzero": values["D_Z"] != 0.0,
                "all_primary_derivative_ports_nonzero": all(
                    values[name] != 0.0 for name in PORT_NAMES if name.endswith("_Z")
                ),
                "independent_Ho_Hf_fields": values["Ho"] != values["Hf"],
                "independent_qyo_qyf_fields": values["qyo"] != values["qyf"],
            }
        )
    else:
        expected_field_difference = values["Rf"] ** 2 * (values["Ho"] - values["Hf"])
        checks.update(
            {
                "D_exactly_zero": values["D"] == 0.0,
                "D_zero_keeps_independent_Ho_minus_Hf": _close(
                    actual.value, expected_field_difference
                ),
                "D_zero_field_difference_nonzero": expected_field_difference != 0.0,
            }
        )
    return checks


def run_fixture(build: Callable[[Mapping[str, Any]], Mapping[str, Any]]) -> dict[str, Any]:
    """Run the independent synthetic dual check against ``build`` once.

    ``build`` is supplied by the caller so this fixture has no import edge to
    the implementation module.  Only opaque token mappings are passed into
    it; the synthetic doubles are consumed solely by this fixture's parser.
    """

    checks: dict[str, bool] = {
        "builder_callable": callable(build),
        "synthetic_inputs_are_not_source_values": True,
        "primary_scenario_checked": False,
        "D_zero_scenario_checked": False,
    }
    if not callable(build):
        return {
            "schema": SCHEMA,
            "status": "rejected",
            "passed": False,
            "checks": checks,
            "source_values_selected": False,
        }

    try:
        built = build(_synthetic_ports())
        checks["builder_result_mapping"] = isinstance(built, Mapping)
        if not isinstance(built, Mapping):
            raise ValueError("builder did not return a mapping")
        checks["builder_declares_source_values_unselected"] = built.get(
            "source_values_selected"
        ) is False
        checks["builder_declares_no_graph_access"] = (
            built.get("source_graph_accessed") is False
            and built.get("source_graph_traversed") is False
        )

        delta_node = _output_node(built, "delta_K")
        delta_z_node = _output_node(built, "delta_K_Z")
        checks.update(
            {
                "primary_scenario_checked": True,
                **{
                    "primary_" + name: passed
                    for name, passed in _scenario_checks(
                        delta_node,
                        delta_z_node,
                        _PRIMARY,
                        require_nonzero_D=True,
                    ).items()
                },
            }
        )
        checks.update(
            {
                "D_zero_scenario_checked": True,
                **{
                    "D_zero_" + name: passed
                    for name, passed in _scenario_checks(
                        delta_node,
                        delta_z_node,
                        _ZERO_D,
                        require_nonzero_D=False,
                    ).items()
                },
            }
        )
    except (ArithmeticError, KeyError, TypeError, ValueError):
        return {
            "schema": SCHEMA,
            "status": "synthetic_comparison_failed",
            "passed": False,
            "checks": checks,
            "source_values_selected": False,
            "source_graph_accessed": False,
            "synthetic_only": True,
        }

    passed = all(checks.values())
    return {
        "schema": SCHEMA,
        "status": "synthetic_comparison_passed" if passed else "synthetic_comparison_failed",
        "passed": passed,
        "checks": checks,
        "synthetic_only": True,
        "source_values_selected": False,
        "source_graph_accessed": False,
        "source_graph_traversed": False,
        "reasonable_float_tolerance": "rel_tol=4e-12, abs_tol=4e-12",
        "limitations": [
            "Synthetic doubles are fixture-local and are not source values.",
            "The expected result uses direct exp dual arithmetic, independently of the DAG factorization.",
            "A pass does not prove owner identity, coordinate semantics, units, error certificates, or seam cancellation.",
        ],
    }


__all__ = ["PORT_NAMES", "SCHEMA", "run_fixture"]
