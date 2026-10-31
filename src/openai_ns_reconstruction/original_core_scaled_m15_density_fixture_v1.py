"""Independent exact-rational fixture for the scaled ``m15`` density leaf.

The fixture owns its interval type and backend and does not import the leaf or
any historical source producer.  ``run_fixture(leaf)`` is intentionally inert
until a caller supplies the already-imported leaf.  Its expected values are
computed from separate rational witnesses, not from a midpoint or from the
leaf's output.
"""

from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
from dataclasses import dataclass
from fractions import Fraction
from typing import Any


SCHEMA = "original_core_scaled_m15_density_fixture_v1"


@dataclass(frozen=True)
class FractionRange:
    lower: Fraction
    upper: Fraction

    def __post_init__(self) -> None:
        if type(self.lower) is not Fraction or type(self.upper) is not Fraction:
            raise TypeError("FractionRange endpoints must be Fraction values")
        if self.lower > self.upper:
            raise ValueError("FractionRange is reversed")


class FractionEndpointBackend:
    """Synthetic exact backend; no floating point or midpoint operation."""

    exact = True

    @staticmethod
    def _fraction(value: Any) -> Fraction:
        if isinstance(value, bool):
            raise TypeError("boolean is not an exact endpoint")
        if isinstance(value, float):
            raise TypeError("float is not an exact endpoint")
        if type(value) is Fraction:
            return value
        if type(value) is int:
            return Fraction(value)
        if isinstance(value, str):
            return Fraction(value)
        raise TypeError("endpoint must be Fraction, int, or rational string")

    def interval(self, lower: Any, upper: Any | None = None) -> FractionRange:
        lo = self._fraction(lower)
        hi = lo if upper is None else self._fraction(upper)
        return FractionRange(lo, hi)

    make_interval = interval

    @staticmethod
    def endpoints(value: FractionRange) -> tuple[Fraction, Fraction]:
        return value.lower, value.upper

    def add(self, left: FractionRange, right: FractionRange) -> FractionRange:
        return FractionRange(left.lower + right.lower, left.upper + right.upper)

    def mul(self, left: FractionRange, right: FractionRange) -> FractionRange:
        values = (
            left.lower * right.lower,
            left.lower * right.upper,
            left.upper * right.lower,
            left.upper * right.upper,
        )
        return FractionRange(min(values), max(values))

    def scale(self, factor: int | Fraction, value: FractionRange) -> FractionRange:
        f = self._fraction(factor)
        if f >= 0:
            return FractionRange(value.lower * f, value.upper * f)
        return FractionRange(value.upper * f, value.lower * f)


def _r(backend: FractionEndpointBackend, lower: int | Fraction, upper: int | Fraction) -> FractionRange:
    return backend.interval(lower, upper)


def _field(
    backend: FractionEndpointBackend,
    name: str,
    z0: tuple[Fraction, Fraction],
    z1: tuple[Fraction, Fraction],
    *,
    scale: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    def row(order: int, bounds: tuple[Fraction, Fraction]) -> dict[str, Any]:
        finite = {
            "interval": _r(backend, *bounds),
            "source_id": f"synthetic:{name}:Z{order}:finite",
            "term_id": f"term:{name}:Z{order}:finite",
            "source_path": f"fixture/{name}/Z{order}/finite",
            "error_class": "synthetic_finite_range",
        }
        # Each tail is independent and straddles zero, so it cannot be hidden
        # in a widened value while the witness remains in the enclosure.
        tail = {
            "interval": _r(backend, Fraction(-1, 100 + order), Fraction(1, 101 + order)),
            "tail_id": f"tail:{name}:Z{order}",
            "source_id": f"synthetic:{name}:Z{order}:tail",
            "term_id": f"term:{name}:Z{order}:tail",
            "source_path": f"fixture/{name}/Z{order}/tail",
            "error_class": "synthetic_additive_tail",
            "uncertainty": True,
        }
        if scale is not None:
            finite["original_scale"] = dict(scale)
            tail["original_scale"] = dict(scale)
            finite["formal_powers"] = [scale.get("epsilon_power", 0), scale.get("F0_power", 0)]
            tail["formal_powers"] = [scale.get("epsilon_power", 0), scale.get("F0_power", 0)]
        return {"finite": finite, "tails": [tail]}

    return {"Z0": row(0, z0), "Z1": row(1, z1), "source_id": f"synthetic:field:{name}"}


def _packet(backend: FractionEndpointBackend) -> dict[str, Any]:
    return {
        "rho": {
            "interval": _r(backend, Fraction(12, 5), Fraction(13, 5)),
            "source_id": "synthetic:rho:fixed-cell",
            "term_id": "term:rho:fixed-cell",
        },
        "rho_Z": 0,
        "epsilon_Z": 0,
        "responses_include_eta": True,
        "eta_Z": 0,
        "include_pressure": False,
        "fields": {
            "U0": _field(
                backend,
                "U0",
                (Fraction(19, 10), Fraction(21, 10)),
                (Fraction(2, 5), Fraction(3, 5)),
            ),
            "Delta": _field(
                backend,
                "Delta",
                (Fraction(9, 10), Fraction(11, 10)),
                (Fraction(-3, 10), Fraction(-1, 5)),
            ),
            "Phi": _field(
                backend,
                "Phi",
                (Fraction(7, 5), Fraction(8, 5)),
                (Fraction(3, 8), Fraction(11, 25)),
            ),
            "F1": _field(
                backend,
                "F1",
                (Fraction(3, 5), Fraction(11, 15)),
                (Fraction(-3, 8), Fraction(-1, 4)),
                scale={"epsilon_power": 0, "F0_power": 1, "eta_included": True, "scale_origin": "synthetic original F1"},
            ),
            "Uz1": _field(
                backend,
                "Uz1",
                (Fraction(5, 4), Fraction(17, 12)),
                (Fraction(1, 6), Fraction(4, 15)),
                scale={"epsilon_power": 0, "F0_power": 0, "eta_included": True, "scale_origin": "synthetic original Uz1"},
            ),
            "L": _field(
                backend,
                "L",
                (Fraction(3, 10), Fraction(2, 5)),
                (Fraction(1, 10), Fraction(3, 20)),
            ),
        },
        "certificate": {
            "whole_cell_verified": False,
            "owner_verified": False,
            "binder_verified": False,
            "scope": "synthetic fixture only",
        },
    }


def _bounds(value: Any) -> tuple[Fraction, Fraction]:
    if isinstance(value, Mapping):
        lower = value.get("lower", value.get("lo", value.get("left")))
        upper = value.get("upper", value.get("hi", value.get("right")))
    else:
        lower, upper = value.lower, value.upper
    return lower, upper


def _contains(terms: list[Mapping[str, Any]], expected: Fraction) -> bool:
    lower = Fraction(0)
    upper = Fraction(0)
    for term in terms:
        lo, hi = _bounds(term["coefficient"])
        lower += lo
        upper += hi
    return lower <= expected <= upper


def _finite_sign(terms: list[Mapping[str, Any]], *, sign: int, piece: str | None = None) -> bool:
    selected = [term for term in terms if not term.get("has_additive_tail") and (piece is None or term.get("piece") == piece)]
    if not selected:
        return False
    for term in selected:
        lower, upper = _bounds(term["coefficient"])
        if sign > 0 and lower <= 0:
            return False
        if sign < 0 and upper >= 0:
            return False
    return True


def _call(leaf: Any, packet: Mapping[str, Any], backend: FractionEndpointBackend) -> Mapping[str, Any]:
    function = getattr(leaf, "density", leaf)
    result = function(packet, backend=backend)
    if not isinstance(result, Mapping):
        raise TypeError("density leaf did not return a mapping")
    return result


def _manual_expected() -> dict[str, Fraction]:
    # Independent rational witnesses, chosen inside the finite intervals above.
    rho = Fraction(5, 2)
    u0, delta, uz1 = Fraction(2), Fraction(1), Fraction(4, 3)
    u0_z, delta_z, uz1_z = Fraction(1, 2), Fraction(-1, 4), Fraction(1, 5)
    phi, f1, log_derivative = Fraction(3, 2), Fraction(2, 3), Fraction(1, 3)
    phi_z, f1_z = Fraction(2, 5), Fraction(-1, 3)
    return {
        "axial_Z0": 2 * (u0 + delta) * uz1,
        "axial_Z1": 2 * ((u0_z + delta_z) * uz1 + (u0 + delta) * uz1_z),
        "swirl_Z0": -2 * rho * phi * f1,
        "swirl_Z1": -2 * rho * (log_derivative * phi * f1 + phi_z * f1 + phi * f1_z),
    }


def run_fixture(leaf: Any) -> dict[str, Any]:
    """Run independent signed-range and negative-guard checks on *leaf*."""

    backend = FractionEndpointBackend()
    packet = _packet(backend)
    checks: dict[str, bool] = {}
    try:
        result = _call(leaf, packet, backend)
    except Exception as exc:
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": {"build_completed": False},
            "check_count": 1,
            "error": f"{type(exc).__name__}: {exc}",
            "synthetic_only": True,
            "implementation_imported": False,
        }

    checks["build_completed"] = isinstance(result, Mapping)
    checks["schema"] = result.get("schema") == "original_core_scaled_m15_density_v1"
    checks["one_fixed_cell"] = result.get("cell_count") == 1 and result.get("cells", [{}])[0].get("rho") is not None
    checks["ordinary_Z_rows"] = result.get("ordinary_Z_orders") == [0, 1] and result.get("ordinary_Z_labels") == ["Z0", "Z1"]
    checks["formal_power_basis"] = result.get("formal_power_basis") == ["epsilon", "F0"]
    checks["velocity_channels_only"] = (
        result.get("pressure", {}).get("included") is False
        and all(term.get("channel") in {"axial", "swirl"} for term in result.get("terms", ()))
        and result.get("W1_axial_swirl_only_by_default") is True
    )
    checks["fixed_zero_premises"] = result.get("rho_Z_is_zero") is True and result.get("epsilon_Z_is_zero") is True and result.get("eta_Z_is_zero_when_eta_included") is True
    checks["eta_marked"] = result.get("responses_include_eta") is True and result.get("eta_status", {}).get("included") is True
    checks["no_materialized_scales"] = result.get("epsilon_materialized") is False and result.get("F0_materialized") is False and result.get("epsilon_divided") is False
    checks["no_integration_or_endpoint_sum"] = result.get("integration_performed") is False and result.get("endpoint_sum_performed") is False
    checks["not_complete_m15"] = result.get("complete_m15_computed") is False and result.get("full_m15_claimed") is False

    axial = result.get("axial", {})
    swirl = result.get("swirl", {})
    checks["manual_Z0_Z1_enclosures"] = all(
        _contains(axial.get(label, []), expected) and _contains(swirl.get(label, []), expected_swirl)
        for label, expected, expected_swirl in (
            ("Z0", _manual_expected()["axial_Z0"], _manual_expected()["swirl_Z0"]),
            ("Z1", _manual_expected()["axial_Z1"], _manual_expected()["swirl_Z1"]),
        )
    )
    checks["signed_Z0_terms"] = _finite_sign(axial.get("Z0", []), sign=1) and _finite_sign(swirl.get("Z0", []), sign=-1)
    checks["signed_swirl_Z1_pieces"] = (
        _finite_sign(swirl.get("Z1", []), sign=-1, piece="rho*L_0*Phi_0*F1_0")
        and _finite_sign(swirl.get("Z1", []), sign=-1, piece="rho*Phi_Z*F1_0")
        and _finite_sign(swirl.get("Z1", []), sign=1, piece="rho*Phi_0*F1_Z")
    )
    all_terms = list(result.get("terms", ()))
    checks["all_terms_have_scale_tags"] = all(
        term.get("formal_powers") in ({"epsilon": 1, "F0": 0}, {"epsilon": 2, "F0": 2})
        and term.get("ordinary_Z_order") in (0, 1)
        for term in all_terms
    )
    # A factored F1=F0*coefficient is a synthetic test of power composition;
    # for an unscaled F1, the only F0 factor is the leading-field factor.
    unscaled = deepcopy(packet)
    for order in ("Z0", "Z1"):
        for component in [unscaled["fields"]["F1"][order]["finite"]] + unscaled["fields"]["F1"][order]["tails"]:
            component["formal_powers"] = [0, 0]
            component["original_scale"]["F0_power"] = 0
    unscaled_result = _call(leaf, unscaled, backend)
    checks["unscaled_F1_has_one_leading_F0"] = all(term["formal_powers"]["F0"] == 1 for term in unscaled_result["swirl"]["Z0"] + unscaled_result["swirl"]["Z1"])
    checks["tails_remain_distributive"] = (
        result.get("tails_distributed") is True
        and result.get("aggregation_performed") is False
        and result.get("tail_term_count", 0) > 0
        and all(term.get("tail_components") for term in all_terms if term.get("has_additive_tail"))
    )
    checks["provenance_and_scale_metadata"] = all(
        term.get("source_ids") and term.get("term_ids")
        for term in all_terms
    ) and any(
        term.get("original_scale_metadata", {}).get("F1")
        or term.get("original_scale_metadata", {}).get("Uz1")
        for term in all_terms
    )
    checks["unverified_certificates"] = all(
        result.get(name, {}).get("verified") is False
        for name in ("whole_cell_certificate", "owner_certificate", "binder_certificate")
    )
    checks["rho_provenance_on_swirl"] = all("synthetic:rho:fixed-cell" in term.get("source_ids", []) for term in all_terms if term.get("channel") == "swirl")
    claimed = deepcopy(packet)
    claimed["certificate"]["whole_cell_verified"] = True
    claimed_result = _call(leaf, claimed, backend)
    checks["caller_coverage_claim_not_promoted"] = claimed_result["whole_cell_certificate"]["verified"] is False

    pressure_packet = deepcopy(packet)
    pressure_packet["include_pressure"] = True
    pressure_packet["N1p"] = _field(backend, "N1p", (Fraction(2), Fraction(2)), (Fraction(3), Fraction(3)), scale={"epsilon_power": -1, "F0_power": 0})
    pressure_result = _call(leaf, pressure_packet, backend)
    checks["pressure_kept_separate_and_scaled"] = (
        pressure_result.get("pressure", {}).get("included") is True
        and _contains(pressure_result["pressure"]["Z0"], Fraction(-5))
        and _contains(pressure_result["pressure"]["Z1"], Fraction(-15, 2))
        and all(term["formal_powers"] == {"epsilon": 1, "F0": 0} for term in pressure_result["pressure"]["Z0"] + pressure_result["pressure"]["Z1"])
    )
    missing_pressure_power = deepcopy(pressure_packet)
    for order in ("Z0", "Z1"):
        del missing_pressure_power["N1p"][order]["finite"]["formal_powers"]
    try:
        _call(leaf, missing_pressure_power, backend)
    except Exception:
        checks["missing_pressure_scale_rejected"] = True
    else:
        checks["missing_pressure_scale_rejected"] = False

    # Missing derivative and nonzero fixed-geometry premises must fail closed.
    missing = deepcopy(packet)
    del missing["fields"]["Uz1"]["Z1"]
    try:
        _call(leaf, missing, backend)
    except Exception:
        checks["missing_derivative_rejected"] = True
    else:
        checks["missing_derivative_rejected"] = False
    missing_tail = deepcopy(packet)
    del missing_tail["fields"]["Phi"]["Z0"]["tails"]
    try:
        _call(leaf, missing_tail, backend)
    except Exception:
        checks["missing_tail_declaration_rejected"] = True
    else:
        checks["missing_tail_declaration_rejected"] = False
    bad_epsilon = deepcopy(packet)
    bad_epsilon["epsilon_Z"] = _r(backend, Fraction(1, 10), Fraction(1, 10))
    try:
        _call(leaf, bad_epsilon, backend)
    except Exception:
        checks["nonzero_epsilon_Z_rejected"] = True
    else:
        checks["nonzero_epsilon_Z_rejected"] = False
    bad_rho = deepcopy(packet)
    bad_rho["rho_Z"] = Fraction(1, 10)
    try:
        _call(leaf, bad_rho, backend)
    except Exception:
        checks["nonzero_rho_Z_rejected"] = True
    else:
        checks["nonzero_rho_Z_rejected"] = False
    eta_free = deepcopy(packet)
    eta_free["responses_include_eta"] = False
    try:
        _call(leaf, eta_free, backend)
    except Exception:
        checks["eta_free_response_rejected"] = True
    else:
        checks["eta_free_response_rejected"] = False
    bad_point = deepcopy(packet)
    bad_point["fields"]["U0"]["Z0"]["midpoint"] = Fraction(2)
    try:
        _call(leaf, bad_point, backend)
    except Exception:
        checks["point_midpoint_claim_rejected"] = True
    else:
        checks["point_midpoint_claim_rejected"] = False

    return {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "implementation_imported": False,
        "old_code_called": False,
        "manual_expected": {key: str(value) for key, value in _manual_expected().items()},
        "scope": "one fixed-rho cell; signed axial/swirl W1 density terms only",
    }


run = run_fixture


__all__ = ["SCHEMA", "FractionRange", "FractionEndpointBackend", "run_fixture", "run"]
