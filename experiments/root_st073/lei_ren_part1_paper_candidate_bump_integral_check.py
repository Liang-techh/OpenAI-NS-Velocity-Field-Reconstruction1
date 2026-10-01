"""Independent physical five-bump moment integration.

The production five-bump map stores precomputed coefficient integrals.  This
module deliberately evaluates the compactly supported fields afresh with an
independent MP Gauss--Legendre rule.  It integrates the *polarized physical
changes* directly, so a large reference moment is never subtracted from a
nearby corrected moment.

For a base angular field ``F`` and axial field ``Uz`` in physical units, the
controls are ordered ``(c1, c2, xi1, xi2, xi3)`` and

``f = xi1*gamma1 + xi2*gamma2 + xi3*gamma3``

``g = c1*gamma1 + c2*gamma3``.

The angular velocity perturbation is
``a = Am*f/sqrt(2*Rm*x)`` and the axial perturbation is ``g``.  The five
physical changes are expanded before integration:

``dtheta   = 2*Rm**2*int(x*a)``
``dz       = Rm*int(g)``
``dtheta_z = 2*Rm**2*int(x*(F*g + Uz*a + a*g))``
``dztheta  = Rm*int(2*Uz*g + g*g) - Rm**2*int(x*(2*F*a + a*a))``
``dp       = Rm*int(2*F*a + a*a)``.

This is the physical form of Section 10.2, equation (10.8), with all
angular--axial and quadratic terms retained.  The default reference fields
are the paper power-law fields ``F=Am*x**(-2/5)/sqrt(2*Rm)`` and
``Uz=4*Z``.  Callers may provide ``reference_fields`` (a mapping, callable,
or object exposing ``evaluate_x``) to use an actual source profile.

The public ``integrate_rows`` wrapper has the short signature used by the
candidate inverse.  ``AxialDual`` inputs are propagated through the same
quadrature, producing analytic first ``Z`` tangents.  This helper does not
solve the inverse problem, fit a profile, or provide a quadrature enclosure.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import hashlib
import json
import operator
from pathlib import Path
from typing import Any

import mpmath as mp

try:
    from .lei_ren_part1_paper_axial_dual import AxialDual
except (ImportError, ValueError):
    from lei_ren_part1_paper_axial_dual import AxialDual


COEFFICIENT_NAMES = ("c1", "c2", "xi1", "xi2", "xi3")
PHYSICAL_MOMENT_NAMES = ("theta", "z", "theta_z", "z_theta", "p")
PAPER_ROW_NAMES = ("z", "z_weighted", "theta", "z_theta", "p")


def _index(value: Any, name: str) -> int:
    try:
        result = operator.index(value)
    except TypeError as exc:
        raise TypeError(f"{name} must be an integer") from exc
    if result < 1:
        raise ValueError(f"{name} must be positive")
    return int(result)


def _mp(value: Any) -> mp.mpf:
    if isinstance(value, mp.mpf):
        return value
    return mp.mpf(str(value))


def _orders(*values: Any) -> tuple[int, int] | None:
    for value in values:
        if isinstance(value, AxialDual):
            return value.orders
        if isinstance(value, Mapping):
            found = _orders(*value.values())
            if found is not None:
                return found
        if isinstance(value, (tuple, list)):
            found = _orders(*value)
            if found is not None:
                return found
    return None


def _dual(value: Any, tangent: Any = 0, orders: tuple[int, int] | None = None) -> AxialDual:
    if isinstance(value, AxialDual):
        if tangent == 0:
            return value
        return AxialDual(
            value.value,
            tangent,
            pressure_order=value.pressure_order,
            width_order=value.width_order,
        )
    kwargs = {} if orders is None else {
        "pressure_order": orders[0],
        "width_order": orders[1],
    }
    return AxialDual(value, tangent, **kwargs)


def _zero_like(value: Any) -> Any:
    try:
        return value * 0
    except Exception:
        return mp.mpf(0)


def _sqrt(value: Any) -> Any:
    method = getattr(value, "sqrt", None)
    if method is not None:
        return method()
    return mp.sqrt(value)


def _sum(values: Sequence[Any], template: Any | None = None) -> Any:
    if template is None:
        if not values:
            return mp.mpf(0)
        result = values[0]
    else:
        result = _zero_like(template)
    for value in values:
        result = result + value
    return result


def _nominal(value: Any) -> Any:
    """Evaluate nested axial/pressure jets only for scalar-only source APIs."""

    if isinstance(value, AxialDual):
        return _nominal(value.value)
    evaluator = getattr(value, "evaluate", None)
    if evaluator is not None:
        return evaluator(pressure=1, width=1)
    return value


def _vector(value: Any) -> tuple[Any, ...]:
    try:
        result = tuple(value)
    except TypeError as exc:
        raise TypeError("controls must be a five-entry sequence") from exc
    if len(result) != 5:
        raise ValueError("controls must contain exactly five entries")
    return result


def _call_with_z(function: Any, x: mp.mpf, z: AxialDual) -> Any:
    """Call a source field while preserving its axial dual when supported."""

    if not callable(function):
        return function
    try:
        return function(x, z)
    except TypeError as first_error:
        try:
            return function(x, _nominal(z))
        except TypeError:
            try:
                return function(_nominal(z), x)
            except TypeError:
                try:
                    return function(x)
                except TypeError:
                    raise first_error


def _reference_mapping(reference_fields: Any, x: mp.mpf, z: AxialDual) -> Mapping[str, Any]:
    if reference_fields is None:
        return {}
    if hasattr(reference_fields, "evaluate_x"):
        try:
            values = reference_fields.evaluate_x(x, z)
        except TypeError:
            # Existing source adapters accept a scalar Z and construct their
            # own AxialDual internally.  Preserve that analytic path rather
            # than forcing a finite-difference source query.
            scalar_z = _nominal(z)
            values = reference_fields.evaluate_x(x, scalar_z)
    elif callable(reference_fields):
        values = _call_with_z(reference_fields, x, z)
    else:
        values = reference_fields
    if isinstance(values, Mapping):
        return values
    if isinstance(values, (tuple, list)) and len(values) >= 2:
        return {"F": values[0], "Uz": values[1]}
    raise TypeError("reference_fields must return a mapping or (F, Uz)")


def _field_from_reference(
    reference_fields: Any,
    name: str,
    x: mp.mpf,
    z: AxialDual,
    orders: tuple[int, int] | None,
    *,
    Rm: Any,
) -> AxialDual:
    values = _reference_mapping(reference_fields, x, z)
    aliases = (name, "Utheta") if name == "F" else (name,)
    selected = None
    selected_name = None
    for alias in aliases:
        if alias in values:
            selected = values[alias]
            selected_name = alias
            break
    if selected is None:
        raise KeyError(f"reference_fields is missing {name!r}")
    selected = _call_with_z(selected, x, z)
    if name == "F" and selected_name == "Utheta":
        selected = selected / (_sqrt(2 * Rm * x))
    if isinstance(selected, AxialDual):
        return selected
    tangent = values.get(f"{selected_name}_Z", values.get(f"{selected_name}Z", 0))
    tangent = _call_with_z(tangent, x, z) if callable(tangent) else tangent
    return _dual(selected, tangent, orders)


def _reference_fields(
    reference_fields: Any,
    x: mp.mpf,
    z: AxialDual,
    *,
    Rm: Any,
    Am: Any,
    orders: tuple[int, int] | None,
) -> tuple[AxialDual, AxialDual]:
    if reference_fields is None:
        F = Am * x ** (-mp.mpf(2) / 5) / (_sqrt(2 * Rm))
        Uz = 4 * z
        return _dual(F, 0, orders), _dual(Uz, 0, orders)
    values = _reference_mapping(reference_fields, x, z)
    F = _field_from_reference(reference_fields, "F", x, z, orders, Rm=Rm)
    if "Uz" not in values:
        raise KeyError("reference_fields is missing 'Uz'")
    Uz_value = _call_with_z(values["Uz"], x, z)
    if isinstance(Uz_value, AxialDual):
        Uz = Uz_value
    else:
        tangent = values.get("Uz_Z", values.get("UzZ", 0))
        tangent = _call_with_z(tangent, x, z) if callable(tangent) else tangent
        Uz = _dual(Uz_value, tangent, orders)
    return F, Uz


def _beta_values(moment_map: Any, x: mp.mpf) -> tuple[mp.mpf, ...]:
    beta = getattr(moment_map, "beta", None)
    centers = getattr(moment_map, "centers", None)
    if beta is None or centers is None or len(centers) != 3:
        raise TypeError("moment_map must expose beta(x, center) and three centers")
    return tuple(_mp(beta(x, center)) for center in centers)


def _breaks(moment_map: Any, left: mp.mpf, right: mp.mpf) -> tuple[mp.mpf, ...]:
    centers = tuple(_mp(value) for value in getattr(moment_map, "centers", ()))
    radius = getattr(moment_map, "radius", None)
    if radius is None:
        return (left, right)
    radius = _mp(radius)
    values = {left, right}
    for center in centers:
        for edge in (center - radius, center, center + radius):
            if left < edge < right:
                values.add(edge)
    return tuple(sorted(values))


def _quadrature(order: int, left: mp.mpf, right: mp.mpf, moment_map: Any):
    nodes, weights = mp.gauss_quadrature(order, "legendre")
    for lower, upper in zip(_breaks(moment_map, left, right), _breaks(moment_map, left, right)[1:]):
        midpoint = (lower + upper) / 2
        half = (upper - lower) / 2
        for node, weight in zip(nodes, weights):
            yield midpoint + half * node, half * weight


def _physical_to_rows(
    changes: Mapping[str, AxialDual],
    *,
    Rm: Any,
    Am: Any,
    Z: AxialDual,
) -> tuple[AxialDual, ...]:
    angular = mp.sqrt(2) * Rm * _sqrt(Rm) * Am
    row_z = changes["z"] / Rm
    row_theta = changes["theta"] / angular
    row_weighted = changes["theta_z"] / angular - 4 * Z * row_theta
    row_ztheta = (changes["z_theta"] - 8 * Z * changes["z"]) / (Rm * Am ** 2)
    row_pressure = changes["p"] / (Am ** 2)
    return row_z, row_weighted, row_theta, row_ztheta, row_pressure


def independent_bump_integral(
    moment_map: Any,
    controls: Sequence[Any],
    *,
    Rm: Any = 1,
    Am: Any,
    Z: Any,
    reference_fields: Any = None,
    left: Any = 1,
    right: Any = 2,
    order: int = 128,
    precision: int | None = None,
) -> dict[str, Any]:
    """Integrate physical bump changes and normalized paper rows independently.

    ``reference_fields`` may be a mapping containing ``F`` and ``Uz``
    callables, an object with ``evaluate_x(x, Z)``, or a callable returning a
    mapping.  Callables accepting the ``AxialDual`` ``Z`` receive analytic
    first derivatives.  Omitting it selects the paper's unperturbed profile.
    """

    h = _vector(controls)
    if not hasattr(moment_map, "beta"):
        raise TypeError("moment_map must expose beta(x, center)")
    order = _index(order, "order")
    left_value, right_value = _mp(left), _mp(right)
    if left_value < 1 or right_value > 2 or right_value <= left_value:
        raise ValueError("require 1 <= left < right <= 2")
    inferred = _orders(h, Rm, Am, Z)
    if isinstance(Z, AxialDual):
        z_dual = Z
    else:
        z_dual = _dual(Z, 1, inferred)
    inferred = _orders(h, Rm, Am, z_dual)
    rm_dual = _dual(Rm, 0, inferred)
    am_dual = _dual(Am, 0, inferred)
    h_dual = tuple(_dual(value, 0, inferred) for value in h)
    precision_value = int(precision or getattr(moment_map, "precision", 120))
    if precision_value < 50:
        raise ValueError("precision must be at least 50 digits")
    with mp.workdps(max(precision_value, 80)):
        # Use fresh accumulators and a fresh Gauss rule; no map matrix or
        # production quadrature weights are consumed here.
        zero = _zero_like(h_dual[0])
        totals = {name: zero for name in PHYSICAL_MOMENT_NAMES}
        for x, weight in _quadrature(order, left_value, right_value, moment_map):
            gamma = _beta_values(moment_map, x)
            f = _sum([h_dual[index + 2] * gamma[index] for index in range(3)], h_dual[0])
            g = h_dual[0] * gamma[0] + h_dual[1] * gamma[2]
            F, Uz = _reference_fields(
                reference_fields,
                x,
                z_dual,
                Rm=rm_dual,
                Am=am_dual,
                orders=inferred,
            )
            angular = am_dual * f / _sqrt(2 * rm_dual * x)
            fg = F * g
            ua = Uz * angular
            ag = angular * g
            fa = F * angular
            aa = angular * angular
            totals["theta"] = totals["theta"] + weight * (2 * rm_dual ** 2 * x * angular)
            totals["z"] = totals["z"] + weight * (rm_dual * g)
            totals["theta_z"] = totals["theta_z"] + weight * (2 * rm_dual ** 2 * x * (fg + ua + ag))
            totals["z_theta"] = totals["z_theta"] + weight * (rm_dual * (2 * Uz * g + g * g) - rm_dual ** 2 * x * (2 * fa + aa))
            totals["p"] = totals["p"] + weight * (rm_dual * (2 * fa + aa))
        rows = _physical_to_rows(totals, Rm=rm_dual, Am=am_dual, Z=z_dual)
        physical_values = {name: value.value for name, value in totals.items()}
        physical_tangents = {name: value.tangent for name, value in totals.items()}
        row_values = tuple(value.value for value in rows)
        row_tangents = tuple(value.tangent for value in rows)
        return {
            "physical_changes": totals,
            "physical_changes_value": physical_values,
            "physical_changes_Z": physical_tangents,
            "paper_rows": rows,
            "paper_rows_value": row_values,
            "paper_rows_Z": row_tangents,
            "paper_row_names": PAPER_ROW_NAMES,
            "physical_moment_names": PHYSICAL_MOMENT_NAMES,
            "controls": h_dual,
            "Rm": rm_dual,
            "Am": am_dual,
            "Z": z_dual,
            "interval": (left_value, right_value),
            "quadrature_order": order,
            "quadrature_rule": "fresh MP Gauss-Legendre split at bump centers and support edges",
            "reference_fields_supplied": reference_fields is not None,
            "polarized_cross_terms_integrated": True,
            "background_subtraction_used": False,
            "production_map_apply_used": False,
            "production_map_matrix_used": False,
            "quadrature_enclosed": False,
            "inverse_solved": False,
            "functional_closure": False,
            "cone_certified": False,
        }


def integrate_rows(
    moment_map: Any,
    h: Sequence[Any],
    Am: Any,
    Z: Any,
    order: int = 128,
    *,
    Rm: Any = 1,
    reference_fields: Any = None,
    left: Any = 1,
    right: Any = 2,
    precision: int | None = None,
) -> tuple[AxialDual, ...]:
    """Short candidate-inverse API returning the five normalized rows."""

    return independent_bump_integral(
        moment_map,
        h,
        Rm=Rm,
        Am=Am,
        Z=Z,
        reference_fields=reference_fields,
        left=left,
        right=right,
        order=order,
        precision=precision,
    )["paper_rows"]


def _scalar_rows(rows: Sequence[AxialDual]) -> tuple[mp.mpf, ...]:
    return tuple(value.value for value in rows)


def _jet_scalar(value: Any) -> Any:
    evaluator = getattr(value, "evaluate", None)
    if evaluator is not None:
        return evaluator(pressure=1, width=1)
    return value


def run_fixture() -> dict[str, Any]:
    """Independent resolved fixture and dual tangent check."""

    try:
        from .lei_ren_part1_paper_five_bump_map import FiveBumpMomentMap
    except (ImportError, ValueError):
        from lei_ren_part1_paper_five_bump_map import FiveBumpMomentMap
    with mp.workdps(120):
        moment_map = FiveBumpMomentMap(precision=110, order=32)
        h0 = tuple(mp.mpf(item) for item in (".013", "-.007", ".011", "-.019", ".005"))
        hZ = tuple(mp.mpf(item) for item in (".002", ".003", "-.004", ".001", ".006"))
        z0 = mp.mpf(".3")
        am0 = mp.mpf("1.4")
        amZ = mp.mpf("-.21")
        z_dual = AxialDual(z0, 1, pressure_order=0, width_order=0)
        am_dual = AxialDual(am0, amZ, pressure_order=0, width_order=0)
        h_dual = tuple(AxialDual(v, t, pressure_order=0, width_order=0) for v, t in zip(h0, hZ))
        result = independent_bump_integral(moment_map, h_dual, Rm=1, Am=am_dual, Z=z_dual, order=96)
        scalar = independent_bump_integral(moment_map, h0, Rm=1, Am=am0, Z=z0, order=96)
        scaled_radius = independent_bump_integral(moment_map, h_dual, Rm=3, Am=am_dual, Z=z_dual, order=96)
        step = mp.mpf("1e-10")
        plus = independent_bump_integral(moment_map, tuple(v + step * t for v, t in zip(h0, hZ)), Rm=1, Am=am0 + step * amZ, Z=z0 + step, order=96)
        minus = independent_bump_integral(moment_map, tuple(v - step * t for v, t in zip(h0, hZ)), Rm=1, Am=am0 - step * amZ, Z=z0 - step, order=96)
        fd = tuple((_jet_scalar(plus["paper_rows_value"][i]) - _jet_scalar(minus["paper_rows_value"][i])) / (2 * step) for i in range(5))
        tangent_errors = tuple(abs(_jet_scalar(result["paper_rows_Z"][i]) - fd[i]) for i in range(5))
        value_errors = tuple(abs(_jet_scalar(result["paper_rows_value"][i]) - _jet_scalar(scalar["paper_rows_value"][i])) for i in range(5))
        radius_scale_errors = tuple(
            abs(_jet_scalar(scaled_radius["paper_rows_value"][i]) - _jet_scalar(result["paper_rows_value"][i]))
            for i in range(5)
        )
        if max(value_errors) > mp.mpf("1e-100"):
            raise AssertionError("dual value disagrees with scalar integration")
        if max(tangent_errors) > mp.mpf("1e-12"):
            raise AssertionError("AxialDual tangent disagrees with independent central difference")
        if max(radius_scale_errors) > mp.mpf("1e-100"):
            raise AssertionError("normalized rows changed under the equivalent Rm = 3 chart")
        report = {
            "branch": "independent_physical_five_bump_integral",
            "precision": 110,
            "quadrature_order": 96,
            "interval": ["1", "2"],
            "paper_row_names": list(PAPER_ROW_NAMES),
            "physical_moment_names": list(PHYSICAL_MOMENT_NAMES),
            "value_max_absolute_error": mp.nstr(max(value_errors), 30),
            "dual_tangent_max_absolute_error": mp.nstr(max(tangent_errors), 30),
            "Rm_three_normalized_scale_max_absolute_error": mp.nstr(max(radius_scale_errors), 30),
            "polarized_cross_terms_integrated": True,
            "background_subtraction_used": False,
            "production_map_apply_used": False,
            "production_map_matrix_used": False,
            "reference_fields_supplied": False,
            "quadrature_enclosed": False,
            "inverse_solved": False,
            "functional_closure": False,
            "cone_certified": False,
            "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        }
        Path(__file__).with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        return report


if __name__ == "__main__":
    print(json.dumps(run_fixture(), indent=2), flush=True)


__all__ = ["independent_bump_integral", "integrate_rows", "run_fixture"]
