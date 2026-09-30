"""Physical five-bump correction in the similarity coordinate ``x=R/Rm``.

This adapter consumes coefficients supplied by another component and applies
the finite Section 10.2 five-bump map on ``1 <= x <= e``.  It deliberately
does not solve for the coefficients.  The coefficient provider, the
reference moments, and the axis pressure datum are caller-owned inputs.

The map rows ``a1,...,a5`` over ``[1,min(x,2)]`` are converted to physical
increments by the exact scales

``(Rm*a1, sqrt(2)*Rm**(3/2)*Am*a3,
  sqrt(2)*Rm**(3/2)*Am*(a2+4*Z*a3),
  Rm*Am**2*a4+8*Z*Rm*a1, Am**2*a5)``.

Values and first ``Z`` derivatives are formed in one :class:`AxialDual`, so
the derivative cannot silently use a separately rounded formula.  The
returned moment and pressure increments remain separate from the inherited
reference fields.  This is a finite coefficient adapter; no coefficient
solve, cone estimate, functional closure, or global field installation is
claimed.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import operator
from typing import Any

import mpmath as mp

try:
    from .lei_ren_part1_paper_axial_dual import AxialDual
    from .lei_ren_part1_paper_five_bump_map import FiveBumpMomentMap
    from .lei_ren_part1_paper_mp_stress import evaluate_mp_stress
    from .lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
except (ImportError, ValueError):
    from lei_ren_part1_paper_axial_dual import AxialDual
    from lei_ren_part1_paper_five_bump_map import FiveBumpMomentMap
    from lei_ren_part1_paper_mp_stress import evaluate_mp_stress
    from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet


MOMENT_KEYS = ("theta", "z", "theta_z", "z_theta", "p")
COEFFICIENT_NAMES = ("c1", "c2", "xi1", "xi2", "xi3")
DOMAIN_LEFT = mp.mpf(1)
DOMAIN_RIGHT = "e"


def _mp(value: Any) -> mp.mpf:
    return value if isinstance(value, mp.mpf) else mp.mpf(str(value))


def _index(value: Any, name: str) -> int:
    try:
        return operator.index(value)
    except TypeError as exc:
        raise TypeError(f"{name} must be an integer") from exc


def _sqrt(value: Any) -> Any:
    method = getattr(value, "sqrt", None)
    if method is not None:
        return method()
    return mp.sqrt(value)


def _pow_three_halves(value: Any) -> Any:
    return value * _sqrt(value)


def _zero_like(value: Any) -> Any:
    try:
        return value * 0
    except Exception:
        return mp.mpf(0)


def _is_zero(value: Any) -> bool:
    if isinstance(value, AxialDual):
        return _is_zero(value.tangent)
    if isinstance(value, PressureWidthJet):
        return not value.atoms
    try:
        return bool(value == 0)
    except Exception:
        return False


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


def _vector(value: Any, name: str, length: int = 5) -> tuple[Any, ...]:
    if hasattr(value, "rows") and hasattr(value, "cols"):
        if int(value.rows) * int(value.cols) != length:
            raise ValueError(f"{name} must have exactly {length} entries")
        return tuple(value[index] for index in range(length))
    try:
        result = tuple(value)
    except TypeError as exc:
        raise TypeError(f"{name} must be an iterable of {length} entries") from exc
    if len(result) != length:
        raise ValueError(f"{name} must have exactly {length} entries")
    return result


def _orders(*values: Any) -> tuple[int, int] | None:
    for value in values:
        if isinstance(value, AxialDual):
            return value.orders
        if isinstance(value, PressureWidthJet):
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
    kwargs = {}
    if orders is not None:
        kwargs = {"pressure_order": orders[0], "width_order": orders[1]}
    return AxialDual(value, tangent, **kwargs)


def _public(value: Any) -> tuple[Any, Any]:
    if isinstance(value, AxialDual):
        return value.value, value.tangent
    return value, _zero_like(value)


def _public_dict(values: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    slots = [_public(value) for value in values.values()]
    return (
        {key: slots[index][0] for index, key in enumerate(values)},
        {key: slots[index][1] for index, key in enumerate(values)},
    )


def _field_tangent(source: Mapping[str, Any], name: str) -> Any:
    for key in (f"{name}_Z", f"{name}Z"):
        if key in source:
            return source[key]
    return 0


def _field_dual(source: Mapping[str, Any], name: str, orders: tuple[int, int] | None) -> AxialDual:
    if name not in source:
        raise KeyError(f"reference_data is missing {name!r}")
    value = source[name]
    if isinstance(value, AxialDual):
        return value
    return _dual(value, _field_tangent(source, name), orders)


def _moment_duals(source: Mapping[str, Any], orders: tuple[int, int] | None) -> dict[str, AxialDual]:
    if "moments" not in source:
        raise KeyError("reference_data is missing 'moments'")
    values = source["moments"]
    tangents = source.get("moments_Z", source.get("momentsZ", {}))
    return {
        key: _dual(values[key], tangents.get(key, 0), orders)
        for key in MOMENT_KEYS
    }


def _raw_duals(source: Mapping[str, Any], orders: tuple[int, int] | None) -> dict[str, AxialDual] | None:
    values = source.get("raw_quadratic_integrals")
    if not isinstance(values, Mapping) or not all(key in values for key in ("axial", "swirl")):
        return None
    tangents = source.get("raw_quadratic_integrals_Z", {})
    return {
        key: _dual(values[key], tangents.get(key, values.get(f"{key}_Z", 0)), orders)
        for key in ("axial", "swirl")
    }


def _constant(value: Any) -> mp.mpf:
    if isinstance(value, AxialDual):
        value = value.value
    if isinstance(value, PressureWidthJet):
        return value.constant
    return _mp(value)


def _gamma_values(moment_map: Any, x: mp.mpf) -> tuple[mp.mpf, ...]:
    centers = getattr(moment_map, "centers", None)
    if centers is None:
        centers = tuple(mp.mpf(value) for value in ("1.25", "1.5", "1.75"))
    if hasattr(moment_map, "beta"):
        return tuple(moment_map.beta(x, center) for center in centers)
    if hasattr(moment_map, "bump"):
        return tuple(moment_map.bump(2 + index, x) for index in range(3))
    raise TypeError("moment_map must expose beta(x, center) or bump(index, x)")


def _gamma_derivatives(moment_map: Any, x: mp.mpf) -> tuple[mp.mpf, ...]:
    if hasattr(moment_map, "bump"):
        try:
            return tuple(
                moment_map.bump(2 + index, x, derivative=1)
                for index in range(3)
            )
        except (TypeError, AttributeError):
            pass
    centers = getattr(moment_map, "centers", None)
    radius = getattr(moment_map, "radius", None)
    normalization = getattr(moment_map, "beta_normalization", None)
    if centers is None or radius is None or normalization is None:
        raise TypeError("moment_map must expose bump derivatives or beta geometry")
    radius = _mp(radius)
    normalization = _mp(normalization)
    derivatives = []
    for center in centers:
        coordinate = (x - _mp(center)) / radius
        if abs(coordinate) >= 1:
            derivatives.append(mp.mpf(0))
            continue
        raw = mp.exp(-1 / (1 - coordinate * coordinate))
        raw_derivative = raw * (-2 * coordinate) / (1 - coordinate * coordinate) ** 2
        derivatives.append(raw_derivative / (radius * radius * normalization))
    return tuple(derivatives)


def _dual_map_rows(moment_map: Any, h: Sequence[Any], Am: Any, x: mp.mpf) -> tuple[AxialDual, ...]:
    rows = moment_map.partial(h, Am, 1, min(x, mp.mpf(2)))
    return tuple(_dual(row) for row in rows)


def _bilinear(left: Sequence[Any], matrix: Any, right: Sequence[Any]) -> Any:
    template = left[0] if left else (right[0] if right else mp.mpf(0))
    return _sum(
        [
            left[row] * matrix[row, column] * right[column]
            for row in range(matrix.rows)
            for column in range(matrix.cols)
        ],
        template,
    )


def _partial_raw_primitives(
    moment_map: Any,
    h: Sequence[Any],
    x: mp.mpf,
) -> tuple[Any, Any, Any] | None:
    """Return ``(integral(g^2), integral(x^.1*f), integral(f^2))``.

    The five-bump map owns the fixed-support weights used here.  Keeping the
    primitives separate preserves the raw axial and swirl receipts instead of
    reconstructing them from a summed moment row.
    """

    if not hasattr(moment_map, "_weights_for_interval"):
        return None
    linear, weights = moment_map._weights_for_interval(1, min(x, mp.mpf(2)))
    c = h[:2]
    xi = h[2:]
    g2 = _bilinear(c, weights["gg"], c)
    f2 = _bilinear(xi, weights["ff"], xi)
    # The linear row four is -integral(x^.1 f).
    xf = _sum(
        [-xi[index] * linear[3, 2 + index] for index in range(3)],
        xi[0],
    )
    return g2, xf, f2


def _correction_duals(
    moment_map: Any,
    h: Any,
    Am: Any,
    Rm: Any,
    Z: Any,
    x: Any,
    *,
    h_Z: Any = None,
    Am_Z: Any = None,
    Rm_Z: Any = None,
    delta: Any = 0,
) -> dict[str, Any]:
    coefficients = _vector(h, "h")
    coefficient_tangents = None if h_Z is None else _vector(h_Z, "h_Z")
    inferred_orders = _orders(coefficients, Am, Rm, Z, coefficient_tangents, Am_Z, Rm_Z)
    h_dual = tuple(
        _dual(
            value,
            0 if coefficient_tangents is None else coefficient_tangents[index],
            inferred_orders,
        )
        for index, value in enumerate(coefficients)
    )
    Am_dual = _dual(Am, 0 if Am_Z is None else Am_Z, inferred_orders)
    Rm_dual = _dual(Rm, 0 if Rm_Z is None else Rm_Z, inferred_orders)
    # Rm is the fixed physical similarity radius in the paper's correction
    # chart.  A nonzero Rm_Z would differentiate at a moving radius and make
    # the supplied moment/radial identity a different problem.
    if Rm_Z is not None and not _is_zero(Rm_Z):
        raise ValueError("Rm_Z must be zero: the physical Rm chart is fixed")
    if isinstance(Rm, AxialDual) and not _is_zero(Rm.tangent):
        raise ValueError("Rm AxialDual tangent must be zero: the physical Rm chart is fixed")
    Z_dual = _dual(Z, 1, inferred_orders)
    x_value = _mp(x)
    if not DOMAIN_LEFT <= x_value <= mp.exp(1):
        raise ValueError("require 1 <= x <= exp(1)")
    rows = _dual_map_rows(moment_map, h_dual, Am_dual, x_value)
    a1, a2, a3, a4, a5 = rows
    raw_primitives = _partial_raw_primitives(moment_map, h_dual, x_value)

    scale_theta = mp.sqrt(2) * _pow_three_halves(Rm_dual) * Am_dual
    scale_ztheta = Rm_dual * Am_dual**2
    increments = {
        "z": Rm_dual * a1,
        "theta": scale_theta * a3,
        "theta_z": scale_theta * (a2 + 4 * Z_dual * a3),
        "z_theta": scale_ztheta * a4 + 8 * Z_dual * Rm_dual * a1,
        "p": Am_dual**2 * a5,
    }
    normalized = {
        "z": a1,
        "theta": a3,
        "theta_z": a2 + 4 * Z_dual * a3,
        "z_theta": increments["z_theta"] / scale_ztheta,
        "p": a5,
    }
    if raw_primitives is None:
        raw_increments = None
    else:
        g2, xf, f2 = raw_primitives
        # Repository convention: swirl = integral(R*F^2) dR = integral(u_theta^2/2) dR.
        raw_increments = {
            "axial": Rm_dual * (8 * Z_dual * a1 + g2),
            "swirl": Rm_dual * Am_dual**2 * (xf + f2 / 2),
        }

    gammas = _gamma_values(moment_map, x_value)
    gamma_x = _gamma_derivatives(moment_map, x_value)
    f = _sum(
        [h_dual[index + 2] * gammas[index] for index in range(3)],
        h_dual[0],
    )
    f_x = _sum(
        [h_dual[index + 2] * gamma_x[index] for index in range(3)],
        h_dual[0],
    )
    g = h_dual[0] * gammas[0] + h_dual[1] * gammas[2]
    g_x = h_dual[0] * gamma_x[0] + h_dual[1] * gamma_x[2]
    x_power = x_value ** mp.mpf(".1")
    u = Am_dual * (x_power + f)
    # y = log R = log Rm + log x, so x*d(x**.1)/dx = .1*x**.1.
    u_y = Am_dual * (mp.mpf(".1") * x_power + x_value * f_x)
    Uz = 4 * Z_dual + g
    Uz_y = x_value * g_x
    R = Rm_dual * x_value
    F = u / (mp.sqrt(2) * _sqrt(R))
    F_y = u_y / (mp.sqrt(2) * _sqrt(R)) - F / 2
    F_R = F_y / R
    delta_value = _mp(delta)
    L = 1 - delta_value * Z_dual * Z_dual
    sqrt_2R = mp.sqrt(2) * _sqrt(R)
    delta_Ur = (
        2 * Z_dual * R * g
        - (1 - delta_value) * Z_dual * increments["z"]
        - (1 - Z_dual * Z_dual) * increments["z"].tangent
    ) / (L * sqrt_2R)
    return {
        "rows": rows,
        "centered_rows": rows,
        "fields": {
            "R": R,
            "F": F,
            "Utheta": u,
            "Utheta_y": u_y,
            "Uz": Uz,
            "Uz_y": Uz_y,
            "F_y": F_y,
            "F_R": F_R,
            "f": f,
            "f_x": f_x,
            "g": g,
            "g_x": g_x,
        },
        "increments": increments,
        "normalized": normalized,
        "raw_increments": raw_increments,
        "delta_Ur": delta_Ur,
        "parameters": {
            "h": h_dual,
            "Am": Am_dual,
            "Rm": Rm_dual,
            "Z": Z_dual,
            "x": x_value,
            "delta": delta_value,
        },
    }


def _dual_public_map(values: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    return _public_dict(values)


def evaluate_correction(
    moment_map: Any,
    h: Any,
    Am: Any,
    Rm: Any,
    Z: Any,
    x: Any,
    delta: Any = 0,
    h_Z: Any = None,
    Am_Z: Any = None,
    Rm_Z: Any = None,
) -> dict[str, Any]:
    """Evaluate the standalone five-bump correction at ``x=R/Rm``.

    ``h_Z``, ``Am_Z``, and ``Rm_Z`` are optional when the corresponding
    values are already :class:`AxialDual` instances.  A plain ``Z`` is given
    tangent one, so direct ``Z`` dependence is differentiated in the same
    dual calculation.
    """

    precision = int(getattr(moment_map, "precision", 100))
    with mp.workdps(max(50, precision)):
        data = _correction_duals(
            moment_map,
            h,
            Am,
            Rm,
            Z,
            x,
            h_Z=h_Z,
            Am_Z=Am_Z,
            Rm_Z=Rm_Z,
            delta=delta,
        )
        rows, rows_Z = _dual_public_map(
            {f"a{index + 1}": value for index, value in enumerate(data["rows"])}
        )
        increments, increments_Z = _dual_public_map(data["increments"])
        normalized, normalized_Z = _dual_public_map(data["normalized"])
        fields, fields_Z = _dual_public_map(data["fields"])
        if data["raw_increments"] is None:
            raw_increments = raw_increments_Z = None
        else:
            raw_increments, raw_increments_Z = _dual_public_map(data["raw_increments"])
        # The value uses the first Z derivative of Delta Mz.  Differentiating
        # that identity once more would require second Z data, so expose only
        # the value and keep the unavailable tangent explicit.
        delta_ur, _ = _public(data["delta_Ur"])
        delta_ur_Z = None
        parameters = data["parameters"]
        return {
            "x": data["parameters"]["x"],
            "R": fields["R"],
            "R_Z": fields_Z["R"],
            "rows": rows,
            "rows_Z": rows_Z,
            "dual_rows": data["rows"],
            "centered_rows": rows,
            "centered_rows_Z": rows_Z,
            "partial_centered_change": rows,
            "partial_centered_change_Z": rows_Z,
            "increments": increments,
            "increments_Z": increments_Z,
            "dual_increments": data["increments"],
            "moment_increments": increments,
            "moment_increments_Z": increments_Z,
            "physical_increments": increments,
            "physical_increments_Z": increments_Z,
            "raw_increments": raw_increments,
            "raw_increments_Z": raw_increments_Z,
            "raw_quadratic_increments": raw_increments,
            "raw_quadratic_increments_Z": raw_increments_Z,
            "normalized_corrections": normalized,
            "normalized_corrections_Z": normalized_Z,
            "normalized_increments": normalized,
            "normalized_increments_Z": normalized_Z,
            "normalized": normalized,
            "normalized_Z": normalized_Z,
            "fields": fields,
            "fields_Z": fields_Z,
            "dual_fields": data["fields"],
            "f": fields["f"],
            "f_Z": fields_Z["f"],
            "g": fields["g"],
            "g_Z": fields_Z["g"],
            "u": fields["Utheta"],
            "u_Z": fields_Z["Utheta"],
            "Utheta": fields["Utheta"],
            "Utheta_Z": fields_Z["Utheta"],
            "Uz": fields["Uz"],
            "Uz_Z": fields_Z["Uz"],
            "F": fields["F"],
            "F_Z": fields_Z["F"],
            "delta_Ur": delta_ur,
            "delta_Ur_Z": delta_ur_Z,
            "parameters": parameters,
            "precision": precision,
            "first_Z_from_same_dual": True,
            "radial_velocity_recovery": "exact moment identity for value using Delta Mz and its first Z derivative",
            "radial_velocity_first_Z_available": False,
            "quadrature_enclosed": False,
            "coefficient_solve_performed": False,
            "functional_closure": False,
            "cone_certified": False,
        }


def _join_reference(
    reference_data: Mapping[str, Any],
    correction: Mapping[str, Any],
    *,
    delta: Any = None,
) -> dict[str, Any]:
    """Merge correction fields with one supplied reference evaluation."""

    if not isinstance(reference_data, Mapping):
        raise TypeError("reference_data must be a mapping from PressureWidthAxialRestore")
    if "P0" in reference_data:
        p0_name = "P0"
    elif "axis_P0" in reference_data:
        p0_name = "axis_P0"
    else:
        raise KeyError("reference_data must provide explicit P0; total P is not an axis datum")
    inferred = _orders(
        reference_data.get("moments", {}),
        reference_data.get("P"),
        correction.get("R"),
    )
    base_moments = _moment_duals(reference_data, inferred)
    increments = {
        key: _dual(value, correction["increments_Z"].get(key, 0), inferred)
        for key, value in correction["increments"].items()
    }
    total_moments = {
        key: base_moments[key] + increments[key] for key in MOMENT_KEYS
    }
    p_base = _field_dual(reference_data, "P", inferred)
    pressure = p_base + _dual(
        correction["increments"]["p"],
        correction["increments_Z"].get("p", 0),
        inferred,
    )
    p0 = _field_dual(reference_data, p0_name, inferred)
    base_raw = _raw_duals(reference_data, inferred)
    raw_increment_values = correction.get("raw_increments")
    raw_increment_tangents = correction.get("raw_increments_Z")
    if base_raw is not None and isinstance(raw_increment_values, Mapping):
        raw_total_duals = {
            key: base_raw[key]
            + _dual(
                raw_increment_values[key],
                raw_increment_tangents.get(key, 0)
                if isinstance(raw_increment_tangents, Mapping)
                else 0,
                inferred,
            )
            for key in ("axial", "swirl")
        }
        raw, raw_Z = _public_dict(raw_total_duals)
        raw_applied = True
    else:
        raw = raw_Z = None
        raw_applied = False
    raw_baseline = reference_data.get("raw_quadratic_integrals")
    raw_baseline_Z = reference_data.get("raw_quadratic_integrals_Z")
    fields = correction["fields"]
    fields_Z = correction["fields_Z"]
    F = _dual(fields["F"], fields_Z["F"], inferred)
    Utheta = _dual(fields["Utheta"], fields_Z["Utheta"], inferred)
    Uz = _dual(fields["Uz"], fields_Z["Uz"], inferred)
    Uz_y = _dual(fields["Uz_y"], fields_Z["Uz_y"], inferred)
    R = _dual(fields["R"], fields_Z["R"], inferred)
    u_y = _dual(fields["Utheta_y"], fields_Z["Utheta_y"], inferred)
    uz_y = _dual(fields["Uz_y"], fields_Z["Uz_y"], inferred)
    z_value = _constant(correction["parameters"]["Z"])
    delta_value = _mp(reference_data.get("delta", 0) if delta is None else delta)
    log_radius = reference_data.get("logR", mp.log(_constant(R)))
    a = 1 - 2 * u_y / Utheta
    shear_theta = -a * F
    shear_z = mp.sqrt(2) * _sqrt(R) * uz_y / R
    moment_values = {key: value.value for key, value in total_moments.items()}
    moment_tangents = {key: value.tangent for key, value in total_moments.items()}
    stress = None
    stress_error = None
    try:
        converter_orders = _orders(
            moment_values,
            moment_tangents,
            pressure.value,
            R.value,
        )

        def converter(value: Any) -> PressureWidthJet:
            if isinstance(value, PressureWidthJet):
                return value
            kwargs = {}
            if converter_orders is not None:
                kwargs = {
                    "pressure_order": converter_orders[0],
                    "width_order": converter_orders[1],
                }
            return PressureWidthJet(value, **kwargs)

        stress = evaluate_mp_stress(
            log_radius,
            z_value,
            delta_value,
            Utheta=Utheta.value,
            Uz=Uz.value,
            Utheta_y=u_y.value,
            Utheta_Z=Utheta.tangent,
            Uz_y=uz_y.value,
            Uz_Z=Uz.tangent,
            moments=moment_values,
            moments_Z=moment_tangents,
            P=pressure.value,
            P_Z=pressure.tangent,
            precision=max(50, int(correction.get("precision", 100))),
            radius_override=R.value,
            scalar_converter=converter,
            axial_override=converter(z_value),
            shear_theta=shear_theta.value,
            shear_z=shear_z.value,
            include_components=True,
        )
    except (KeyError, TypeError, ValueError, ZeroDivisionError) as exc:
        stress_error = str(exc)

    moments, moments_Z = _public_dict(total_moments)
    result = dict(reference_data)
    a_value, a_tangent = _public(a)
    b = 2 * uz_y / Utheta
    b_value, b_tangent = _public(b)
    Am_parameter = correction["parameters"]["Am"]
    Z_parameter = correction["parameters"]["Z"]
    x_power = _mp(correction["x"]) ** mp.mpf(".1")
    unperturbed_u = Am_parameter * x_power
    unperturbed_uz = 4 * Z_parameter
    unperturbed_f = unperturbed_u / (mp.sqrt(2) * _sqrt(R))
    unperturbed_u_value, _ = _public(unperturbed_u)
    unperturbed_uz_value, _ = _public(unperturbed_uz)
    unperturbed_f_value, _ = _public(unperturbed_f)
    baseline_compatibility: dict[str, Any] = {}
    if "moments" in reference_data and "p" in reference_data["moments"]:
        baseline_compatibility["P_minus_P0_minus_Mp"] = (
            reference_data["P"] - reference_data[p0_name] - reference_data["moments"]["p"]
        )
    if "R" in reference_data:
        baseline_compatibility["R_minus_correction_R"] = reference_data["R"] - R.value
    if "Utheta" in reference_data:
        baseline_compatibility["Utheta_minus_unperturbed"] = reference_data["Utheta"] - unperturbed_u_value
    if "Uz" in reference_data:
        baseline_compatibility["Uz_minus_unperturbed"] = reference_data["Uz"] - unperturbed_uz_value
    if "F" in reference_data:
        baseline_compatibility["F_minus_unperturbed"] = reference_data["F"] - unperturbed_f_value
    uz_r = uz_y / R
    uz_r_value, uz_r_tangent = _public(uz_r)
    result.update(
        {
            "R": R.value,
            "R_Z": R.tangent,
            "logR": log_radius,
            "x": correction["x"],
            "a": a_value,
            "a_Z": a_tangent,
            "b": b_value,
            "b_Z": b_tangent,
            "F": F.value,
            "FZ": F.tangent,
            "F_Z": F.tangent,
            "Utheta": Utheta.value,
            "Utheta_Z": Utheta.tangent,
            "Uz": Uz.value,
            "UZ": Uz.tangent,
            "Uz_Z": Uz.tangent,
            "P": pressure.value,
            "PZ": pressure.tangent,
            "P_Z": pressure.tangent,
            "P0": p0.value,
            "P0_Z": p0.tangent,
            "axis_P0": p0.value,
            "axis_P0_Z": p0.tangent,
            "moments": moments,
            "momentsZ": moments_Z,
            "moments_Z": moments_Z,
            "moment_increments": correction["increments"],
            "moment_increments_Z": correction["increments_Z"],
            "physical_increments": correction["increments"],
            "physical_increments_Z": correction["increments_Z"],
            "normalized_corrections": correction["normalized_corrections"],
            "normalized_corrections_Z": correction["normalized_corrections_Z"],
            "normalized": correction["normalized_corrections"],
            "normalized_Z": correction["normalized_corrections_Z"],
            "partial_centered_change": correction["partial_centered_change"],
            "partial_centered_change_Z": correction["partial_centered_change_Z"],
            "f": correction["f"],
            "f_Z": correction["f_Z"],
            "g": correction["g"],
            "g_Z": correction["g_Z"],
            "u": correction["u"],
            "u_Z": correction["u_Z"],
            "u_y": correction["fields"]["Utheta_y"],
            "u_y_Z": correction["fields_Z"]["Utheta_y"],
            "Utheta_y": correction["fields"]["Utheta_y"],
            "Utheta_y_Z": correction["fields_Z"]["Utheta_y"],
            "Uz_y": correction["fields"]["Uz_y"],
            "Uz_y_Z": correction["fields_Z"]["Uz_y"],
            "g_y": correction["fields"]["Uz_y"],
            "g_y_Z": correction["fields_Z"]["Uz_y"],
            "F_R": correction["fields"]["F_R"],
            "F_R_Z": correction["fields_Z"]["F_R"],
            "Uz_R": uz_r_value,
            "Uz_R_Z": uz_r_tangent,
            "Ur": None if stress is None else stress["U_r"],
            "stress": stress,
            "stress_error": stress_error,
            "raw_quadratic_integrals": raw,
            "raw_quadratic_integrals_Z": raw_Z,
            "raw_quadratic_increments": correction.get("raw_increments"),
            "raw_quadratic_increments_Z": correction.get("raw_increments_Z"),
            "raw_quadratic_integrals_baseline_only": raw_baseline,
            "raw_quadratic_integrals_Z_baseline_only": raw_baseline_Z,
            "raw_quadratic_baseline_only": not raw_applied,
            "P0_preserved": True,
            "raw_quadratic_correction_applied": raw_applied,
            "baseline_compatibility": baseline_compatibility,
            "baseline_compatibility_checked": True,
            "baseline_compatibility_certified": False,
            "first_Z_from_same_dual": True,
            "radial_velocity_recovery": correction["radial_velocity_recovery"],
            "radial_velocity_first_Z_available": False,
            "quadrature_enclosed": False,
            "coefficient_solve_performed": False,
            "functional_closure": False,
            "cone_certified": False,
            "source_constants_certified": False,
            "finite_energy_certified": False,
            "global_moment_repair_certified": False,
            "independent_Cartesian_divergence_verified": False,
            "region": "five_bump_physical_similarity_correction",
        }
    )
    if stress is not None:
        result["Ur"] = stress["U_r"]
    result["delta_Ur"] = correction["delta_Ur"]
    result["delta_Ur_Z"] = correction["delta_Ur_Z"]
    # Keep corrected interval contributions visible to downstream moment consumers.
    result['moment_parts'] = dict(reference_data.get('moment_parts', {}))
    result['moment_parts_Z'] = dict(reference_data.get('moment_parts_Z', {}))
    result['moment_parts']['five_bump_correction'] = dict(correction['increments'])
    result['moment_parts_Z']['five_bump_correction'] = dict(correction['increments_Z'])
    if correction.get('raw_increments') is not None:
        result['raw_quadratic_parts'] = dict(reference_data.get('raw_quadratic_parts') or {})
        result['raw_quadratic_parts_Z'] = dict(reference_data.get('raw_quadratic_parts_Z') or {})
        result['raw_quadratic_parts']['five_bump_correction'] = dict(correction['raw_increments'])
        result['raw_quadratic_parts_Z']['five_bump_correction'] = dict(correction['raw_increments_Z'])
        result['raw_quadratic_parts_complete'] = reference_data.get('raw_quadratic_integrals') is not None
    return result


def join(
    reference_data: Mapping[str, Any],
    correction: Mapping[str, Any] | None = None,
    *args: Any,
    moment_map: Any = None,
    h: Any = None,
    Am: Any = None,
    Rm: Any = None,
    Z: Any = None,
    x: Any = None,
    delta: Any = None,
    h_Z: Any = None,
    Am_Z: Any = None,
    Rm_Z: Any = None,
) -> dict[str, Any]:
    """Join a correction to one ``PressureWidthAxialRestore`` evaluation.

    The preferred form is ``join(reference_data, correction)`` after
    :func:`evaluate_correction`.  For convenience, callers may pass the map
    and physical arguments as keyword arguments, or use positional
    ``join(reference_data, map, h, Am, Rm, Z, x)``.
    """

    if correction is not None and hasattr(correction, "partial"):
        if args:
            values = list(args)
        else:
            values = [h, Am, Rm, Z, x]
        if len(values) != 5:
            raise TypeError("map form requires h, Am, Rm, Z, and x")
        moment_map = correction
        delta_for_eval = reference_data.get("delta", 0) if delta is None else delta
        correction = evaluate_correction(
            moment_map,
            values[0],
            values[1],
            values[2],
            values[3],
            values[4],
            delta_for_eval,
            h_Z=h_Z,
            Am_Z=Am_Z,
            Rm_Z=Rm_Z,
        )
    elif correction is None:
        if moment_map is None or h is None or Am is None or Rm is None or Z is None or x is None:
            raise TypeError("join requires a correction or map, h, Am, Rm, Z, and x")
        delta_for_eval = reference_data.get("delta", 0) if delta is None else delta
        correction = evaluate_correction(
            moment_map,
            h,
            Am,
            Rm,
            Z,
            x,
            delta_for_eval,
            h_Z=h_Z,
            Am_Z=Am_Z,
            Rm_Z=Rm_Z,
        )
    return _join_reference(reference_data, correction, delta=delta)


class FiveBumpField:
    """Reference-backed adapter using a caller-supplied coefficient provider."""

    branch = "five_bump_physical_similarity_correction"

    def __init__(
        self,
        reference: Any,
        moment_map: Any,
        coefficients: Any,
        *,
        Am: Any = None,
        Rm: Any = None,
        delta: Any = None,
    ) -> None:
        if not hasattr(reference, "evaluate_phase"):
            raise TypeError("reference must expose evaluate_phase")
        if not hasattr(moment_map, "partial"):
            raise TypeError("moment_map must expose partial")
        self.reference = reference
        self.moment_map = moment_map
        self.coefficients = coefficients
        self.Am = Am
        self.Rm = Rm
        self.delta = delta
        self.precision = int(getattr(reference, "work_precision", getattr(moment_map, "precision", 100)))

    def _provider_values(self, Z: Any, baseline: Mapping[str, Any] | None = None) -> dict[str, Any]:
        values = self.coefficients(Z) if callable(self.coefficients) else self.coefficients
        if isinstance(values, Mapping):
            result = dict(values)
        else:
            try:
                sequence = tuple(values)
            except TypeError as exc:
                raise TypeError("coefficient provider must return a mapping or sequence") from exc
            if len(sequence) < 3:
                raise ValueError("coefficient provider sequence requires h, Am, and Rm")
            result = {"h": sequence[0], "Am": sequence[1], "Rm": sequence[2]}
            if len(sequence) > 3:
                result["h_Z"] = sequence[3]
            if len(sequence) > 4:
                result["Am_Z"] = sequence[4]
            if len(sequence) > 5:
                result["Rm_Z"] = sequence[5]
        for key in ("h", "Am", "Rm"):
            if key not in result and key == "Am" and self.Am is not None:
                result[key] = self.Am
            if key not in result and key == "Rm" and self.Rm is not None:
                result[key] = self.Rm
        if "Rm" not in result or "Am" not in result:
            if baseline is None:
                raise KeyError("coefficient provider requires Am and Rm or a reference endpoint")
            endpoint = self.reference.evaluate_phase(2, Z)
            result.setdefault("Rm", endpoint["R"])
            result.setdefault("Am", mp.sqrt(2) * _sqrt(endpoint["R"]) * endpoint["F"])
        if "h" not in result:
            raise KeyError("coefficient provider is missing 'h'")
        return result

    def evaluate_x(self, x: Any, Z: Any) -> dict[str, Any]:
        with mp.workdps(self.precision):
            x_value = _mp(x)
            if not DOMAIN_LEFT <= x_value <= mp.exp(1):
                raise ValueError("require 1 <= x <= exp(1)")
            baseline = self.reference.evaluate_phase(2 + mp.log(x_value), Z)
            values = self._provider_values(Z, baseline)
            return join(
                baseline,
                moment_map=self.moment_map,
                h=values["h"],
                Am=values["Am"],
                Rm=values["Rm"],
                Z=Z,
                x=x_value,
                delta=self.delta,
                h_Z=values.get("h_Z"),
                Am_Z=values.get("Am_Z"),
                Rm_Z=values.get("Rm_Z"),
            )

    evaluate = evaluate_x

    def metadata(self) -> dict[str, Any]:
        return {
            "branch": self.branch,
            "coordinate": "x=R/Rm",
            "domain": "1 <= x <= exp(1)",
            "coefficients_external": True,
            "reference_pressure_width_axial_restore": True,
            "P0_preserved": True,
            "raw_quadratic_updated_when_baseline_available": True,
            "baseline_compatibility_certified": False,
            "first_Z_from_same_dual": True,
            "radial_velocity_first_Z_available": False,
            "quadrature_enclosed": False,
            "coefficient_solve_performed": False,
            "functional_closure": False,
            "cone_certified": False,
        }


FiveBumpCorrectionField = FiveBumpField


__all__ = [
    "FiveBumpField",
    "FiveBumpCorrectionField",
    "evaluate_correction",
    "join",
    "MOMENT_KEYS",
    "COEFFICIENT_NAMES",
]
