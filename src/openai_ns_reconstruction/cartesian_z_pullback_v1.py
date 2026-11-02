"""Exact scalar-profile pullback from ordinary ``Z`` to physical ``z``.

This leaf is intentionally narrower than a velocity, pressure, or
Navier--Stokes operator. It evaluates a scalar profile ``f(rho, Z)`` whose
ordinary jet is supplied by the caller. The caller also supplies the
already-bound ``rho`` and ``lambda`` values; the relations

    lambda**2 = tau/(1 - Z**2),
    R = (X**2 + Y**2) / (2*lambda**2),
    rho = R/epsilon,
    z_phys = lambda**(1-delta)*Z

are recorded as geometry premises. This leaf never silently replaces a
caller-owned ``rho`` or ``lambda`` with values from another source sector.
All arithmetic is performed through the directed context ``c``: no midpoint,
float conversion, or sampled ``Z`` is introduced.

The ordinary jet names are explicit: ``value, rho, Z, rhoZ, rhorho, ZZ``
mean ``f``, ``f_rho``, ``f_Z``, ``f_rhoZ``, ``f_rhorho``, and ``f_ZZ``.
Missing radial jets are rejected. The optional constant ``weight_power``
evaluates the full product ``lambda**p * f`` and never treats the weight as a
constant under ``Z``. Results are conditional scalar-profile data only; they
do not claim a physical velocity, pressure, NS validation, or m15 term.
Source-bound identifiers and indices can be passed under
``params['source_provenance']`` or ``provenance`` and are retained.
"""

from __future__ import annotations

from collections.abc import Mapping
from fractions import Fraction
from typing import Any


SCHEMA = "cartesian_z_pullback_v1"
REQUIRED_PARAMS = ("Z", "rho", "lambda", "delta")
REQUIRED_JETS = ("value", "rho", "Z", "rhoZ", "rhorho", "ZZ")
FIXED_COORDINATES = ("X", "Y", "tau", "epsilon", "delta")


class PullbackInputError(ValueError):
    """Raised when a caller-owned scalar or coordinate contract is absent."""


def _is_interval(value: Any) -> bool:
    return hasattr(value, "a") and hasattr(value, "b")


def _bounds(value: Any) -> tuple[Any, Any]:
    if _is_interval(value):
        return getattr(value, "a"), getattr(value, "b")
    return value, value


def _as_value(c: Any, value: Any, name: str) -> Any:
    """Convert exact caller data through ``c`` without accepting binary floats."""

    if _is_interval(value):
        return value
    if isinstance(value, Mapping):
        endpoint_keys = ("lower_exact_mpf_tuple", "upper_exact_mpf_tuple")
        present = [key in value for key in endpoint_keys]
        if not any(present):
            raise TypeError(
                f"{name} mapping must contain lower_exact_mpf_tuple and upper_exact_mpf_tuple"
            )
        if not all(present):
            raise PullbackInputError(f"{name} has an incomplete directed endpoint mapping")
        if value.get("directed") is not True:
            raise PullbackInputError(f"{name} directed endpoint mapping requires directed=True")
        lower = _exact_mpf_tuple(value[endpoint_keys[0]], name + ".lower")
        upper = _exact_mpf_tuple(value[endpoint_keys[1]], name + ".upper")
        maker = getattr(c, "make_mpf", None)
        if callable(maker):
            try:
                result = maker((lower, upper))
            except (TypeError, ValueError, ArithmeticError) as exc:
                raise PullbackInputError(f"{name} directed endpoint mapping is invalid") from exc
        else:
            converter = getattr(c, "mpf", None)
            if not callable(converter):
                raise TypeError("directed context c must expose make_mpf or mpf")
            try:
                result = converter((lower, upper))
            except (TypeError, ValueError, ArithmeticError) as exc:
                raise PullbackInputError(f"{name} directed endpoint mapping is invalid") from exc
        if not _is_interval(result):
            raise TypeError(f"{name} directed endpoint mapping did not produce an interval")
        result_lower, result_upper = _bounds(result)
        try:
            if not bool(result_lower <= result_upper):
                raise PullbackInputError(f"{name} directed endpoint mapping is reversed")
        except PullbackInputError:
            raise
        except (TypeError, ValueError, ArithmeticError) as exc:
            raise PullbackInputError(f"{name} directed endpoint ordering is unavailable") from exc
        return result
    if isinstance(value, bool):
        raise TypeError(f"{name} must be a scalar, not bool")
    if isinstance(value, float):
        raise TypeError(f"{name} must be exact (use int, Fraction, str, or c.mpf)")
    if isinstance(value, Fraction):
        return c.mpf(value.numerator) / value.denominator
    if isinstance(value, (int, str)):
        return c.mpf(value)
    try:
        return c.mpf(value)
    except (TypeError, ValueError, ArithmeticError) as exc:
        raise TypeError(f"{name} is not an exact caller scalar") from exc


def _exact_mpf_tuple(value: Any, name: str) -> tuple[int, int, int, int]:
    """Validate one finite four-integer mpmath ``_mpf_`` tuple."""

    if not isinstance(value, (list, tuple)) or len(value) != 4:
        raise PullbackInputError(f"{name} must be an integer 4-tuple")
    if any(isinstance(item, bool) or not isinstance(item, int) for item in value):
        raise PullbackInputError(f"{name} must contain integers only")
    sign, man, exponent, bit_count = tuple(value)
    if sign not in (0, 1) or man < 0 or bit_count < 0:
        raise PullbackInputError(f"{name} has invalid mpmath tuple fields")
    # mpmath reserves man=0, exponent!=0 for infinities/NaN.  A finite zero
    # is exactly (sign, 0, 0, 0); a nonzero finite mantissa carries its bit
    # length in the fourth field.
    if man == 0:
        if exponent != 0 or bit_count != 0:
            raise PullbackInputError(f"{name} must be finite")
    elif bit_count != man.bit_length():
        raise PullbackInputError(f"{name} has a noncanonical mantissa bit count")
    return sign, man, exponent, bit_count


def _copy_mapping(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _copy_mapping(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_copy_mapping(item) for item in value]
    return value


def _fixed_contract(params: Mapping[str, Any]) -> dict[str, Any]:
    """Require an explicit fixed-physical-coordinate declaration."""

    raw = None
    raw_name = None
    for key in ("fixed_conditions", "coordinate_fixed", "fixed_coordinates"):
        if key in params:
            raw = params[key]
            raw_name = key
            break
    if raw is None:
        raise PullbackInputError(
            "explicit fixed_conditions must declare X,Y,tau,epsilon,delta as True"
        )
    if not isinstance(raw, Mapping):
        raise PullbackInputError("fixed coordinate conditions must be a mapping")
    aliases = {
        "X": ("X", "physical_X"),
        "Y": ("Y", "physical_Y"),
        "tau": ("tau",),
        "epsilon": ("epsilon",),
        "delta": ("delta",),
    }
    result: dict[str, Any] = {"provided": True, "required": True, "source_key": raw_name}
    for label, choices in aliases.items():
        present = [raw[name] for name in choices if name in raw]
        if not present:
            raise PullbackInputError(f"explicit fixed conditions omit {label}")
        if any(item is not True for item in present):
            raise PullbackInputError(f"fixed condition {label} must be True")
        result[label] = True
    return result


def _validate_params(params: Any) -> dict[str, Any]:
    if not isinstance(params, Mapping):
        raise TypeError("params must be a mapping")
    missing = [name for name in REQUIRED_PARAMS if name not in params]
    if missing:
        raise PullbackInputError("missing geometry parameters: " + ", ".join(missing))
    return _fixed_contract(params)


def _jet(c: Any, jet: Any) -> dict[str, Any]:
    if not isinstance(jet, Mapping):
        raise TypeError("jet must be a mapping with explicit scalar profile fields")
    missing = [name for name in REQUIRED_JETS if name not in jet]
    if missing:
        raise PullbackInputError("missing required profile jets: " + ", ".join(missing))
    return {name: _as_value(c, jet[name], "jet." + name) for name in REQUIRED_JETS}


def _contains_zero(value: Any) -> bool:
    lower, upper = _bounds(value)
    return bool(lower <= 0 and upper >= 0)


def _require_domain(z: Any, rho: Any, lam: Any, denominator: Any) -> None:
    z_lower, z_upper = _bounds(z)
    rho_lower, _ = _bounds(rho)
    lam_lower, _ = _bounds(lam)
    if not bool(z_lower > -1 and z_upper < 1):
        raise PullbackInputError("domain requires |Z| < 1 on the full directed interval")
    if not bool(rho_lower >= 0):
        raise PullbackInputError("domain requires rho >= 0 on the full directed interval")
    if not bool(lam_lower > 0):
        raise PullbackInputError("domain requires lambda > 0 on the full directed interval")
    if _contains_zero(denominator):
        raise PullbackInputError("domain excludes zero of 1 - delta*Z**2")


def _power(c: Any, base: Any, exponent: Any) -> Any:
    function = getattr(c, "power", None)
    if callable(function):
        return function(base, exponent)
    return base ** exponent


def _geometry(c: Any, params: Mapping[str, Any]) -> dict[str, Any]:
    z = _as_value(c, params["Z"], "params.Z")
    rho = _as_value(c, params["rho"], "params.rho")
    lam = _as_value(c, params["lambda"], "params.lambda")
    delta = _as_value(c, params["delta"], "params.delta")
    one = _as_value(c, 1, "one")
    denominator = one - delta * z * z
    _require_domain(z, rho, lam, denominator)
    one_minus_z2 = one - z * z
    rho_Z = -2 * z * rho / one_minus_z2
    rho_ZZ = -2 * rho / one_minus_z2
    lambda_power = _power(c, lam, one - delta)
    B = one_minus_z2 / (lambda_power * denominator)
    connection = (-(3 - delta) * z / one_minus_z2
                  + 2 * delta * z / denominator)
    z_phys = lambda_power * z
    sqrt = getattr(c, "sqrt", None)
    if not callable(sqrt):
        raise TypeError("directed context c must expose sqrt")
    x_normalized = sqrt(rho)
    return {
        "Z": z,
        "rho": rho,
        "lambda": lam,
        "delta": delta,
        "one_minus_Z2": one_minus_z2,
        "one_minus_delta_Z2": denominator,
        "rho_Z": rho_Z,
        "rho_ZZ": rho_ZZ,
        "lambda_power_1_minus_delta": lambda_power,
        "B": B,
        "connection": connection,
        "z_phys": z_phys,
        "x_normalized": x_normalized,
    }


def _ordinary(geometry: Mapping[str, Any], jet: Mapping[str, Any]) -> dict[str, Any]:
    rho_Z = geometry["rho_Z"]
    rho_ZZ = geometry["rho_ZZ"]
    f = jet["value"]
    f_rho = jet["rho"]
    f_Z = jet["Z"]
    f_rhoZ = jet["rhoZ"]
    f_rhorho = jet["rhorho"]
    f_ZZ = jet["ZZ"]
    D = f_Z + rho_Z * f_rho
    D2 = (f_ZZ + 2 * rho_Z * f_rhoZ
          + rho_Z * rho_Z * f_rhorho + rho_ZZ * f_rho)
    B = geometry["B"]
    connection = geometry["connection"]
    return {
        "value": f,
        "partial_rho": f_rho,
        "partial_Z": f_Z,
        "partial_rhoZ": f_rhoZ,
        "partial_rhorho": f_rhorho,
        "partial_ZZ": f_ZZ,
        "D": D,
        "D2": D2,
        "dz_phys": B * D,
        "dzz_phys": B * B * (D2 + connection * D),
    }


def _zero(value: Any) -> bool:
    lower, upper = _bounds(value)
    return bool(lower == 0 and upper == 0)


def _weighted(
    c: Any,
    geometry: Mapping[str, Any],
    ordinary: Mapping[str, Any],
    weight_power: Any,
) -> dict[str, Any]:
    p = _as_value(c, 0 if weight_power is None else weight_power, "weight_power")
    z = geometry["Z"]
    one_minus_z2 = geometry["one_minus_Z2"]
    scale = _power(c, geometry["lambda"], p)
    q = p * z / one_minus_z2
    q_Z = p * (1 + z * z) / (one_minus_z2 * one_minus_z2)
    f = ordinary["value"]
    D = scale * (ordinary["D"] + q * f)
    D2 = scale * (
        ordinary["D2"] + 2 * q * ordinary["D"] + (q_Z + q * q) * f
    )
    B = geometry["B"]
    connection = geometry["connection"]
    return {
        "value": scale * f,
        "D": D,
        "D2": D2,
        "dz_phys": B * D,
        "dzz_phys": B * B * (D2 + connection * D),
        "lambda_power": scale,
        "weight_power": p,
        "lambda_log_derivative": z / one_minus_z2,
        "weight_log_derivative": q,
        "weight_log_second_derivative": q_Z,
        "product_rule": "D(lambda**p*f)=lambda**p*(Df+p*Z/(1-Z**2)*f)",
    }


def _provenance(params: Mapping[str, Any], supplied: Any) -> dict[str, Any]:
    candidates = []
    for key in ("source_provenance", "provenance"):
        value = params.get(key)
        if isinstance(value, Mapping):
            candidates.append(value)
    if isinstance(supplied, Mapping):
        candidates.append(supplied)
    record: dict[str, Any] = {
        "inputs": [_copy_mapping(item) for item in candidates],
        "source_owner_proved": False,
    }
    for key in ("source_ids", "source_indices", "source_paths", "error_classes"):
        values = []
        for item in candidates:
            if key not in item:
                continue
            raw = item[key]
            values.extend(list(raw) if isinstance(raw, (list, tuple)) else [raw])
        if values:
            record[key] = values
    return record


def _stencil_metadata(fixed: Mapping[str, Any], provenance: Any) -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "status": "conditional_scalar_profile_pullback",
        "profile_kind": "scalar_f(rho,Z)",
        "required_params": list(REQUIRED_PARAMS),
        "required_jet_fields": list(REQUIRED_JETS),
        "ordinary_Z_orders": [0, 1, 2],
        "fixed_coordinate_contract": dict(fixed),
        "formula": {
            "rho_Z": "-2*Z*rho/(1-Z**2)",
            "rho_ZZ": "-2*rho/(1-Z**2)",
            "D": "f_Z + rho_Z*f_rho",
            "D2": "f_ZZ + 2*rho_Z*f_rhoZ + rho_Z**2*f_rhorho + rho_ZZ*f_rho",
            "B": "(1-Z**2)/(lambda**(1-delta)*(1-delta*Z**2))",
            "dz_phys": "B*D",
            "dzz_phys": "B**2*(D2+(-(3-delta)*Z/(1-Z**2)+2*delta*Z/(1-delta*Z**2))*D)",
            "x_normalized": "sqrt(rho)",
            "lambda_relation": "lambda**2=tau/(1-Z**2)",
            "rho_relation": "rho=((X**2+Y**2)/(2*lambda**2))/epsilon",
        },
        "weight_product_rule": {
            "p_is_Z_constant": True,
            "formula_D": "D(lambda**p*f)=lambda**p*(D+p*Z/(1-Z**2)*f)",
            "formula_D2": "D2(lambda**p*f)=lambda**p*(D2+2*q*D+(q_Z+q**2)*f)",
            "q": "p*Z/(1-Z**2)",
            "q_Z": "p*(1+Z**2)/(1-Z**2)**2",
        },
        "input_provenance": _provenance({}, provenance),
        "physical_layer_complete": False,
        "pressure": False,
        "m15": False,
        "ns_validation": False,
    }


def stencil(params: Mapping[str, Any], *, provenance: Any = None) -> dict[str, Any]:
    """Return the unevaluated scalar stencil and its caller contracts."""

    fixed = _validate_params(params)
    result = _stencil_metadata(fixed, provenance)
    result["input_provenance"] = _provenance(params, provenance)
    return result


def build(params: Mapping[str, Any], *, provenance: Any = None) -> dict[str, Any]:
    """Alias for :func:`stencil` kept for callers that use a build phase."""

    return stencil(params, provenance=provenance)


def evaluate(
    c: Any,
    params: Mapping[str, Any],
    jet: Mapping[str, Any],
    *,
    weight_power: Any = 0,
    provenance: Any = None,
) -> dict[str, Any]:
    """Evaluate the scalar pullback with directed interval arithmetic."""

    fixed = _validate_params(params)
    values = _jet(c, jet)
    geometry = _geometry(c, params)
    ordinary = _ordinary(geometry, values)
    weighted = _weighted(c, geometry, ordinary, weight_power)
    selected = ordinary if _zero(weighted["weight_power"]) else weighted
    source = _provenance(params, provenance)
    result = {
        "schema": SCHEMA,
        "status": "conditional_scalar_profile_pullback",
        "profile_kind": "scalar_f(rho,Z)",
        "coordinate_fixed_conditions": fixed,
        "geometry_premise": {
            "lambda_relation": "lambda**2=tau/(1-Z**2)",
            "rho_relation": "rho=((X**2+Y**2)/(2*lambda**2))/epsilon",
            "z_phys_relation": "z_phys=lambda**(1-delta)*Z",
            "caller_supplied_lambda_and_rho": True,
            "physical_values_auto_verified": False,
            "physical_values_consistency_is_caller_premise": True,
            "source_owner_proved": False,
        },
        "geometry": geometry,
        "ordinary_Z": {
            "Z0": ordinary["value"],
            "Z1": ordinary["D"],
            "Z2": ordinary["D2"],
            "partial_Z": ordinary["partial_Z"],
            "partial_ZZ": ordinary["partial_ZZ"],
        },
        "unweighted": ordinary,
        "weighted": weighted,
        "profile": selected,
        "value": selected["value"],
        "D": selected["D"],
        "D2": selected["D2"],
        "dz_phys": selected["dz_phys"],
        "dzz_phys": selected["dzz_phys"],
        "weight_power": weighted["weight_power"],
        "weight_applied": not _zero(weighted["weight_power"]),
        "jet_fields": list(REQUIRED_JETS),
        "input_provenance": source,
        "physical_layer_complete": False,
        "velocity_weight_is_scalar_profile_only": True,
        "pressure": False,
        "m15": False,
        "ns_validation": False,
    }
    return result


__all__ = [
    "SCHEMA",
    "REQUIRED_PARAMS",
    "REQUIRED_JETS",
    "FIXED_COORDINATES",
    "PullbackInputError",
    "stencil",
    "build",
    "evaluate",
]
