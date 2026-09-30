"""Independent resolved check for the Section 9.30 reference extension.

The production reference-extension adapter is supplied by the neighbouring
``pressure_width_reference_extension`` module.  This fixture starts from the
actual resolved endpoint of ``PressureWidthLongReshape`` at ``Rsh`` and
replays the constant-``V`` reference branch with independent scalar
``mp.quad`` integrals.  The first ``Z`` derivatives are checked with a
fourth-order centred stencil of those scalar integrals.

The receipt checks deliberately keep the original endpoint seeds, the long
reshape increments, and the reference-extension increments as three separate
objects.  A summed moment is useful for stress evaluation, but is not a
replacement for any of those receipts.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

import mpmath as mp

try:
    from .lei_ren_part1_paper_axial_dual import AxialDual
    from .lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
    from .lei_ren_part1_paper_pressure_width_long_reshape import (
        PressureWidthLongReshape,
    )
    from .lei_ren_part1_paper_pressure_width_long_reshape_fixture import (
        _ResolvedR110Provider,
        _build_comparison,
    )
    from .lei_ren_part1_paper_pressure_width_reference_extension import (
        PressureWidthReferenceExtension,
    )
except (ImportError, ValueError):
    from lei_ren_part1_paper_axial_dual import AxialDual
    from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
    from lei_ren_part1_paper_pressure_width_long_reshape import (
        PressureWidthLongReshape,
    )
    from lei_ren_part1_paper_pressure_width_long_reshape_fixture import (
        _ResolvedR110Provider,
        _build_comparison,
    )
    from lei_ren_part1_paper_pressure_width_reference_extension import (
        PressureWidthReferenceExtension,
    )


PRECISION = 100
PRESSURE_ORDER = 3
WIDTH_ORDER = 2
RESHAPE_A = mp.mpf(1)
RESHAPE_T = mp.mpf(400)
LOG_C = mp.mpf(-4)
QUADRATURE_ORDER = 24
TAIL_DIGITS = 50
Z0 = mp.mpf(".3")
LOG_RREF = mp.log(110) + RESHAPE_T + 8 + 5
EXTENSION_LENGTH = mp.mpf(5)
EXTENSION_POINTS = (mp.mpf(1), mp.mpf(5))
FINITE_DIFFERENCE_STEP = mp.mpf("1e-4")
ERROR_THRESHOLD = mp.mpf("1e-8")
RECEIPT_THRESHOLD = mp.mpf("1e-30")

MOMENT_NAMES = ("theta", "z", "theta_z", "z_theta", "p")
RAW_NAMES = ("axial", "swirl")


def _scalar(value: Any) -> mp.mpf:
    """Materialize the resolved value of an axial dual or pressure jet."""

    if isinstance(value, AxialDual):
        value = value.value
    if isinstance(value, PressureWidthJet):
        return value.evaluate(pressure=1, width=1)
    return mp.mpf(str(value))


def _scaled(actual: Any, expected: Any) -> mp.mpf:
    a = _scalar(actual)
    b = _scalar(expected)
    return abs(a - b) / max(mp.mpf(1), abs(a), abs(b))


def _field(source: Mapping[str, Any], *names: str) -> Any:
    for name in names:
        if name in source:
            return source[name]
    raise KeyError(f"missing field; tried {names}")


def _moments(source: Mapping[str, Any], *, derivative: bool = False) -> Mapping[str, Any]:
    names = ("moments_Z", "momentsZ") if derivative else ("moments",)
    for name in names:
        value = source.get(name)
        if isinstance(value, Mapping):
            return value
    raise KeyError(f"missing {'first-Z ' if derivative else ''}moment mapping")


def _receipt(source: Mapping[str, Any], *names: str) -> Mapping[str, Any]:
    for name in names:
        value = source.get(name)
        if isinstance(value, Mapping):
            return value
    raise KeyError(f"missing receipt; tried {names}")


def _fourth_difference(values: Mapping[int, mp.mpf], step: mp.mpf) -> mp.mpf:
    return (-values[2] + 8 * values[1] - 8 * values[-1] + values[-2]) / (12 * step)


def _reference_integrals(
    x: mp.mpf,
    start: Mapping[str, Any],
) -> dict[str, mp.mpf]:
    """Independent physical moment increments on ``Rsh <= R <= Rsh exp(x)``."""

    if x == 0:
        return {name: mp.mpf(0) for name in ("theta", "z", "theta_z", "p", "swirl", "axial", "z_theta")}
    Rs = _scalar(_field(start, "R"))
    u_start = _field(start, "Utheta", "u")
    if u_start is None:
        u_start = mp.sqrt(2 * Rs) * _scalar(_field(start, "F"))
    u_start = _scalar(u_start)
    V = _scalar(_field(start, "Uz", "V"))

    def radius(t: mp.mpf) -> mp.mpf:
        return Rs * mp.exp(t)

    def swirl(t: mp.mpf) -> mp.mpf:
        R = radius(t)
        u = u_start * mp.exp(t / 10)
        return u

    def quad(function: Any) -> mp.mpf:
        # Splitting at the midpoint keeps the resolved check independent of
        # the production Gauss rule while remaining stable for x = 5.
        return mp.quad(function, [0, x / 2, x])

    dtheta = lambda t: mp.sqrt(2) * radius(t) ** mp.mpf("1.5") * swirl(t)
    dp = lambda t: swirl(t) ** 2 / 2
    dswirl = lambda t: radius(t) * swirl(t) ** 2 / 2
    dz = lambda t: V * radius(t)
    daxial = lambda t: V * V * radius(t)
    return {
        "theta": quad(dtheta),
        "z": quad(dz),
        "theta_z": V * quad(dtheta),
        "p": quad(dp),
        "swirl": quad(dswirl),
        "axial": quad(daxial),
        "z_theta": quad(daxial) - quad(dswirl),
    }


def _scalar_start(
    reshape: PressureWidthLongReshape,
    z: mp.mpf,
) -> dict[str, Any]:
    return reshape.evaluate_log_offset(reshape.T, z)


def _expected_at(
    reshape: PressureWidthLongReshape,
    x: mp.mpf,
    z: mp.mpf,
) -> dict[str, Any]:
    start = _scalar_start(reshape, z)
    Rs = _scalar(_field(start, "R"))
    u_start = _scalar(_field(start, "Utheta", "u"))
    V = _scalar(_field(start, "Uz", "V"))
    increments = _reference_integrals(x, start)
    R = Rs * mp.exp(x)
    u = u_start * mp.exp(x / 10)
    F = u / mp.sqrt(2 * R)
    moments = {
        name: _scalar(_moments(start)[name])
        for name in MOMENT_NAMES
    }
    moments["theta"] += increments["theta"]
    moments["z"] += increments["z"]
    moments["theta_z"] += increments["theta_z"]
    moments["p"] += increments["p"]
    moments["z_theta"] += increments["z_theta"]
    return {
        "R": R,
        "Utheta": u,
        "F": F,
        "Uz": V,
        "P": _scalar(_field(start, "P")) + increments["p"],
        "moments": moments,
        "increments": increments,
    }


def _extension_increments(output: Mapping[str, Any]) -> Mapping[str, Any]:
    """Read the explicit reference receipt from the production adapter."""

    return _receipt(
        output,
        "moment_increments",
        "reference_moment_increments",
        "extension_moment_increments",
        "reference_increments",
        "extension_increments",
    )


def _extension_raw_increments(output: Mapping[str, Any]) -> Mapping[str, Any]:
    parts = output.get("raw_quadratic_parts")
    if isinstance(parts, Mapping) and isinstance(parts.get("reference"), Mapping):
        return parts["reference"]
    return _receipt(
        output,
        "reference_raw_quadratic_increments",
        "extension_raw_quadratic_increments",
        "raw_reference_increments",
        "raw_extension_increments",
    )


def _extension_increment_tangents(output: Mapping[str, Any]) -> Mapping[str, Any]:
    parts = output.get("moment_parts_Z")
    if isinstance(parts, Mapping) and isinstance(parts.get("reference"), Mapping):
        return parts["reference"]
    return _receipt(
        output,
        "moment_increments_Z",
        "reference_moment_increments_Z",
        "extension_moment_increments_Z",
        "reference_increments_Z",
        "extension_increments_Z",
    )


def _extension_receipt_names(output: Mapping[str, Any]) -> dict[str, str]:
    moment_name = next(
        name
        for name in (
            "moment_increments",
            "reference_moment_increments",
            "extension_moment_increments",
            "reference_increments",
            "extension_increments",
        )
        if isinstance(output.get(name), Mapping)
    )
    raw_name = next(
        name
        for name in (
            "raw_quadratic_parts",
            "reference_raw_quadratic_increments",
            "extension_raw_quadratic_increments",
            "raw_reference_increments",
            "raw_extension_increments",
        )
        if isinstance(output.get(name), Mapping)
    )
    if raw_name == "raw_quadratic_parts":
        raw_name = "raw_quadratic_parts.reference"
    return {"moment_increments": moment_name, "raw_increments": raw_name}


def _build_fixture_objects(precision: int) -> tuple[Any, PressureWidthLongReshape, Any]:
    comparison = _build_comparison(precision)
    switches = _ResolvedR110Provider(comparison)
    reshape = PressureWidthLongReshape(
        switches,
        A=RESHAPE_A,
        logC=LOG_C,
        order=QUADRATURE_ORDER,
        tail_digits=TAIL_DIGITS,
    )
    extension = PressureWidthReferenceExtension(reshape, logRref=LOG_RREF)
    return switches, reshape, extension


def run_fixture(
    *,
    precision: int = PRECISION,
    points: tuple[Any, ...] = EXTENSION_POINTS,
    finite_difference_step: Any = FINITE_DIFFERENCE_STEP,
) -> dict[str, Any]:
    precision = max(80, int(precision))
    step = mp.mpf(str(finite_difference_step))
    if step <= 0:
        raise ValueError("finite_difference_step must be positive")

    with mp.workdps(precision):
        switches, reshape, extension = _build_fixture_objects(precision)
        points = tuple(mp.mpf(str(point)) for point in points)
        if not points or any(point < 0 or point > EXTENSION_LENGTH for point in points):
            raise ValueError("resolved points must lie in [0, 5]")

        source_endpoint = switches.evaluate_R(110, Z0)
        reshape_start = _scalar_start(reshape, Z0)
        zero = extension.evaluate_log_offset(0, Z0)

        continuity_errors = {
            name: _scaled(_field(zero, name), _field(reshape_start, name))
            for name in ("R", "F", "Utheta", "Uz", "P", "P0")
        }
        continuity_errors.update(
            {
                f"field_Z:{name}": _scaled(
                    _field(zero, f"{name}_Z", f"{name}Z"),
                    _field(reshape_start, f"{name}_Z", f"{name}Z"),
                )
                for name in ("F", "Utheta", "Uz", "P", "P0")
            }
        )
        continuity_errors.update(
            {
                f"moment:{name}": _scaled(_moments(zero)[name], _moments(reshape_start)[name])
                for name in MOMENT_NAMES
            }
        )

        # Confirm that the extension exposes all three receipts at x = 0.
        # The maps may have different descriptive names, but each must be
        # present and numerically equal to the corresponding long-reshape
        # receipt or zero extension increment.
        zero_extension_moments = _extension_increments(zero)
        zero_extension_raw = _extension_raw_increments(zero)
        for name in MOMENT_NAMES:
            continuity_errors[f"receipt:extension_zero:{name}"] = abs(
                _scalar(zero_extension_moments[name])
            )
        for name in RAW_NAMES:
            continuity_errors[f"receipt:raw_extension_zero:{name}"] = abs(
                _scalar(zero_extension_raw[name])
            )

        integral_errors: dict[str, mp.mpf] = {}
        tangent_errors: dict[str, mp.mpf] = {}
        receipt_errors: dict[str, mp.mpf] = {}
        rows: list[dict[str, Any]] = []
        receipt_names: dict[str, str] | None = None

        for x in points:
            output = extension.evaluate_log_offset(x, Z0)
            expected = _expected_at(reshape, x, Z0)
            receipt_names = receipt_names or _extension_receipt_names(output)

            for name in ("R", "F", "Utheta", "Uz", "P", "P0"):
                if name == "P0":
                    expected_value = _field(reshape_start, "P0", "axis_P0")
                elif name == "R":
                    expected_value = expected["R"]
                else:
                    expected_value = expected[name]
                integral_errors[f"field:{name}@x={mp.nstr(x, 6)}"] = _scaled(
                    _field(output, name), expected_value
                )
            for name in MOMENT_NAMES:
                integral_errors[f"moment:{name}@x={mp.nstr(x, 6)}"] = _scaled(
                    _moments(output)[name], expected["moments"][name]
                )

            extension_moments = _extension_increments(output)
            for name in MOMENT_NAMES:
                receipt_errors[f"moment:{name}@x={mp.nstr(x, 6)}"] = _scaled(
                    extension_moments[name], expected["increments"][name]
                )
            extension_raw = _extension_raw_increments(output)
            for name in RAW_NAMES:
                receipt_errors[f"raw:{name}@x={mp.nstr(x, 6)}"] = _scaled(
                    extension_raw[name], expected["increments"][name]
                )

            # The reference branch keeps the inherited axial slope exactly;
            # this is checked independently of the moment primitives below.
            integral_errors[f"inherited_slope:Uz_Z@x={mp.nstr(x, 6)}"] = _scaled(
                _field(output, "Uz_Z", "UZ"),
                _field(reshape_start, "Uz_Z", "UZ"),
            )

            fd_values: dict[str, dict[int, mp.mpf]] = {
                name: {} for name in ("F", "Utheta", "Uz", "P", "P0")
            }
            fd_moments: dict[str, dict[int, mp.mpf]] = {
                name: {} for name in MOMENT_NAMES
            }
            fd_increments: dict[str, dict[int, mp.mpf]] = {
                name: {} for name in MOMENT_NAMES
            }
            for index in (-2, -1, 1, 2):
                shifted = extension.evaluate_log_offset(x, Z0 + index * step)
                for name in fd_values:
                    if name == "P0":
                        fd_values[name][index] = _scalar(_field(shifted, "P0", "axis_P0"))
                    else:
                        fd_values[name][index] = _scalar(_field(shifted, name))
                for name in MOMENT_NAMES:
                    fd_moments[name][index] = _scalar(_moments(shifted, derivative=False)[name])
                    fd_increments[name][index] = _scalar(
                        _extension_increments(shifted)[name]
                    )
            for name in fd_values:
                expected_tangent = _fourth_difference(fd_values[name], step)
                actual_tangent = _scalar(
                    _field(
                        output,
                        f"{name}_Z",
                        f"{name}Z",
                    )
                )
                tangent_errors[f"field:{name}@x={mp.nstr(x, 6)}"] = abs(
                    actual_tangent - expected_tangent
                ) / max(mp.mpf(1), abs(actual_tangent), abs(expected_tangent))
            for name in MOMENT_NAMES:
                expected_tangent = _fourth_difference(fd_moments[name], step)
                actual_tangent = _scalar(_moments(output, derivative=True)[name])
                tangent_errors[f"moment:{name}@x={mp.nstr(x, 6)}"] = abs(
                    actual_tangent - expected_tangent
                ) / max(mp.mpf(1), abs(actual_tangent), abs(expected_tangent))
                expected_tangent = _fourth_difference(fd_increments[name], step)
                actual_tangent = _scalar(_extension_increment_tangents(output)[name])
                tangent_errors[f"receipt:moment:{name}@x={mp.nstr(x, 6)}"] = abs(
                    actual_tangent - expected_tangent
                ) / max(mp.mpf(1), abs(actual_tangent), abs(expected_tangent))

            # Raw quadratic extension receipts carry the same first axial
            # tangent convention as the moment receipt.
            for name in RAW_NAMES:
                fd_raw: dict[int, mp.mpf] = {}
                for index in (-2, -1, 1, 2):
                    shifted = extension.evaluate_log_offset(x, Z0 + index * step)
                    fd_raw[index] = _scalar(_extension_raw_increments(shifted)[name])
                expected_tangent = _fourth_difference(fd_raw, step)
                raw_parts_Z = output.get("raw_quadratic_parts_Z")
                if isinstance(raw_parts_Z, Mapping) and isinstance(raw_parts_Z.get("reference"), Mapping):
                    actual_tangent = _scalar(raw_parts_Z["reference"][name])
                else:
                    actual_tangent = _scalar(
                        _receipt(
                            output,
                            "reference_raw_quadratic_increments_Z",
                            "extension_raw_quadratic_increments_Z",
                        )[name]
                    )
                tangent_errors[f"receipt:raw:{name}@x={mp.nstr(x, 6)}"] = abs(
                    actual_tangent - expected_tangent
                ) / max(mp.mpf(1), abs(actual_tangent), abs(expected_tangent))

            rows.append(
                {
                    "x": mp.nstr(x, 12),
                    "R": mp.nstr(_scalar(_field(output, "R")), 24),
                    "maximum_field_and_moment_error": mp.nstr(
                        max(
                            (
                                integral_errors[key]
                                for key in integral_errors
                                if f"@x={mp.nstr(x, 6)}" in key
                            ),
                            default=mp.mpf(0),
                        ),
                        24,
                    ),
                    "extension_increments_nonzero": all(
                        abs(_scalar(extension_moments[name])) > RECEIPT_THRESHOLD
                        for name in MOMENT_NAMES
                    ),
                    "raw_extension_increments_nonzero": all(
                        abs(_scalar(extension_raw[name])) > RECEIPT_THRESHOLD
                        for name in RAW_NAMES
                    ),
                }
            )

        # The original endpoint and reshape receipts must remain available in
        # the extension output.  Check the first positive point, where all
        # three receipts are nontrivial.
        sample = extension.evaluate_log_offset(points[0], Z0)
        parts = sample.get("moment_parts")
        if isinstance(parts, Mapping):
            source_receipt = _receipt(parts, "inner_seed", "source", "original")
            reshape_receipt = _receipt(parts, "reshape", "long_reshape")
        else:
            source_receipt = _receipt(
                sample,
                "source_moment_seeds",
                "original_moment_seeds",
                "endpoint_moment_seeds",
                "original_moments",
                "source_seeds",
            )
            reshape_receipt = _receipt(
                sample,
                "reshape_moment_increments",
                "long_reshape_moment_increments",
                "reshape_increments",
                "long_reshape_increments",
            )
        for name in MOMENT_NAMES:
            receipt_errors[f"source_receipt:{name}"] = _scaled(
                source_receipt[name], _moments(source_endpoint)[name]
            )
            receipt_errors[f"reshape_receipt:{name}"] = _scaled(
                reshape_receipt[name],
                _receipt(reshape_start, "moment_increments", "increment_moments")[name],
            )

        maximum_integral = max(integral_errors.values(), default=mp.mpf(0))
        maximum_tangent = max(tangent_errors.values(), default=mp.mpf(0))
        maximum_continuity = max(continuity_errors.values(), default=mp.mpf(0))
        maximum_receipt = max(receipt_errors.values(), default=mp.mpf(0))
        if maximum_integral >= ERROR_THRESHOLD:
            key = max(integral_errors, key=integral_errors.get)
            raise AssertionError(f"reference extension integral mismatch {key}: {maximum_integral}")
        if maximum_tangent >= ERROR_THRESHOLD:
            key = max(tangent_errors, key=tangent_errors.get)
            raise AssertionError(f"reference extension Z tangent mismatch {key}: {maximum_tangent}")
        if maximum_continuity >= ERROR_THRESHOLD:
            key = max(continuity_errors, key=continuity_errors.get)
            raise AssertionError(f"reference extension x=0 mismatch {key}: {maximum_continuity}")
        if maximum_receipt >= ERROR_THRESHOLD:
            key = max(receipt_errors, key=receipt_errors.get)
            raise AssertionError(f"reference extension receipt mismatch {key}: {maximum_receipt}")

        report = {
            "passed": True,
            "precision": precision,
            "pressure_order": PRESSURE_ORDER,
            "width_order": WIDTH_ORDER,
            "reshape_A": mp.nstr(RESHAPE_A, 30),
            "reshape_T": mp.nstr(RESHAPE_T, 30),
            "logC": mp.nstr(LOG_C, 30),
            "logRref": mp.nstr(LOG_RREF, 30),
            "extension_length": mp.nstr(EXTENSION_LENGTH, 30),
            "points": [mp.nstr(point, 20) for point in points],
            "finite_difference_step": mp.nstr(step, 20),
            "maximum_scaled_integral_error": mp.nstr(maximum_integral, 30),
            "maximum_scaled_Z_tangent_error": mp.nstr(maximum_tangent, 30),
            "maximum_scaled_x0_continuity_error": mp.nstr(maximum_continuity, 30),
            "maximum_scaled_receipt_error": mp.nstr(maximum_receipt, 30),
            "rows": rows,
            "integral_errors": {key: mp.nstr(value, 24) for key, value in sorted(integral_errors.items())},
            "Z_tangent_errors": {key: mp.nstr(value, 24) for key, value in sorted(tangent_errors.items())},
            "continuity_errors": {key: mp.nstr(value, 24) for key, value in sorted(continuity_errors.items())},
            "receipt_errors": {key: mp.nstr(value, 24) for key, value in sorted(receipt_errors.items())},
            "receipt_names": receipt_names,
            "independent_integrals": "mp.quad of physical reference moment derivatives in x=log(R/Rsh)",
            "independent_Z_tangents": "fourth-order centered finite differences of scalar extension evaluations",
            "original_reshape_extension_receipts_separate": True,
            "quadrature_error_enclosed": False,
            "pressure_width_truncation_remainder_enclosed": False,
            "C2_claim": False,
            "cone_certified": False,
            "global_field_installed": False,
            "metadata": extension.metadata() if hasattr(extension, "metadata") else {},
        }
        Path(__file__).with_suffix(".json").write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return report


def fixture(**kwargs: Any) -> dict[str, Any]:
    return run_fixture(**kwargs)


def main() -> None:
    print(json.dumps(run_fixture(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()


__all__ = ["fixture", "run_fixture"]
