"""Independent resolved check for the pressure/width axial restoration.

The production adapter starts at the reference-extension endpoint ``Rz`` and
uses ``t = log(R/Rz)`` on ``0 <= t <= 3``.  This fixture replays the physical
moment derivatives with scalar ``mp.quad``: the first logarithmic unit uses
the flat interpolation from ``V0`` to ``4 Z`` and the remaining two units use
``V = 4 Z``.  First ``Z`` derivatives are checked by fourth-order centred
differences of independent scalar evaluations.

The inherited inner, reshape, and reference receipts are checked separately
from the new axial-restore and post-restore-reference increments.
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
    from .lei_ren_part1_paper_pressure_width_axial_restore import (
        PressureWidthAxialRestore,
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
    from lei_ren_part1_paper_pressure_width_axial_restore import (
        PressureWidthAxialRestore,
    )


PRECISION = 100
PRESSURE_ORDER = 3
WIDTH_ORDER = 2
RESHAPE_A = mp.mpf(1)
RESHAPE_T = mp.mpf(400)
LOG_C = mp.mpf(-4)
LOG_RREF = mp.log(110) + RESHAPE_T + 8 + 5
QUADRATURE_ORDER = 24
TAIL_DIGITS = 50
Z0 = mp.mpf(".3")
RESTORE_POINTS = (mp.mpf(".5"), mp.mpf(1), mp.mpf(3))
FINITE_DIFFERENCE_STEP = mp.mpf("1e-4")
ERROR_THRESHOLD = mp.mpf("1e-8")
RECEIPT_THRESHOLD = mp.mpf("1e-30")

MOMENT_NAMES = ("theta", "z", "theta_z", "z_theta", "p")
RAW_NAMES = ("axial", "swirl")


def _scalar(value: Any) -> mp.mpf:
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
    raise KeyError("missing moment mapping")


def _parts(source: Mapping[str, Any], *, derivative: bool = False) -> Mapping[str, Any]:
    names = ("moment_parts_Z", "moment_partsZ") if derivative else ("moment_parts",)
    for name in names:
        value = source.get(name)
        if isinstance(value, Mapping):
            return value
    raise KeyError("missing moment parts")


def _raw_parts(source: Mapping[str, Any], *, derivative: bool = False) -> Mapping[str, Any]:
    names = (
        ("raw_quadratic_parts_Z", "raw_quadratic_partsZ")
        if derivative
        else ("raw_quadratic_parts",)
    )
    for name in names:
        value = source.get(name)
        if isinstance(value, Mapping):
            return value
    raise KeyError("missing raw quadratic parts")


def _fourth_difference(values: Mapping[int, mp.mpf], step: mp.mpf) -> mp.mpf:
    return (-values[2] + 8 * values[1] - 8 * values[-1] + values[-2]) / (12 * step)


def _sigma(s: mp.mpf) -> mp.mpf:
    if s <= 0:
        return mp.mpf(0)
    if s >= 1:
        return mp.mpf(1)
    f = lambda value: mp.exp(-1 / value**2) if value > 0 else mp.mpf(0)
    return f(s) / (f(s) + f(1 - s))


def _sigma_prime(s: mp.mpf) -> mp.mpf:
    sig = _sigma(s)
    if 0 < s < 1:
        return sig * (1 - sig) * (2 / s**3 + 2 / (1 - s) ** 3)
    return mp.mpf(0)


def _quad(function: Any, endpoint: mp.mpf) -> mp.mpf:
    if endpoint == 0:
        return mp.mpf(0)
    if endpoint <= 1:
        return mp.quad(function, [0, endpoint / 2, endpoint])
    return mp.quad(function, [0, mp.mpf(1) / 2, 1, (endpoint + 1) / 2, endpoint])


def _start(reference: PressureWidthReferenceExtension, z: mp.mpf) -> dict[str, Any]:
    return reference.evaluate_log_offset(reference.extension_length, z)


def _reference_parts(reference_start: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    parts = _parts(reference_start)
    value = parts.get(name)
    if not isinstance(value, Mapping):
        raise KeyError(f"reference endpoint lacks moment part {name!r}")
    return value


def _independent(
    reference: PressureWidthReferenceExtension,
    t: mp.mpf,
    z: mp.mpf,
) -> dict[str, Any]:
    """Scalar physical replay from Rz through a restore phase ``t``."""

    start = _start(reference, z)
    Rz = _scalar(_field(start, "R"))
    u0 = _scalar(_field(start, "Utheta"))
    V0 = _scalar(_field(start, "Uz"))
    V0_Z = _scalar(_field(start, "Uz_Z", "UZ"))
    delta_v = 4 * z - V0
    delta_v_Z = 4 - V0_Z

    def radius(s: mp.mpf) -> mp.mpf:
        return Rz * mp.exp(s)

    def theta(s: mp.mpf) -> mp.mpf:
        return u0 * mp.exp(s / 10)

    def velocity(s: mp.mpf) -> mp.mpf:
        return V0 + delta_v * _sigma(s) if s <= 1 else 4 * z

    def dtheta(s: mp.mpf) -> mp.mpf:
        R = radius(s)
        return mp.sqrt(2) * R**mp.mpf("1.5") * theta(s)

    def dpressure(s: mp.mpf) -> mp.mpf:
        return theta(s) ** 2 / 2

    def dswirl(s: mp.mpf) -> mp.mpf:
        return radius(s) * theta(s) ** 2 / 2

    def dz(s: mp.mpf) -> mp.mpf:
        return radius(s) * velocity(s)

    def daxial(s: mp.mpf) -> mp.mpf:
        return radius(s) * velocity(s) ** 2

    def dmixed(s: mp.mpf) -> mp.mpf:
        return velocity(s) * dtheta(s)

    increments = {
        "theta": _quad(dtheta, t),
        "z": _quad(dz, t),
        "theta_z": _quad(dmixed, t),
        "z_theta": _quad(daxial, t) - _quad(dswirl, t),
        "p": _quad(dpressure, t),
    }
    raw = {
        "axial": _quad(daxial, t),
        "swirl": _quad(dswirl, t),
    }
    R = radius(t)
    u = theta(t)
    V = velocity(t)
    moments = {name: _scalar(_moments(start)[name]) for name in MOMENT_NAMES}
    moments.update({name: moments[name] + increments[name] for name in MOMENT_NAMES})
    return {
        "R": R,
        "F": u / mp.sqrt(2 * R),
        "Utheta": u,
        "Uz": V,
        "Uz_Z": V0_Z + delta_v_Z * _sigma(t) if t <= 1 else mp.mpf(4),
        "Uz_y": delta_v * _sigma_prime(t) if t < 1 else mp.mpf(0),
        "P": _scalar(_field(start, "P")) + increments["p"],
        "P0": _scalar(_field(start, "P0", "axis_P0")),
        "moments": moments,
        "increments": increments,
        "raw": raw,
        "start": start,
    }


def _new_increments(output: Mapping[str, Any]) -> Mapping[str, Any]:
    axial = output.get("axial_restore_increments")
    restored = output.get("restored_reference_increments")
    if isinstance(axial, Mapping) and isinstance(restored, Mapping):
        return {
            key: _scalar(axial[key]) + _scalar(restored[key])
            for key in MOMENT_NAMES
        }
    for name in (
        "restore_moment_increments",
        "axial_restore_moment_increments",
        "new_moment_increments",
        "moment_increments",
    ):
        value = output.get(name)
        if isinstance(value, Mapping) and all(key in value for key in MOMENT_NAMES):
            return value
    parts = _parts(output)
    # Some implementations expose only stage receipts.  Their sum is the
    # accumulated increment from Rz and is the quantity checked here.
    stage_values = [parts[name] for name in ("axial_restore", "restored_reference") if isinstance(parts.get(name), Mapping)]
    if stage_values:
        return {
            key: mp.fsum((_scalar(stage[key]) for stage in stage_values))
            for key in MOMENT_NAMES
        }
    raise KeyError("missing new axial-restore moment increment mapping")


def _new_raw_increments(output: Mapping[str, Any]) -> Mapping[str, Any]:
    axial = output.get("axial_restore_raw_increments")
    restored = output.get("restored_reference_raw_increments")
    if isinstance(axial, Mapping) and isinstance(restored, Mapping):
        return {
            key: _scalar(axial[key]) + _scalar(restored[key])
            for key in RAW_NAMES
        }
    for name in (
        "restore_raw_quadratic_increments",
        "axial_restore_raw_quadratic_increments",
        "new_raw_quadratic_increments",
        "raw_quadratic_increments",
    ):
        value = output.get(name)
        if isinstance(value, Mapping) and all(key in value for key in RAW_NAMES):
            return value
    parts = _raw_parts(output)
    stage_values = [parts[name] for name in ("axial_restore", "restored_reference") if isinstance(parts.get(name), Mapping)]
    if stage_values:
        return {
            key: mp.fsum((_scalar(stage[key]) for stage in stage_values))
            for key in RAW_NAMES
        }
    raise KeyError("missing new axial-restore raw increment mapping")


def _new_increment_tangents(output: Mapping[str, Any]) -> Mapping[str, Any]:
    axial = output.get("axial_restore_increments_Z")
    restored = output.get("restored_reference_increments_Z")
    if isinstance(axial, Mapping) and isinstance(restored, Mapping):
        return {
            key: _scalar(axial[key]) + _scalar(restored[key])
            for key in MOMENT_NAMES
        }
    for name in (
        "restore_moment_increments_Z",
        "axial_restore_moment_increments_Z",
        "new_moment_increments_Z",
        "moment_increments_Z",
    ):
        value = output.get(name)
        if isinstance(value, Mapping) and all(key in value for key in MOMENT_NAMES):
            return value
    parts = _parts(output, derivative=True)
    stage_values = [parts[name] for name in ("axial_restore", "restored_reference") if isinstance(parts.get(name), Mapping)]
    if stage_values:
        return {
            key: mp.fsum((_scalar(stage[key]) for stage in stage_values))
            for key in MOMENT_NAMES
        }
    raise KeyError("missing new axial-restore moment tangent mapping")


def _new_raw_increment_tangents(output: Mapping[str, Any]) -> Mapping[str, Any]:
    axial = output.get("axial_restore_raw_increments_Z")
    restored = output.get("restored_reference_raw_increments_Z")
    if isinstance(axial, Mapping) and isinstance(restored, Mapping):
        return {
            key: _scalar(axial[key]) + _scalar(restored[key])
            for key in RAW_NAMES
        }
    for name in (
        "restore_raw_quadratic_increments_Z",
        "axial_restore_raw_quadratic_increments_Z",
        "new_raw_quadratic_increments_Z",
        "raw_quadratic_increments_Z",
    ):
        value = output.get(name)
        if isinstance(value, Mapping) and all(key in value for key in RAW_NAMES):
            return value
    parts = _raw_parts(output, derivative=True)
    stage_values = [parts[name] for name in ("axial_restore", "restored_reference") if isinstance(parts.get(name), Mapping)]
    if stage_values:
        return {
            key: mp.fsum((_scalar(stage[key]) for stage in stage_values))
            for key in RAW_NAMES
        }
    raise KeyError("missing new axial-restore raw tangent mapping")


def _difference_maps(
    right: Mapping[str, Any],
    left: Mapping[str, Any],
    names: tuple[str, ...],
) -> dict[str, mp.mpf]:
    return {
        name: _scalar(right[name]) - _scalar(left[name])
        for name in names
    }


def _build_objects(precision: int) -> tuple[Any, PressureWidthLongReshape, PressureWidthReferenceExtension, Any]:
    comparison = _build_comparison(precision)
    switches = _ResolvedR110Provider(comparison)
    reshape = PressureWidthLongReshape(
        switches,
        A=RESHAPE_A,
        logC=LOG_C,
        order=QUADRATURE_ORDER,
        tail_digits=TAIL_DIGITS,
    )
    reference = PressureWidthReferenceExtension(reshape, logRref=LOG_RREF)
    restore = PressureWidthAxialRestore(reference, order=48)
    return switches, reshape, reference, restore


def run_fixture(
    *,
    precision: int = PRECISION,
    points: tuple[Any, ...] = RESTORE_POINTS,
    finite_difference_step: Any = FINITE_DIFFERENCE_STEP,
) -> dict[str, Any]:
    precision = max(80, int(precision))
    step = mp.mpf(str(finite_difference_step))
    if step <= 0:
        raise ValueError("finite_difference_step must be positive")
    with mp.workdps(precision):
        switches, reshape, reference, restore = _build_objects(precision)
        points = tuple(mp.mpf(str(point)) for point in points)
        if not points or any(point < 0 or point > 3 for point in points):
            raise ValueError("restore points must lie in [0, 3]")

        reference_end = _start(reference, Z0)
        zero = restore.evaluate_phase(0, Z0)
        continuity_errors = {
            name: _scaled(_field(zero, name), _field(reference_end, name))
            for name in ("R", "F", "Utheta", "Uz", "P", "P0")
        }
        continuity_errors.update(
            {
                f"field_Z:{name}": _scaled(
                    _field(zero, f"{name}_Z", f"{name}Z"),
                    _field(reference_end, f"{name}_Z", f"{name}Z"),
                )
                for name in ("F", "Utheta", "Uz", "P", "P0")
            }
        )
        continuity_errors.update(
            {
                f"moment:{name}": _scaled(_moments(zero)[name], _moments(reference_end)[name])
                for name in MOMENT_NAMES
            }
        )
        continuity_errors.update(
            {
                f"prior_part:{part}:{name}": _scaled(
                    _parts(zero)[part][name], _parts(reference_end)[part][name]
                )
                for part in ("inner_seed", "reshape", "reference")
                for name in MOMENT_NAMES
            }
        )
        continuity_errors.update(
            {
                f"prior_raw_part:{part}:{name}": _scaled(
                    _raw_parts(zero)[part][name], _raw_parts(reference_end)[part][name]
                )
                for part in ("inner_seed", "reshape", "reference")
                for name in RAW_NAMES
            }
        )

        integral_errors: dict[str, mp.mpf] = {}
        tangent_errors: dict[str, mp.mpf] = {}
        receipt_errors: dict[str, mp.mpf] = {}
        rows: list[dict[str, Any]] = []
        for t in points:
            output = restore.evaluate_phase(t, Z0)
            expected = _independent(reference, t, Z0)
            suffix = f"@t={mp.nstr(t, 6)}"
            for name in ("R", "F", "Utheta", "Uz", "P", "P0"):
                integral_errors[f"field:{name}{suffix}"] = _scaled(
                    _field(output, name), expected[name]
                )
            for name in MOMENT_NAMES:
                integral_errors[f"moment:{name}{suffix}"] = _scaled(
                    _moments(output)[name], expected["moments"][name]
                )
            for name in MOMENT_NAMES:
                receipt_errors[f"new_moment:{name}{suffix}"] = _scaled(
                    _new_increments(output)[name], expected["increments"][name]
                )
            for name in RAW_NAMES:
                receipt_errors[f"new_raw:{name}{suffix}"] = _scaled(
                    _new_raw_increments(output)[name], expected["raw"][name]
                )

            # Validate the two new receipt labels independently, then verify
            # that their sum is the accumulated Rz increment used above.
            stage_parts = _parts(output)
            stage_raw_parts = _raw_parts(output)
            if all(isinstance(stage_parts.get(name), Mapping) for name in ("axial_restore", "restored_reference")):
                stage_one = _independent(reference, min(t, mp.mpf(1)), Z0)
                stage_two = (
                    _difference_maps(expected["increments"], stage_one["increments"], MOMENT_NAMES)
                    if t > 1
                    else {name: mp.mpf(0) for name in MOMENT_NAMES}
                )
                for name in MOMENT_NAMES:
                    receipt_errors[f"stage:axial_restore:{name}{suffix}"] = _scaled(
                        stage_parts["axial_restore"][name], stage_one["increments"][name]
                    )
                    receipt_errors[f"stage:restored_reference:{name}{suffix}"] = _scaled(
                        stage_parts["restored_reference"][name], stage_two[name]
                    )
            if all(isinstance(stage_raw_parts.get(name), Mapping) for name in ("axial_restore", "restored_reference")):
                stage_one = _independent(reference, min(t, mp.mpf(1)), Z0)
                stage_two = (
                    _difference_maps(expected["raw"], stage_one["raw"], RAW_NAMES)
                    if t > 1
                    else {name: mp.mpf(0) for name in RAW_NAMES}
                )
                for name in RAW_NAMES:
                    receipt_errors[f"stage_raw:axial_restore:{name}{suffix}"] = _scaled(
                        stage_raw_parts["axial_restore"][name], stage_one["raw"][name]
                    )
                    receipt_errors[f"stage_raw:restored_reference:{name}{suffix}"] = _scaled(
                        stage_raw_parts["restored_reference"][name], stage_two[name]
                    )

            if t in (1, 3):
                integral_errors[f"join:Uz={suffix}"] = _scaled(_field(output, "Uz"), 4 * Z0)
                integral_errors[f"join:Uz_Z={suffix}"] = _scaled(
                    _field(output, "Uz_Z", "UZ"), 4
                )
                integral_errors[f"join:Uz_y={suffix}"] = _scaled(
                    _field(output, "Uz_y", "u_y"), 0
                )

            fd_fields = {name: {} for name in ("R", "F", "Utheta", "Uz", "P", "P0")}
            fd_moments = {name: {} for name in MOMENT_NAMES}
            fd_increment = {name: {} for name in MOMENT_NAMES}
            fd_raw = {name: {} for name in RAW_NAMES}
            for index in (-2, -1, 1, 2):
                shifted = _independent(reference, t, Z0 + index * step)
                for name in fd_fields:
                    fd_fields[name][index] = _scalar(_field(shifted, name))
                for name in MOMENT_NAMES:
                    fd_moments[name][index] = _scalar(_moments(shifted)[name])
                    fd_increment[name][index] = _scalar(shifted["increments"][name])
                for name in RAW_NAMES:
                    fd_raw[name][index] = _scalar(shifted["raw"][name])
            for name in fd_fields:
                actual = _scalar(_field(output, f"{name}_Z", f"{name}Z"))
                # Rz is fixed by the shared logRref input and has exactly
                # zero Z tangent.  A finite difference of its ~1e178 scalar
                # value would measure only decimal cancellation noise.
                expected_tangent = (
                    mp.mpf(0)
                    if name == "R"
                    else _fourth_difference(fd_fields[name], step)
                )
                tangent_errors[f"field:{name}{suffix}"] = abs(actual - expected_tangent) / max(
                    mp.mpf(1), abs(actual), abs(expected_tangent)
                )
            for name in MOMENT_NAMES:
                actual = _scalar(_moments(output, derivative=True)[name])
                expected_tangent = _fourth_difference(fd_moments[name], step)
                tangent_errors[f"moment:{name}{suffix}"] = abs(actual - expected_tangent) / max(
                    mp.mpf(1), abs(actual), abs(expected_tangent)
                )
                actual = _scalar(_new_increment_tangents(output)[name])
                expected_tangent = _fourth_difference(fd_increment[name], step)
                tangent_errors[f"new_moment:{name}{suffix}"] = abs(actual - expected_tangent) / max(
                    mp.mpf(1), abs(actual), abs(expected_tangent)
                )
            for name in RAW_NAMES:
                actual = _scalar(_new_raw_increment_tangents(output)[name])
                expected_tangent = _fourth_difference(fd_raw[name], step)
                tangent_errors[f"new_raw:{name}{suffix}"] = abs(actual - expected_tangent) / max(
                    mp.mpf(1), abs(actual), abs(expected_tangent)
                )

            rows.append(
                {
                    "t": mp.nstr(t, 12),
                    "R": mp.nstr(_scalar(_field(output, "R")), 24),
                    "new_moment_increments_nonzero": all(
                        abs(_scalar(_new_increments(output)[name])) > RECEIPT_THRESHOLD
                        for name in MOMENT_NAMES
                    ),
                    "new_raw_increments_nonzero": all(
                        abs(_scalar(_new_raw_increments(output)[name])) > RECEIPT_THRESHOLD
                        for name in RAW_NAMES
                    ),
                }
            )

        max_integral = max(integral_errors.values(), default=mp.mpf(0))
        max_tangent = max(tangent_errors.values(), default=mp.mpf(0))
        max_continuity = max(continuity_errors.values(), default=mp.mpf(0))
        max_receipt = max(receipt_errors.values(), default=mp.mpf(0))
        if max_integral >= ERROR_THRESHOLD:
            key = max(integral_errors, key=integral_errors.get)
            raise AssertionError(f"axial restore integral mismatch {key}: {max_integral}")
        if max_tangent >= ERROR_THRESHOLD:
            key = max(tangent_errors, key=tangent_errors.get)
            raise AssertionError(f"axial restore Z tangent mismatch {key}: {max_tangent}")
        if max_continuity >= ERROR_THRESHOLD:
            key = max(continuity_errors, key=continuity_errors.get)
            raise AssertionError(f"axial restore join mismatch {key}: {max_continuity}")
        if max_receipt >= ERROR_THRESHOLD:
            key = max(receipt_errors, key=receipt_errors.get)
            raise AssertionError(f"axial restore receipt mismatch {key}: {max_receipt}")

        report = {
            "passed": True,
            "precision": precision,
            "pressure_order": PRESSURE_ORDER,
            "width_order": WIDTH_ORDER,
            "reshape_A": mp.nstr(RESHAPE_A, 30),
            "reshape_T": mp.nstr(RESHAPE_T, 30),
            "logC": mp.nstr(LOG_C, 30),
            "logRref": mp.nstr(LOG_RREF, 30),
            "points": [mp.nstr(point, 20) for point in points],
            "finite_difference_step": mp.nstr(step, 20),
            "maximum_scaled_integral_error": mp.nstr(max_integral, 30),
            "maximum_scaled_Z_tangent_error": mp.nstr(max_tangent, 30),
            "maximum_scaled_join_error": mp.nstr(max_continuity, 30),
            "maximum_scaled_receipt_error": mp.nstr(max_receipt, 30),
            "rows": rows,
            "integral_errors": {key: mp.nstr(value, 24) for key, value in sorted(integral_errors.items())},
            "Z_tangent_errors": {key: mp.nstr(value, 24) for key, value in sorted(tangent_errors.items())},
            "join_errors": {key: mp.nstr(value, 24) for key, value in sorted(continuity_errors.items())},
            "receipt_errors": {key: mp.nstr(value, 24) for key, value in sorted(receipt_errors.items())},
            "independent_integrals": "mp.quad of physical angular, axial, mixed, pressure, and quadratic derivatives",
            "independent_Z_tangents": "fourth-order centered finite differences of independent scalar mp.quad evaluations",
            "prior_three_parts_preserved": True,
            "axis_P0_preserved": True,
            "inherited_Uz_slope_preserved": True,
            "quadrature_error_enclosed": False,
            "pressure_width_truncation_remainder_enclosed": False,
            "C2_claim": False,
            "cone_certified": False,
            "global_field_installed": False,
            "metadata": restore.metadata() if hasattr(restore, "metadata") else {},
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
