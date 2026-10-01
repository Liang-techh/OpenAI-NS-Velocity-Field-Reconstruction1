"""Directed finite-core interval evaluation at the collar inlet.

This module consumes the trusted finite radial/Z core pickle at ``Z=.3``.  It
does not rebuild the collar or its transition ODE.  Every pressure atom is
therefore a directed interval enclosure of the stored finite coefficient;
source, radial-series, and continuation remainders remain outside its scope.

The cached radial rows are Taylor rows in the axial coordinate.  For each
radial coefficient we retain the three available axial slots
``row[k]``, ``(k+1) row[k+1]``, and ``(k+1)(k+2) row[k+2]``.  The first two
slots are supplied to the shared component-core evaluator, while the third
slot makes the returned value jet's second axial derivative explicit.
"""

from __future__ import annotations

import hashlib
import json
import pickle
from pathlib import Path
from typing import Any, Mapping

import mpmath as mp

try:  # package import
    from .lei_ren_part1_paper_component_pressure_core import (
        evaluate_component_core_coefficients,
    )
    from .lei_ren_part1_paper_interval_axial_second_jet import (
        IntervalAxialSecondJet,
    )
    from .lei_ren_part1_paper_interval_difference import IntervalDifference
    from .lei_ren_part1_paper_interval_pressure_width_jet import (
        IntervalPressureWidthJet,
    )
except (ImportError, ValueError):  # direct experiment execution
    from lei_ren_part1_paper_component_pressure_core import (
        evaluate_component_core_coefficients,
    )
    from lei_ren_part1_paper_interval_axial_second_jet import (
        IntervalAxialSecondJet,
    )
    from lei_ren_part1_paper_interval_difference import IntervalDifference
    from lei_ren_part1_paper_interval_pressure_width_jet import (
        IntervalPressureWidthJet,
    )


MINIMUM_PRECISION = 473
PRESSURE_ORDER = 9
WIDTH_ORDER = 2
LOCAL_PICKLE_CACHE = (
    Path(__file__).resolve().parents[2]
    / "work_paper_cache"
    / "second_axial_core_Z03.pkl"
)
DEFAULT_CACHE = Path(__file__).with_name('lei_ren_part1_paper_interval_collar_core_inlet_input.json')
if not DEFAULT_CACHE.exists():
    DEFAULT_CACHE = LOCAL_PICKLE_CACHE
MOMENT_KEYS = ("theta", "z", "theta_z", "z_theta", "p")


def _endpoints(value: Any) -> tuple[mp.mpf, mp.mpf]:
    return tuple(mp.make_mpf(item) for item in value._mpi_)


def _interval_text(value: Any, digits: int = 30) -> dict[str, Any]:
    lower, upper = _endpoints(value)
    return {
        "lower": mp.nstr(lower, digits),
        "upper": mp.nstr(upper, digits),
        # Display strings are human-readable; exact MP tuples make the
        # receipt replayable even for the extremely small stored F atoms.
        "lower_exact_mpf_tuple": list(value._mpi_[0]),
        "upper_exact_mpf_tuple": list(value._mpi_[1]),
    }


def _pair(ctx: Any, value: Any) -> IntervalDifference:
    """Convert one stored MP coefficient without global-context coercion."""

    # ``ctx.mpf(mp.mpf)`` is intentional: the source coefficient is already
    # an arbitrary-precision stored value, and the interval context performs
    # the only rounding used by this evaluator.
    nominal = ctx.mpf(value)
    return IntervalDifference(ctx, nominal, ctx.mpf(0))


def _ring(ctx: Any, value: Any) -> IntervalPressureWidthJet:
    return IntervalPressureWidthJet(
        ctx,
        value,
        pressure_order=PRESSURE_ORDER,
        width_order=WIDTH_ORDER,
    )


def _zero_ring(ctx: Any) -> IntervalPressureWidthJet:
    return _ring(ctx, 0)


def _axial_constant(ctx: Any, value: Any) -> IntervalAxialSecondJet:
    ring = _ring(ctx, value)
    zero = _zero_ring(ctx)
    return IntervalAxialSecondJet._from_slots(ctx, ring, zero, zero)


def _row_atom(ctx: Any, polynomial: Any) -> IntervalPressureWidthJet:
    """Convert one cached ``PressurePolynomial`` to pressure atoms."""

    atoms = {
        (int(power), 0): _pair(ctx, coefficient)
        for power, coefficient in polynomial.atoms.items()
    }
    return IntervalPressureWidthJet._from_atoms(
        ctx,
        atoms,
        PRESSURE_ORDER,
        WIDTH_ORDER,
    )


def _axial_row(ctx: Any, row: list[Any]) -> list[IntervalAxialSecondJet]:
    """Build value, first, and second axial slots for one radial row."""

    values = [_row_atom(ctx, item) for item in row]
    zero = _zero_ring(ctx)
    output: list[IntervalAxialSecondJet] = []
    for k in range(len(values)):
        first = (k + 1) * values[k + 1] if k + 1 < len(values) else zero
        second = (
            (k + 1) * (k + 2) * values[k + 2]
            if k + 2 < len(values)
            else zero
        )
        output.append(IntervalAxialSecondJet._from_slots(
            ctx,
            values[k],
            first,
            second,
        ))
    return output


def _interval_coefficients(
    ctx: Any, raw: Mapping[str, Any]
) -> dict[str, Any]:
    converted = dict(raw)
    for name in ("F", "Uz", "P"):
        converted[name] = [_axial_row(ctx, list(row)) for row in raw[name]]
    # evaluate_core_jets uses this only to select a global work precision;
    # all coefficient arithmetic remains in the supplied interval context.
    converted["precision"] = int(raw.get("precision", MINIMUM_PRECISION))
    return converted


def _inertial_and_radial_derivatives(
    snapshot: Mapping[str, Any],
    radius: IntervalAxialSecondJet,
    z: IntervalAxialSecondJet,
    delta: IntervalAxialSecondJet,
) -> dict[str, IntervalAxialSecondJet]:
    """Evaluate the algebraic inlet inertial formulas and their radial slope."""

    one = _axial_constant(radius.ctx, 1)
    two = _axial_constant(radius.ctx, 2)
    d = one - z * z
    L = one - delta * z * z
    root = (two * radius).sqrt()
    f = snapshot["F"]
    f_r = snapshot["F_R"]
    f_z = snapshot["F_Z"]
    u = snapshot["Uz"]
    u_r = snapshot["Uz_R"]
    u_z = snapshot["Uz_Z"]
    moments = snapshot["moments"]
    moments_z = snapshot["moments_Z"]
    pressure = snapshot["P"]
    pressure_z = snapshot["P_Z"]

    # This is the same algebra as evaluate_mp_stress; no MP scalar conversion
    # or alternate source is used here.
    transport = -radius + (one - delta) * z * moments["z"] + d * moments_z["z"]
    utheta = root * f
    i_theta = utheta * transport / (L * root)
    i_theta += (
        (one - delta / 2) * moments["theta"]
        - (one - delta) * z * moments_z["theta"] / 2
        - d * moments_z["theta_z"]
        + (2 * delta - one) * z * moments["theta_z"]
    ) / (2 * L * radius)

    i_z = (
        transport * u
        + (one - delta) * (moments["z"] - z * moments_z["z"]) / 2
        + 2 * delta * z * moments["z_theta"]
        - d * moments_z["z_theta"]
        + radius * (2 * (one + delta) * z * pressure - d * pressure_z)
    ) / (L * root)
    d_value = i_theta / f

    # The radial derivatives are the explicit differentiations needed by the
    # dynamic second-width inlet transfer.  They use only the same snapshot.
    v_r = -one + (one - delta) * z * u + d * u_z
    c = (
        (one - delta / 2) * moments["theta"]
        - (one - delta) * z * moments_z["theta"] / 2
        - d * moments_z["theta_z"]
        + (2 * delta - one) * z * moments["theta_z"]
    )
    c_r = (
        (one - delta / 2) * 2 * radius * f
        - (one - delta) * z * radius * f_z
        - d * 2 * radius * (f_z * u + f * u_z)
        + (2 * delta - one) * z * 2 * radius * f * u
    )
    # utheta/root=F in the inertial transport term; its radial derivative
    # is (F_R*transport + F*transport_R)/L.
    i_theta_r = (f_r * transport + f * v_r) / L
    i_theta_r += c_r / (2 * L * radius) - c / (2 * L * radius**2)
    d_r = i_theta_r / f - i_theta * f_r / f**2

    q = 2 * (one + delta) * z * pressure - d * pressure_z
    n_r = (
        v_r * u
        + transport * u_r
        + (one - delta) * (u - z * u_z) / 2
        + 2 * delta * z * (u**2 - radius * f**2)
        - d * (2 * u * u_z - 2 * radius * f * f_z)
        + q
        + radius * (2 * (one + delta) * z * f**2 - 2 * d * f * f_z)
    )
    i_z_r = n_r / (L * root) - i_z / (2 * radius)

    return {
        "I_theta": i_theta,
        "I_z": i_z,
        "D": d_value,
        "I_theta_R": i_theta_r,
        "I_z_R": i_z_r,
        "D_R": d_r,
    }


def _read_cache(cache_path: Any) -> tuple[dict[str, Any], Path, str]:
    path = Path(cache_path).resolve()
    if path.suffix=='.json':
        from lei_ren_part1_paper_component_pressure_core import PressurePolynomial
        stored=json.loads(path.read_text())
        if stored.get('format')!='exact-finite-core-atoms-v1':
            raise ValueError('Unsupported finite-core JSON format')
        raw=dict(version=1,Z=stored['Z'],Z_depth=stored['Z_depth'],
                 Lambda=mp.make_mpf(tuple(stored['Lambda_exact_mpf_tuple'])),
                 delta=mp.make_mpf(tuple(stored['delta_exact_mpf_tuple'])),
                 source_pickle_sha256=stored['source_pickle_sha256'],
                 coefficients=dict(precision=stored['precision']))
        # PressurePolynomial's constructor uses the global MP context. Keep
        # enough bits even when this public loader is called at default dps.
        with mp.workdps(max(mp.mp.dps,int(stored['precision'])+20)):
            for name in ('F','Uz','P'):
                raw['coefficients'][name]=[
                    [PressurePolynomial({int(p):mp.make_mpf(tuple(t)) for p,t in atom.items()}) for atom in row]
                    for row in stored['coefficients'][name]]
    else:
        # Only trusted locally generated pickles are supported as the optional
        # legacy input. The committed JSON route performs no deserialization.
        with path.open("rb") as handle:
            raw = pickle.load(handle)
    if not isinstance(raw, Mapping):
        raise TypeError("finite core cache must contain a mapping")
    if int(raw.get("version", 0)) != 1:
        raise ValueError("unsupported finite core cache version")
    if str(raw.get("Z")) != ".3":
        raise ValueError("inlet evaluator requires the stored Z=.3 core")
    coefficients = raw.get("coefficients")
    if not isinstance(coefficients, Mapping):
        raise ValueError("finite core cache has no coefficient mapping")
    if any(name not in coefficients for name in ("F", "Uz", "P")):
        raise ValueError("finite core cache lacks F, Uz, or P rows")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return dict(raw), path, digest


def build_inlet(
    cache_path: Any = DEFAULT_CACHE,
    *,
    precision: int = MINIMUM_PRECISION,
) -> dict[str, Any]:
    """Return the directed finite-core snapshot at ``R=4/Lambda, Z=.3``.

    The return value is a mapping so callers can directly consume ``Uz``,
    ``R``, ``I_z``, ``D``, and their radial derivatives.  ``ctx`` is included
    for downstream interval operations; it is the same context used for every
    field and coefficient atom.
    """

    precision = int(precision)
    if precision < MINIMUM_PRECISION:
        raise ValueError(f"precision must be at least {MINIMUM_PRECISION}")
    raw, path, digest = _read_cache(cache_path)
    ctx = mp.ctx_iv.MPIntervalContext()
    ctx.dps = precision
    coefficients = _interval_coefficients(ctx, raw["coefficients"])
    lambda_interval = ctx.mpf(raw["Lambda"])
    delta_interval = ctx.mpf(raw["delta"])
    radius = _axial_constant(ctx, _ring(ctx, 4) / _ring(ctx, lambda_interval))
    z = IntervalAxialSecondJet._from_slots(
        ctx,
        _ring(ctx, ctx.mpf(".3")),
        _ring(ctx, 1),
        _zero_ring(ctx),
    )
    delta = _axial_constant(ctx, delta_interval)
    fields = evaluate_component_core_coefficients(
        coefficients,
        radius,
        z,
        delta,
        square_root=lambda value: value.sqrt(),
    )
    # evaluate_core_jets returns the radial derivatives for each profile.  The
    # component evaluator also restores the cached P0 exactly before moments.
    snapshot = {
        key: fields[key]
        for key in (
            "F",
            "Uz",
            "F_R",
            "Uz_R",
            "F_Z",
            "Uz_Z",
            "moments",
            "moments_Z",
            "P",
            "P_Z",
        )
    }
    snapshot.update(_inertial_and_radial_derivatives(fields, radius, z, delta))
    snapshot.update(
        {
            "R": radius,
            "z": z,
            "P0": coefficients["P"][0][0],
            "ctx": ctx,
            "Lambda": lambda_interval,
            "delta": delta_interval,
            "cache_path": str(path),
            "cache_sha256": digest,
            "source_pickle_sha256": raw.get('source_pickle_sha256',digest),
            "pressure_order": PRESSURE_ORDER,
            "width_order": WIDTH_ORDER,
            "precision": precision,
            "scope": "stored finite radial/Z core only at Z=.3; no collar/transition RK rebuild",
            "stored_finite_core_only": True,
            "source_remainder_enclosed": False,
            "core_remainder_enclosed": False,
            "collar_RK_rebuilt": False,
            "temporal_recursion": False,
        }
    )
    return snapshot


def _encode_ring(ring: IntervalPressureWidthJet) -> dict[str, Any]:
    return {
        f"{p},{w}": {
            "nominal": _interval_text(coefficient.nominal),
            "difference": _interval_text(coefficient.difference),
        }
        for (p, w), coefficient in sorted(ring.atoms.items())
    }


def encode_jet(jet: IntervalAxialSecondJet) -> dict[str, Any]:
    """JSON-safe receipt for every nominal/perturbation atom."""

    return {
        "value": _encode_ring(jet.value),
        "tangent": _encode_ring(jet.tangent),
        "second": _encode_ring(jet.second),
    }


def encode_snapshot(snapshot: Mapping[str, Any]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in snapshot.items():
        if isinstance(value, IntervalAxialSecondJet):
            result[key] = encode_jet(value)
        elif isinstance(value, IntervalPressureWidthJet):
            result[key] = _encode_ring(value)
        elif isinstance(value, Mapping) and value and all(
            isinstance(item, IntervalAxialSecondJet) for item in value.values()
        ):
            result[key] = {name: encode_jet(item) for name, item in value.items()}
        elif isinstance(value, Mapping) and value and all(
            isinstance(item, IntervalPressureWidthJet) for item in value.values()
        ):
            result[key] = {name: _encode_ring(item) for name, item in value.items()}
        elif hasattr(value, "_mpi_"):
            result[key] = _interval_text(value)
        elif key != "ctx":
            result[key] = value
    return result


__all__ = [
    "MINIMUM_PRECISION",
    "PRESSURE_ORDER",
    "WIDTH_ORDER",
    "DEFAULT_CACHE",
    "build_inlet",
    "encode_jet",
    "encode_snapshot",
]
