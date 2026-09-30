"""Profile-derived Part I stress formulas from Lei--Ren Section 3.4.

The independent profile data are the swirl ``F``, axial velocity ``U`` and
their radial and axial derivatives, the five cumulative moments and their
axial derivatives, and ``P, P_Z``.  This module evaluates the formulas in
(3.9), (3.14)--(3.18) pointwise.  It never fits an inertial stress to sampled
stress values and never invents a moment boundary value.

The canonical profile protocol is::

    F(R, Z), U(R, Z), F_R(R, Z), F_eta(R, Z), U_R(R, Z), U_eta(R, Z)
    moments(R, Z) -> {"theta", "z", "theta_z", "z_theta", "p"}
    moments_Z(R, Z) -> the same five keys
    P(R, Z), P_Z(R, Z)

``F_Z`` and ``U_Z`` are accepted as aliases for the axial derivatives.  The
five moment values are supplied by the caller, so a boundary-target or
axis-origin convention remains explicit in the profile rather than hidden in
this stress evaluator.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
import math
from pathlib import Path
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_VERSION = "2609.35406v1"
SOURCE_SECTIONS = "Section 3.1--3.4, equations (3.5), (3.9), (3.12)--(3.18)"
MOMENT_KEYS = ("theta", "z", "theta_z", "z_theta", "p")


def _scalar(value: Any, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be finite") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _coordinates(R: Any, Z: Any, delta: Any) -> tuple[float, float, float, float]:
    radius = _scalar(R, "R")
    z = _scalar(Z, "Z")
    delta_value = _scalar(delta, "delta")
    if radius <= 0.0:
        raise ValueError("R must be positive")
    if abs(z) > 1.0:
        raise ValueError("|Z| must be at most one")
    if not 0.0 <= delta_value < 1.0:
        raise ValueError("delta must lie in [0,1)")
    d = 1.0 - z * z
    L = 1.0 - delta_value * z * z
    if L <= 0.0:
        raise ValueError("L=1-delta*Z^2 must be positive")
    return radius, z, d, L


def _moment_values(values: Mapping[str, Any] | tuple[Any, ...] | list[Any], name: str) -> dict[str, float]:
    if isinstance(values, Mapping):
        result: dict[str, float] = {}
        for key in MOMENT_KEYS:
            if key not in values:
                raise KeyError(f"{name} is missing moment key {key!r}")
            result[key] = _scalar(values[key], f"{name}[{key}]")
        return result
    if isinstance(values, (tuple, list)) and len(values) == len(MOMENT_KEYS):
        return {
            key: _scalar(value, f"{name}[{key}]")
            for key, value in zip(MOMENT_KEYS, values)
        }
    raise TypeError(
        f"{name} must be a mapping with keys {MOMENT_KEYS} or a five-entry tuple"
    )


def _call_or_value(value: Any, R: float, Z: float, name: str) -> float:
    result = value(R, Z) if callable(value) else value
    return _scalar(result, name)


def evaluate_profile_stress(
    R: Any,
    Z: Any,
    delta: Any,
    *,
    F: Any,
    U: Any,
    F_R: Any,
    F_eta: Any,
    U_R: Any,
    U_eta: Any,
    moments: Mapping[str, Any] | tuple[Any, ...] | list[Any],
    moments_Z: Mapping[str, Any] | tuple[Any, ...] | list[Any],
    P: Any,
    P_Z: Any,
) -> dict[str, float | dict[str, float]]:
    """Evaluate profile-derived ``U^r``, ``N``, ``I``, ``S`` and ``T``.

    All profile quantities may be scalar values or two-argument callables.
    The formulas are the paper's moment representation (3.16)--(3.18), with
    ``U^r`` from (3.9).  ``N_theta`` and ``N_z`` are also returned from (3.12)
    for an independent radial-equation check.
    """

    radius, z, d, L = _coordinates(R, Z, delta)
    delta_value = _scalar(delta, "delta")
    values = {
        "F": _call_or_value(F, radius, z, "F"),
        "U": _call_or_value(U, radius, z, "U"),
        "F_R": _call_or_value(F_R, radius, z, "F_R"),
        "F_eta": _call_or_value(F_eta, radius, z, "F_eta"),
        "U_R": _call_or_value(U_R, radius, z, "U_R"),
        "U_eta": _call_or_value(U_eta, radius, z, "U_eta"),
        "P": _call_or_value(P, radius, z, "P"),
        "P_Z": _call_or_value(P_Z, radius, z, "P_Z"),
    }
    m = _moment_values(moments, "moments")
    mz = _moment_values(moments_Z, "moments_Z")
    root_2R = math.sqrt(2.0 * radius)
    root_R_over_2 = math.sqrt(radius / 2.0)

    U_theta = root_2R * values["F"]
    U_theta_R = root_2R * (values["F_R"] + values["F"] / (2.0 * radius))
    U_theta_Z = root_2R * values["F_eta"]
    U_z = values["U"]
    U_z_R = values["U_R"]
    U_z_Z = values["U_eta"]

    # (3.9), followed by its radial derivative using M^z_R=U^z and
    # (M^z_Z)_R=U^z_Z.
    radial_numerator = (
        2.0 * z * radius * U_z
        - (1.0 - delta_value) * z * m["z"]
        - d * mz["z"]
    )
    U_r = radial_numerator / (L * root_2R)
    radial_numerator_R = (
        (1.0 + delta_value) * z * U_z
        + 2.0 * z * radius * U_z_R
        - d * U_z_Z
    )
    U_r_R = radial_numerator_R / (L * root_2R) - U_r / (2.0 * radius)

    # Moment representation (3.16).
    moment_transport = (
        -radius + (1.0 - delta_value) * z * m["z"] + d * mz["z"]
    )
    theta_moment_bracket = (
        (1.0 - delta_value / 2.0) * m["theta"]
        - (1.0 - delta_value) * z * mz["theta"] / 2.0
        - d * mz["theta_z"]
        + (2.0 * delta_value - 1.0) * z * m["theta_z"]
    )
    I_theta = (
        U_theta * moment_transport / (L * root_2R)
        + theta_moment_bracket / (2.0 * L * radius)
    )

    # Moment representation (3.17), with M^p-M^p(infinity)=P and its
    # axial derivative P_Z.
    z_moment_bracket = (
        moment_transport * U_z
        + (1.0 - delta_value) * (m["z"] - z * mz["z"]) / 2.0
        + 2.0 * delta_value * z * m["z_theta"]
        - d * mz["z_theta"]
        + radius * (2.0 * (1.0 + delta_value) * z * values["P"] - d * values["P_Z"])
    )
    I_z = z_moment_bracket / (L * root_2R)

    # (3.18).
    S_theta = 2.0 * radius * values["F_R"]
    S_z = root_2R * values["U_R"]

    # Source formulas (3.12).  These are returned to make the radial ODE
    # independently testable; they are not used to fit I.
    theta_product_Z = U_z_Z * U_theta + U_z * U_theta_Z
    theta_product_R = U_z_R * U_theta + U_z * U_theta_R
    radial_theta_product_R = U_r_R * U_theta + U_r * U_theta_R
    N_theta = (
        -root_R_over_2
        / L
        * (
            (1.0 + delta_value) * U_theta / 2.0
            + (1.0 - delta_value) * z * U_theta_Z / 2.0
            + radius * U_theta_R
        )
        - root_R_over_2
        / L
        * (
            -2.0 * (1.0 + delta_value) * z * U_z * U_theta
            + d * theta_product_Z
            - 2.0 * z * radius * theta_product_R
        )
        - (radius * radial_theta_product_R + U_r * U_theta)
    )
    axial_product_Z = 2.0 * U_z * U_z_Z + values["P_Z"]
    axial_product_R = 2.0 * U_z * U_z_R + values["F"] ** 2
    radial_axial_product_R = U_r_R * U_z + U_r * U_z_R
    N_z = (
        -root_R_over_2
        / L
        * (
            (1.0 + delta_value) * U_z / 2.0
            + (1.0 - delta_value) * z * U_z_Z / 2.0
            + radius * U_z_R
        )
        - root_R_over_2
        / L
        * (
            -2.0 * (1.0 + delta_value) * z * (U_z**2 + values["P"])
            + d * axial_product_Z
            - 2.0 * z * radius * axial_product_R
        )
        - (radius * radial_axial_product_R + 0.5 * U_r * U_z)
    )

    return {
        "R": radius,
        "Z": z,
        "delta": delta_value,
        "d": d,
        "L": L,
        **values,
        "moments": m,
        "moments_Z": mz,
        "U_theta": U_theta,
        "U_theta_R": U_theta_R,
        "U_theta_Z": U_theta_Z,
        "U_r": U_r,
        "U_r_R": U_r_R,
        "I_theta": I_theta,
        "I_z": I_z,
        "S_theta": S_theta,
        "S_z": S_z,
        "T_theta": I_theta + S_theta,
        "T_z": I_z + S_z,
        "N_theta": N_theta,
        "N_z": N_z,
    }


def _find_method(profile: Any, names: tuple[str, ...], label: str) -> Callable[..., Any]:
    for name in names:
        value = getattr(profile, name, None)
        if callable(value):
            return value
    raise AttributeError(f"profile must provide callable {label}; tried {names}")


def _call_profile_method(profile: Any, names: tuple[str, ...], R: float, Z: float, label: str) -> Any:
    return _find_method(profile, names, label)(R, Z)


class ProfileStress:
    """Evaluate Section 3 stress from one profile's supplied data."""

    def __init__(self, profile: Any, *, delta: float | None = None) -> None:
        self.profile = profile
        if delta is None:
            candidate = getattr(profile, "delta", None)
            if candidate is None and hasattr(profile, "heat"):
                candidate = 2.0 * getattr(profile.heat, "h")
            if candidate is None:
                raise ValueError("delta must be supplied when profile has no delta")
            delta = candidate
        self.delta = _scalar(delta, "delta")
        _coordinates(1.0, 0.0, self.delta)

    def _value(self, names: tuple[str, ...], R: float, Z: float, label: str) -> float:
        return _scalar(_call_profile_method(self.profile, names, R, Z, label), label)

    def F(self, R: float, Z: float) -> float:
        return self._value(("F", "f"), R, Z, "F")

    def U(self, R: float, Z: float) -> float:
        return self._value(("U", "U_z", "Uz", "u_z"), R, Z, "U")

    def F_R(self, R: float, Z: float) -> float:
        return self._value(("F_R", "f_R"), R, Z, "F_R")

    def F_eta(self, R: float, Z: float) -> float:
        return self._value(("F_eta", "F_Z", "f_eta", "f_Z"), R, Z, "F_eta")

    def U_R(self, R: float, Z: float) -> float:
        return self._value(("U_R", "U_r", "Uz_R", "U_z_R"), R, Z, "U_R")

    def U_eta(self, R: float, Z: float) -> float:
        return self._value(("U_eta", "U_Z", "Uz_eta", "U_z_Z"), R, Z, "U_eta")

    def moments(self, R: float, Z: float) -> dict[str, float]:
        method = getattr(self.profile, "moments", None)
        if callable(method):
            return _moment_values(method(R, Z), "moments")
        aliases = {
            "theta": ("M_theta", "Mtheta"),
            "z": ("M_z", "Mz"),
            "theta_z": ("M_theta_z", "Mtheta_z", "Mthetaz"),
            "z_theta": ("M_z_theta", "Mz_theta", "Mztheta"),
            "p": ("M_p", "Mp"),
        }
        values: dict[str, float] = {}
        for key, names in aliases.items():
            values[key] = self._value(names, R, Z, f"M_{key}")
        return values

    def moments_Z(self, R: float, Z: float) -> dict[str, float]:
        method = getattr(self.profile, "moments_Z", None)
        if callable(method):
            return _moment_values(method(R, Z), "moments_Z")
        aliases = {
            "theta": ("M_theta_Z", "Mtheta_Z", "Mtheta_eta"),
            "z": ("M_z_Z", "Mz_Z", "Mz_eta"),
            "theta_z": ("M_theta_z_Z", "Mtheta_z_Z", "Mthetaz_Z"),
            "z_theta": ("M_z_theta_Z", "Mz_theta_Z", "Mztheta_Z"),
            "p": ("M_p_Z", "Mp_Z", "Mp_eta"),
        }
        values: dict[str, float] = {}
        for key, names in aliases.items():
            values[key] = self._value(names, R, Z, f"M_{key}_Z")
        return values

    def P(self, R: float, Z: float) -> float:
        return self._value(("P", "pressure", "p"), R, Z, "P")

    def P_Z(self, R: float, Z: float) -> float:
        return self._value(("P_Z", "P_eta", "pressure_Z", "p_Z"), R, Z, "P_Z")

    def evaluate(self, R: float, Z: float) -> dict[str, float | dict[str, float]]:
        return evaluate_profile_stress(
            R,
            Z,
            self.delta,
            F=self.F(R, Z),
            U=self.U(R, Z),
            F_R=self.F_R(R, Z),
            F_eta=self.F_eta(R, Z),
            U_R=self.U_R(R, Z),
            U_eta=self.U_eta(R, Z),
            moments=self.moments(R, Z),
            moments_Z=self.moments_Z(R, Z),
            P=self.P(R, Z),
            P_Z=self.P_Z(R, Z),
        )

    def inertial(self, R: float, Z: float) -> tuple[float, float]:
        data = self.evaluate(R, Z)
        return float(data["I_theta"]), float(data["I_z"])

    def shear(self, R: float, Z: float) -> tuple[float, float]:
        data = self.evaluate(R, Z)
        return float(data["S_theta"]), float(data["S_z"])

    def stress(self, R: float, Z: float) -> tuple[float, float]:
        data = self.evaluate(R, Z)
        return float(data["T_theta"]), float(data["T_z"])

    def N(self, R: float, Z: float) -> tuple[float, float]:
        data = self.evaluate(R, Z)
        return float(data["N_theta"]), float(data["N_z"])

    def metadata(self) -> dict[str, Any]:
        return {
            "source": SOURCE,
            "source_version": SOURCE_VERSION,
            "source_sections": SOURCE_SECTIONS,
            "delta": self.delta,
            "moment_keys": MOMENT_KEYS,
            "stress_is_profile_derived": True,
            "inertial_formula": "(3.16)-(3.17)",
            "shear_formula": "(3.18)",
            "radial_velocity_formula": "(3.9)",
            "global_closure_claim": False,
            "pde_claim": False,
        }


__all__ = [
    "MOMENT_KEYS",
    "ProfileStress",
    "SOURCE",
    "SOURCE_SECTIONS",
    "SOURCE_VERSION",
    "evaluate_profile_stress",
]
