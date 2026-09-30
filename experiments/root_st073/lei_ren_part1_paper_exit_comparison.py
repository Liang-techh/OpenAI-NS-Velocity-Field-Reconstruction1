"""Section 9.23 flat comparison profile on the regular-core exit collar.

The comparison is the paper's auxiliary profile, not the prescribed shear
bridge.  It starts from the actual ``CorePolynomial`` at ``R_a`` and, with

    alpha(y) = 1 - sigma((y-h_b)/h_b),  h_b = .005,

integrates

    d_y log(bar F) = alpha R F_c,R/F_c,
    d_y bar U^z   = alpha R U_c,R.

The same analytic coefficient rows supply the Z and mixed RZ derivatives.
Moments and pressure are accumulated from these comparison fields, then the
source Section 3 stress evaluator supplies ``D=I_theta/bar F`` and
``E=I_z/bar F``.  No finite Z differences are used.

The default builder calls ``core_adapter.build_source_core`` with the actual
dominant AxisPressureJets setup, Lambda = 1e36, and radial degree 18.  The
module is intentionally local to 0 <= y <= 2 h_b (with frozen-value support
after 2 h_b); it makes no cone, matching, or global pressure claim.
"""

from __future__ import annotations

from collections.abc import Mapping
import json
from pathlib import Path
import sys
from typing import Any

import mpmath as mp


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from lei_ren_part1_paper_core_adapter import build_source_core  # noqa: E402
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress  # noqa: E402


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_VERSION = "2609.35406v1"
MOMENT_KEYS = ("theta", "z", "theta_z", "z_theta", "p")


def _mp(value: Any) -> mp.mpf:
    if isinstance(value, mp.mpf):
        return mp.mpf(value)
    return mp.mpf(str(value))


def _nstr(value: Any, digits: int = 40) -> str:
    if value is None:
        return "null"
    return mp.nstr(value if isinstance(value, mp.mpf) else _mp(value), digits)


def _signed_log(value: Any, digits: int = 50) -> dict[str, Any]:
    with mp.workdps(max(80, digits + 20)):
        x = _mp(value)
        if x == 0:
            return {"sign": 0, "log_abs": None, "value": "0"}
        return {
            "sign": 1 if x > 0 else -1,
            "log_abs": _nstr(mp.log(abs(x)), digits),
            "value": _nstr(x, digits),
        }


def _relative(value: mp.mpf, reference: mp.mpf) -> mp.mpf:
    if reference == 0:
        return abs(value)
    return abs(value - reference) / abs(reference)


def _sigma(x: mp.mpf) -> mp.mpf:
    """The source smooth step edge(1,x)/(edge(1,x)+edge(1,1-x))."""

    if x <= 0:
        return mp.mpf(0)
    if x >= 1:
        return mp.mpf(1)
    left = mp.exp(-1 / (x * x))
    right = mp.exp(-1 / ((1 - x) * (1 - x)))
    return left / (left + right)


def _alpha(y: mp.mpf, h_b: mp.mpf) -> mp.mpf:
    return 1 - _sigma((y - h_b) / h_b)


def _state_copy(state: dict[str, Any]) -> dict[str, Any]:
    result = {key: value for key, value in state.items() if key not in ("moments", "moments_Z")}
    result["moments"] = dict(state["moments"])
    result["moments_Z"] = dict(state["moments_Z"])
    return result


def _state_add(*states: dict[str, Any]) -> dict[str, Any]:
    out = {
        key: sum(state[key] for state in states)
        for key in ("F", "F_Z", "Uz", "Uz_Z", "P", "P_Z")
    }
    out["moments"] = {
        key: sum(state["moments"][key] for state in states)
        for key in MOMENT_KEYS
    }
    out["moments_Z"] = {
        key: sum(state["moments_Z"][key] for state in states)
        for key in MOMENT_KEYS
    }
    return out


def _state_scale(state: dict[str, Any], scalar: mp.mpf) -> dict[str, Any]:
    out = {
        key: scalar * state[key]
        for key in ("F", "F_Z", "Uz", "Uz_Z", "P", "P_Z")
    }
    out["moments"] = {key: scalar * state["moments"][key] for key in MOMENT_KEYS}
    out["moments_Z"] = {key: scalar * state["moments_Z"][key] for key in MOMENT_KEYS}
    return out


class Section923Comparison:
    """Cached Section 9.23 comparison around one finite source core."""

    def __init__(
        self,
        bundle: Mapping[str, Any],
        *,
        h_b: Any = ".005",
        transition_steps: int = 32,
    ) -> None:
        self.bundle = bundle
        self.core = bundle["core"]
        self.precision = int(bundle.get("precision", getattr(self.core, "precision", 160)))
        self.transition_steps = int(transition_steps)
        with mp.workdps(self.precision):
            self.axis = bundle["axis"]
            self.Lambda = _mp(bundle["Lambda"])
            self.delta = _mp(getattr(self.core, "delta", "1e-32"))
            self.h_b = _mp(h_b)
            # Public short alias used by the ExitBridge adapter.
            self.hb = self.h_b
            if self.h_b <= 0 or self.transition_steps < 4:
                raise ValueError("h_b must be positive and transition_steps >= 4")
            self.R_a = mp.mpf(4) / self.Lambda
            self.R_limit = mp.mpf("4.1") / self.Lambda
        self._state_cache: dict[tuple[str, str], dict[str, Any]] = {}
        self._core_cache: dict[tuple[str, str], dict[str, Any]] = {}

    # ---------------------------------------------------------------
    # Core analytic jets and state ODE.
    # ---------------------------------------------------------------
    def _z_key(self, Z: Any) -> str:
        return mp.nstr(_mp(Z), self.precision)

    def _core_snapshot(self, R: mp.mpf, Z: mp.mpf) -> dict[str, Any]:
        key = (_nstr(R, self.precision), self._z_key(Z))
        if key in self._core_cache:
            return self._core_cache[key]
        with mp.workdps(self.precision):
            z = _mp(Z)
            coefficients = self.core.coefficients(self._z_key(z))
            value = self.core.evaluate(R, self._z_key(z))
            f_rz = mp.mpf(0)
            uz_rz = mp.mpf(0)
            for n, row in enumerate(coefficients["F"]):
                if n:
                    f_rz += n * row[1] * R ** (n - 1)
            for n, row in enumerate(coefficients["Uz"]):
                if n:
                    uz_rz += n * row[1] * R ** (n - 1)
            value["F_RZ"] = f_rz
            value["Uz_RZ"] = uz_rz
            self._core_cache[key] = value
            return value

    def _state_from_core(self, R: mp.mpf, Z: mp.mpf) -> dict[str, Any]:
        core = self._core_snapshot(R, Z)
        return {
            "F": core["F"],
            "F_Z": core["F_Z"],
            "Uz": core["Uz"],
            "Uz_Z": core["Uz_Z"],
            "P": core["P"],
            "P_Z": core["P_Z"],
            "moments": dict(core["moments"]),
            "moments_Z": dict(core["moments_Z"]),
        }

    def _state_derivative(self, y: mp.mpf, state: dict[str, Any], Z: mp.mpf) -> dict[str, Any]:
        with mp.workdps(self.precision):
            R = self.R_a * mp.exp(y)
            core = self._core_snapshot(R, Z)
            alpha = _alpha(y, self.h_b)
            F_c = core["F"]
            F_c_R = core["F_R"]
            F_c_Z = core["F_Z"]
            U_c_R = core["Uz_R"]
            U_c_Z = core["Uz_Z"]
            log_slope = R * F_c_R / F_c
            log_slope_Z = R * (core["F_RZ"] * F_c - F_c_R * F_c_Z) / (F_c * F_c)
            out = {
                "F": alpha * log_slope * state["F"],
                "F_Z": alpha * (log_slope_Z * state["F"] + log_slope * state["F_Z"]),
                "Uz": alpha * R * U_c_R,
                "Uz_Z": alpha * R * core["Uz_RZ"],
                "P": R * state["F"] ** 2,
                "P_Z": 2 * R * state["F"] * state["F_Z"],
            }
            F = state["F"]
            F_Z = state["F_Z"]
            U = state["Uz"]
            U_Z = state["Uz_Z"]
            out["moments"] = {
                "theta": 2 * R * R * F,
                "z": R * U,
                "theta_z": 2 * R * R * F * U,
                "z_theta": R * (U * U - R * F * F),
                "p": R * F * F,
            }
            out["moments_Z"] = {
                "theta": 2 * R * R * F_Z,
                "z": R * U_Z,
                "theta_z": 2 * R * R * (F_Z * U + F * U_Z),
                "z_theta": R * (2 * U * U_Z - 2 * R * F * F_Z),
                "p": 2 * R * F * F_Z,
            }
            return out

    def _transition_state(self, y: mp.mpf, Z: mp.mpf) -> dict[str, Any]:
        """Integrate the comparison state from h_b to y in MP arithmetic."""

        key = ("transition", _nstr(y, self.precision), self._z_key(Z))
        if key in self._state_cache:
            return _state_copy(self._state_cache[key])
        if y <= self.h_b:
            state = self._state_from_core(self.R_a * mp.exp(y), Z)
            self._state_cache[key] = _state_copy(state)
            return state
        start_y = self.h_b
        state = self._state_from_core(self.R_a * mp.exp(start_y), Z)
        length = y - start_y
        steps = max(1, int(mp.ceil(length / self.h_b * self.transition_steps)))
        step = length / steps
        t = start_y
        for _ in range(steps):
            k1 = self._state_derivative(t, state, Z)
            k2 = self._state_derivative(t + step / 2, _state_add(state, _state_scale(k1, step / 2)), Z)
            k3 = self._state_derivative(t + step / 2, _state_add(state, _state_scale(k2, step / 2)), Z)
            k4 = self._state_derivative(t + step, _state_add(state, _state_scale(k3, step)), Z)
            increment = _state_scale(
                _state_add(k1, _state_scale(k2, 2), _state_scale(k3, 2), k4),
                step / 6,
            )
            state = _state_add(state, increment)
            t += step
        self._state_cache[key] = _state_copy(state)
        return state

    def _state(self, y: mp.mpf, Z: mp.mpf) -> tuple[dict[str, Any], mp.mpf]:
        if y < 0:
            raise ValueError("y must be nonnegative")
        if self.R_a * mp.exp(y) > self.R_limit:
            raise ValueError("comparison query lies outside the finite core radius")
        if y <= 2 * self.h_b:
            return self._transition_state(y, Z), _alpha(y, self.h_b)
        endpoint_y = 2 * self.h_b
        endpoint = self._transition_state(endpoint_y, Z)
        endpoint_R = self.R_a * mp.exp(endpoint_y)
        R = self.R_a * mp.exp(y)
        delta_R = R - endpoint_R
        frozen = _state_copy(endpoint)
        frozen["P"] += delta_R * endpoint["F"] ** 2
        frozen["P_Z"] += delta_R * 2 * endpoint["F"] * endpoint["F_Z"]
        delta_R2 = R * R - endpoint_R * endpoint_R
        frozen["moments"]["theta"] += endpoint["F"] * delta_R2
        frozen["moments_Z"]["theta"] += endpoint["F_Z"] * delta_R2
        frozen["moments"]["z"] += endpoint["Uz"] * delta_R
        frozen["moments_Z"]["z"] += endpoint["Uz_Z"] * delta_R
        frozen["moments"]["theta_z"] += endpoint["F"] * endpoint["Uz"] * delta_R2
        frozen["moments_Z"]["theta_z"] += (
            endpoint["F_Z"] * endpoint["Uz"] + endpoint["F"] * endpoint["Uz_Z"]
        ) * delta_R2
        frozen["moments"]["z_theta"] += (
            endpoint["Uz"] ** 2 * delta_R - endpoint["F"] ** 2 * delta_R2 / 2
        )
        frozen["moments_Z"]["z_theta"] += (
            2 * endpoint["Uz"] * endpoint["Uz_Z"] * delta_R
            - endpoint["F"] * endpoint["F_Z"] * delta_R2
        )
        frozen["moments"]["p"] += endpoint["F"] ** 2 * delta_R
        frozen["moments_Z"]["p"] += 2 * endpoint["F"] * endpoint["F_Z"] * delta_R
        return frozen, mp.mpf(0)

    # ---------------------------------------------------------------
    # Public comparison evaluation.
    # ---------------------------------------------------------------
    def evaluate(self, y: Any, Z: Any) -> dict[str, Any]:
        with mp.workdps(self.precision):
            yy = _mp(y)
            zz = _mp(Z)
            if abs(zz) >= 1:
                raise ValueError("require |Z| < 1")
            state, alpha = self._state(yy, zz)
            R = self.R_a * mp.exp(yy)
            core = self._core_snapshot(R, zz)
            if yy <= 2 * self.h_b:
                F_R = alpha * core["F_R"] * state["F"] / core["F"]
                Uz_R = alpha * core["Uz_R"]
            else:
                F_R = mp.mpf(0)
                Uz_R = mp.mpf(0)
            root = mp.sqrt(2 * R)
            stress = evaluate_mp_stress(
                mp.log(R),
                zz,
                self.delta,
                Utheta=root * state["F"],
                Uz=state["Uz"],
                Utheta_y=root * (R * F_R + state["F"] / 2),
                Utheta_Z=root * state["F_Z"],
                Uz_y=R * Uz_R,
                Uz_Z=state["Uz_Z"],
                moments=state["moments"],
                moments_Z=state["moments_Z"],
                P=state["P"],
                P_Z=state["P_Z"],
                precision=self.precision,
            )
            return {
                "y": yy,
                "Z": zz,
                "R": R,
                "alpha": alpha,
                "F": state["F"],
                "Uz": state["Uz"],
                "FR": F_R,
                "UR": Uz_R,
                "FZ": state["F_Z"],
                "UZ": state["Uz_Z"],
                "moments": state["moments"],
                "momentsZ": state["moments_Z"],
                "P": state["P"],
                "PZ": state["P_Z"],
                "D": stress["I_theta"] / state["F"],
                "E": stress["I_z"] / state["F"],
                "I_theta": stress["I_theta"],
                "I_z": stress["I_z"],
                "S_theta": stress["S_theta"],
                "S_z": stress["S_z"],
                "T_theta": stress["T_theta"],
                "T_z": stress["T_z"],
                "stress": stress,
                "core": {
                    "F": core["F"],
                    "Uz": core["Uz"],
                    "FR": core["F_R"],
                    "UR": core["Uz_R"],
                    "FZ": core["F_Z"],
                    "UZ": core["Uz_Z"],
                },
            }

    def metadata(self) -> dict[str, Any]:
        return {
            "source": SOURCE,
            "source_version": SOURCE_VERSION,
            "equation": "Section 9.23 / source (9.23)",
            "Lambda": _nstr(self.Lambda, self.precision),
            "radial_degree": self.bundle.get("radial_degree"),
            "R_a": _nstr(self.R_a, self.precision),
            "h_b": _nstr(self.h_b, self.precision),
            "comparison_domain": "0 <= y <= 2 h_b; frozen support after 2 h_b while within the finite core radius",
            "transition_steps": self.transition_steps,
            "Z_derivatives": "analytic coefficient rows including mixed RZ rows; no finite Z differences",
            "pressure_source": "core_adapter.build_source_core dominant AxisPressureJets factory",
            "global_cone_certified": False,
            "outer_connection_complete": False,
        }


def build_comparison(*, precision: int = 160, degree: int = 18, Lambda: Any = "1e36", h_b: Any = ".005") -> Section923Comparison:
    bundle = build_source_core(precision=precision, degree=degree, Lambda=Lambda)
    return Section923Comparison(bundle, h_b=h_b)


def _json_value(value: Any) -> Any:
    if isinstance(value, mp.mpf):
        return _nstr(value, 50)
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    return value


def _point_receipt(comparison: Section923Comparison, y: Any, Z: Any) -> dict[str, Any]:
    value = comparison.evaluate(y, Z)
    return {
        "y": _nstr(value["y"]),
        "Z": _nstr(value["Z"]),
        "alpha": _nstr(value["alpha"]),
        "F": _signed_log(value["F"]),
        "Uz": _signed_log(value["Uz"]),
        "FR": _signed_log(value["FR"]),
        "UR": _signed_log(value["UR"]),
        "FZ": _signed_log(value["FZ"]),
        "UZ": _signed_log(value["UZ"]),
        "P": _signed_log(value["P"]),
        "PZ": _signed_log(value["PZ"]),
        "D": _signed_log(value["D"]),
        "E": _signed_log(value["E"]),
        "stress": {key: _signed_log(item) for key, item in value["stress"].items()},
    }


def run_focused_checks() -> dict[str, Any]:
    comparison = build_comparison(precision=160, degree=18, Lambda="1e36", h_b=".005")
    h_b = comparison.h_b
    points = []
    for Z in (".3", _nstr(comparison.axis.Z0, comparison.precision)):
        for y in ("0", ".0025", ".005", ".0075", ".01", ".012"):
            points.append(_point_receipt(comparison, y, Z))
    # Before h_b the comparison is exactly the core.  Check the shared fields
    # directly at two interior points; this is a replay check, not a theorem.
    equality_checks = []
    for Z in (".3", _nstr(comparison.axis.Z0, comparison.precision)):
        for y in ("0", ".0025", ".005"):
            value = comparison.evaluate(y, Z)
            max_rel = max(
                _relative(value[name], value["core"][name])
                for name in ("F", "Uz", "FR", "UR", "FZ", "UZ")
            )
            equality_checks.append({"Z": Z, "y": y, "max_core_equality_relative": _nstr(max_rel, 30)})
    endpoint = comparison.evaluate(".01", ".3")
    frozen = comparison.evaluate(".012", ".3")
    R_endpoint = endpoint["R"]
    R_frozen = frozen["R"]
    delta_R = R_frozen - R_endpoint
    delta_R2 = R_frozen * R_frozen - R_endpoint * R_endpoint
    expected_frozen_moments = {
        "theta": endpoint["moments"]["theta"] + endpoint["F"] * delta_R2,
        "z": endpoint["moments"]["z"] + endpoint["Uz"] * delta_R,
        "theta_z": endpoint["moments"]["theta_z"] + endpoint["F"] * endpoint["Uz"] * delta_R2,
        "z_theta": endpoint["moments"]["z_theta"] + endpoint["Uz"] ** 2 * delta_R - endpoint["F"] ** 2 * delta_R2 / 2,
        "p": endpoint["moments"]["p"] + endpoint["F"] ** 2 * delta_R,
    }
    expected_frozen_moments_Z = {
        "theta": endpoint["momentsZ"]["theta"] + endpoint["FZ"] * delta_R2,
        "z": endpoint["momentsZ"]["z"] + endpoint["UZ"] * delta_R,
        "theta_z": endpoint["momentsZ"]["theta_z"] + (endpoint["FZ"] * endpoint["Uz"] + endpoint["F"] * endpoint["UZ"]) * delta_R2,
        "z_theta": endpoint["momentsZ"]["z_theta"] + 2 * endpoint["Uz"] * endpoint["UZ"] * delta_R - endpoint["F"] * endpoint["FZ"] * delta_R2,
        "p": endpoint["momentsZ"]["p"] + 2 * endpoint["F"] * endpoint["FZ"] * delta_R,
    }
    max_frozen_moment_error = max(
        _relative(frozen["moments"][key], expected_frozen_moments[key])
        for key in MOMENT_KEYS
    )
    max_frozen_moment_Z_error = max(
        _relative(frozen["momentsZ"][key], expected_frozen_moments_Z[key])
        for key in MOMENT_KEYS
    )
    freeze_checks = {
        "F_relative": _nstr(_relative(frozen["F"], endpoint["F"]), 30),
        "Uz_relative": _nstr(_relative(frozen["Uz"], endpoint["Uz"]), 30),
        "FR_zero": frozen["FR"] == 0,
        "UR_zero": frozen["UR"] == 0,
        "max_frozen_moment_relative_error": _nstr(max_frozen_moment_error, 30),
        "max_frozen_moment_Z_relative_error": _nstr(max_frozen_moment_Z_error, 30),
    }
    return {
        "metadata": comparison.metadata(),
        "points": points,
        "core_equality_checks": equality_checks,
        "freeze_checks": freeze_checks,
        "scope": "Finite Section 9.23 comparison collar with actual dominant axis-pressure core factory; no prescribed shear bridge or cone claim.",
        "analytic_Z_derivatives": True,
        "finite_Z_difference_used": False,
        "global_cone_certified": False,
        "outer_connection_complete": False,
    }


if __name__ == "__main__":
    receipt = run_focused_checks()
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(_json_value(receipt), indent=2) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "point_count": len(receipt["points"]),
                "finite_Z_difference_used": receipt["finite_Z_difference_used"],
                "global_cone_certified": receipt["global_cone_certified"],
            }
        ),
        flush=True,
    )


__all__ = ["Section923Comparison", "build_comparison", "run_focused_checks"]
