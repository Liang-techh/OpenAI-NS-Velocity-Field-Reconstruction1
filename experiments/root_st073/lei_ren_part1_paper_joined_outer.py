"""Bounded join of the corrected inner field to its shared source profile.

The inner adapter owns the field through ``Rh``.  Beyond that radius this
prototype calls the *same* ``CorrectedSourceProfile`` held by
``inner.provider.reshape.switches.comparison.bundle['profile']``.  The
terminal axial mass and pressure differences are retained as explicit join
offsets.  They are not rounded to zero and they are not used to claim finite
energy or a completed five-moment construction.

The public coordinate is ``y = log(R / inner.r)``.  ``evaluate(y, Z)`` has
the subset of the profile record needed by ``LocalCoreExitField``.  The
outer record deliberately marks the full five moments unavailable: only the
axial mass/mean and its Z derivative are propagated for the radial transport
prototype.
"""

from __future__ import annotations

from functools import lru_cache
import json
from pathlib import Path
import sys
from typing import Any, Callable

import mpmath as mp


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from lei_ren_part1_paper_exit_field import LocalCoreExitField, physical_chart  # noqa: E402
from lei_ren_part1_paper_inner_corrected_field import (  # noqa: E402
    CorrectedInnerField,
    build_candidate,
)


MOMENT_KEYS = ("theta", "z", "theta_z", "z_theta", "p")


def _mp(value: Any) -> mp.mpf:
    return value if isinstance(value, mp.mpf) else mp.mpf(str(value))


def _signed_log(value: Any, digits: int = 60) -> dict[str, Any]:
    """Keep sign, log magnitude and an arbitrary-exponent MP value."""

    x = _mp(value)
    if x == 0:
        return {"sign": 0, "log_abs": None, "arbitrary_exponent_value": "0"}
    return {
        "sign": 1 if x > 0 else -1,
        "log_abs": mp.nstr(mp.log(abs(x)), digits),
        "arbitrary_exponent_value": mp.nstr(x, digits),
    }


def _from_signed_log(row: Any) -> mp.mpf:
    """Read the full MP value when present, with a legacy log fallback."""

    if not isinstance(row, dict):
        return _mp(row)
    exact = row.get("arbitrary_exponent_value")
    if exact is not None:
        value = _mp(exact)
        sign = int(row.get("sign", 0))
        if sign and value and (1 if value > 0 else -1) != sign:
            raise ValueError("signed-log receipt sign disagrees with arbitrary value")
        return value
    sign = int(row.get("sign", 0))
    log_abs = row.get("log_abs")
    if sign == 0 or log_abs is None:
        return mp.mpf(0)
    return sign * mp.exp(_mp(log_abs))


def _relative(left: Any, right: Any) -> mp.mpf:
    a = _mp(left)
    b = _mp(right)
    scale = max(abs(a), abs(b))
    return abs(a - b) / scale if scale else abs(a - b)


def _n(value: Any, digits: int = 60) -> str:
    return mp.nstr(_mp(value), digits)


def _finite_difference(
    function: Callable[[mp.mpf], mp.mpf], z: mp.mpf, *, step: mp.mpf
) -> mp.mpf:
    """Fourth-order centered Z derivative on the profile's float-Z API."""

    h = _mp(step)
    if abs(z) + 2 * h >= 1:
        raise ValueError("Z stencil leaves the open profile domain")
    return (
        function(z - 2 * h)
        - 8 * function(z - h)
        + 8 * function(z + h)
        - function(z + 2 * h)
    ) / (12 * h)


class JoinedLocalCoreExitField(LocalCoreExitField):
    """LocalCoreExitField with an honest joined-field scope label."""

    def evaluate(self, x: Any, y: Any, z: Any, t: Any) -> dict[str, Any]:
        result = super().evaluate(x, y, z, t)
        result["scope"] = (
            "Joined corrected inner through shared source outer/heat; "
            "full five moments, finite energy and NS closure remain open"
        )
        result["outer_matching_complete"] = False
        return result


class JoinedOuterField:
    """Inner field plus a same-profile outer/heat continuation."""

    def __init__(self, inner: CorrectedInnerField, *, derivative_step: Any = "1e-6"):
        self.inner = inner
        self.provider = inner.provider
        self.core = inner.core
        self.precision = int(inner.precision)
        self.work_precision = self.precision
        self.r = _mp(inner.r)
        self.derivative_step = _mp(derivative_step)
        self.outer = self.provider.reshape.switches.comparison.bundle["profile"]
        self.schedule = self.outer.schedule
        self.delta = _mp(self.core.delta)
        with mp.workdps(self.precision):
            # The terminal construction is defined by provider phase 3 and
            # is independent of Z in its physical radius.  The corrected
            # field's x=e evaluation is used for values and moments below;
            # provider phase 3 alone would omit the five-bump repair.
            terminal = self.provider.evaluate_phase(3, mp.mpf(0))
            self.Rh = _mp(terminal["R"])
            self.logRh = _mp(terminal["logR"])
            self.join_y = self.logRh - mp.log(self.r)
            self.logRref = _mp(str(self.schedule.logRref))
            self.outer_reference_y = self.logRref - mp.log(self.r)
            self.outer_join_offset = self.logRh - self.logRref
        self._terminal_cache: dict[str, dict[str, mp.mpf]] = {}

    def _z_key(self, z: Any) -> str:
        return mp.nstr(_mp(z), self.precision)

    @lru_cache(maxsize=32)
    def terminal_offsets(self, z_key: str) -> dict[str, mp.mpf]:
        """Actual inner mass and pressure offsets at ``Rh``.

        The reference mass is exactly ``4 Z Rh`` and its Z derivative is
        ``4 Rh``.  Pressure is compared with the same outer profile at Rh.
        """

        with mp.workdps(self.precision):
            z = _mp(z_key)
            value = self.inner.evaluate_x(mp.e, z)
            inner_mass = _mp(value["moments"]["z"])
            inner_mass_z = _mp(value["momentsZ"]["z"])
            reference_mass = 4 * z * self.Rh
            reference_mass_z = 4 * self.Rh
            base = self._outer_pressure(self.logRh, z)
            base_z = self._outer_pressure_Z(self.logRh, z)
            inner_pressure = _mp(value["P"])
            inner_pressure_z = _mp(value["PZ"])
            return {
                "mass_offset": inner_mass - reference_mass,
                "mass_offset_Z": inner_mass_z - reference_mass_z,
                "pressure_offset": inner_pressure - base,
                "pressure_offset_Z": inner_pressure_z - base_z,
                "inner_mass": inner_mass,
                "inner_mass_Z": inner_mass_z,
                "inner_pressure": inner_pressure,
                "inner_pressure_Z": inner_pressure_z,
                "reference_mass": reference_mass,
                "reference_mass_Z": reference_mass_z,
                "outer_base_pressure": base,
                "outer_base_pressure_Z": base_z,
            }

    def _outer_values(self, log_radius: mp.mpf, z: mp.mpf) -> dict[str, mp.mpf]:
        """Read source values and analytic schedule jets in MP arithmetic."""

        # CorrectedSourceProfile.values keeps Utheta as an MP value, while
        # the schedule supplies the shared analytic Z and log-R jets.  The
        # sampled heat continuation is outside the angular bump support, so
        # these are the actual jets of the selected outer values there.
        row = self.outer.values(mp.nstr(log_radius, self.precision), z)
        base = self.schedule.at_log_radius(mp.nstr(log_radius, self.precision), z)
        Utheta = _mp(row["Utheta"])
        Uz = _mp(row["Uz"])
        # Reconstruct the Z jet from the Decimal log jet.  The schedule's
        # float ``Utheta_Z`` field can underflow on the exact heat tail.
        Utheta_Z = Utheta * _mp(base["dlogU_dZ"])
        Uz_Z = _mp(base["Uz_Z"])
        R = mp.exp(log_radius)
        root = mp.sqrt(2 * R)
        return {
            "Utheta": Utheta,
            "Utheta_Z": Utheta_Z,
            "Uz": Uz,
            "Uz_Z": Uz_Z,
            "R": R,
            "F": Utheta / root,
            "FZ": Utheta_Z / root,
            "logF_slope": _mp(base["logF_slope"]),
            "logarithmic_slope": _mp(base["logarithmic_slope"]),
            "heat_method": row.get("heat_method", base.get("heat_method")),
        }

    def _outer_average(self, log_radius: mp.mpf, z: mp.mpf) -> tuple[mp.mpf, mp.mpf]:
        """Actual profile axial average and sampled Z derivative."""

        average = _mp(self.outer.axial_average(mp.nstr(log_radius, self.precision), z))
        if log_radius - self.logRref <= 0:
            return average, mp.mpf(4)
        h = self.derivative_step
        derivative = _finite_difference(
            lambda zz: _mp(self.outer.axial_average(mp.nstr(log_radius, self.precision), float(zz))),
            z,
            step=h,
        )
        return average, derivative

    def _outer_pressure(self, log_radius: mp.mpf, z: mp.mpf) -> mp.mpf:
        row = self.outer.pressure_at_log_radius(mp.nstr(log_radius, self.precision), float(z))
        return _from_signed_log(row["corrected_pressure_nominal"])

    def _outer_pressure_Z(self, log_radius: mp.mpf, z: mp.mpf) -> mp.mpf:
        h = self.derivative_step
        return _finite_difference(
            lambda zz: self._outer_pressure(log_radius, zz), z, step=h
        )

    def _outer_state(self, log_radius: mp.mpf, z: mp.mpf, *, carry: bool = True) -> dict[str, Any]:
        with mp.workdps(self.work_precision):
            values = self._outer_values(log_radius, z)
            R = values["R"]
            average, average_Z = self._outer_average(log_radius, z)
            base_mass = R * average
            base_mass_Z = R * average_Z
            offsets = self.terminal_offsets(self._z_key(z))
            mean_seeded = getattr(self, 'axial_mean_already_seeded', False)
            if carry:
                mass = base_mass if mean_seeded else base_mass + offsets["mass_offset"]
                mass_Z = base_mass_Z if mean_seeded else base_mass_Z + offsets["mass_offset_Z"]
                pressure_offset = offsets["pressure_offset"]
                pressure_offset_Z = offsets["pressure_offset_Z"]
            else:
                mass = base_mass
                mass_Z = base_mass_Z
                pressure_offset = mp.mpf(0)
                pressure_offset_Z = mp.mpf(0)
            pressure_base = self._outer_pressure(log_radius, z)
            pressure_base_Z = self._outer_pressure_Z(log_radius, z)
            pressure = pressure_base + pressure_offset
            pressure_Z = pressure_base_Z + pressure_offset_Z
            d = 1 - z * z
            L = 1 - self.delta * z * z
            root = mp.sqrt(2 * R)
            # This is the same radial transport identity used by the inner
            # stress evaluator, with the actual propagated Mz/Mz_Z.
            Ur = (
                2 * z * R * values["Uz"]
                - (1 - self.delta) * z * mass
                - d * mass_Z
            ) / (L * root)
            return {
                "R": R,
                "logR": log_radius,
                "Z": z,
                "F": values["F"],
                "FZ": values["FZ"],
                "Uz": values["Uz"],
                "UZ": values["Uz_Z"],
                "Ur": Ur,
                "P": pressure,
                "PZ": pressure_Z,
                "Utheta": values["Utheta"],
                "Utheta_Z": values["Utheta_Z"],
                "axial_average": average + (offsets["mass_offset"] / R if carry and not mean_seeded else 0),
                "axial_average_Z": average_Z + (offsets["mass_offset_Z"] / R if carry and not mean_seeded else 0),
                "axial_mean_already_seeded": mean_seeded,
                "base_axial_average": average,
                "base_axial_average_Z": average_Z,
                "mass": mass,
                "mass_Z": mass_Z,
                "base_mass": base_mass,
                "base_mass_Z": base_mass_Z,
                "mass_offset": offsets["mass_offset"] if carry else mp.mpf(0),
                "mass_offset_Z": offsets["mass_offset_Z"] if carry else mp.mpf(0),
                "pressure_base": pressure_base,
                "pressure_base_Z": pressure_base_Z,
                "pressure_offset": pressure_offset,
                "pressure_offset_Z": pressure_offset_Z,
                "heat_method": values["heat_method"],
                "moments": None,
                "momentsZ": None,
                "full_five_moments_available": False,
                "outer_velocity_jets_complete": False,
                "velocity_jet_scope": "FZ/UZ use schedule jets; corrected angular bumps and axial pulse jets are not implemented.",
                "outer_moment_scope": (
                    "Only axial mass/mean is propagated from Rh; full five moments "
                    "beyond the early exterior are unavailable."
                ),
                "region": "shared_corrected_outer_profile",
                "outer_matching_complete": False,
                "finite_energy_certified": False,
                "scale_recursion_established": False,
            }

    def evaluate(self, y: Any, Z: Any) -> dict[str, Any]:
        """Global profile dispatch in ``y=log(R/inner.r)`` coordinates."""

        with mp.workdps(self.work_precision):
            yy = _mp(y)
            z = _mp(Z)
            if not abs(z) < 1:
                raise ValueError("Require |Z|<1")
            log_radius = mp.log(self.r) + yy
            if log_radius <= self.logRh:
                return self.inner.evaluate(yy, z)
            return self._outer_state(log_radius, z, carry=True)

    def physical_field(self, *, nu: Any = ".01", T: Any = 0) -> LocalCoreExitField:
        """Return the requested LocalCoreExitField wrapper over this join."""

        return JoinedLocalCoreExitField(self, nu=nu, T=T)

    def _tail_state(self, z: mp.mpf) -> dict[str, Any]:
        """Evaluate one exact-heat point and expose the 1/R transport term."""

        with mp.workdps(self.work_precision):
            log_radius = self.logRref + _mp(str(self.schedule.y_b)) + 1
            state = self._outer_state(log_radius, z, carry=True)
            M = state["mass"]
            MZ = state["mass_Z"]
            d = 1 - z * z
            L = 1 - self.delta * z * z
            coefficient = -((1 - self.delta) * z * M + d * MZ) / L
            state["tail_log_radius"] = log_radius
            state["tail_mass"] = M
            state["tail_mass_Z"] = MZ
            state["tail_radial_transport_coefficient"] = coefficient
            state["tail_radial_transport_coefficient_over_sqrt2"] = coefficient / mp.sqrt(2)
            state["tail_radial_energy_integrand"] = (
                "coefficient^2 / R (logarithmic divergence when coefficient != 0)"
            )
            state["tail_log_energy_divergence_detected"] = bool(coefficient != 0)
            return state


_DEFAULT: JoinedOuterField | None = None


def build_joined_field() -> JoinedOuterField:
    global _DEFAULT
    if _DEFAULT is None:
        _DEFAULT = JoinedOuterField(CorrectedInnerField(build_candidate()))
    return _DEFAULT


def evaluate(y: Any, Z: Any) -> dict[str, Any]:
    """Evaluate the global inner-through-heat profile."""

    return build_joined_field().evaluate(y, Z)


def physical_field(*, nu: Any = ".01", T: Any = 0) -> LocalCoreExitField:
    """Expose the same joined field through ``LocalCoreExitField``."""

    return build_joined_field().physical_field(nu=nu, T=T)


def _receipt_row(field: JoinedOuterField, name: str, log_radius: mp.mpf, z: mp.mpf) -> dict[str, Any]:
    with mp.workdps(field.work_precision):
        y = log_radius - mp.log(field.r)
        state = field.evaluate(y, z)
        values = {
            key: _signed_log(state[key])
            for key in ("F", "FZ", "Uz", "UZ", "Ur", "P", "PZ")
        }
        return {
            "name": name,
            "y": _n(y, 80),
            "schedule_offset_from_logRref": _n(log_radius - field.logRref, 80),
            "region": state.get("region"),
            "heat_method": state.get("heat_method"),
            "values": values,
            "axial_average": _signed_log(state.get("axial_average", 0)),
            "axial_average_Z": _signed_log(state.get("axial_average_Z", 0)),
            "mass_offset": _signed_log(state.get("mass_offset", 0)),
            "mass_offset_Z": _signed_log(state.get("mass_offset_Z", 0)),
            "pressure_offset": _signed_log(state.get("pressure_offset", 0)),
            "pressure_offset_Z": _signed_log(state.get("pressure_offset_Z", 0)),
            "full_five_moments_available": state.get("full_five_moments_available", True),
            "outer_velocity_jets_complete": state.get("outer_velocity_jets_complete"),
            "outer_moment_scope": state.get("outer_moment_scope"),
        }


def run() -> dict[str, Any]:
    """Run a small Z=.3 join receipt through the source heat connection."""

    field = build_joined_field()
    with mp.workdps(field.work_precision):
        z = mp.mpf(".3")
        terminal = field.inner.evaluate_x(mp.e, z)
        log_rh = _mp(terminal["logR"])
        offsets = field.terminal_offsets(field._z_key(z))
        unjoined = field._outer_state(log_rh, z, carry=False)
        joined = field._outer_state(log_rh, z, carry=True)
        # The inner names are F,FZ,Uz,UZ,Ur,P,PZ; keeping this explicit makes
        # the receipt readable and avoids any target-moment substitution.
        join_relative_errors = {
            "F": _n(_relative(unjoined["F"], terminal["F"]), 60),
            "FZ": _n(_relative(unjoined["FZ"], terminal["FZ"]), 60),
            "Uz": _n(_relative(unjoined["Uz"], terminal["Uz"]), 60),
            "UZ": _n(_relative(unjoined["UZ"], terminal["UZ"]), 60),
            "Ur": _n(_relative(unjoined["Ur"], terminal["Ur"]), 60),
            "P": _n(_relative(unjoined["P"], terminal["P"]), 60),
            "PZ": _n(_relative(unjoined["PZ"], terminal["PZ"]), 60),
        }
        carried_join_errors = {
            "mass": _n(_relative(joined["mass"], terminal["moments"]["z"]), 60),
            "mass_Z": _n(_relative(joined["mass_Z"], terminal["momentsZ"]["z"]), 60),
            "pressure": _n(_relative(joined["P"], terminal["P"]), 60),
            "pressure_Z": _n(_relative(joined["PZ"], terminal["PZ"]), 60),
        }
        rows = [
            _receipt_row(field, "Rh_join", log_rh, z),
            _receipt_row(field, "Rref", field.logRref, z),
            _receipt_row(field, "heat_tail_start", field.logRref + _mp(str(field.schedule.y_tail)), z),
            _receipt_row(field, "heat_connection", field.logRref + _mp(str(field.schedule.y_b)), z),
            _receipt_row(field, "exact_heat", field.logRref + _mp(str(field.schedule.y_b)) + 1, z),
        ]
        tail = field._tail_state(z)
        # One global-callable roundtrip through the same LocalCoreExitField at
        # the exact heat point.  The chart itself uses arbitrary MP radii.
        physical = field.physical_field(nu=".01", T=0)
        heat_log_radius = field.logRref + _mp(str(field.schedule.y_b)) + 1
        heat_state = field._outer_state(heat_log_radius, z, carry=True)
        chart = physical_chart(
            heat_state,
            z,
            "-4",
            delta=field.delta,
            precision=field.work_precision,
        )
        back = physical.evaluate(*chart["xyz"], -chart["tau"])
        physical_roundtrip = {
            "region": back["region"],
            "Z_relative_error": _n(_relative(back["Z"], z), 60),
            "R_relative_error": _n(_relative(back["R"], heat_state["R"]), 60),
            "cylindrical_velocity_relative_errors": [
                _n(_relative(a, b), 60)
                for a, b in zip(back["cylindrical_velocity"], chart["cylindrical_velocity"])
            ],
            "wrapper": "LocalCoreExitField",
            "wrapper_subclass": "JoinedLocalCoreExitField",
            "scope": "Joined inner through shared source heat; no finite-energy certificate",
        }
        report = {
            "Z": ".3",
            "precision": field.work_precision,
            "files": {
                "python": str(Path(__file__).resolve()),
                "json": str(Path(__file__).with_suffix(".json").resolve()),
                "markdown": str(Path(__file__).with_suffix(".md").resolve()),
            },
            "construction": {
                "inner": "CorrectedInnerField(build_candidate())",
                "outer_profile": "inner.provider.reshape.switches.comparison.bundle['profile']",
                "coordinate": "y=log(R/inner.r)",
                "Rh_logR": _n(log_rh, 100),
                "logRref": _n(field.logRref, 100),
                "outer_join_offset_from_logRref": _n(field.outer_join_offset, 100),
            },
            "terminal_offsets_relative_to_reference_4ZR": {
                "mass_offset": _signed_log(offsets["mass_offset"]),
                "mass_offset_Z": _signed_log(offsets["mass_offset_Z"]),
                "mass_offset_over_Rh": _signed_log(offsets["mass_offset"] / field.Rh),
                "mass_offset_Z_over_Rh": _signed_log(offsets["mass_offset_Z"] / field.Rh),
                "pressure_offset": _signed_log(offsets["pressure_offset"]),
                "pressure_offset_Z": _signed_log(offsets["pressure_offset_Z"]),
                "small_nonzero_offsets_were_preserved": bool(
                    offsets["mass_offset"] != 0 or offsets["mass_offset_Z"] != 0
                ),
                "pressure_base_at_Rh": _signed_log(offsets["outer_base_pressure"]),
                "inner_pressure_at_Rh": _signed_log(offsets["inner_pressure"]),
            },
            "join_relative_errors_before_carrying_offsets": join_relative_errors,
            "join_relative_errors_after_carrying_offsets": carried_join_errors,
            "rows": rows,
            "tail_radial_transport": {
                "sample": "exact_heat y_b+1",
                "tail_logR": _n(tail["tail_log_radius"], 100),
                "tail_mass": _signed_log(tail["tail_mass"]),
                "tail_mass_Z": _signed_log(tail["tail_mass_Z"]),
                "tail_mass_over_Rh": _signed_log(tail["tail_mass"] / field.Rh),
                "tail_mass_Z_over_Rh": _signed_log(tail["tail_mass_Z"] / field.Rh),
                "coefficient_for_V_over_R": _signed_log(
                    tail["tail_radial_transport_coefficient"]
                ),
                "coefficient_over_Rh": _signed_log(
                    tail["tail_radial_transport_coefficient"] / field.Rh
                ),
                "coefficient_for_Ur_asymptotic": _signed_log(
                    tail["tail_radial_transport_coefficient_over_sqrt2"]
                ),
                "energy_integrand": tail["tail_radial_energy_integrand"],
                "nonzero_coefficient_detected": tail["tail_log_energy_divergence_detected"],
                "axial_average_Z_method": (
                    "Fourth-order centered finite difference with h=1e-6 through "
                    "the profile API, whose Z argument is float-backed."
                ),
                "coefficient_status": (
                    "Sampled arbitrary-MP transport diagnostic; the float-backed "
                    "outer Z derivative and cancellation are not a continuous "
                    "source coefficient proof."
                ),
                "interpretation": (
                    "A nonzero carried axial mass leaves V/R ~ C/R, so the "
                    "corresponding radial transport energy has C^2/R logarithmic "
                    "divergence. This is a defect receipt, not a finite-energy claim."
                ),
            },
            "physical_callable_roundtrip": physical_roundtrip,
            "availability": {
                "full_five_moments_beyond_early_exterior": False,
                "outer_velocity_jets_complete": False,
                "moments_were_fabricated": False,
                "actual_axial_average_used": True,
                "actual_axial_average_Z_used": True,
                "actual_axial_average_Z_method": "sampled fourth-order float-Z stencil h=1e-6",
                "outer_corrected_Z_jet_limitation": (
                    "Outer values use the shared schedule log-Z jet; a complete "
                    "angular-correction Z jet is unavailable in the profile API."
                ),
                "pressure_offset_carried": True,
                "heat_profile_reached": True,
                "finite_energy_certified": False,
                "scale_recursion_established": False,
                "global_NS_or_stress_certificate": False,
            },
        }
    path = Path(__file__).with_suffix(".json")
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "json": str(path),
                "rows": len(rows),
                "mass_offset": report["terminal_offsets_relative_to_reference_4ZR"]["mass_offset"],
                "tail_transport": report["tail_radial_transport"]["coefficient_for_V_over_R"],
                "nonzero_tail_coefficient": report["tail_radial_transport"]["nonzero_coefficient_detected"],
            }
        ),
        flush=True,
    )
    return report


if __name__ == "__main__":
    run()
