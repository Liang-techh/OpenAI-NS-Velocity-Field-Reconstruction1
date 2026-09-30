"""Cumulative source axial moments through the Section 7.5 pulse.

This adapter evaluates the two linear cumulative moments using the actual
incoming rows and the coefficients returned by ``CorrectedSourceProfile``.
For the pulse it uses the dedicated normalized primitives; for the two end
bumps it integrates the incomplete translated bump with the same canonical
quadrature order as the coefficient solve.  The dimensional moments retain
the unflattened angular amplitude powers

``Mz = R E * m1`` and
``Mtheta_z = sqrt(2) R**(3/2) E**2 * m2``.

The terminal values are reported as computed residuals.  They are not
replaced by zero, since the source inputs and all quadratures are numerical.
This module does not claim exact continuous moment closure, pressure closure,
PDE validity, cone admissibility, or global energy control.
"""

from __future__ import annotations

from decimal import Decimal, localcontext
from functools import lru_cache
import json
from pathlib import Path
import sys
from typing import Any

import mpmath as mp
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SRC = ROOT / "src"
for path in (SRC, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from lei_ren_part1_paper_axial_correction import (  # noqa: E402
    CONSISTENT_QUADRATURE_ORDER,
    from_signed_log,
    independent_end_bump_replay,
)
from lei_ren_part1_paper_axial_primitive import (  # noqa: E402
    normalized_pulse_primitive,
)
from lei_ren_part1_paper_corrected_profile import CorrectedSourceProfile  # noqa: E402
from lei_ren_part1_paper_outer_closure import _paper_raw_bump  # noqa: E402


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_VERSION = "2609.35406v1"
SOURCE_EQUATIONS = "(2.21)-(2.22), (7.29), (7.31), (7.36)"
ELL = mp.mpf(".15")
BUMP_CENTERS = (mp.mpf("-3"), mp.mpf("-1"))


def _signed_log(value: Any, *, precision: int) -> dict[str, Any]:
    value_mp = mp.mpf(value)
    if value_mp == 0:
        return {
            "sign": 0,
            "log_abs": None,
            "arbitrary_exponent_value": "0",
        }
    return {
        "sign": 1 if value_mp > 0 else -1,
        "log_abs": mp.nstr(mp.log(abs(value_mp)), precision),
        "arbitrary_exponent_value": mp.nstr(value_mp, precision),
    }


def _signed_from_log(log_abs: Any, sign: Any, *, precision: int) -> dict[str, Any]:
    sign_int = int(sign)
    if sign_int == 0 or log_abs is None:
        return _signed_log(0, precision=precision)
    return _signed_log(sign_int * mp.exp(mp.mpf(str(log_abs))), precision=precision)


def _primitive_value(receipt: dict[str, Any]) -> mp.mpf:
    sign = int(receipt["sign"])
    if sign == 0 or receipt.get("log_J") is None:
        return mp.mpf("0")
    return sign * mp.exp(mp.mpf(str(receipt["log_J"])))


@lru_cache(maxsize=8)
def _bump_nodes(order: int) -> tuple[tuple[mp.mpf, ...], tuple[mp.mpf, ...], mp.mpf]:
    """Return canonical nodes, weights, and the source bump normalization."""

    nodes, weights = np.polynomial.legendre.leggauss(int(order))
    nodes_mp = tuple(mp.mpf(str(float(node))) for node in nodes)
    weights_mp = tuple(mp.mpf(str(float(weight))) for weight in weights)
    raw = tuple(mp.mpf(str(float(_paper_raw_bump(float(node))))) for node in nodes)
    normalization = sum(weight * value for weight, value in zip(weights_mp, raw))
    return nodes_mp, weights_mp, normalization


def _bump_primitive(
    end: Any,
    lam: mp.mpf,
    center: mp.mpf,
    *,
    order: int,
) -> mp.mpf:
    """Return ``int_{-inf}^{end} exp(lam*s) beta(s-center) ds``."""

    end_mp = mp.mpf(str(end))
    left = center - ELL
    right = center + ELL
    if end_mp <= left:
        return mp.mpf("0")
    nodes, weights, normalization = _bump_nodes(int(order))
    if end_mp >= right:
        # Use the exact canonical nodes used to build the coefficient rows.
        total = mp.mpf("0")
        for node, weight in zip(nodes, weights):
            raw = mp.mpf(str(float(_paper_raw_bump(float(node)))))
            s = ELL * node
            beta = raw / (ELL * normalization)
            total += weight * mp.exp(lam * (center + s)) * beta * ELL
        return total

    half = (end_mp - left) / 2
    midpoint = (end_mp + left) / 2
    total = mp.mpf("0")
    for node, weight in zip(nodes, weights):
        s = (midpoint + half * node - center) / ELL
        raw = mp.mpf(str(float(_paper_raw_bump(float(s)))))
        beta = raw / (ELL * normalization)
        total += weight * mp.exp(lam * (midpoint + half * node)) * beta * half
    return total


def _decimal_log_at_xi(schedule: Any, xi: Any) -> Decimal:
    with localcontext() as context:
        context.prec = max(int(schedule.decimal_precision), 180)
        return schedule.logR_p + Decimal(str(xi)) / Decimal(str(schedule.mu))


def _decimal_log_at_end(schedule: Any, end: Any) -> Decimal:
    with localcontext() as context:
        context.prec = max(int(schedule.decimal_precision), 180)
        return schedule.logR_v + Decimal(str(end))


def _mp_log_radius(log_radius: Any) -> mp.mpf:
    return mp.mpf(str(log_radius))


class PaperPulseCumulative:
    """Evaluate actual normalized cumulative axial moments.

    The supported radial domain starts at ``R_p`` and continues through the
    pulse, the two translated end bumps, and the post-``R_v`` constant-moment
    tail.  Values before ``R_p`` require a separate incoming primitive and are
    rejected rather than inferred from the endpoint receipt.
    """

    def __init__(
        self,
        profile: Any | None = None,
        *,
        precision: int = 160,
        order: int = 128,
    ) -> None:
        if precision < 80:
            raise ValueError("precision must be at least 80")
        if int(order) != order or order < 16:
            raise ValueError("order must be an integer at least 16")
        self.profile = profile or CorrectedSourceProfile(precision=precision, order=order)
        self.precision = int(getattr(self.profile, "precision", precision))
        self.requested_order = int(order)

    def _data(self, Z: Any) -> tuple[float, dict[str, Any], dict[str, Any], mp.mpf, mp.mpf]:
        z = float(Z)
        if not np.isfinite(z) or abs(z) > 1:
            raise ValueError("Z must lie in [-1,1]")
        data = self.profile.coefficients(z)
        incoming = data["incoming"]
        axial = data["axial"]
        m1 = from_signed_log(incoming["m1_Mz_over_RpEp"])
        m2 = from_signed_log(incoming["m2_Mtheta_z_over_Rp_sqrt2RpEp2"])
        return z, incoming, axial, m1, m2

    def _row_values(
        self,
        log_radius: Any,
        Z: Any,
        *,
        incoming: dict[str, Any],
        axial: dict[str, Any],
        m1: mp.mpf,
        m2: mp.mpf,
    ) -> dict[str, Any]:
        schedule = self.profile.schedule
        log_radius_decimal = schedule._coerce_log_radius(log_radius)
        if log_radius_decimal < schedule.logR_p:
            raise ValueError("cumulative pulse moments require log_radius >= logR_p")
        end_decimal = log_radius_decimal - schedule.logR_v
        end = mp.mpf(str(end_decimal))
        start_decimal = log_radius_decimal - schedule.logR_p
        start = mp.mpf(str(start_decimal))
        mu = mp.mpf(str(schedule.mu))
        lam1 = mp.mpf(".5") - mu
        lam2 = mp.mpf(".5") - 2 * mu
        order = int(axial["quadrature_order"])
        ap = mp.mpf(str(axial["a_p"]))
        c = [from_signed_log(row) for row in axial["c"]]
        row_pulses = []
        row_bumps = []
        row_norms = []

        if end < 0:
            for row, lam, base in ((1, lam1, m1), (2, lam2, m2)):
                xi = mu * start
                pulse = _primitive_value(
                    normalized_pulse_primitive(
                        mu,
                        row,
                        mp.nstr(xi, self.precision),
                        precision=self.precision,
                        order=order,
                    )
                )
                bump = mp.mpf("0")
                if end >= -mp.mpf("3.15"):
                    bump = sum(
                        coefficient * _bump_primitive(
                            end,
                            lam,
                            center,
                            order=order,
                        )
                        for coefficient, center in zip(c, BUMP_CENTERS)
                    ) * mp.exp(-lam * end)
                row_pulses.append(ap * pulse)
                row_bumps.append(bump)
                row_norms.append(base * mp.exp(-lam * start) + ap * pulse + bump)
            terminal = False
        else:
            # The dimensional moments are constant after R_v.  First evaluate
            # the endpoint with the same normalized pulse and full bump rows,
            # then convert that endpoint to the current normalization using
            # the actual unflattened E power.
            endpoint = self._row_values(
                schedule.logR_v,
                Z,
                incoming=incoming,
                axial=axial,
                m1=m1,
                m2=m2,
            ) if end != 0 else None
            if endpoint is None:
                endpoint = self._row_values_at_endpoint(
                    Z,
                    incoming=incoming,
                    axial=axial,
                    m1=m1,
                    m2=m2,
                )
            endpoint_log_mz = endpoint["moments"]["Mz"]["log_abs"]
            endpoint_log_mtheta = endpoint["moments"]["Mtheta_z"]["log_abs"]
            endpoint_sign_mz = endpoint["moments"]["Mz"]["sign"]
            endpoint_sign_mtheta = endpoint["moments"]["Mtheta_z"]["sign"]
            row = schedule.at_log_radius(log_radius_decimal, float(Z))
            log_E = mp.mpf(str(row["log_angular_amplitude"]))
            log_R = _mp_log_radius(log_radius_decimal)
            current_mz = _signed_from_log(
                endpoint_log_mz,
                endpoint_sign_mz,
                precision=self.precision,
            )
            current_mtheta = _signed_from_log(
                endpoint_log_mtheta,
                endpoint_sign_mtheta,
                precision=self.precision,
            )
            norm1 = _signed_from_log(
                (mp.mpf(endpoint_log_mz) - log_R - log_E)
                if endpoint_log_mz is not None
                else None,
                endpoint_sign_mz,
                precision=self.precision,
            )
            norm2 = _signed_from_log(
                (
                    mp.mpf(endpoint_log_mtheta)
                    - mp.log(mp.sqrt(2))
                    - mp.mpf("1.5") * log_R
                    - 2 * log_E
                )
                if endpoint_log_mtheta is not None
                else None,
                endpoint_sign_mtheta,
                precision=self.precision,
            )
            return {
                "log_radius": str(log_radius_decimal),
                "stage_local": schedule.stage_local_coordinate(log_radius_decimal),
                "start": mp.nstr(start, self.precision),
                "xi": mp.nstr(mu * start, self.precision),
                "end_offset": mp.nstr(end, self.precision),
                "terminal": True,
                "normalization": {
                    "Mz_over_R_E": norm1,
                    "Mtheta_z_over_sqrt2_R32_E2": norm2,
                },
                "moments": {"Mz": current_mz, "Mtheta_z": current_mtheta},
                "decomposition": endpoint["decomposition"],
                "log_E": mp.nstr(log_E, self.precision),
                "log_R": mp.nstr(log_R, self.precision),
                "terminal_endpoint": endpoint,
            }

        row = schedule.at_log_radius(log_radius_decimal, float(Z))
        log_E = mp.mpf(str(row["log_angular_amplitude"]))
        log_R = _mp_log_radius(log_radius_decimal)
        norm1, norm2 = row_norms
        mz = mp.exp(log_R + log_E) * norm1
        mtheta = mp.sqrt(2) * mp.exp(mp.mpf("1.5") * log_R + 2 * log_E) * norm2
        return {
            "log_radius": str(log_radius_decimal),
            "stage_local": schedule.stage_local_coordinate(log_radius_decimal),
            "start": mp.nstr(start, self.precision),
            "xi": mp.nstr(mu * start, self.precision),
            "end_offset": mp.nstr(end, self.precision),
            "terminal": False,
            "normalization": {
                "Mz_over_R_E": _signed_log(norm1, precision=self.precision),
                "Mtheta_z_over_sqrt2_R32_E2": _signed_log(norm2, precision=self.precision),
            },
            "moments": {
                "Mz": _signed_log(mz, precision=self.precision),
                "Mtheta_z": _signed_log(mtheta, precision=self.precision),
            },
            "decomposition": {
                "Mz_over_R_E": {
                    "incoming": _signed_log(m1 * mp.exp(-lam1 * start), precision=self.precision),
                    "pulse": _signed_log(row_pulses[0], precision=self.precision),
                    "end_bumps": _signed_log(row_bumps[0], precision=self.precision),
                },
                "Mtheta_z_over_sqrt2_R32_E2": {
                    "incoming": _signed_log(m2 * mp.exp(-lam2 * start), precision=self.precision),
                    "pulse": _signed_log(row_pulses[1], precision=self.precision),
                    "end_bumps": _signed_log(row_bumps[1], precision=self.precision),
                },
            },
            "log_E": mp.nstr(log_E, self.precision),
            "log_R": mp.nstr(log_R, self.precision),
        }

    def _row_values_at_endpoint(
        self,
        Z: Any,
        *,
        incoming: dict[str, Any],
        axial: dict[str, Any],
        m1: mp.mpf,
        m2: mp.mpf,
    ) -> dict[str, Any]:
        # Keep endpoint construction separate so post-terminal evaluation does
        # not recurse through the ``end >= 0`` branch.
        schedule = self.profile.schedule
        mu = mp.mpf(str(schedule.mu))
        lam1 = mp.mpf(".5") - mu
        lam2 = mp.mpf(".5") - 2 * mu
        order = int(axial["quadrature_order"])
        ap = mp.mpf(str(axial["a_p"]))
        c = [from_signed_log(row) for row in axial["c"]]
        start = mp.mpf("13") / mu
        xi = mp.mpf("13")
        rows = []
        pulse_terms = []
        bumps = []
        for row, lam, base in ((1, lam1, m1), (2, lam2, m2)):
            pulse = _primitive_value(
                normalized_pulse_primitive(
                    mu,
                    row,
                    xi,
                    precision=self.precision,
                    order=order,
                )
            )
            bump = sum(
                coefficient * _bump_primitive(
                    0,
                    lam,
                    center,
                    order=order,
                )
                for coefficient, center in zip(c, BUMP_CENTERS)
            ) * mp.exp(-lam * 0)
            pulse_term = ap * pulse
            pulse_terms.append(pulse_term)
            bumps.append(bump)
            rows.append(base * mp.exp(-lam * start) + pulse_term + bump)
        row = schedule.at_log_radius(schedule.logR_v, float(Z))
        log_E = mp.mpf(str(row["log_angular_amplitude"]))
        log_R = _mp_log_radius(schedule.logR_v)
        mz = mp.exp(log_R + log_E) * rows[0]
        mtheta = mp.sqrt(2) * mp.exp(mp.mpf("1.5") * log_R + 2 * log_E) * rows[1]
        return {
            "log_radius": str(schedule.logR_v),
            "stage_local": schedule.stage_local_coordinate(schedule.logR_v),
            "start": mp.nstr(start, self.precision),
            "xi": mp.nstr(xi, self.precision),
            "end_offset": "0",
            "terminal": True,
            "normalization": {
                "Mz_over_R_E": _signed_log(rows[0], precision=self.precision),
                "Mtheta_z_over_sqrt2_R32_E2": _signed_log(rows[1], precision=self.precision),
            },
            "moments": {
                "Mz": _signed_log(mz, precision=self.precision),
                "Mtheta_z": _signed_log(mtheta, precision=self.precision),
            },
            "decomposition": {
                "Mz_over_R_E": {
                    "incoming": _signed_log(m1 * mp.exp(-lam1 * start), precision=self.precision),
                    "pulse": _signed_log(pulse_terms[0], precision=self.precision),
                    "end_bumps": _signed_log(bumps[0], precision=self.precision),
                },
                "Mtheta_z_over_sqrt2_R32_E2": {
                    "incoming": _signed_log(m2 * mp.exp(-lam2 * start), precision=self.precision),
                    "pulse": _signed_log(pulse_terms[1], precision=self.precision),
                    "end_bumps": _signed_log(bumps[1], precision=self.precision),
                },
            },
            "log_E": mp.nstr(log_E, self.precision),
            "log_R": mp.nstr(log_R, self.precision),
        }

    def evaluate(self, log_radius: Any, Z: Any) -> dict[str, Any]:
        with mp.workdps(self.precision):
            z, incoming, axial, m1, m2 = self._data(Z)
            value = self._row_values(
                log_radius,
                z,
                incoming=incoming,
                axial=axial,
                m1=m1,
                m2=m2,
            )
            value["Z"] = z
            value["quadrature_order"] = int(axial["quadrature_order"])
            value["coefficient_requested_order"] = int(axial.get("requested_quadrature_order", self.requested_order))
            return value

    def moments(self, log_radius: Any, Z: Any) -> dict[str, Any]:
        """Return raw and normalized moments, with the full receipt attached."""

        receipt = self.evaluate(log_radius, Z)
        return {
            "Mz": mp.mpf(receipt["moments"]["Mz"]["arbitrary_exponent_value"]),
            "Mtheta_z": mp.mpf(
                receipt["moments"]["Mtheta_z"]["arbitrary_exponent_value"]
            ),
            "Mz_over_R_E": mp.mpf(
                receipt["normalization"]["Mz_over_R_E"]["arbitrary_exponent_value"]
            ),
            "Mtheta_z_over_sqrt2_R32_E2": mp.mpf(
                receipt["normalization"]["Mtheta_z_over_sqrt2_R32_E2"][
                    "arbitrary_exponent_value"
                ]
            ),
            "receipt": receipt,
        }

    def terminal_replay(self, Z: Any) -> dict[str, Any]:
        z, incoming, axial, _, _ = self._data(Z)
        base_rows = [
            incoming["row_normalization"][key]
            for key in ("scaled_base_m1", "scaled_base_m2")
        ]
        replay = independent_end_bump_replay(
            self.profile.schedule.mu,
            axial,
            base_rows,
            self.profile.pulse,
            order=int(axial["quadrature_order"]),
            precision=self.precision,
        )
        endpoint = self.evaluate(self.profile.schedule.logR_v, z)
        replay["terminal_normalization"] = endpoint["normalization"]
        replay["terminal_moments"] = endpoint["moments"]
        replay["Z"] = z
        replay["scope"] = (
            "Actual finite-order row replay and cumulative endpoint; terminal "
            "residual is retained and not replaced by zero."
        )
        return replay


def _relative(defect: mp.mpf, rhs: mp.mpf) -> dict[str, Any]:
    if rhs == 0:
        return {"absolute": _signed_log(defect, precision=80), "relative": None}
    return {
        "absolute": _signed_log(defect, precision=80),
        "relative": mp.nstr(abs(defect / rhs), 50),
    }


def _derivative_check(
    adapter: PaperPulseCumulative,
    log_radius: Any,
    Z: Any,
    label: str,
    *,
    h: str = "1e-6",
) -> dict[str, Any]:
    with mp.workdps(adapter.precision):
        h_mp = mp.mpf(h)
        log0 = adapter.profile.schedule._coerce_log_radius(log_radius)
        # Decimal's ambient default precision is too small to preserve a
        # 1e-6 local increment beside the 1e28-scale stage coordinate.
        plus_log = adapter.profile.log_at(log0, h)
        minus_log = adapter.profile.log_at(log0, "-" + h)
        plus = adapter.evaluate(plus_log, Z)
        minus = adapter.evaluate(minus_log, Z)
        center = adapter.evaluate(log0, Z)
        mz_plus = mp.mpf(plus["moments"]["Mz"]["arbitrary_exponent_value"])
        mz_minus = mp.mpf(minus["moments"]["Mz"]["arbitrary_exponent_value"])
        mt_plus = mp.mpf(plus["moments"]["Mtheta_z"]["arbitrary_exponent_value"])
        mt_minus = mp.mpf(minus["moments"]["Mtheta_z"]["arbitrary_exponent_value"])
        d_mz = (mz_plus - mz_minus) / (2 * h_mp)
        d_mt = (mt_plus - mt_minus) / (2 * h_mp)
        values = adapter.profile.values(log0, float(Z))
        uz = mp.mpf(str(values["Uz"]))
        utheta = mp.mpf(str(values["Utheta"]))
        log_R = _mp_log_radius(log0)
        radius = mp.exp(log_R)
        rhs_mz = radius * uz
        rhs_mt = mp.sqrt(2) * radius ** mp.mpf("1.5") * utheta * uz
        return {
            "label": label,
            "log_radius": str(log0),
            "Z": float(Z),
            "Mz_dlogR": _signed_log(d_mz, precision=adapter.precision),
            "Mz_rhs_R_Uz": _signed_log(rhs_mz, precision=adapter.precision),
            "Mz_error": _relative(d_mz - rhs_mz, rhs_mz),
            "Mtheta_z_dlogR": _signed_log(d_mt, precision=adapter.precision),
            "Mtheta_z_rhs": _signed_log(rhs_mt, precision=adapter.precision),
            "Mtheta_z_error": _relative(d_mt - rhs_mt, rhs_mt),
            "finite_difference_step": h,
            "center": center["stage_local"],
        }


def _sample_log_radii(profile: Any) -> dict[str, Decimal]:
    schedule = profile.schedule
    return {
        "bulk_xi_5": _decimal_log_at_xi(schedule, "5"),
        "end_bump_1_offset_-3": _decimal_log_at_end(schedule, "-3"),
        "end_bump_2_offset_-1": _decimal_log_at_end(schedule, "-1"),
        "terminal_Rv": schedule.logR_v,
        "post_terminal_offset_1": _decimal_log_at_end(schedule, "1"),
    }


def run() -> dict[str, Any]:
    profile = CorrectedSourceProfile(precision=160, order=128)
    adapter = PaperPulseCumulative(profile)
    samples = _sample_log_radii(profile)
    Z = 0.5
    evaluations = {
        name: adapter.evaluate(log_radius, Z)
        for name, log_radius in samples.items()
    }
    derivative_checks = [
        _derivative_check(adapter, samples["bulk_xi_5"], Z, "bulk_xi_5"),
        _derivative_check(adapter, samples["end_bump_1_offset_-3"], Z, "end_bump_1_offset_-3"),
        _derivative_check(adapter, samples["end_bump_2_offset_-1"], Z, "end_bump_2_offset_-1"),
    ]
    replay = adapter.terminal_replay(Z)
    replay_relative = [mp.mpf(value) for value in replay["row_relative_differences"]]
    derivative_relative = [
        mp.mpf(check[key]["relative"])
        for check in derivative_checks
        for key in ("Mz_error", "Mtheta_z_error")
    ]
    checks_passed = bool(
        all(value <= mp.mpf("1e-30") for value in replay_relative)
        and all(value <= mp.mpf("1e-8") for value in derivative_relative)
    )
    report = {
        "source": SOURCE,
        "source_version": SOURCE_VERSION,
        "source_equations": SOURCE_EQUATIONS,
        "api": {
            "class": "PaperPulseCumulative",
            "evaluate": "evaluate(log_radius, Z)",
            "terminal_replay": "terminal_replay(Z)",
        },
        "profile": {
            "precision": profile.precision,
            "requested_order": profile.order,
            "canonical_order": CONSISTENT_QUADRATURE_ORDER,
            "schedule": profile.schedule.metadata(),
        },
        "samples": {name: str(value) for name, value in samples.items()},
        "Z": Z,
        "evaluations": evaluations,
        "terminal_replay": replay,
        "radial_derivative_checks": derivative_checks,
        "checks": {
            "terminal_row_relative_threshold": "1e-30",
            "radial_identity_relative_threshold": "1e-8",
            "terminal_rows_passed": all(value <= mp.mpf("1e-30") for value in replay_relative),
            "radial_identities_passed": all(value <= mp.mpf("1e-8") for value in derivative_relative),
        },
        "checks_passed": checks_passed,
        "normalization": {
            "Mz": "Mz/(R*E)",
            "Mtheta_z": "Mtheta_z/(sqrt(2)*R^(3/2)*E^2)",
            "pulse_E_power": "unflattened E and E^2 from schedule at query radius",
        },
        "scope": (
            "Actual corrected-profile incoming rows, pulse primitives, and end "
            "bumps. Terminal mass residuals are retained. This is a cumulative "
            "moment diagnostic, not exact continuous closure or a PDE/cone/"
            "pressure/global certificate."
        ),
    }
    path = Path(__file__).with_suffix(".json")
    path.write_text(json.dumps(report, indent=2, default=str) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "samples": list(evaluations),
                "checks_passed": checks_passed,
                "terminal_Mz_norm": evaluations["terminal_Rv"]["normalization"]["Mz_over_R_E"],
                "terminal_Mtheta_norm": evaluations["terminal_Rv"]["normalization"]["Mtheta_z_over_sqrt2_R32_E2"],
                "replay_rows": replay["row_relative_differences"],
                "derivative_relative": [
                    (check["label"], check["Mz_error"]["relative"], check["Mtheta_z_error"]["relative"])
                    for check in derivative_checks
                ],
            },
            default=str,
        )
    )
    return report


if __name__ == "__main__":
    run()
