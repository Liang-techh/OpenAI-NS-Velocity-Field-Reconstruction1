"""Actual angular tail constant and matching diagnosis.

The returned ``Ctheta`` is the constant left by the actual angular primitive
after subtracting the growing reference primitive.  Reference, heat collar,
exterior polynomial and actual-anchor pieces stay separate.  This module
only extracts a target; it does not overwrite the angular field or claim
terminal closure. The exterior atoms below describe the implemented quadratic
heat approximation; they are not a certified exact Lei-Ren heat-tail integral.
"""

from __future__ import annotations

import json
from typing import Any

import mpmath as mp


def _n(value: Any) -> mp.mpf:
    return value if isinstance(value, mp.mpf) else mp.mpf(str(value))


def _tail_terms(
    *,
    delta: Any,
    log_r_tail: Any,
    log_c_inf: Any,
    xi: Any,
    xi_z: Any,
    collar_ref_theta: Any,
    collar_heat_correction_theta: Any,
    collar_heat_correction_theta_z: Any,
):
    """Assemble the normalized infinite heat defect from shared heat atoms."""

    delta = _n(delta)
    h = delta / 2
    if not (0 < h < 1):
        raise ValueError("the tail formula requires 0 < delta/2 < 1")
    xi = _n(xi)
    xi_z = _n(xi_z)
    log_r_tail = _n(log_r_tail)
    e_tail = _n(log_c_inf) - (1 + delta) * log_r_tail / 2
    scale = mp.sqrt(2) * mp.exp(mp.mpf("1.5") * log_r_tail + e_tail)

    c1 = h * (1 + h)
    c2 = h * (h + 1) ** 2 * (h + 2)
    b = {1: -c1 * xi, 2: c2 * xi * xi / 2}
    b_z = {1: -c1 * xi_z, 2: c2 * xi * xi_z}

    # The reference term is the finite collar reference atom minus the
    # reference primitive accumulated from t=0 to the collar endpoint t=3.
    reference_endpoint = -mp.exp(3 * (1 - h)) / (1 - h)
    normalized_reference_terms = {
        "collar_ref_theta": _n(collar_ref_theta),
        "reference_endpoint_subtraction": reference_endpoint,
        "reference_endpoint_subtraction_Z": mp.mpf(0),
    }
    normalized_reference = mp.fsum(
        (normalized_reference_terms["collar_ref_theta"], reference_endpoint)
    )

    # For t >= 3, the two heat Taylor atoms have exponents 1-h-j.  Their
    # integrals converge because j >= 1, while retaining the 1/h factor in
    # the j=1 term when h is very small.
    exterior = {
        j: b[j] * mp.exp(3 * (1 - h - j)) / (j - 1 + h) for j in (1, 2)
    }
    exterior_z = {
        j: b_z[j] * mp.exp(3 * (1 - h - j)) / (j - 1 + h) for j in (1, 2)
    }
    normalized_heat_terms = {
        "collar_heat_correction_theta": _n(collar_heat_correction_theta),
        "exterior_b1": exterior[1],
        "exterior_b2": exterior[2],
        "collar_heat_correction_theta_Z": _n(collar_heat_correction_theta_z),
        "exterior_b1_Z": exterior_z[1],
        "exterior_b2_Z": exterior_z[2],
    }
    normalized_heat = mp.fsum(
        (
            normalized_heat_terms["collar_heat_correction_theta"],
            exterior[1],
            exterior[2],
        )
    )
    normalized_heat_z = mp.fsum(
        (
            normalized_heat_terms["collar_heat_correction_theta_Z"],
            exterior_z[1],
            exterior_z[2],
        )
    )
    normalized_total = mp.fsum((normalized_reference, normalized_heat))
    normalized_total_z = normalized_heat_z
    return {
        "delta": delta,
        "h": h,
        "xi": xi,
        "xi_Z": xi_z,
        "E_tail": e_tail,
        "scale": scale,
        "b": b,
        "b_Z": b_z,
        "normalized_reference_terms": normalized_reference_terms,
        "normalized_heat_terms": normalized_heat_terms,
        "normalized_reference": normalized_reference,
        "normalized_reference_Z": mp.mpf(0),
        "normalized_heat_correction": normalized_heat,
        "normalized_heat_correction_Z": normalized_heat_z,
        "normalized_total": normalized_total,
        "normalized_total_Z": normalized_total_z,
        "reference_increment": scale * normalized_reference,
        "reference_increment_Z": mp.mpf(0),
        "heat_correction_increment": scale * normalized_heat,
        "heat_correction_increment_Z": scale * normalized_heat_z,
        "tail_increment": scale * normalized_total,
        "tail_increment_Z": scale * normalized_total_z,
    }


def _anchor_components(anchor: dict[str, Any]) -> dict[str, Any]:
    """Keep the actual Rtail anchor atoms visible in the target receipt."""

    return {
        "theta": _n(anchor["theta"]),
        "theta_Z": _n(anchor["theta_Z"]),
        "baseline_normalized": list(anchor.get("baseline_normalized", ())),
        "bump_normalized": list(anchor.get("bump_normalized", ())),
        "inner_offsets": list(anchor.get("inner_offsets", ())),
        "inner_offsets_reapplied_once": anchor.get(
            "inner_offsets_reapplied_once", False
        ),
    }


class ContinuousAngularTailTarget:
    """Extract ``Ctheta`` and its jet from an installed continuous profile."""

    def __init__(self, profile):
        if hasattr(profile, "outer") and not hasattr(profile, "angular_moment_provider"):
            profile = profile.outer
        self.profile = profile
        self.schedule = profile.schedule
        self.precision = int(profile.precision)
        self.angular = profile.angular_moment_provider
        self.heat = self.angular.heat_provider
        # The correction object is shared by the actual angular provider.  It
        # is retained as provenance even though the complete anchor API owns
        # its two compact bump atoms.
        self.correction = getattr(profile, "angular_correction_provider", None)

    def _terminal_mixed(self, z: mp.mpf) -> dict[str, mp.mpf]:
        provider = getattr(self.profile, "mixed_moment_provider", None)
        if provider is not None:
            row = provider.moments_jet(
                mp.nstr(_n(self.schedule.logR_v), self.precision), z
            )
        else:
            method = getattr(self.profile, "linear_axial_moments_jet", None)
            if method is None:
                raise AttributeError(
                    "actual mixed moment provider is required for angular coupling"
                )
            row = method(mp.nstr(_n(self.schedule.logR_v), self.precision), z)
        return {
            "theta_z": _n(row.get("theta_z", row.get("Mtheta_z"))),
            "theta_z_Z": _n(row.get("theta_z_Z", row.get("Mtheta_z_Z"))),
            "method": row.get("method", "actual_mixed_provider_Rv"),
        }

    def evaluate(self, Z: Any, *, include_inner: bool = True) -> dict[str, Any]:
        with mp.workdps(self.precision):
            z = _n(Z)
            if not mp.isfinite(z) or abs(z) >= 1:
                raise ValueError("actual angular target requires finite |Z| < 1")
            s = self.schedule
            rtail = _n(s.logR_tail)
            tail_key = mp.nstr(rtail, self.precision)
            anchor = self.angular.moments_jet(
                tail_key, z, include_inner=include_inner
            )
            collar = self.heat.increments(mp.mpf(3), z)
            heat = self.profile.angular_schedule_provider.heat_jet(s.logR_tail, z)
            terms = _tail_terms(
                delta=s.delta,
                log_r_tail=rtail,
                log_c_inf=s._log_c_inf,
                xi=heat["xi"],
                xi_z=heat["xi_Z"],
                collar_ref_theta=collar[0],
                collar_heat_correction_theta=collar[2],
                collar_heat_correction_theta_z=collar[4],
            )
            anchor_parts = _anchor_components(anchor)
            ctheta = mp.fsum((anchor_parts["theta"], terms["tail_increment"]))
            ctheta_z = mp.fsum(
                (anchor_parts["theta_Z"], terms["tail_increment_Z"])
            )

            mixed = self._terminal_mixed(z)
            delta = terms["delta"]
            h = terms["h"]
            coupling_terms = {
                "(1-h)*Ctheta": (1 - h) * ctheta,
                "-(1-delta)*Z*Ctheta_Z/2": -(1 - delta) * z * ctheta_z / 2,
                "-(1-Z^2)*terminal_Mtheta_z_Z": -(1 - z * z) * mixed["theta_z_Z"],
                "+(2*delta-1)*Z*terminal_Mtheta_z": (2 * delta - 1) * z * mixed["theta_z"],
            }
            coupling = mp.fsum(coupling_terms.values())
            return {
                "Z": z,
                "Ctheta": ctheta,
                "Ctheta_Z": ctheta_z,
                "theta_tail_constant": ctheta,
                "theta_tail_constant_Z": ctheta_z,
                "actual_anchor": anchor_parts,
                "normalization": {
                    "reference_Utheta": "Cinf*R^(-(1+delta)/2)",
                    "angular_reference_integral": "sqrt(2)*Cinf*R^(1-h)/(1-h)",
                    "E_tail": terms["E_tail"],
                    "scale": terms["scale"],
                },
                "complete_heat_angular_defect": terms,
                "coupling": {
                    "terminal_Mtheta_z": mixed["theta_z"],
                    "terminal_Mtheta_z_Z": mixed["theta_z_Z"],
                    "terms": coupling_terms,
                    "value": coupling,
                    "formula": "(1-h)Ctheta-(1-delta)Z Ctheta_Z/2-(1-Z^2) terminal_Mtheta_z_Z+(2delta-1)Z terminal_Mtheta_z",
                    "provider_method": mixed["method"],
                },
                "continuous_angular_correction_shared": self.correction is not None,
                "matching_target_only": True,
                "angular_field_overwritten": False,
                "terminal_angular_closure_certified": False,
                "quadrature_error_enclosed": False,
                "exterior_model": "implemented quadratic heat atoms; exact heat remainder not enclosed",
            }

    __call__ = evaluate
    tail_constant_jet = evaluate


def complete_heat_angular_defect(
    profile, Z: Any, *, include_inner: bool = True
) -> dict[str, Any]:
    """Return the actual angular tail constant and matching coupling."""

    return ContinuousAngularTailTarget(profile).evaluate(
        Z, include_inner=include_inner
    )


angular_tail_target = complete_heat_angular_defect
actual_angular_tail_target = complete_heat_angular_defect
actual_angular_tail_constant = complete_heat_angular_defect
angular_tail_constant = complete_heat_angular_defect


def analytic_exponential_fixture() -> dict[str, Any]:
    """Independent exponential-tail check without constructing the full field."""

    with mp.workdps(90):
        delta = mp.mpf(".5")
        h = delta / 2
        xi = mp.mpf(".03")
        xi_z = mp.mpf("-.02")
        collar_ref = mp.mpf(".7")
        collar_heat = mp.mpf("-.0012")
        collar_heat_z = mp.mpf(".00017")
        pieces = _tail_terms(
            delta=delta,
            log_r_tail=mp.mpf(".3"),
            log_c_inf=mp.mpf("-.2"),
            xi=xi,
            xi_z=xi_z,
            collar_ref_theta=collar_ref,
            collar_heat_correction_theta=collar_heat,
            collar_heat_correction_theta_z=collar_heat_z,
        )
        b1, b2 = pieces["b"][1], pieces["b"][2]
        b1z, b2z = pieces["b_Z"][1], pieces["b_Z"][2]
        q = 1 - h
        # Integrate the decaying variable u=t-3 rather than subtracting two
        # growing primitives.  This is an independent check of the reference
        # endpoint exp(3*q)/q used by the constant-defect formula.
        direct_reference = mp.exp(3 * q) * mp.quad(
            lambda u: mp.exp(-q * u), [0, mp.inf]
        )
        direct_correction = mp.quad(
            lambda t: b1 * mp.exp((q - 1) * t) + b2 * mp.exp((q - 2) * t),
            [3, mp.inf],
        )
        direct_correction_z = mp.quad(
            lambda t: b1z * mp.exp((q - 1) * t) + b2z * mp.exp((q - 2) * t),
            [3, mp.inf],
        )
        expected = collar_ref - direct_reference + collar_heat + direct_correction
        expected_z = collar_heat_z + direct_correction_z
        err = abs(expected - pieces["normalized_total"])
        err_z = abs(expected_z - pieces["normalized_total_Z"])
        scale = max(abs(expected), abs(pieces["normalized_total"]), mp.mpf("1e-80"))
        scale_z = max(abs(expected_z), abs(pieces["normalized_total_Z"]), mp.mpf("1e-80"))
        return {
            "normalized_total_relative_error": mp.nstr(err / scale, 30),
            "normalized_total_Z_relative_error": mp.nstr(err_z / scale_z, 30),
            "passed": err / scale < mp.mpf("1e-70") and err_z / scale_z < mp.mpf("1e-70"),
            "independent_integral": "mp.quad exponential exterior atoms",
        }


if __name__ == "__main__":
    print(json.dumps(analytic_exponential_fixture(), indent=2))
