"""Leading formal-width signed integrals on the frozen comparison macro.

This module is intentionally narrower than the actual bridge provider.  It reads
the accepted coefficientwise comparison-history receipt, reconstructs the exact
zero-width exponential moment modes after the 2 h_b smoothing chart, and
integrates the signed (9.13) direction and axial-drive modes from R_a to 100.
The returned coefficients are the formal leading terms

    log(F_100/F_2) = h_b * J_F + O(h_b^2),
    V_100 - V_2     = h_b * J_V + O(h_b^2).

No numerical cap is used as h_b and no h_b is materialized.  Higher width
orders, the first switch, and the complete actual bridge remain unresolved.
"""

import hashlib
import json
import math
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_compliant_inner_bridge_profiles import (
    CompliantInnerBridgeProfiles,
    IntervalTaylor,
    MTH,
    MZ,
    MTHZ,
    MZT,
    MP,
    direction,
)
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_candidate_exact_amplitude import _encode


HERE = Path(__file__).parent
PREFIX = "lei_ren_part1_paper_compliant_"
COMPARISON_NAME = PREFIX + "comparison_point_integrals.json"
COMPARISON_CHECK_NAME = PREFIX + "comparison_point_integrals_check.json"


def _sha256(path):
    return hashlib.sha256((HERE / path).read_bytes()).hexdigest()


def _verify_hashes(receipt):
    for path, digest in receipt.get("input_hashes", {}).items():
        if _sha256(path) != digest:
            raise ValueError("Accepted source changed: " + path)


def _jet(c, row):
    return IntervalTaylor(c, [read_interval(c, value) for value in row])


def _jet_list(value):
    return list(value.coefficients)


def _zero(value):
    return value * 0


def _zero_moments(phi):
    z = _zero(phi)
    return {
        MTH: z,
        MZ: z,
        MTHZ: z,
        MZT: {"axial": z, "swirl": z},
        MP: z,
    }


def _copy_moments(value):
    return {
        MTH: value[MTH],
        MZ: value[MZ],
        MTHZ: value[MTHZ],
        MZT: {"axial": value[MZT]["axial"], "swirl": value[MZT]["swirl"]},
        MP: value[MP],
    }


def _add_moments(left, right):
    return {
        MTH: left[MTH] + right[MTH],
        MZ: left[MZ] + right[MZ],
        MTHZ: left[MTHZ] + right[MTHZ],
        MZT: {
            "axial": left[MZT]["axial"] + right[MZT]["axial"],
            "swirl": left[MZT]["swirl"] + right[MZT]["swirl"],
        },
        MP: left[MP] + right[MP],
    }


def _sub_moments(left, right):
    return {
        MTH: left[MTH] - right[MTH],
        MZ: left[MZ] - right[MZ],
        MTHZ: left[MTHZ] - right[MTHZ],
        MZT: {
            "axial": left[MZT]["axial"] - right[MZT]["axial"],
            "swirl": left[MZT]["swirl"] - right[MZT]["swirl"],
        },
        MP: left[MP] - right[MP],
    }


def _scale_moments(value, scalar):
    return {
        MTH: value[MTH] * scalar,
        MZ: value[MZ] * scalar,
        MTHZ: value[MTHZ] * scalar,
        MZT: {
            "axial": value[MZT]["axial"] * scalar,
            "swirl": value[MZT]["swirl"] * scalar,
        },
        MP: value[MP] * scalar,
    }


def _packet_for_coordinate(c, packets, coordinate):
    for packet in packets:
        value = read_interval(c, packet["coordinate"])
        lo, hi = endpoints(value)
        if lo == coordinate and hi == coordinate:
            return packet
    raise ValueError("Comparison macro packet coordinate not found")


def _packet_phi(c, packet):
    return _jet(c, packet["comparison_phi"]["signed_hb_power_axial_coefficients"][0])


def _packet_v(c, packet):
    return _jet(c, packet["comparison_raw_V"]["signed_hb_power_axial_coefficients"][0])


def _packet_moments_correct(c, packet):
    """Map receipt H,M,K,A,B,C rows to direction's moment names."""
    raw = packet["comparison_own_six_moments"]
    # comparison_point_integrals uses the canonical six scalar names H,M,K,A,B,C:
    # H=angular, M=axial, K=angular_axial, A=axial_quadratic,
    # B=swirl, C=pressure.
    return {
        MTH: _jet(c, raw["H"]["signed_hb_power_axial_coefficients"][0]),
        MZ: _jet(c, raw["M"]["signed_hb_power_axial_coefficients"][0]),
        MTHZ: _jet(c, raw["K"]["signed_hb_power_axial_coefficients"][0]),
        MZT: {
            "axial": _jet(c, raw["A"]["signed_hb_power_axial_coefficients"][0]),
            "swirl": _jet(c, raw["B"]["signed_hb_power_axial_coefficients"][0]),
        },
        MP: _jet(c, raw["C"]["signed_hb_power_axial_coefficients"][0]),
    }


def _target_moments(phi, value):
    return {
        MTH: phi,
        MZ: value,
        MTHZ: phi * value,
        MZT: {"axial": value * value, "swirl": phi * phi / 2},
        MP: phi * phi,
    }


def _mode_direction(bridge, z, inputs, phi, value, moments):
    return direction(
        bridge.ctx,
        z,
        bridge.delta,
        phi,
        value,
        moments,
        inputs["p0"],
        inputs["F0_ratios"],
        inputs["F0_squared_ratios"],
    )


def _mode_rows(bridge, z, inputs, phi, value, inlet):
    """Return zero-width direction modes in e^{-Delta}, e^{-2 Delta}."""
    target = _target_moments(phi, value)
    delta = _sub_moments(inlet, target)
    delta_one = {
        MTH: _zero(phi), MZ: delta[MZ], MTHZ: _zero(phi),
        MZT: {"axial": delta[MZT]["axial"], "swirl": _zero(phi)},
        MP: delta[MP],
    }
    delta_two = {
        MTH: delta[MTH], MZ: _zero(phi), MTHZ: delta[MTHZ],
        MZT: {"axial": _zero(phi), "swirl": delta[MZT]["swirl"]},
        MP: _zero(phi),
    }
    base = _mode_direction(bridge, z, inputs, phi, value, target)
    one = _mode_direction(bridge, z, inputs, phi, value, _add_moments(target, delta_one))
    two = _mode_direction(bridge, z, inputs, phi, value, _add_moments(target, delta_two))
    return {
        "D_over_R": [base["D_over_R"], one["D_over_R"] - base["D_over_R"],
                     two["D_over_R"] - base["D_over_R"]],
        "drive_hydro": [base["drive_hydro"], one["drive_hydro"] - base["drive_hydro"],
                        _zero(base["drive_hydro"])],
        "drive_pressure": [base["drive_pressure"], _zero(base["drive_pressure"]),
                           _zero(base["drive_pressure"])],
        "drive_swirl": [base["drive_swirl"], one["drive_swirl"] - base["drive_swirl"],
                        two["drive_swirl"] - base["drive_swirl"]],
    }


def _radius_kernels(c, radius, Y, radius_power):
    """Integrals of R^p exp(-j Delta), simplified with exp(Y)=100/Ra."""
    if radius_power == 1:
        return [c.mpf(100) - radius, radius * Y, radius - radius * radius / 100]
    if radius_power == 2:
        return [(c.mpf(10000) - radius * radius) / 2,
                100 * radius - radius * radius, radius * radius * Y]
    raise ValueError("Only the original R and R^2 source weights are supported")


def _weighted_modes(c, modes, kernels):
    result = modes[0] * 0
    for row, kernel in zip(modes, kernels):
        result += row * kernel
    return result


def _factored_component(value, scale_log):
    return dict(positive_scale_log=scale_log,
                normalized_axial_coefficients=_as_output(value),
                exact_representation="exp(positive_scale_log)*normalized_axial_coefficients")


def _as_output(value):
    return _jet_list(value)


class CompliantMacroSignedIntegrals:
    def __init__(self):
        print("Macro signed integrals: constructing admitted bridge source", flush=True)
        self.bridge = CompliantInnerBridgeProfiles()
        print("Macro signed integrals: bridge source ready", flush=True)
        self.ctx = c = self.bridge.ctx
        self.comparison = json.loads((HERE / COMPARISON_NAME).read_bytes())
        check = json.loads((HERE / COMPARISON_CHECK_NAME).read_bytes())
        _verify_hashes(check)
        _verify_hashes(self.comparison)
        if not check.get("all_passed"):
            raise ValueError("Accepted comparison-point receipt required")
        if self.comparison["actual_five_defect_family_sha256"] != self.bridge.family:
            raise ValueError("Comparison family does not match current bridge")
        if self.comparison["implicit_source_sha256"] != self.bridge.source:
            raise ValueError("Comparison source does not match current bridge")
        if self.comparison["datum_enclosure_sha256"] != self.bridge.core.datum.datum_sha:
            raise ValueError("Comparison datum does not match current bridge")
        self.Y = c.ln(100 / self.bridge.r)
        # Keep the extraordinary F0^2 scale logarithmic.  Materializing
        # exp(-2 log Cstar - 2 Lambda Gbar) can require an exponent outside
        # the numeric format; later factored algebra consumes this source.
        self.F02_log = -2 * self.bridge.core.logC - 2 * self.bridge.core.Lambda * self.bridge.core.Gbar
        self.hashes = dict(self.comparison.get("input_hashes", {}))
        self.hashes[COMPARISON_NAME] = _sha256(COMPARISON_NAME)
        self.hashes[COMPARISON_CHECK_NAME] = _sha256(COMPARISON_CHECK_NAME)
        for path in (
            "lei_ren_part1_paper_compliant_inner_bridge_profiles.py",
            "lei_ren_part1_paper_compliant_bridge_mixed_C4.py",
            "lei_ren_part1_paper_interval_taylor.py",
        ):
            self.hashes[path] = _sha256(path)

    def packet(self, label, z_text=None):
        c = self.ctx
        print("Macro signed integrals: reconstructing Z=" + label, flush=True)
        entries = self.comparison["comparison_point_packets"][label]["macro"]
        at_start = _packet_for_coordinate(c, entries, mp.mpf(0))
        at_end = _packet_for_coordinate(c, entries, mp.mpf(1))
        # The exact shared-root packet carries its certified nonzero root in
        # the receipt.  Never replace that source coordinate by the nominal
        # label ``0``; the root displacement is part of the common source.
        z = read_interval(c, at_start["Z"])
        phi = _packet_phi(c, at_start)
        value = _packet_v(c, at_start)
        inlet = _packet_moments_correct(c, at_start)
        endpoint = _packet_moments_correct(c, at_end)
        inputs = self.bridge.inputs(z)
        print("Macro signed integrals: point inputs ready for Z=" + label, flush=True)
        modes = _mode_rows(self.bridge, z, inputs, phi, value, inlet)
        radius = self.bridge.r
        # The radius is astronomically small and Y=log(100/Ra) is very large.
        # Simplify Ra*exp(Y) exactly before interval arithmetic; never form
        # exp(Y), whose intermediate exponent is needlessly enormous.
        D_modes = modes["D_over_R"]
        J_F = _weighted_modes(c, D_modes, _radius_kernels(c, radius, self.Y, 1)) * (-c.mpf("0.5"))

        # R * hydro, R * P0^2 * pressure, and R^2 * F0^2 * swirl.
        hydro_kernels = _radius_kernels(c, radius, self.Y, 1)
        swirl_kernels = _radius_kernels(c, radius, self.Y, 2)
        JV_hydro = _weighted_modes(c, modes["drive_hydro"], hydro_kernels) * -1
        JV_pressure_normalized = _weighted_modes(c, modes["drive_pressure"], hydro_kernels) * -1
        JV_swirl_normalized = _weighted_modes(c, modes["drive_swirl"], swirl_kernels) * -1
        pressure_log = 2 * self.bridge.core.logP
        hydro_term = _factored_component(JV_hydro, c.mpf(0))
        pressure_term = _factored_component(JV_pressure_normalized, pressure_log)
        swirl_term = _factored_component(JV_swirl_normalized, self.F02_log)
        JV = dict(exact_representation="sum(exp(scale_log_i)*normalized_i(Z))",
                  terms=dict(hydro=hydro_term, pressure=pressure_term, swirl=swirl_term),
                  no_amplitude_scale_materialized=True)

        target = _target_moments(phi, value)
        endpoint_expected = {}
        for name, rate, target_value in (
            ("H", 2, target[MTH]), ("M", 1, target[MZ]),
            ("K", 2, target[MTHZ]), ("A", 1, target[MZT]["axial"]),
            ("B", 2, target[MZT]["swirl"]), ("C", 1, target[MP]),
        ):
            source = {"H": inlet[MTH], "M": inlet[MZ], "K": inlet[MTHZ],
                      "A": inlet[MZT]["axial"], "B": inlet[MZT]["swirl"],
                      "C": inlet[MP]}[name]
            endpoint_expected[name] = target_value + (source - target_value) * c.exp(-rate * self.Y)

        return {
            "Z_label": label,
            "Z": z,
            "macro_coordinate": "Delta=log(R/Ra)-2hb",
            "Y": self.Y,
            "Ra": radius,
            "comparison_phi_two_hb": _as_output(phi),
            "comparison_V_two_hb": _as_output(value),
            "comparison_moments_two_hb": {
                MTH: _as_output(inlet[MTH]), MZ: _as_output(inlet[MZ]),
                MTHZ: _as_output(inlet[MTHZ]),
                MZT: {"axial": _as_output(inlet[MZT]["axial"]), "swirl": _as_output(inlet[MZT]["swirl"])},
                MP: _as_output(inlet[MP]),
            },
            "comparison_moments_R100_zero_width": {
                MTH: _as_output(endpoint[MTH]), MZ: _as_output(endpoint[MZ]),
                MTHZ: _as_output(endpoint[MTHZ]),
                MZT: {"axial": _as_output(endpoint[MZT]["axial"]), "swirl": _as_output(endpoint[MZT]["swirl"])},
                MP: _as_output(endpoint[MP]),
            },
            "comparison_target_moments": {
                MTH: _as_output(target[MTH]), MZ: _as_output(target[MZ]),
                MTHZ: _as_output(target[MTHZ]),
                MZT: {"axial": _as_output(target[MZT]["axial"]), "swirl": _as_output(target[MZT]["swirl"])},
                MP: _as_output(target[MP]),
            },
            "comparison_R100_endpoint_identity": {
                name: _as_output(expected)
                for name, expected in endpoint_expected.items()
            },
            "direction_mode_exponents": {
                "D_over_R": [0, -1, -2],
                "drive_hydro": [0, -1, -2],
                "drive_pressure": [0, -1, -2],
                "drive_swirl": [0, -1, -2],
            },
            "D_over_R_modes": [_as_output(row) for row in D_modes],
            "drive_hydro_modes": [_as_output(row) for row in modes["drive_hydro"]],
            "drive_pressure_modes": [_as_output(row) for row in modes["drive_pressure"]],
            "drive_swirl_modes": [_as_output(row) for row in modes["drive_swirl"]],
            "J_F": _as_output(J_F),
            "J_V_hydro": hydro_term,
            "J_V_pressure": pressure_term,
            "J_V_swirl": swirl_term,
            "J_V": JV,
            "formal_changes": {
                "log_F100_over_F2": "hb*J_F + O(hb^2)",
                "V100_minus_V2": "hb*J_V + O(hb^2)",
                "hb_source": "hb=cstar*K^-100; retained formally",
                "F0_squared_source_log": self.F02_log,
            },
            "signed_source_formulas": {
                "J_F": "-1/2*integral_0^Y Dbar(Delta)dDelta",
                "J_V": "-integral_0^Y [R*hydro+R*Pstar^2*pressure+R^2*F0^2*swirl]dDelta",
                "R": "Ra*exp(Delta)",
                "Dbar": "R*D_over_R from original (9.13)",
            },
            "comparison_receipt_modes_used": True,
            "old_constant_cover_not_used": True,
            "higher_hb_powers_resolved": False,
            "first_switch_resolved": False,
            "actual_signed_bridge_completed": False,
        }

    def report(self):
        packets = {
            ".5": self.packet(".5"),
            "0": self.packet("0"),
            "exact_shared_root": self.packet("exact_shared_root"),
        }
        return {
            "actual_five_defect_family_sha256": self.comparison["actual_five_defect_family_sha256"],
            "implicit_source_sha256": self.comparison["implicit_source_sha256"],
            "datum_enclosure_sha256": self.comparison["datum_enclosure_sha256"],
            "comparison_source_namespace": self.comparison["original_comparison_source"]["source_namespace"],
            "Y_exact_log_100_over_Ra": self.Y,
            "formal_width_source": "hb=cstar*K^-100; no cap substituted or inverted",
            "macro_signed_integral_packets": packets,
            "leading_formal_hb_coefficient_only": True,
            "higher_hb_powers_resolved": False,
            "first_switch_resolved": False,
            "actual_signed_bridge_completed": False,
            "old_inner_bridge_cumulative_cover_used": False,
            "input_hashes": self.hashes,
        }


def run():
    with mp.workdps(300):
        producer = CompliantMacroSignedIntegrals()
        result = producer.report()
    result["input_hashes"][Path(__file__).name] = _sha256(Path(__file__).name)
    path = HERE / (Path(__file__).stem + ".json")
    path.write_text(json.dumps(_encode(result), indent=2) + "\n", encoding="utf8")
    print("Leading signed macro comparison integrals generated", flush=True)
    return result


if __name__ == "__main__":
    run()
