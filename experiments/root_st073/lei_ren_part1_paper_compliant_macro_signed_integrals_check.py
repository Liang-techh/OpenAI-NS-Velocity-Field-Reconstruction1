"""Independent receipt checker for the leading formal-width macro integrals.

The checker deliberately does not import the producer's modal algebra.  It
rebuilds the six moment modes, the signed direction and drive integrals, and
one high-precision quadrature directly from the accepted comparison receipt.
The width h_b remains a formal symbol throughout.
"""

import hashlib
import json
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
PRODUCER_NAME = PREFIX + "macro_signed_integrals.json"
PRODUCER_SOURCE_NAME = PREFIX + "macro_signed_integrals.py"
CHECK_SOURCE_NAME = PREFIX + "macro_signed_integrals_check.py"


def _sha256(path):
    return hashlib.sha256((HERE / path).read_bytes()).hexdigest()


def _verify_hashes(receipt):
    for path, digest in receipt.get("input_hashes", {}).items():
        if _sha256(path) != digest:
            raise ValueError("Accepted source changed: " + path)


def _jet(c, row):
    return IntervalTaylor(c, [read_interval(c, item) for item in row])


def _zero(value):
    return value * 0


def _add(left, right):
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


def _sub(left, right):
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


def _packet_for_coordinate(c, packets, coordinate):
    for packet in packets:
        lo, hi = endpoints(read_interval(c, packet["coordinate"]))
        if lo == coordinate and hi == coordinate:
            return packet
    raise ValueError("Comparison macro packet coordinate not found")


def _moments(c, packet):
    raw = packet["comparison_own_six_moments"]
    # The source receipt names the six scalar histories H,M,K,A,B,C.  The
    # direction API uses angular, axial, angular_axial, axial_quadratic and
    # pressure, with A and B occupying the two axial_quadratic slots.
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


def _phi(c, packet):
    return _jet(c, packet["comparison_phi"]["signed_hb_power_axial_coefficients"][0])


def _value(c, packet):
    return _jet(c, packet["comparison_raw_V"]["signed_hb_power_axial_coefficients"][0])


def _target(phi, value):
    return {
        MTH: phi,
        MZ: value,
        MTHZ: phi * value,
        MZT: {"axial": value * value, "swirl": phi * phi / 2},
        MP: phi * phi,
    }


def _direction(bridge, z, inputs, phi, value, moments):
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


def _modes(bridge, z, inputs, phi, value, inlet):
    target = _target(phi, value)
    difference = _sub(inlet, target)
    one_difference = {
        MTH: _zero(phi),
        MZ: difference[MZ],
        MTHZ: _zero(phi),
        MZT: {"axial": difference[MZT]["axial"], "swirl": _zero(phi)},
        MP: difference[MP],
    }
    two_difference = {
        MTH: difference[MTH],
        MZ: _zero(phi),
        MTHZ: difference[MTHZ],
        MZT: {"axial": _zero(phi), "swirl": difference[MZT]["swirl"]},
        MP: _zero(phi),
    }
    base = _direction(bridge, z, inputs, phi, value, target)
    one = _direction(bridge, z, inputs, phi, value, _add(target, one_difference))
    two = _direction(bridge, z, inputs, phi, value, _add(target, two_difference))
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


def _radius_kernels(c, radius, length, power):
    # Algebraically simplify exp(log(100/Ra)) before interval arithmetic.
    if power == 1:
        return [c.mpf(100) - radius, radius * length, radius - radius * radius / 100]
    if power == 2:
        return [(c.mpf(10000) - radius * radius) / 2,
                100 * radius - radius * radius, radius * radius * length]
    raise ValueError("Unexpected original bridge radius power")


def _weighted(rows, kernels):
    result = rows[0] * 0
    for row, kernel in zip(rows, kernels):
        result += row * kernel
    return result


def _output_jet(c, row):
    return _jet(c, row)


def _contains(stored, expected):
    """Whether every expected interval is enclosed by stored interval."""
    if len(stored.coefficients) != len(expected.coefficients):
        return False
    for lhs, rhs in zip(stored.coefficients, expected.coefficients):
        slo, shi = endpoints(lhs)
        elo, ehi = endpoints(rhs)
        if slo > elo or shi < ehi:
            return False
    return True


def _overlap(left, right):
    for lhs, rhs in zip(left.coefficients, right.coefficients):
        llo, lhi = endpoints(lhs)
        rlo, rhi = endpoints(rhs)
        if lhi < rlo or rhi < llo:
            return False
    return True


def _midpoint(value):
    lo, hi = endpoints(value)
    return (lo + hi) / 2


def _midpoint_jet(value):
    c = value.ctx
    return [_midpoint(row) for row in value.coefficients]


def _sigma_mp(x):
    """Independent scalar pulse evaluator used only in the micro quadrature."""
    if x <= 0:
        return mp.mpf(0)
    if x >= 1:
        return mp.mpf(1)
    if x > mp.mpf("0.5"):
        return 1 - _sigma_mp(1 - x)
    exponent = 1 / (1 - x) ** 2 - 1 / x ** 2
    value = mp.exp(exponent)
    return value / (1 + value)


class CompliantMacroSignedIntegralCheck:
    def __init__(self):
        self.bridge = CompliantInnerBridgeProfiles()
        self.ctx = c = self.bridge.ctx
        self.producer = json.loads((HERE / PRODUCER_NAME).read_bytes())
        self.comparison = json.loads((HERE / COMPARISON_NAME).read_bytes())
        comparison_check = json.loads((HERE / COMPARISON_CHECK_NAME).read_bytes())
        _verify_hashes(comparison_check)
        _verify_hashes(self.comparison)
        _verify_hashes(self.producer)
        if not comparison_check.get("all_passed"):
            raise ValueError("Comparison-point receipt is not accepted")
        if not self.producer.get("leading_formal_hb_coefficient_only"):
            raise ValueError("Producer must be leading formal coefficient only")
        for flag in ("higher_hb_powers_resolved", "first_switch_resolved",
                     "actual_signed_bridge_completed"):
            if self.producer.get(flag):
                raise ValueError("Unexpected resolved-scope flag: " + flag)
        if self.producer["actual_five_defect_family_sha256"] != self.bridge.family:
            raise ValueError("Producer family hash mismatch")
        if self.producer["implicit_source_sha256"] != self.bridge.source:
            raise ValueError("Producer source hash mismatch")
        if self.producer["datum_enclosure_sha256"] != self.bridge.core.datum.datum_sha:
            raise ValueError("Producer datum hash mismatch")
        self.Y = c.ln(100 / self.bridge.r)
        self.Pstar2_log = 2 * self.bridge.core.logP
        self.F02_log = c.mpf([endpoints(-2*self.bridge.core.logC-2*self.bridge.core.Lambda*self.bridge.core.Gbar)[0],
                              endpoints(-2*self.bridge.core.logC)[1]])
        self.failures = []
        self.hashes = {
            COMPARISON_NAME: _sha256(COMPARISON_NAME),
            COMPARISON_CHECK_NAME: _sha256(COMPARISON_CHECK_NAME),
            PRODUCER_NAME: _sha256(PRODUCER_NAME),
            PRODUCER_SOURCE_NAME: _sha256(PRODUCER_SOURCE_NAME),
            CHECK_SOURCE_NAME: _sha256(CHECK_SOURCE_NAME),
        }

    def require(self, condition, message):
        if not condition:
            self.failures.append(message)

    def rebuild(self, label):
        c = self.ctx
        entries = self.comparison["comparison_point_packets"][label]["macro"]
        start = _packet_for_coordinate(c, entries, mp.mpf(0))
        end = _packet_for_coordinate(c, entries, mp.mpf(1))
        z = read_interval(c, start["Z"])
        phi = _phi(c, start)
        value = _value(c, start)
        inlet = _moments(c, start)
        endpoint = _moments(c, end)
        target = _target(phi, value)
        inputs = self.bridge.inputs(z)
        rows = _modes(self.bridge, z, inputs, phi, value, inlet)
        core_direction = _direction(self.bridge, z, inputs, phi, value, inlet)
        micro_weight = c.mpf("0.5")
        radius = self.bridge.r
        d_rows = rows["D_over_R"]
        jf_macro = _weighted(d_rows, _radius_kernels(c, radius, self.Y, 1)) * (-c.mpf("0.5"))
        jf_micro = core_direction["D_over_R"] * (-c.mpf("0.5") * micro_weight * radius)
        jf = jf_micro + jf_macro
        radial_kernels = _radius_kernels(c, radius, self.Y, 1)
        swirl_kernels = _radius_kernels(c, radius, self.Y, 2)
        hydro_rows = rows["drive_hydro"]
        pressure_rows = rows["drive_pressure"]
        swirl_rows = rows["drive_swirl"]
        jvh_macro = -_weighted(hydro_rows, radial_kernels)
        jvp_macro = -_weighted(pressure_rows, radial_kernels)
        jvs_macro = -_weighted(swirl_rows, swirl_kernels)
        jvh_micro = core_direction["drive_hydro"] * (-micro_weight * radius)
        jvp_micro = core_direction["drive_pressure"] * (-micro_weight * radius)
        jvs_micro = core_direction["drive_swirl"] * (-micro_weight * radius * radius)
        jvh = jvh_micro + jvh_macro
        jvp = jvp_micro + jvp_macro
        jvs = jvs_micro + jvs_macro
        expected_endpoint = {}
        for name, rate, source, target_value in (
            ("H", 2, inlet[MTH], target[MTH]),
            ("M", 1, inlet[MZ], target[MZ]),
            ("K", 2, inlet[MTHZ], target[MTHZ]),
            ("A", 1, inlet[MZT]["axial"], target[MZT]["axial"]),
            ("B", 2, inlet[MZT]["swirl"], target[MZT]["swirl"]),
            ("C", 1, inlet[MP], target[MP]),
        ):
            expected_endpoint[name] = target_value + (source - target_value) * c.exp(-rate * self.Y)
        return dict(label=label, z=z, start=start, end=end, endpoint=endpoint,
                    inlet=inlet, target=target, rows=rows, d_rows=d_rows,
                    radius=radius, jf=jf, jf_micro=jf_micro, jf_macro=jf_macro,
                    core_direction=core_direction, hydro_rows=hydro_rows, pressure_rows=pressure_rows,
                    swirl_rows=swirl_rows, jvh=jvh, jvp=jvp, jvs=jvs,
                    components={"hydro": (jvh, jvh_micro, jvh_macro),
                                "pressure": (jvp, jvp_micro, jvp_macro),
                                "swirl": (jvs, jvs_micro, jvs_macro)},
                    expected_endpoint=expected_endpoint)

    def check_packet(self, data):
        label = data["label"]
        stored = self.producer["leading_signed_bridge_packets"][label]
        c = self.ctx
        stored_jf = _output_jet(c, stored["J_F"])
        self.require(_contains(stored_jf, data["jf"]), label + ": J_F enclosure")
        self.require(_contains(_output_jet(c, stored["J_F_micro"]), data["jf_micro"]),
                     label + ": first microscopic J_F")
        self.require(_contains(_output_jet(c, stored["J_F_macro"]), data["jf_macro"]),
                     label + ": macro J_F")
        self.require(stored.get("second_micro_leading_coefficient_zero") is True,
                     label + ": second microscopic chart leading order")
        for component, expected, scale_log in (
            ("hydro", data["jvh"], c.mpf(0)),
            ("pressure", data["jvp"], self.Pstar2_log),
            ("swirl", data["jvs"], self.F02_log),
        ):
            item = stored["J_V_" + component]
            self.require(_contains(_output_jet(c, item["normalized_axial_coefficients"]), expected),
                         label + ": J_V " + component + " factored coefficient")
            actual_log = read_interval(c, item["positive_scale_log"])
            elo, ehi = endpoints(scale_log); alo, ahi = endpoints(actual_log)
            self.require(alo <= elo and ahi >= ehi, label + ": J_V " + component + " full source scale log")
            self.require(item.get("positive_scale_log_is_enclosure") is True,
                         label + ": source amplitude log is an enclosure")
            self.require(stored["J_V"]["terms"][component] == item,
                         label + ": J_V signed factored sum binding")
            _, micro_expected, macro_expected = data["components"][component]
            self.require(_contains(_output_jet(c, item["micro_normalized_axial_coefficients"]), micro_expected),
                         label + ": J_V " + component + " micro contribution")
            self.require(_contains(_output_jet(c, item["macro_normalized_axial_coefficients"]), macro_expected),
                         label + ": J_V " + component + " macro contribution")
        self.require(stored["J_V"].get("no_amplitude_scale_materialized") is True,
                     label + ": amplitude scale must remain factored")
        stored_modes = stored["D_over_R_modes"]
        self.require(len(stored_modes) == 3, label + ": Dbar mode count")
        for index, row in enumerate(data["d_rows"]):
            self.require(_contains(_output_jet(c, stored_modes[index]), row),
                         label + ": Dbar mode " + str(index))
        for name, key in (("drive_hydro", "hydro_rows"),
                          ("drive_pressure", "pressure_rows"),
                          ("drive_swirl", "swirl_rows")):
            stored_rows = stored["drive_" + name.split("drive_")[-1] + "_modes"]
            for index, row in enumerate(data[key]):
                self.require(_contains(_output_jet(c, stored_rows[index]), row),
                             label + ": " + name + " mode " + str(index))
        # Check the actual receipt endpoint against the independently derived
        # zero-width modal endpoint, coefficient by coefficient.
        for name, key in (("H", MTH), ("M", MZ), ("K", MTHZ),
                          ("A", (MZT, "axial")), ("B", (MZT, "swirl")),
                          ("C", MP)):
            actual = data["endpoint"][key] if isinstance(key, str) else data["endpoint"][key[0]][key[1]]
            expected = data["expected_endpoint"][name]
            self.require(_overlap(actual, expected), label + ": endpoint identity " + name)
            saved_expected = _output_jet(c, stored["comparison_R100_endpoint_identity"][name])
            self.require(_contains(saved_expected, expected), label + ": saved endpoint identity " + name)

    def independent_quadrature(self, data):
        """Use midpoint coefficients and mp.quad, independent of interval algebra."""
        c = self.ctx
        stored = self.producer["leading_signed_bridge_packets"][data["label"]]
        modes = [_output_jet(c, row) for row in stored["D_over_R_modes"]]
        radius = _midpoint(data["radius"])
        alpha = radius / 100
        y = _midpoint(read_interval(c, stored["Y"]))
        coeff = [_midpoint(row[0]) for row in modes]
        with mp.workdps(180):
            # Independent regularized substitutions in x=R/100.  They avoid
            # integrating across the ~10^18-wide logarithmic interval.
            q0 = 100 * coeff[0] * mp.quad(lambda x: mp.mpf(1), [alpha, 1])
            q1 = 100 * coeff[1] * alpha * y * mp.quad(lambda t: mp.mpf(1), [0, 1])
            q2 = 100 * coeff[2] * alpha * mp.quad(lambda s: mp.mpf(1), [0, 1 - alpha])
            macro_quad = -mp.mpf("0.5") * (q0 + q1 + q2)
            micro_weight = mp.quad(lambda x: 1 - _sigma_mp(x), [0, mp.mpf("0.5"), 1])
            core_d = _midpoint(data["core_direction"]["D_over_R"].coefficients[0])
            micro_quad = -mp.mpf("0.5") * micro_weight * radius * core_d
            quad = macro_quad + micro_quad
            self.require(abs(micro_weight - mp.mpf("0.5")) < mp.mpf("1e-150"),
                         data["label"] + ": independent pulse micro weight quadrature")
        lo, hi = endpoints(_output_jet(c, stored["J_F"]).coefficients[0])
        tol = max(mp.mpf("1e-150"), abs(quad) * mp.mpf("1e-140"))
        self.require(lo - tol <= quad <= hi + tol,
                     data["label"] + ": independent high-precision J_F quadrature")

    def run(self):
        c = self.ctx
        producer_y = read_interval(c, self.producer["Y_exact_log_100_over_Ra"])
        pylo, pyhi = endpoints(producer_y)
        eylo, eyhi = endpoints(self.Y)
        self.require(max(pylo, eylo) <= min(pyhi, eyhi), "exact Y=log(100/Ra)")
        self.require(self.producer.get("formal_width_source") == "hb=cstar*K^-100; no cap substituted or inverted",
                     "formal width source")
        self.require(self.producer.get("old_inner_bridge_cumulative_cover_used") is False,
                     "old constant cover flag")
        self.require(self.producer.get("leading_formal_hb_full_Ra_to_R100") is True,
                     "leading formal coefficient includes core-to-R100 first chart")
        results = {}
        for label in (".5", "0", "exact_shared_root"):
            rebuilt = self.rebuild(label)
            self.check_packet(rebuilt)
            results[label] = {
                "endpoint_identities_checked": 6,
                "mode_integrals_checked": 12,
                "label": label,
            }
        self.independent_quadrature(self.rebuild(".5"))
        result = {
            "all_passed": not self.failures,
            "failures": self.failures,
            "producer_source": PRODUCER_SOURCE_NAME,
            "producer_json": PRODUCER_NAME,
            "comparison_source_namespace": self.comparison["original_comparison_source"]["source_namespace"],
            "actual_five_defect_family_sha256": self.bridge.family,
            "implicit_source_sha256": self.bridge.source,
            "datum_enclosure_sha256": self.bridge.core.datum.datum_sha,
            "formal_width_not_materialized": True,
            "higher_hb_powers_resolved": False,
            "first_switch_resolved": False,
            "actual_signed_bridge_completed": False,
            "independent_modal_rebuild": results,
            "independent_high_precision_quadrature": True,
            "input_hashes": self.hashes,
        }
        result["input_hashes"][CHECK_SOURCE_NAME] = _sha256(CHECK_SOURCE_NAME)
        return result


def run():
    with mp.workdps(320):
        checker = CompliantMacroSignedIntegralCheck()
        result = checker.run()
    path = HERE / (Path(__file__).stem + ".json")
    path.write_text(json.dumps(_encode(result), indent=2) + "\n", encoding="utf8")
    print("Leading signed macro integral checker:", "PASS" if result["all_passed"] else "FAIL", flush=True)
    if not result["all_passed"]:
        raise RuntimeError("; ".join(result["failures"]))
    return result


if __name__ == "__main__":
    run()
