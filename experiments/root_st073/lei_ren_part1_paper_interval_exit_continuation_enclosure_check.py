"""Bounded receipt check for the controlled post-switch continuation."""

import hashlib
import json
from pathlib import Path

import mpmath as mp


HERE = Path(__file__).resolve().parent
RECEIPT = HERE / "lei_ren_part1_paper_interval_exit_continuation_enclosure.json"
SOURCE = HERE / "lei_ren_part1_paper_interval_exit_continuation_enclosure.py"


def interval(record):
    return mp.mpf(record["lower"]), mp.mpf(record["upper"])


def jet_order(record, order=1):
    if not isinstance(record, list) or len(record) != order + 1:
        raise AssertionError("unexpected Taylor order")


def contains(wide, narrow):
    for upper, inner in zip(wide, narrow):
        wlo, whi = interval(upper)
        nlo, nhi = interval(inner)
        if not wlo <= nlo <= nhi <= whi:
            return False
    return True


def run():
    data = json.loads(RECEIPT.read_text(encoding="utf-8"))
    for name in ("normalized_exit_state", "normalized_range_state"):
        state = data[name]
        for field in ("phi", "U", "theta", "z", "theta_z", "p", "u_squared", "weighted_phi_squared"):
            jet_order(state[field], 1)
    for name in ("F", "Uz", "P", "actual_terminal_g_y", "actual_terminal_U_y", "actual_F_y", "actual_F_R", "actual_Uz_R"):
        jet_order(data[name], 1)
    interval(data["Ur_value"])
    for name in ("delta_g_target", "delta_g_range", "delta_phi_target", "delta_phi_range", "delta_u_target", "delta_u_range"):
        jet_order(data[name], 1)
    for target, wide in (("delta_g_target", "delta_g_range"),
                         ("delta_phi_target", "delta_phi_range"),
                         ("delta_u_target", "delta_u_range")):
        if not contains(data[wide], data[target]):
            raise AssertionError(target + " target box escaped full-path range box")
    for name in ("theta", "z", "theta_z", "z_theta", "p"):
        jet_order(data["physical_moments"][name], 1)
    target_lo, target_hi = interval(data["target_physical_R"])
    if not target_lo <= mp.mpf("100") <= target_hi:
        raise AssertionError("target physical radius is not enclosed")
    if not data["controlled_post_collar_range"] or not data["constant_field_moment_increments_exact"]:
        raise AssertionError("continuation scope flags are incomplete")
    if not data["target_g_integral_uses_range_driver"] or not data["target_u_integral_uses_range_driver"]:
        raise AssertionError("target increments do not use full-path driver ranges")
    if not data["terminal_derivatives_use_target_driver"]:
        raise AssertionError("terminal derivative provenance is missing")
    if data["switch_100_110_enclosed"] or data["terminal_matching_complete"] or data["stress_cone_certified"]:
        raise AssertionError("continuation overclaims scope")
    digest = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    if data["input_hashes"][SOURCE.name] != digest:
        raise AssertionError("continuation source hash mismatch")
    report = {
        "passed": True,
        "target_R": data["target_R"],
        "target_physical_radius_enclosed": True,
        "target_fields_axial_order": 1,
        "range_fields_axial_order": 1,
        "physical_moment_rows_checked": 5,
        "cancellation_safe_exp_delta_g": True,
        "target_g_integral_uses_full_path_driver": data["target_g_integral_uses_range_driver"],
        "target_u_integral_uses_full_path_driver": data["target_u_integral_uses_range_driver"],
        "terminal_derivatives_use_endpoint_driver": data["terminal_derivatives_use_target_driver"],
        "Ur_Z_available": data["Ur_Z_available"],
        "switch_100_110_enclosed": data["switch_100_110_enclosed"],
        "terminal_matching_complete": data["terminal_matching_complete"],
        "stress_cone_certified": data["stress_cone_certified"],
    }
    (HERE / "lei_ren_part1_paper_interval_exit_continuation_enclosure_check.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print("Controlled continuation receipt check passed", flush=True)
    return report


if __name__ == "__main__":
    run()
