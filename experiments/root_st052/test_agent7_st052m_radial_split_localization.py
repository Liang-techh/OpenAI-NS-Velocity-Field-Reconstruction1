from __future__ import annotations

import numpy as np

import agent7_st052m_radial_split_localization as audit


def _summary(n: int, inward: bool):
    count = n if inward else 0
    return {
        "path_count": n,
        "inward_path_count": count,
        "outward_or_zero_path_count": n - count,
        "mean_delta_r": -0.1 if inward else 0.1,
        "min_delta_r": -0.2 if inward else 0.05,
        "max_delta_r": -0.05 if inward else 0.2,
        "mean_start_radial_velocity": -0.2 if inward else 0.2,
        "min_start_radial_velocity": -0.3 if inward else 0.1,
        "max_start_radial_velocity": -0.1 if inward else 0.3,
        "uniform_inward_label": True,
        "all_inward": inward,
        "all_outward_or_zero": not inward,
    }


def _records(sign_pattern=(True, False), mixed_band=True):
    keys = ("0.250-0.375", "0.375-0.500", "0.500-0.625", "0.625-0.750")
    result = {}
    for field in ("linear", "child"):
        segments = {}
        for key in keys:
            sign = {"-1": _summary(24, sign_pattern[0]), "1": _summary(24, sign_pattern[1])}
            if mixed_band:
                shoulder = _summary(24, False)
                shoulder["uniform_inward_label"] = False
                shoulder["inward_path_count"] = 12
                shoulder["outward_or_zero_path_count"] = 12
                shoulder["all_outward_or_zero"] = False
                tip = dict(shoulder)
                band = {"shoulder": shoulder, "tip": tip}
            else:
                band = {"shoulder": _summary(24, True), "tip": _summary(24, False)}
            radius = {}
            for r in ("0.6", "0.9", "1.2"):
                row = _summary(16, False)
                row["uniform_inward_label"] = False
                row["inward_path_count"] = 8
                row["outward_or_zero_path_count"] = 8
                row["all_outward_or_zero"] = False
                radius[r] = row
            segments[key] = {"groups": {"sign": sign, "band": band, "radius": radius}}
        result[field] = {"segments": segments}
    return result


def test_group_summary_localizes_sign_and_preserves_counts():
    delta_r = np.array([-0.2, -0.1, 0.1, 0.2])
    start_ur = np.array([-0.3, -0.2, 0.2, 0.3])
    metadata = [
        {"sign": -1, "band": "shoulder", "radius": 0.6, "angle_index": 0},
        {"sign": -1, "band": "tip", "radius": 0.9, "angle_index": 1},
        {"sign": 1, "band": "shoulder", "radius": 0.6, "angle_index": 2},
        {"sign": 1, "band": "tip", "radius": 0.9, "angle_index": 3},
    ]
    grouped = audit._group_summary(delta_r, start_ur, metadata, "sign")
    assert grouped["-1"]["path_count"] == 2
    assert grouped["-1"]["all_inward"] is True
    assert grouped["1"]["all_outward_or_zero"] is True
    assert grouped["-1"]["mean_start_radial_velocity"] < 0.0
    assert grouped["1"]["mean_start_radial_velocity"] > 0.0


def test_preregistered_decision_identifies_exact_z_sign_parity_split():
    decision = audit.localization_decision(_records())
    assert decision["z_sign_exact_split"] is True
    assert decision["inward_z_sign"] == "-1"
    assert decision["classification"] == "z_sign_parity"
    assert decision["basis_growth_justified_by_this_audit"] is False


def test_decision_rejects_sign_pattern_drift_and_falls_back_to_mixed():
    records = _records()
    records["child"]["segments"]["0.500-0.625"]["groups"]["sign"] = {
        "-1": _summary(24, False),
        "1": _summary(24, True),
    }
    decision = audit.localization_decision(records)
    assert decision["z_sign_exact_split"] is False
    assert decision["classification"] == "mixed"


def test_decision_can_route_exact_band_local_pattern_without_claiming_z_parity():
    records = _records(sign_pattern=(True, True), mixed_band=False)
    # Make sign groups nontriviality fail while bands remain exact and stable.
    for field in ("linear", "child"):
        for row in records[field]["segments"].values():
            row["groups"]["sign"] = {"-1": _summary(24, True), "1": _summary(24, True)}
    decision = audit.localization_decision(records)
    assert decision["z_sign_exact_split"] is False
    assert decision["classification"] == "band_local"


def test_truth_boundary_is_nonpromotional_and_parent_is_frozen():
    assert audit.PREREG_ISSUE == 691
    assert audit.SOURCE_PARENT_PR == 683
    assert audit.SOURCE_PARENT_HEAD == "c8046f42814a45ecaa9dbe44d5a0521ee0289948"
    assert audit.parent.witness.SWIRL_A == 0.065899695471146
    assert abs(audit.parent.witness.SHOULDER_LAMBDA + 0.02) <= 1.0e-15
    assert audit.TRUTH["new_spatial_basis_added"] is False
    assert audit.TRUTH["new_temporal_basis_added"] is False
    assert audit.TRUTH["witness_retuned"] is False
    assert audit.TRUTH["visual_correspondence_verified"] is False
    assert audit.TRUTH["pde_validated"] is False
    assert audit.TRUTH["openai_field_identified"] is False
