from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
MODULE_PATH = HERE / "agent7_st052m_tip_return_concentration_audit.py"
SPEC = importlib.util.spec_from_file_location("agent7_tip_return_concentration_audit", MODULE_PATH)
assert SPEC and SPEC.loader
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


def test_frozen_lineage_and_scope() -> None:
    report = audit.evaluate()
    assert report["task_id"] == "CR003-ST052M-TIP-RETURN-CONCENTRATION-AUDIT-104"
    assert report["issue"] == 749
    assert report["source_nonlinear_pr"] == {
        "pr": 740,
        "head": "21019df839557e64b8d8689f68cdf8fc685a3a36",
    }
    assert report["frozen_comparison"]["baseline"] == {"p": 4, "m": 4}
    assert report["frozen_comparison"]["reshape"] == {"p": 9, "m": 4}
    assert report["basis_dimension_changed"] is False
    assert report["candidate_velocity_changed"] is False
    assert report["fresh_714_path_data_used"] is False
    assert report["held_out_pde_residual_evaluated"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["pde_validated"] is False


def test_refined_metrics_are_stable_and_compact_balance_holds() -> None:
    report = audit.evaluate()
    assert report["stability"]["passes"] is True
    fine = report["finest_comparison"]
    for key in ("baseline_p4", "reshaped_p9"):
        assert abs(fine[key]["signed_integral_balance"]) <= audit.SIGNED_BALANCE_TOL


def test_p9_thins_geometry_but_concentrates_return_response() -> None:
    report = audit.evaluate()
    d = report["decision"]
    c = report["finest_comparison"]
    assert d["geometric_return_collar_thinner"] is True
    assert d["outward_squared_response_fraction_increased"] is True
    assert d["normalized_derivative_sharpness_increased"] is True
    assert d["normalized_peak_response_increased"] is True
    assert d["return_flow_concentrated_by_reshape"] is True
    assert d["geometric_collar_thinning_is_monotonic_free_gain"] is False
    assert d["increase_inner_exponent_again_justified_now"] is False
    assert c["profile_cosine"] < 0.5
    assert c["return_collar_reduction_fraction"] > 0.35


def test_axial_response_binds_live_reshape_implementation() -> None:
    z = audit.np.asarray([1.25, 1.55, 1.70], dtype=float)
    g = audit.axial_response(z, audit.RESHAPED_P)
    # C_r has the opposite sign of G for positive coefficient/radial bump.
    assert g[0] > 0.0
    assert g[1] > 0.0
    assert g[2] < 0.0
    _, z_star = audit.relocation.physical_sign_change(audit.RESHAPED_P, audit.M)
    assert z_star == pytest.approx(1.6020544561, abs=2e-10)
