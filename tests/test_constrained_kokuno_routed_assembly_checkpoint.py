import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_core_composite_checkpoint import (
    KokunoCoreCompositeCandidate,
)
from openai_ns_reconstruction.kokuno_heat_exterior_profile import (
    KokunoHeatExteriorProfile,
)
from openai_ns_reconstruction.kokuno_routed_assembly_checkpoint import (
    PDE_THRESHOLD,
    global_velocity,
    load_routed_assembly_checkpoint,
    write_routed_assembly_checkpoint,
)


def test_routed_assembly_bundle_round_trips_and_keeps_truth_boundary(tmp_path):
    receipt = write_routed_assembly_checkpoint(tmp_path)
    loaded = load_routed_assembly_checkpoint(tmp_path / "assembly_checkpoint.json")
    assert loaded["receipt_sha256"] == receipt["receipt_sha256"]

    routing = receipt["routing"]
    assert routing["oscillatory_amplitude"] == pytest.approx(0.125)
    assert routing["oscillatory_phase"] == pytest.approx(0.0)
    assert routing["phase0_retained_after_agent2_holdout"] is True
    assert routing["agent3_mean_correction_applied"] is False
    assert routing["heat_exterior_combined_with_core"] is False

    a2 = receipt["upstream"]["agent2_phase_screen"]
    assert a2["selected_generalizes_vs_phase0"] is False
    assert a2["accepted_for_next_cycle"] is False
    a3 = receipt["upstream"]["agent3_core_mean_cycle"]
    assert a3["accepted_for_next_cycle"] is False
    assert a3["held_out_total_mean_defect_ratio"] > 1.0
    assert a3["held_out_theta_ratio"] > 1.0

    a4 = receipt["upstream"]["agent4_core_independent"]
    assert a4["formal_full_domain_gate_assessed"] is False
    assert a4["local_reference_threshold_met"] is False
    assert a4["pde_validated"] is False
    assert a4["stage_ratios"]["leading_plus_oscillatory_over_leading"] < 1.0
    assert (
        a4["stage_ratios"]["after_rejected_mean_correction_over_leading_plus_oscillatory"]
        < 1.0
    )

    state = receipt["stage_state"]
    assert state["leading_core_ready"] is True
    assert state["heat_exterior_component_ready"] is True
    assert state["leading_ready"] is False
    assert state["correction_ready"] is False
    assert state["core_velocity_export_ready"] is True
    assert state["velocity_export_ready"] is False
    assert state["pde_validated"] is False
    assert receipt["registered_gate"]["held_out_normalized_ns_residual_threshold"] == PDE_THRESHOLD
    assert receipt["registered_gate"]["assessed"] is False

    core = KokunoCoreCompositeCandidate.load_json(tmp_path / "core_candidate.json")
    heat = KokunoHeatExteriorProfile.load_json(tmp_path / "heat_exterior_component.json")
    assert core.sha256 == receipt["artifacts"]["core_candidate"]["sha256"]
    assert heat.sha256 == receipt["artifacts"]["heat_exterior_component"]["sha256"]

    core_u = core.velocity(0.1, 0.05, 0.0, 0.5)
    heat_u = heat.velocity(2.0, 0.75, 0.0, 0.4)
    assert np.asarray(core_u).shape == (3,)
    assert np.asarray(heat_u).shape == (3,)
    assert np.isfinite(core_u).all() and np.linalg.norm(core_u) > 0.0
    assert np.isfinite(heat_u).all() and np.linalg.norm(heat_u) > 0.0

    with pytest.raises(RuntimeError, match="matching are not complete"):
        global_velocity(0.1, 0.0, 0.0, 0.5)


def test_routed_assembly_receipt_fails_closed_on_promotion(tmp_path):
    write_routed_assembly_checkpoint(tmp_path)
    path = tmp_path / "assembly_checkpoint.json"
    payload = json.loads(path.read_text())
    payload["stage_state"]["pde_validated"] = True
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="SHA mismatch"):
        load_routed_assembly_checkpoint(path)
