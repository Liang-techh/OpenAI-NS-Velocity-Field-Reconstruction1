import json

from openai_ns_reconstruction.constrained_bipolar_capped_fresh import run
from openai_ns_reconstruction.eq45_supported_delivery import Eq45SupportedDeliveryField


def test_fresh_refinement_reloads_public_candidate_and_catches_mutation(tmp_path):
    output = tmp_path / "fresh_refinement.json"
    report = run(
        output=output,
        seed=9172631,
        points=64,
        times=(0.3125, 0.6875),
        steps=(0.02, 0.01, 0.005),
    )

    loaded = Eq45SupportedDeliveryField.load_candidate(
        "artifacts/bipolar_joint_capped/candidate.json"
    )
    assert report["candidate_sha256"] == loaded.sha256
    assert report["steps"] == [0.02, 0.01, 0.005]
    assert len(report["rows"]) == 6
    assert report["preregistered_thresholds_applied_without_change"] is True
    assert report["formal_registered_divergence_gate_assessed"] is False
    assert report["formal_full_momentum_gate_assessed"] is False
    assert report["necessary_axisymmetric_pressure_PDE_condition_assessed"] is True
    assert report["necessary_axisymmetric_pressure_PDE_condition_passed"] is False
    assert report["mutation"]["caught"] is True
    assert report["classification"]["visualization_candidate_only"] is True
    assert report["classification"]["pde_validated"] is False

    persisted = json.loads(output.read_text())
    assert persisted["candidate_sha256"] == report["candidate_sha256"]
    assert persisted["seed"] == 9172631
    assert persisted["points"] == 64
