import copy
import json

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_supported_candidate import (
    Eq45SupportedVelocityCandidate,
)
from openai_ns_reconstruction.constrained_eq45_supported_phi10_compact_quartic_blend_delivery_capsule import (
    DEPENDENCY_HEAD,
    EXPECTED_BASE_SHA256,
    REPRESENTATIVE_SHA256,
    build_candidate,
    delivery_capsule,
    load_recipe,
    write_delivery_bundle,
)
from openai_ns_reconstruction.constrained_eq45_supported_phi10_compact_quartic_blend_temporal_mode import (
    Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate,
)


def _probes():
    return np.asarray(
        [
            [0.37, 0.11, -0.22],
            [0.95, -0.35, 0.48],
            [1.45, 0.15, -0.62],
            [1.72, 0.00, 0.30],
            [0.44, -0.52, 1.68],
            [1.70, 0.17, -1.72],
        ],
        dtype=float,
    )


def test_blend_capsule_binds_family_identities_and_temporal_derivative_evidence():
    recipe = load_recipe()
    capsule = delivery_capsule(0.5)

    assert recipe["base_dependency_head"] == DEPENDENCY_HEAD
    assert capsule["base_supported_sha256"] == EXPECTED_BASE_SHA256
    assert capsule["declared_blend_weight"] == 0.5
    assert capsule["blend_weight_selected_by_capsule"] is False

    for weight, key in ((0.0, "0.000000"), (0.5, "0.500000"), (1.0, "1.000000")):
        assert build_candidate(weight).sha256 == REPRESENTATIVE_SHA256[key]

    evidence = capsule["temporal_derivative_evidence"]
    assert evidence["audit_pr"] == 187
    assert evidence["actions_run"] == 35234995835
    assert evidence["probe_seed"] == 914617
    assert evidence["finest_normalized_rms"]["0.000000"] == pytest.approx(
        0.00019000133460634884, rel=0.0, abs=1.0e-18
    )
    assert evidence["finest_normalized_rms"]["0.500000"] == pytest.approx(
        0.0006986742153697458, rel=0.0, abs=1.0e-18
    )
    assert evidence["finest_normalized_rms"]["1.000000"] == pytest.approx(
        0.0012568302260969515, rel=0.0, abs=1.0e-18
    )
    assert evidence["formal_pde_gate_assessed"] is False
    assert evidence["pde_validated"] is False

    siblings = capsule["sibling_evidence"]
    assert siblings["restricted_force_pr"] == 188
    assert siblings["restricted_force_status"] == "open_unconsumed"
    assert siblings["spectral_fingerprint_pr"] == 189
    assert siblings["spectral_fingerprint_status"] == "open_unconsumed"
    assert siblings["numeric_summary_bound_here"] is False

    status = capsule["status"]
    assert status["velocity_export_ready"] is True
    assert status["caller_declared_blend_weight"] is True
    assert status["blend_weight_selected"] is False
    assert status["candidate_selection_resolved"] is False
    assert status["visualization_ready"] is False
    assert status["pde_validated"] is False
    assert status["openai_field_identified"] is False
    assert status["blowup_proved"] is False


def test_blend_capsule_replays_public_uvws_static_anchors_affinity_and_grid():
    points = _probes()
    quartic = build_candidate(0.0)
    midpoint = build_candidate(0.5)
    compact = build_candidate(1.0)
    static = Eq45SupportedVelocityCandidate(parent=Eq45VelocityCandidate.seed())

    for time in (0.25, 0.3125, 0.375, 0.4375, 0.5, 0.75):
        q_value = quartic.at_points(points, time)
        m_value = midpoint.at_points(points, time)
        c_value = compact.at_points(points, time)
        assert np.allclose(m_value, 0.5 * (q_value + c_value), rtol=0.0, atol=2.0e-16)

    for candidate in (quartic, midpoint, compact):
        for time in (0.5, 0.625, 0.75):
            assert np.array_equal(
                candidate.at_points(points, time), static.at_points(points, time)
            )

    arbitrary = build_candidate(0.25)
    grid = arbitrary.grid(
        np.asarray([-0.4, 0.4]),
        np.asarray([-0.3, 0.3]),
        np.asarray([-0.2, 0.2]),
        np.asarray([0.25, 0.5, 0.75]),
    )
    assert grid.shape == (3, 2, 2, 2, 3)
    assert np.all(np.isfinite(grid))


def test_blend_delivery_bundle_roundtrips_arbitrary_caller_weight(tmp_path):
    original = build_candidate(0.25)
    paths = write_delivery_bundle(tmp_path, 0.25)
    loaded = Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate.load_json(
        paths["candidate"]
    )
    capsule = json.loads(paths["capsule"].read_text(encoding="utf-8"))

    assert loaded.sha256 == original.sha256 == capsule["candidate_sha256"]
    assert loaded.blend_weight == capsule["declared_blend_weight"] == 0.25
    assert capsule["blend_weight_selected_by_capsule"] is False

    points = _probes()
    times = np.asarray([0.25, 0.3125, 0.375, 0.4375, 0.5, 0.75])
    assert np.array_equal(
        loaded.at_points(points, times), original.at_points(points, times)
    )


def test_blend_recipe_fails_closed_on_identity_evidence_sibling_or_claim_promotion(tmp_path):
    recipe = load_recipe()
    cases = []

    mutated = copy.deepcopy(recipe)
    mutated["representation"]["representative_candidate_sha256"]["0.500000"] = "0" * 64
    cases.append(mutated)

    mutated = copy.deepcopy(recipe)
    mutated["representation"]["blend_weight_selected"] = True
    cases.append(mutated)

    mutated = copy.deepcopy(recipe)
    mutated["temporal_derivative_evidence"]["audit_exact_head"] = "0" * 40
    cases.append(mutated)

    mutated = copy.deepcopy(recipe)
    mutated["temporal_derivative_evidence"]["probe_seed"] = 914618
    cases.append(mutated)

    mutated = copy.deepcopy(recipe)
    mutated["temporal_derivative_evidence"]["formal_pde_gate_assessed"] = True
    cases.append(mutated)

    mutated = copy.deepcopy(recipe)
    mutated["sibling_evidence"]["restricted_force_status"] = "consumed"
    cases.append(mutated)

    mutated = copy.deepcopy(recipe)
    mutated["sibling_evidence"]["numeric_summary_bound_here"] = True
    cases.append(mutated)

    mutated = copy.deepcopy(recipe)
    mutated["status"]["blend_weight_selected"] = True
    cases.append(mutated)

    mutated = copy.deepcopy(recipe)
    mutated["status"]["visualization_ready"] = True
    cases.append(mutated)

    mutated = copy.deepcopy(recipe)
    mutated["status"]["pde_validated"] = True
    cases.append(mutated)

    for index, payload in enumerate(cases):
        path = tmp_path / f"mutated-{index}.json"
        path.write_text(json.dumps(payload), encoding="utf-8")
        with pytest.raises(ValueError):
            load_recipe(path)


@pytest.mark.parametrize("bad_weight", [-0.01, 1.01, np.nan, np.inf])
def test_blend_capsule_rejects_invalid_caller_weights(bad_weight):
    with pytest.raises(ValueError):
        build_candidate(bad_weight)
