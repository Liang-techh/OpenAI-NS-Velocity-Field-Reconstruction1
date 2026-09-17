import copy
import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_closure_frontier_checkpoint import (
    REFERENCE_COMPOSITE_SHA256,
    REGISTERED_PDE_THRESHOLD,
    build_checkpoint,
    load_checkpoint,
    validate_checkpoint,
    write_bundle,
)
from openai_ns_reconstruction.kokuno_five_moment_repair import KokunoFiveMomentRepair
from openai_ns_reconstruction.kokuno_oscillatory_phase import KokunoOscillatoryPhaseContract
from openai_ns_reconstruction.kokuno_reference_composite_checkpoint import (
    KokunoReferenceCompositeCandidate,
)


def _resign(payload):
    from openai_ns_reconstruction.kokuno_closure_frontier_checkpoint import _sha

    unsigned = dict(payload)
    unsigned.pop("checkpoint_sha256", None)
    payload["checkpoint_sha256"] = _sha(unsigned)
    return payload


def test_frontier_binds_latest_agent_routes_without_promotion():
    report = build_checkpoint()
    assert report["reference_stage"]["candidate_sha256"] == REFERENCE_COMPOSITE_SHA256
    assert report["fixed_gates"]["held_out_normalized_full_momentum_residual"] == REGISTERED_PDE_THRESHOLD
    assert report["agent1_five_moment_interface"]["probe"]["linearized_rank"] == 5
    assert report["agent1_five_moment_interface"]["probe"]["real_target_used"] is False
    assert np.isfinite(report["agent1_five_moment_interface"]["probe"]["condition_number"])
    assert len(report["agent2_source_phase_interface"]["manufactured_interface_probe"]["n_Phi"]) == 3
    assert report["upstream"]["agent3_signed_amplitude_curl_cycle"]["public_velocity_correction_materialized"] is True
    assert report["upstream"]["agent3_signed_amplitude_curl_cycle"]["accepted_for_next_cycle"] is False
    assert report["routing"]["apply_agent3_signed_curl_correction"] is False
    assert report["routing"]["compose_source_phase_into_agent2_complete_curl"] is False
    assert report["routing"]["run_formal_full_domain_gate"] is False
    assert report["states"]["leading_ready"] is False
    assert report["states"]["correction_ready"] is False
    assert report["states"]["velocity_export_ready"] is False
    assert report["states"]["pde_validated"] is False


def test_bundle_roundtrips_public_artifacts(tmp_path):
    output = tmp_path / "bundle"
    report = write_bundle(output)
    loaded = load_checkpoint(output / "closure_frontier_checkpoint.json")
    assert loaded["checkpoint_sha256"] == report["checkpoint_sha256"]
    assert KokunoReferenceCompositeCandidate.load_json(
        output / "reference_composite_candidate.json"
    ).sha256 == REFERENCE_COMPOSITE_SHA256
    assert KokunoFiveMomentRepair.load_json(
        output / "five_moment_repair_primitive.json"
    ).sha256 == report["agent1_five_moment_interface"]["primitive_sha256"]
    assert KokunoOscillatoryPhaseContract.load_json(
        output / "source_phase_contract.json"
    ).sha256 == report["agent2_source_phase_interface"]["contract_sha256"]
    uvw = np.asarray(report["reference_stage"]["velocity_probe"]["uvw"], dtype=float)
    assert uvw.shape == (3,)
    assert np.all(np.isfinite(uvw))
    assert np.linalg.norm(uvw) > 0.0
    with pytest.raises(FileExistsError):
        write_bundle(output)


def test_checkpoint_fails_closed_on_scientific_laundering(tmp_path):
    base = build_checkpoint()
    mutations = []

    relaxed = copy.deepcopy(base)
    relaxed["fixed_gates"]["held_out_normalized_full_momentum_residual"] = 1.0e-2
    mutations.append(_resign(relaxed))

    promoted = copy.deepcopy(base)
    promoted["states"]["pde_validated"] = True
    mutations.append(_resign(promoted))

    routed = copy.deepcopy(base)
    routed["routing"]["apply_agent3_signed_curl_correction"] = True
    mutations.append(_resign(routed))

    composed = copy.deepcopy(base)
    composed["routing"]["compose_source_phase_into_agent2_complete_curl"] = True
    mutations.append(_resign(composed))

    for index, payload in enumerate(mutations):
        path = tmp_path / f"mutated_{index}.json"
        path.write_text(json.dumps(payload), encoding="utf-8")
        with pytest.raises(ValueError):
            load_checkpoint(path)

    assert validate_checkpoint(base)["checkpoint_sha256"] == base["checkpoint_sha256"]
