import copy
import json

import pytest

from openai_ns_reconstruction.kokuno_independent_routing_checkpoint import (
    AGENT4_UPSTREAM,
    REGISTERED_PDE_THRESHOLD,
    STATES,
    _sha,
    build_checkpoint,
    summarize_agent4_report,
    validate_checkpoint,
)
from openai_ns_reconstruction.kokuno_independent_signed_curl_cycle import TASK


def _agent4_fixture():
    return {
        "task": TASK,
        "validator": {
            "agent3_cr006_residual_operator_used": False,
            "training_loss_or_tensor_used": False,
            "fixed_residual_reference": 1.0e-3,
            "fixed_divergence_reference": 1.0e-5,
        },
        "sampling": {"annulus_count": 32, "axis_count": 16},
        "aggregates": {
            "0.001": {
                "leading_only": {
                    "sample_rms": 10.4529319480,
                    "sample_max": 23.0848371365,
                    "divergence_rms": 3.60520e-5,
                    "divergence_max": 9.56164e-5,
                },
                "leading_plus_oscillatory": {
                    "sample_rms": 10.3067008431,
                    "sample_max": 22.7787616638,
                    "divergence_rms": 3.58091e-5,
                    "divergence_max": 9.48993e-5,
                },
                "after_signed_correction": {
                    "sample_rms": 10.3075981486,
                    "sample_max": 22.7787616638,
                    "divergence_rms": 3.66924e-5,
                    "divergence_max": 1.02282e-4,
                },
            }
        },
        "finest_stage_ratios": {
            "oscillatory_over_leading_rms": 0.98601051785,
            "corrected_over_oscillatory_rms": 1.00008706040,
            "corrected_over_leading_rms": 0.98609636032,
        },
        "signed_correction_support": {
            "outside_probe_exact_zero": True,
            "annulus_nontrivial_rms": 0.0042096510,
        },
        "local_reference_checks": {
            "finest_corrected_normalized_sample_rms_le_1e_minus_3": False,
            "finest_corrected_normalized_sample_max_le_1e_minus_3": False,
            "finest_corrected_divergence_rms_le_1e_minus_5": False,
            "finest_corrected_divergence_max_le_1e_minus_5": False,
        },
        "truth_boundary": {
            "formal_full_domain_gate_assessed": False,
            "pde_validated": False,
        },
    }


def _resign(payload):
    unsigned = dict(payload)
    unsigned.pop("checkpoint_sha256", None)
    payload["checkpoint_sha256"] = _sha(unsigned)
    return payload


def test_independent_report_summary_preserves_negative_routing():
    summary = summarize_agent4_report(_agent4_fixture())
    assert summary["oscillatory_over_leading_rms"] < 1.0
    assert summary["corrected_over_oscillatory_rms"] > 1.0
    assert summary["corrected_sample_rms"] > REGISTERED_PDE_THRESHOLD
    assert summary["correction_outside_probe_exact_zero"] is True
    assert summary["correction_annulus_nontrivial_rms"] > 0.0
    assert summary["formal_full_domain_gate_assessed"] is False
    assert summary["pde_validated"] is False


def test_checkpoint_binds_agent4_without_promoting_signed_correction():
    checkpoint = build_checkpoint(_agent4_fixture())
    assert checkpoint["agent4_independent_signed_curl_audit"]["provenance"] == AGENT4_UPSTREAM
    assert checkpoint["agent4_independent_signed_curl_audit"]["replay_matches_exact_head_reference"] is True
    assert checkpoint["agent4_independent_signed_curl_audit"]["used_for_parameter_selection"] is False
    assert checkpoint["routing"]["apply_agent3_signed_curl_correction"] is False
    assert checkpoint["routing"]["agent3_rejection_independently_supported"] is True
    assert checkpoint["routing"]["retune_single_signed_column_now"] is False
    assert checkpoint["routing"]["run_formal_full_domain_gate"] is False
    assert checkpoint["states"] == STATES
    assert checkpoint["states"]["independent_signed_curl_validation_ready"] is True
    assert checkpoint["states"]["correction_ready"] is False
    assert checkpoint["states"]["velocity_export_ready"] is False
    assert checkpoint["states"]["pde_validated"] is False
    assert validate_checkpoint(checkpoint)["checkpoint_sha256"] == checkpoint["checkpoint_sha256"]


def test_checkpoint_fails_closed_on_validation_laundering(tmp_path):
    base = build_checkpoint(_agent4_fixture())
    mutations = []

    relaxed = copy.deepcopy(base)
    relaxed["fixed_gates"]["held_out_normalized_full_momentum_residual"] = 1.0e-2
    mutations.append(_resign(relaxed))

    routed = copy.deepcopy(base)
    routed["routing"]["apply_agent3_signed_curl_correction"] = True
    mutations.append(_resign(routed))

    selected = copy.deepcopy(base)
    selected["agent4_independent_signed_curl_audit"]["used_for_parameter_selection"] = True
    mutations.append(_resign(selected))

    promoted = copy.deepcopy(base)
    promoted["states"]["pde_validated"] = True
    mutations.append(_resign(promoted))

    for index, payload in enumerate(mutations):
        path = tmp_path / f"mutation_{index}.json"
        path.write_text(json.dumps(payload), encoding="utf-8")
        with pytest.raises(ValueError):
            validate_checkpoint(json.loads(path.read_text(encoding="utf-8")))
