import copy
import importlib.util
import json
import tempfile
from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction import kokuno_current_main_pulse_logx_leading_oscillatory_identity as m


def _has_a1():
    try:
        return importlib.util.find_spec(m.A1_MODULE) is not None
    except ModuleNotFoundError:
        return False


def test_public_contract_and_truth_boundary():
    contract = m.public_contract()
    assert contract["public_inputs"] == ["x", "y", "z", "t"]
    assert contract["forbidden_velocity_inputs_present"] == []
    assert contract["full_source_xi_11_current_cartesian_leading_materialized"] is True
    for key in (
        "source_exact_main_pulse_amplitude_materialized",
        "terminal_global_leading_velocity_materialized",
        "velocity_export_ready",
        "paper_exact",
        "pde_validated",
    ):
        assert contract[key] is False
    truth = m._truth()
    assert truth["current_logX_main_pulse_leading_plus_frozen_complete_curl_oscillation_materialized"] is True
    assert truth["full_source_xi_11_current_cartesian_leading_materialized"] is True
    assert truth["agent3_correction_velocity_composed"] is False
    assert truth["complete_ns_residual_assessed"] is False


def test_enforcer_rejects_truth_and_scope_promotion():
    receipt = {
        "schema": m.RECEIPT_SCHEMA,
        "task": m.TASK,
        "source_provenance": m.SOURCE,
        "truth_boundary": m._truth(),
        "domain": {
            "finite_parent_xi_max": 9.0,
            "source_main_xi_end": 11.0,
            "full_source_xi_11_current_cartesian_leading_materialized": True,
            "terminal_global_leading_velocity_materialized": False,
        },
        "checks": {
            "inner_additive_max_abs": 0.0,
            "tail_additive_max_abs": 0.0,
            "inner_oscillatory_vector_rms": 1.0e-4,
            "axis_oscillatory_max_abs": 0.0,
            "outside_support_oscillatory_max_abs": 0.0,
            "save_load_velocity_replay_max_abs": 0.0,
            "xi_11_endpoint_velocity_finite": True,
            "post_xi_11_boundary_fail_closed": True,
        },
    }
    m.enforce_receipt(receipt)

    promoted = copy.deepcopy(receipt)
    promoted["truth_boundary"]["pde_validated"] = True
    with pytest.raises(AssertionError):
        m.enforce_receipt(promoted)

    promoted = copy.deepcopy(receipt)
    promoted["domain"]["terminal_global_leading_velocity_materialized"] = True
    with pytest.raises(AssertionError):
        m.enforce_receipt(promoted)

    malformed = copy.deepcopy(receipt)
    malformed["domain"]["finite_parent_xi_max"] = 11.0
    with pytest.raises(AssertionError):
        m.enforce_receipt(malformed)


@pytest.mark.skipif(not _has_a1(), reason="exact A1 #1107 sibling supplied by dedicated CI")
def test_exact_logx_tail_composition_receipt():
    field = m.default_field()
    assert 0.0 < field.finite_xi_max < m.SOURCE_MAIN_XI_END
    x, y, z, t, log_x, xi, radius = m._tail_points(field)
    assert np.all(xi > field.finite_xi_max)
    assert np.all(xi < m.SOURCE_MAIN_XI_END)
    assert np.all(np.isfinite(log_x))
    assert np.all(np.isfinite(radius))

    evaluated = m._comp._evaluate_with_backend(field.leading_backend, x, y, z, t)
    np.testing.assert_allclose(
        field.velocity(x, y, z, t),
        evaluated.leading + evaluated.oscillatory,
        rtol=0.0,
        atol=m.COMP_ATOL,
    )
    assert np.max(np.abs(evaluated.oscillatory)) == 0.0

    receipt = m.materialize_receipt()
    m.enforce_receipt(receipt)
    assert receipt["checks"]["xi_11_endpoint_velocity_finite"] is True
    assert receipt["checks"]["post_xi_11_boundary_fail_closed"] is True
    assert min(receipt["probes"]["tail_xi"]) > field.finite_xi_max


@pytest.mark.skipif(not _has_a1(), reason="exact A1 #1107 sibling supplied by dedicated CI")
def test_inner_signal_axis_support_and_batch_shapes():
    field = m.default_field()
    total = field.velocity(m.INNER[:, 0], m.INNER[:, 1], m.INNER[:, 2], m.INNER_T)
    assert total.shape == (len(m.INNER), 3)
    assert np.all(np.isfinite(total))

    inner = m._comp._evaluate_with_backend(
        field.leading_backend, m.INNER[:, 0], m.INNER[:, 1], m.INNER[:, 2], m.INNER_T
    )
    rms = float(np.sqrt(np.mean(np.sum(inner.oscillatory ** 2, axis=-1))))
    assert rms >= m.SIGNAL_FLOOR

    axis = m._comp._evaluate_with_backend(
        field.leading_backend, m.AXIS[:, 0], m.AXIS[:, 1], m.AXIS[:, 2], m.AXIS_T
    )
    assert np.max(np.abs(axis.oscillatory)) == 0.0
    outside = m._comp.velocity_osc_batch(m.OUTSIDE, m.OUTSIDE_T)
    assert np.max(np.abs(outside)) == 0.0


@pytest.mark.skipif(not _has_a1(), reason="exact A1 #1107 sibling supplied by dedicated CI")
def test_save_load_rejects_rehashed_runtime_source_and_a1_mutations():
    field = m.default_field()
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "candidate.json"
        field.save_candidate(path)
        raw = json.loads(path.read_text())

        changed = copy.deepcopy(raw)
        osc = changed["configuration"]["oscillatory_runtime"]
        osc["payload"]["mutation_probe"] = "x"
        osc["payload_sha256"] = m._id._sha256(osc["payload"])
        changed["semantic_sha256"] = m._id._sha256(changed["configuration"])
        path.write_text(json.dumps(changed))
        with pytest.raises(ValueError, match="oscillatory runtime"):
            m.CurrentMainPulseLogXLeadingOscillatoryField.load_candidate(path)

        changed = copy.deepcopy(raw)
        changed["configuration"]["source_provenance"]["corrected_release_date"] = "2099-01-01"
        changed["semantic_sha256"] = m._id._sha256(changed["configuration"])
        path.write_text(json.dumps(changed))
        with pytest.raises(ValueError, match="configuration/provenance/truth"):
            m.CurrentMainPulseLogXLeadingOscillatoryField.load_candidate(path)

        changed = copy.deepcopy(raw)
        changed["configuration"]["agent1"]["head"] = "0" * 40
        changed["semantic_sha256"] = m._id._sha256(changed["configuration"])
        path.write_text(json.dumps(changed))
        with pytest.raises(ValueError, match="A1 #1107 provenance"):
            m.CurrentMainPulseLogXLeadingOscillatoryField.load_candidate(path)
