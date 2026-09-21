import copy

import numpy as np

from openai_ns_reconstruction.kokuno_a4_actual_xi_prefix_moments_independent_audit import (
    A1_HEAD,
    A1_SOURCE_BLOB,
    FROZEN_FINAL_GATES,
    PROTOCOL,
    TRUTH_BOUNDARY,
    _checksum_valid_C_drift_rejected,
    _metric_gate,
    _row_metrics,
    _simpson_segment,
    public_contract,
)
from openai_ns_reconstruction.kokuno_pa10_actual_xi_prefix_moments import (
    KokunoPA10ActualXiPrefixMoments,
)


class _ManufacturedProfile:
    def profile_values(self, X, eta):
        x = np.asarray(X, dtype=float)
        return {
            "F_actual_prefix": np.full_like(x, 2.0),
            "U_actual_prefix": 1.0 + x,
            "E_actual_prefix": 2.0 * np.sqrt(2.0 * x),
        }


def test_public_contract_freezes_exact_a1_identity_and_final_gates():
    contract = public_contract()
    assert contract["agent1_exact_head"] == A1_HEAD
    assert contract["agent1_source_blob"] == A1_SOURCE_BLOB
    assert contract["protocol"]["seed"] == 9173631
    assert contract["protocol"]["simpson_panels_per_segment"] == (48, 96, 192)
    assert FROZEN_FINAL_GATES["normalized_momentum_max"] == 1.0e-3
    assert FROZEN_FINAL_GATES["divergence_volume_l2"] == 1.0e-5
    assert TRUTH_BOUNDARY["heldout_complete_ns_residual_assessed"] is False
    assert TRUTH_BOUNDARY["pde_validated"] is False


def test_composite_simpson_locks_public_moment_density_signs_and_weights():
    got = _simpson_segment(_ManufacturedProfile(), 0.0, 0.0, 2.0, 2)
    expected = np.asarray(
        [
            4.0,
            8.0,
            56.0 / 3.0,
            2.0 / 3.0,
            8.0,
        ]
    )
    np.testing.assert_allclose(got, expected, rtol=0.0, atol=2.0e-14)


def test_metric_firewall_rejects_one_percent_row_drift_and_sign_flip():
    independent = np.asarray(
        [
            [1.0, 2.0, -3.0, 4.0, 5.0],
            [1.5, 2.5, -2.0, 3.5, 4.5],
            [2.0, 3.0, -1.0, 3.0, 4.0],
        ]
    )
    assert _metric_gate(_row_metrics(independent, independent))

    one_percent = independent.copy()
    one_percent[:, 0] *= 0.99
    assert not _metric_gate(_row_metrics(one_percent, independent))

    sign_flip = independent.copy()
    sign_flip[:, 2] *= -1.0
    assert not _metric_gate(_row_metrics(sign_flip, independent))


def test_checksum_valid_outer_C_drift_fails_closed():
    field = KokunoPA10ActualXiPrefixMoments()
    configuration = copy.deepcopy(field.configuration())
    assert configuration["outer_schedule"]["parameters"]["C"] == field.C
    assert _checksum_valid_C_drift_rejected(configuration)


def test_protocol_keeps_final_thresholds_out_of_scoped_moment_gate():
    assert PROTOCOL.relative_rms_gate == 5.0e-3
    assert PROTOCOL.relative_sampled_max_gate == 2.0e-2
    assert PROTOCOL.refinement_ratio_gate == 6.0
    assert FROZEN_FINAL_GATES == {
        "normalized_momentum_max": 1.0e-3,
        "normalized_momentum_volume_l2": 1.0e-3,
        "divergence_max": 1.0e-5,
        "divergence_volume_l2": 1.0e-5,
    }
