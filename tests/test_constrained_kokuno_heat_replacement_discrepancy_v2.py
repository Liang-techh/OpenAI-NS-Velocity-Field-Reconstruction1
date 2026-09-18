import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_heat_replacement_discrepancy import (
    KokunoHeatReplacementDiscrepancy,
)
from openai_ns_reconstruction.kokuno_heat_replacement_discrepancy_independent import (
    IndependentHeatDiscrepancyAudit,
)
from openai_ns_reconstruction.kokuno_heat_replacement_discrepancy_v2 import (
    KokunoHeatReplacementDiscrepancyV2,
    SCHEMA_V2,
)


def test_v2_restores_exactly_the_missing_first_order_cp_tail_channel():
    old = KokunoHeatReplacementDiscrepancy()
    corrected = KokunoHeatReplacementDiscrepancyV2()
    audit = IndependentHeatDiscrepancyAudit()
    etas = np.asarray([-0.8, -0.37, 0.0, 0.2, 0.8])

    before = old.discrepancy(etas)
    after = corrected.discrepancy(etas)
    change = after - before

    # The correction is local to Delta C_p.  Delta S and Delta I_sub retain
    # the independently audited v1 values bit-for-bit through the same path.
    np.testing.assert_allclose(change[:, 1:], 0.0, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(
        change[:, 0], corrected.missing_first_order_cp_tail(etas), rtol=2.0e-12, atol=2.0e-18
    )

    # Agent-4 isolated the v1 deficit as one half of the source first-order
    # tail.  V2's added copy must equal that half through an independent
    # integration path rather than by replaying the production quadrature.
    independent_full_first_order = np.asarray(
        [audit.first_order_tail_cp(corrected, float(eta)) for eta in etas]
    )
    np.testing.assert_allclose(
        change[:, 0],
        0.5 * independent_full_first_order,
        rtol=3.0e-5,
        atol=2.0e-18,
    )


def test_v2_matches_independent_discrepancy_across_eta():
    corrected = KokunoHeatReplacementDiscrepancyV2()
    audit = IndependentHeatDiscrepancyAudit()
    etas = np.asarray([-0.8, -0.37, 0.0, 0.2, 0.8])
    public = corrected.discrepancy(etas)
    reference = np.stack(
        [audit.discrepancy(corrected, float(eta)) for eta in etas], axis=0
    )
    relative = np.abs(public - reference) / np.maximum(np.abs(reference), 1.0e-300)

    assert float(np.max(relative[:, 0])) < 2.0e-6
    assert float(np.max(relative[:, 1])) < 2.0e-8
    assert float(np.max(relative[:, 2])) < 2.0e-6
    np.testing.assert_allclose(
        corrected.discrepancy(0.2),
        np.asarray([-2.21476226138e-09, 6.25423009649e-06, -2.63159129494]),
        rtol=2.0e-6,
        atol=2.0e-13,
    )


def test_v2_normalized_target_uses_corrected_discrepancy():
    model = KokunoHeatReplacementDiscrepancyV2()
    eta = np.asarray([-0.4, 0.2])
    X_star = 17.0
    e_star = 0.3
    target = model.normalized_repair_target(eta, X_star=X_star, e_star=e_star)
    scales = np.asarray([e_star**2, X_star * e_star**2, X_star**1.5 * e_star])
    np.testing.assert_allclose(
        target * scales + model.discrepancy(eta), 0.0, rtol=0.0, atol=2.0e-14
    )


def test_v2_serialization_has_distinct_fail_closed_identity(tmp_path):
    old = KokunoHeatReplacementDiscrepancy()
    corrected = KokunoHeatReplacementDiscrepancyV2()
    assert corrected.to_payload()["schema"] == SCHEMA_V2
    assert corrected.sha256 != old.sha256

    destination = corrected.save_json(tmp_path / "heat_discrepancy_v2.json")
    loaded = KokunoHeatReplacementDiscrepancyV2.load_json(destination)
    assert loaded.sha256 == corrected.sha256
    assert loaded.to_payload() == corrected.to_payload()

    payload = json.loads(destination.read_text(encoding="utf-8"))
    payload["implementation_revision"]["delta_s_changed"] = True
    with pytest.raises(ValueError, match="payload hash or content mismatch"):
        KokunoHeatReplacementDiscrepancyV2.from_payload(payload)
