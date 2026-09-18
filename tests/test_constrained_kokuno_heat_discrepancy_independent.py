import numpy as np

from openai_ns_reconstruction.kokuno_heat_replacement_discrepancy import (
    KokunoHeatReplacementDiscrepancy,
)
from openai_ns_reconstruction.kokuno_heat_replacement_discrepancy_independent import (
    IndependentHeatDiscrepancyAudit,
)


def test_independent_log_space_reference_matches_two_channels_and_exposes_cp_tail_factor():
    model = KokunoHeatReplacementDiscrepancy()
    audit = IndependentHeatDiscrepancyAudit()

    reference = audit.discrepancy(model, 0.2)
    np.testing.assert_allclose(
        reference,
        np.asarray([-2.21476226138e-09, 6.25423009649e-06, -2.63159129494]),
        rtol=4.0e-9,
        atol=2.0e-13,
    )

    public = model.discrepancy(0.2)
    relative = np.abs((public - reference) / reference)
    assert 2.3e-3 < relative[0] < 2.6e-3
    assert relative[1] < 2.0e-8
    assert relative[2] < 2.0e-6

    first_tail = audit.first_order_tail_cp(model, 0.2)
    # The public-minus-independent discrepancy is exactly the missing half of
    # the first-order y>=3 Delta C_p tail term to numerical accuracy.
    ratio = (public[0] - reference[0]) / (-0.5 * first_tail)
    assert abs(ratio - 1.0) < 2.0e-5
    np.testing.assert_allclose(
        public[0] + 0.5 * first_tail,
        reference[0],
        rtol=6.0e-8,
        atol=2.0e-16,
    )


def test_independent_audit_is_stable_across_eta_and_does_not_promote_pde_state():
    model = KokunoHeatReplacementDiscrepancy()
    report = IndependentHeatDiscrepancyAudit().compare_public(model)

    rel = np.asarray(report["max_relative_by_component"])
    assert 2.3e-3 < rel[0] < 2.6e-3
    assert rel[1] < 2.0e-8
    assert rel[2] < 2.0e-6

    ratios = np.asarray(report["cp_missing_half_factor_ratio"])
    np.testing.assert_allclose(ratios, 1.0, rtol=0.0, atol=3.0e-5)
    assert report["formal_full_domain_pde_gate_assessed"] is False
    assert report["pde_validated"] is False


def test_independent_heat_series_refines_without_using_public_heat_quadrature():
    model = KokunoHeatReplacementDiscrepancy()
    coarse = IndependentHeatDiscrepancyAudit(series_order=10).discrepancy(model, 0.37)
    fine = IndependentHeatDiscrepancyAudit(series_order=14).discrepancy(model, 0.37)
    np.testing.assert_allclose(coarse, fine, rtol=2.0e-11, atol=2.0e-14)
