from __future__ import annotations

import copy
import math
from typing import Callable

import pytest

from openai_ns_reconstruction.audit_kokuno_first_cell_extrapolation_observability_scope import (
    audit_contract,
    hidden_cell_mechanics_witness,
    load_contract,
)


def test_registered_scope_audits_fail_closed() -> None:
    report = audit_contract()
    assert report["scope_ok"] is True
    assert report["truth_boundary"] == {
        "affine_quadratic_agreement_is_error_certificate": False,
        "formal_first_cell_accuracy_verified": False,
        "formal_axis_based_inverse_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }


def test_hidden_cell_witness_has_zero_estimator_disagreement_but_nonzero_formal_integral() -> None:
    witness = hidden_cell_mechanics_witness(a=0.2, amplitude=1.0)
    assert witness["mechanics_only"] is True
    assert witness["real_source_claim"] is False
    assert witness["public_source_fact"] is False
    assert witness["sample_values"] == [0.0, 0.0, 0.0, 0.0]

    e1 = witness["channels"]["e1"]
    e2 = witness["channels"]["e2"]
    for channel in (e1, e2):
        assert channel["degree_1_estimate"] == 0.0
        assert channel["degree_2_estimate"] == 0.0
        assert channel["degree_1_degree_2_disagreement"] == 0.0
        assert channel["formal_weighted_integral"] > 0.0
        assert channel["formal_minus_degree_2"] == channel["formal_weighted_integral"]

    assert math.isclose(e1["formal_weighted_integral"], 0.2**6 / 30.0, rel_tol=2e-15)
    assert math.isclose(e2["formal_weighted_integral"], 0.2**7 / 105.0, rel_tol=2e-15)


def test_hidden_cell_amplitude_changes_formal_integral_without_changing_observed_nodes() -> None:
    unit = hidden_cell_mechanics_witness(a=0.2, amplitude=1.0)
    scaled = hidden_cell_mechanics_witness(a=0.2, amplitude=7.0)
    assert unit["sample_values"] == scaled["sample_values"]
    for exponent in (1, 2):
        key = f"e{exponent}"
        assert unit["channels"][key]["degree_2_estimate"] == scaled["channels"][key]["degree_2_estimate"] == 0.0
        assert math.isclose(
            scaled["channels"][key]["formal_weighted_integral"],
            7.0 * unit["channels"][key]["formal_weighted_integral"],
            rel_tol=2e-15,
        )


def _set_realization(contract: dict, key: str, value: object) -> None:
    contract["audited_realization"][key] = value


def _set_status(contract: dict, key: str, value: object) -> None:
    contract["status_firewall"][key] = value


def _set_cr001(contract: dict, key: str, value: object) -> None:
    contract["canonical_cr001"][key] = value


Mutation = Callable[[dict], None]


def _mutations() -> list[tuple[str, Mutation]]:
    return [
        (
            "affine-quadratic disagreement laundered into formal error bound",
            lambda c: _set_realization(c, "linear_quadratic_disagreement_is_formal_error_bound", True),
        ),
        (
            "quadratic mechanics exactness laundered into real-source error bound",
            lambda c: _set_realization(c, "polynomial_exactness_through_degree_2_is_real_source_error_bound", True),
        ),
        (
            "real-source convergence promoted without refinement evidence",
            lambda c: _set_realization(c, "real_source_first_cell_convergence_verified", True),
        ),
        (
            "formal first-cell accuracy promoted",
            lambda c: _set_realization(c, "formal_first_cell_accuracy_verified", True),
        ),
        (
            "formal axis inverse promoted",
            lambda c: _set_realization(c, "formal_axis_based_inverse_verified", True),
        ),
        (
            "full-domain moment promoted",
            lambda c: _set_realization(c, "formal_full_domain_moment_verified", True),
        ),
        (
            "axis regularity promoted",
            lambda c: _set_realization(c, "axis_regularity_verified", True),
        ),
        (
            "contract falsely claims sub-first-cell observations",
            lambda c: _set_realization(c, "samples_inside_open_interval_0_rmin_used", True),
        ),
        (
            "hidden-cell mechanics promoted to real-source claim",
            lambda c: c["observability_boundary"].__setitem__("hidden_cell_witness_is_real_source_claim", True),
        ),
        (
            "hidden-cell mechanics promoted to public-source fact",
            lambda c: c["observability_boundary"].__setitem__("hidden_cell_witness_is_public_source_fact", True),
        ),
        (
            "autonomous witness copied into public-source provenance",
            lambda c: c["provenance"]["public_source_fact"].append(c["provenance"]["autonomous_design"][3]),
        ),
        (
            "momentum max gate relaxed",
            lambda c: _set_cr001(c, "pde_residual_max", 0.002),
        ),
        (
            "momentum L2 gate relaxed",
            lambda c: _set_cr001(c, "pde_residual_L2", 0.002),
        ),
        (
            "divergence gate relaxed",
            lambda c: _set_cr001(c, "divergence_max", 2e-5),
        ),
        (
            "free residual-defined force enabled",
            lambda c: _set_cr001(c, "residual_defined_pointwise_free_force_allowed", True),
        ),
        (
            "amplitude collapse enabled",
            lambda c: _set_cr001(c, "amplitude_collapse_allowed", True),
        ),
        (
            "post-hoc threshold relaxation enabled",
            lambda c: _set_cr001(c, "post_hoc_threshold_relaxation_allowed", True),
        ),
        (
            "scoped audit promotes PDE validation",
            lambda c: _set_status(c, "pde_validated", True),
        ),
        (
            "scoped audit promotes visual correspondence",
            lambda c: _set_status(c, "visual_correspondence_verified", True),
        ),
        (
            "scoped audit promotes paper exactness",
            lambda c: _set_status(c, "paper_exact", True),
        ),
        (
            "scoped audit promotes OpenAI-field identity",
            lambda c: _set_status(c, "openai_field_identified", True),
        ),
        (
            "PDE pending blocks callable velocity delivery",
            lambda c: _set_status(c, "pde_pending_blocks_callable_velocity_delivery", True),
        ),
        (
            "audited source identity drift",
            lambda c: c.__setitem__("audited_source_blob", "0" * 40),
        ),
    ]


@pytest.mark.parametrize("label,mutate", _mutations(), ids=lambda item: item if isinstance(item, str) else None)
def test_scope_mutations_fail_closed(label: str, mutate: Mutation) -> None:
    contract = copy.deepcopy(load_contract())
    mutate(contract)
    with pytest.raises(AssertionError):
        audit_contract(contract)
