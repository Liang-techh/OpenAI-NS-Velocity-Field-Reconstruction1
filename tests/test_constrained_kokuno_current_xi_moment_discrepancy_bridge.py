from __future__ import annotations

import copy
import hashlib
import inspect
import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_current_xi_moment_discrepancy_bridge import (
    CurrentXiMomentBridgeError,
    UPSTREAM_AGENT1_HEAD,
    UPSTREAM_AGENT1_PR,
    UPSTREAM_AGENT1_SCHEMA,
    UPSTREAM_AGENT1_SOURCE_BLOB,
    UPSTREAM_AGENT1_SOURCE_PATH,
    UpstreamA1Identity,
    materialize_current_xi_moment_bridge,
    truth_boundary,
)


def _canonical_sha(value: object) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _valid_report() -> dict[str, object]:
    source = {
        "repository": "KokunoYumeto/yang-mills-interacting-workbench",
        "commit": "143f6773feb424ad9ed3a8d116653200f20346b7",
        "path": "navier-stokes/navier_stokes_workbench.tex",
        "release": "zenodo:22678406",
        "release_date": "2026-09-09",
    }
    formulas = {
        "five_prefix_moments": "fixture-public-five-moment-map",
        "PA15_scaling": "fixture-public-PA15-scaling",
    }
    realization = {
        "outer_schedule": "fixture-autonomous-source-hierarchy-schedule",
        "prefix_quadrature": "fixture-fixed-quadrature",
    }
    configuration = {
        "schema": UPSTREAM_AGENT1_SCHEMA,
        "bridge_semantic_sha256": "1" * 64,
        "outer_schedule": {"schema": "fixture-schedule"},
        "quadrature_order": 48,
    }
    semantic = _canonical_sha(
        {
            "schema": UPSTREAM_AGENT1_SCHEMA,
            "source_commit": source["commit"],
            "source_formulas": formulas,
            "numerical_realization": realization,
            "configuration": configuration,
        }
    )
    eta = [-0.5, 0.0, 0.5]
    rows = [
        [0.10, -0.20, 0.30, -0.40, 0.50],
        [0.02, 0.03, -0.04, 0.05, -0.06],
        [-0.07, 0.08, -0.09, 0.10, 0.11],
    ]
    schedule_sha = "a" * 64
    pa16 = []
    for i, value in enumerate(eta):
        pa16.append(
            {
                "eta": value,
                "ell_i": 1.0 + 0.1 * i,
                "G_i": -0.2 + 0.2 * i,
                "incoming_scaled_discrepancy": rows[i],
                "outer_schedule_sha256": schedule_sha,
                "source_T_sh_lower_bound_verified": False,
            }
        )
    return {
        "schema": UPSTREAM_AGENT1_SCHEMA,
        "source": source,
        "source_formulas": formulas,
        "numerical_realization": realization,
        "configuration": configuration,
        "geometry": {"X_i": 110.0},
        "eta_probe": eta,
        "incoming_PA15_discrepancy": rows,
        "pa16_inputs": pa16,
        "semantic_sha256": semantic,
        "truth_boundary": {
            "current_lineage_profile_0_to_Xi_executable": True,
            "public_five_prefix_moment_map_executable": True,
            "public_PA15_scaling_executable": True,
            "candidate_side_upstream_five_moment_discrepancy_at_Xi_materialized": True,
            "candidate_side_PA16_input_tuple_materialized": True,
            "source_prepared_appendixA_pressure_stress_materialized": False,
            "source_prepared_upstream_five_moment_discrepancy_materialized": False,
            "source_T_sh_lower_bound_verified": False,
            "five_moment_repair_applied": False,
            "inner_to_outer_join_completed": False,
            "outer_global_leading_velocity_materialized": False,
            "matched_global_pressure_materialized": False,
            "restricted_forcing_materialized": False,
            "heldout_ns_residual_assessed": False,
            "same_protocol_comparable_to_st006": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }


class _Backend:
    def __init__(
        self,
        report: dict[str, object] | None = None,
        identity: UpstreamA1Identity | None = None,
    ) -> None:
        self._report = _valid_report() if report is None else report
        self._identity = identity or UpstreamA1Identity(
            pr_number=UPSTREAM_AGENT1_PR,
            exact_head=UPSTREAM_AGENT1_HEAD,
            source_path=UPSTREAM_AGENT1_SOURCE_PATH,
            source_blob=UPSTREAM_AGENT1_SOURCE_BLOB,
            report_schema=UPSTREAM_AGENT1_SCHEMA,
        )

    def upstream_identity(self) -> UpstreamA1Identity:
        return self._identity

    def current_xi_moment_report(self) -> dict[str, object]:
        return copy.deepcopy(self._report)


def test_materializes_public_pa16_rows_without_ns_promotion() -> None:
    receipt = materialize_current_xi_moment_bridge(_Backend())
    assert len(receipt.samples) == 3
    assert receipt.any_correction_needed
    assert receipt.max_raw_l2_norm > 0.0
    assert len(receipt.handoff_sha256) == 64

    first = receipt.samples[0]
    raw = np.asarray(first.incoming_scaled_discrepancy)
    eta = first.eta
    expected = np.asarray(
        [raw[0], raw[2] - 4.0 * eta * raw[1], raw[1], raw[3] - 8.0 * eta * raw[0], raw[4]]
    )
    np.testing.assert_array_equal(first.pa16_row_transformed_discrepancy, expected)

    truth = receipt.to_dict()["truth_boundary"]
    assert truth["current_candidate_side_xi_source_moment_discrepancy_ingested"] is True
    assert truth["current_candidate_side_pa16_rows_materialized"] is True
    assert truth["historical_pa16_solver_duplicated"] is False
    assert truth["discrepancy_from_complete_ns_defect"] is False
    assert truth["authorized_for_gain_gated_ns_stage"] is False
    assert truth["current_real_ns_correction_velocity_materialized"] is False
    assert truth["heldout_normalized_ns_residual_assessed"] is False
    assert truth["residual_reduction_claimed"] is False
    assert truth["pde_validated"] is False


def test_exact_upstream_identity_is_required() -> None:
    bad = UpstreamA1Identity(
        pr_number=UPSTREAM_AGENT1_PR,
        exact_head="0" * 40,
        source_path=UPSTREAM_AGENT1_SOURCE_PATH,
        source_blob=UPSTREAM_AGENT1_SOURCE_BLOB,
        report_schema=UPSTREAM_AGENT1_SCHEMA,
    )
    with pytest.raises(CurrentXiMomentBridgeError, match="pinned Agent-1"):
        materialize_current_xi_moment_bridge(_Backend(identity=bad))


def test_semantic_sha_tamper_is_rejected() -> None:
    report = _valid_report()
    report["configuration"]["quadrature_order"] = 49  # type: ignore[index]
    with pytest.raises(CurrentXiMomentBridgeError, match="semantic SHA"):
        materialize_current_xi_moment_bridge(_Backend(report=report))


def test_top_level_and_pa16_discrepancy_must_be_identical() -> None:
    report = _valid_report()
    report["pa16_inputs"][1]["incoming_scaled_discrepancy"][2] += 1.0e-3  # type: ignore[index]
    with pytest.raises(CurrentXiMomentBridgeError, match="top-level incoming discrepancy"):
        materialize_current_xi_moment_bridge(_Backend(report=report))


def test_truth_boundary_laundering_is_rejected() -> None:
    report = _valid_report()
    report["truth_boundary"]["heldout_ns_residual_assessed"] = True  # type: ignore[index]
    with pytest.raises(CurrentXiMomentBridgeError, match="attempted promotion"):
        materialize_current_xi_moment_bridge(_Backend(report=report))


def test_unverified_tsh_may_not_be_promoted() -> None:
    report = _valid_report()
    report["pa16_inputs"][0]["source_T_sh_lower_bound_verified"] = True  # type: ignore[index]
    with pytest.raises(CurrentXiMomentBridgeError, match="T_sh"):
        materialize_current_xi_moment_bridge(_Backend(report=report))


def test_public_api_exposes_no_scientific_tuning_surface() -> None:
    parameters = set(inspect.signature(materialize_current_xi_moment_bridge).parameters)
    assert parameters == {"backend"}
    truth = truth_boundary()
    assert truth["forbidden_public_parameters_absent"] is True
    assert truth["final_normalized_momentum_gate"] == 1.0e-3
    assert truth["final_normalized_divergence_gate"] == 1.0e-5
