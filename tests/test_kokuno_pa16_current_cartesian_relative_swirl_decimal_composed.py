from __future__ import annotations

import copy
import inspect
import json
import math
from decimal import Decimal

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa16_current_cartesian_relative_swirl_decimal_composed import (
    DECIMAL_DIGITS,
    PARENT_EXACT_HEAD,
    KokunoPA16CurrentCartesianRelativeSwirlDecimalComposed,
)


def _candidate() -> KokunoPA16CurrentCartesianRelativeSwirlDecimalComposed:
    return KokunoPA16CurrentCartesianRelativeSwirlDecimalComposed()


def _active_radius(c: KokunoPA16CurrentCartesianRelativeSwirlDecimalComposed) -> float:
    log_X = c.log_X_flatten_end + c.compensator.y1
    return math.exp(0.5 * (math.log(2.0) + log_X))


def test_public_decimal_velocity_retains_parent_delta_that_binary64_erases():
    c = _candidate()
    radius = _active_radius(c)
    base, delta = c.split.velocity_split(radius, 0.0, 0.0, 0.0)
    base_d, delta_d = c.velocity_decimal_split(radius, 0.0, 0.0, 0.0)
    total = c.velocity(radius, 0.0, 0.0, 0.0)

    active = [idx for idx in np.ndindex(delta.shape) if float(delta[idx]) != 0.0]
    assert active
    np.testing.assert_array_equal(base + delta, base)
    for idx in active:
        assert isinstance(total[idx], Decimal)
        assert total[idx] != base_d[idx]
        assert total[idx] - base_d[idx] == delta_d[idx]
        assert delta_d[idx] != Decimal(0)


def test_vectorized_shape_axis_regularity_and_finite_decimal_entries():
    c = _candidate()
    radius = _active_radius(c)
    total = c.velocity(
        np.array([radius, 0.0]),
        np.array([0.0, 0.0]),
        np.array([0.0, 0.0]),
        np.array([0.0, 0.0]),
    )
    assert total.shape == (2, 3)
    assert total.dtype == object
    for value in total.reshape(-1):
        assert isinstance(value, Decimal)
        assert value.is_finite()
    assert total[1, 0] == Decimal(0)
    assert total[1, 1] == Decimal(0)


def test_low_dimensional_total_profile_and_analytic_logx_derivatives_retain_delta():
    c = _candidate()
    log_X = c.log_X_flatten_end + c.compensator.y1 + 0.027
    eta = np.array([-0.17, 0.11])
    parent = c.split.split_profile_logX(np.full(eta.shape, log_X), eta)
    total = c.profile_total_decimal_logX(np.full(eta.shape, log_X), eta)
    for key in ("E", "F", "U", "DlogX_E", "DlogX_F", "DlogX_U"):
        assert total[key].shape == eta.shape
        assert all(v.is_finite() for v in total[key].reshape(-1))

    for i in range(eta.size):
        base_E = Decimal.from_float(float(parent["E_base"][i]))
        # The public total must not collapse to the base at an active bump point.
        assert total["E"][i] != base_E
        assert total["F"][i] != Decimal.from_float(float(parent["F_base"][i]))


def test_representation_report_machine_locks_precision_qualified_composition_only():
    c = _candidate()
    report = c.representation_report()
    assert report["decimal_digits"] == DECIMAL_DIGITS
    assert report["precision_margin_available"] >= 0
    assert report["active_cartesian_components"] >= 1
    assert report["binary64_erases_all_active_components"] is True
    assert report["decimal_total_retains_all_active_components"] is True
    assert report["decimal_total_minus_base_replays_delta"] is True
    truth = report["truth_boundary"]
    assert truth["current_cartesian_relative_swirl_composed"] is True
    assert truth["precision_qualified_public_velocity_materialized"] is True
    assert truth["binary64_total_relative_swirl_sum_is_resolved"] is False
    assert truth["unified_global_cartesian_velocity_export_ready"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False


def test_decimal_json_export_preserves_exact_public_total_after_save_load(tmp_path):
    c = _candidate()
    config_path = tmp_path / "candidate.json"
    payload = c.save_configuration(config_path)
    loaded = KokunoPA16CurrentCartesianRelativeSwirlDecimalComposed.load_configuration(
        config_path
    )
    assert loaded.configuration() == payload
    assert loaded.semantic_sha256 == c.semantic_sha256

    radius = _active_radius(loaded)
    x = np.array([radius, 0.0])
    y = np.zeros(2)
    z = np.zeros(2)
    t = np.zeros(2)
    expected = loaded.velocity(x, y, z, t)
    export_path = tmp_path / "velocity_decimal.json"
    exported = loaded.export_velocity_decimal_json(export_path, x, y, z, t)
    disk = json.loads(export_path.read_text())
    assert disk == exported
    assert disk["candidate_semantic_sha256"] == loaded.semantic_sha256
    assert disk["shape"] == list(expected.shape)
    replay = np.asarray([Decimal(v) for v in disk["values_row_major"]], dtype=object)
    replay = replay.reshape(expected.shape)
    for idx in np.ndindex(expected.shape):
        assert replay[idx] == expected[idx]


def test_configuration_provenance_and_precision_fail_closed(tmp_path):
    c = _candidate()
    payload = c.configuration()
    assert payload["parent_exact_head"] == PARENT_EXACT_HEAD
    assert payload["decimal_digits"] == DECIMAL_DIGITS

    mutated = copy.deepcopy(payload)
    mutated["decimal_digits"] = DECIMAL_DIGITS - 1
    with pytest.raises(ValueError, match="configuration/provenance drift"):
        KokunoPA16CurrentCartesianRelativeSwirlDecimalComposed.from_configuration(mutated)

    mutated = copy.deepcopy(payload)
    mutated["truth_boundary"]["binary64_total_relative_swirl_sum_is_resolved"] = True
    with pytest.raises(ValueError, match="configuration/provenance drift"):
        KokunoPA16CurrentCartesianRelativeSwirlDecimalComposed.from_configuration(mutated)

    mutated = copy.deepcopy(payload)
    mutated["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError, match="configuration/provenance drift"):
        KokunoPA16CurrentCartesianRelativeSwirlDecimalComposed.from_configuration(mutated)


def test_public_velocity_has_no_precision_or_scientific_tuning_knobs():
    parameters = set(
        inspect.signature(KokunoPA16CurrentCartesianRelativeSwirlDecimalComposed.velocity).parameters
    )
    forbidden = {
        "precision",
        "decimal_digits",
        "residual",
        "forcing",
        "force",
        "pressure",
        "viscosity",
        "threshold",
        "gain",
        "optimizer",
        "tolerance",
    }
    assert not (parameters & forbidden)
    assert parameters == {"self", "x", "y", "z", "t"}
