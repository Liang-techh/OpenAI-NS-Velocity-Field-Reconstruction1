from __future__ import annotations

import inspect

import numpy as np

import openai_ns_reconstruction.kokuno_pa10_leading_axis_profile_independent_audit as audit_module
from openai_ns_reconstruction.kokuno_pa10_leading_axis_profile_contract import (
    KokunoPA10LeadingAxisProfileContract,
)
from openai_ns_reconstruction.kokuno_pa10_leading_axis_profile_independent_audit import (
    audit_leading_axis_profile,
)


def test_independent_profile_audit_passes_fresh_offgrid_and_fd4_checks() -> None:
    report = audit_leading_axis_profile(seed=9173441, sample_count=512)
    assert report["independent_profile_contract_passed"] is True
    assert report["failed_guards"] == []
    assert all(report["guards"].values())
    assert report["fresh_sample_count"] == 512
    assert report["fd_steps"] == [0.004, 0.002, 0.001]
    assert report["minimum_fd4_refinement_ratio"] >= 12.0
    assert report["maximum_finest_fd4_relative_rms"] <= 2.0e-8
    assert report["max_value_relative_error"] <= 5.0e-12
    assert report["max_derivative_relative_error"] <= 5.0e-11

    gates = report["gates"]
    assert gates["final_momentum_max_L2"] == 1.0e-3
    assert gates["final_divergence_max_L2"] == 1.0e-5

    truth = report["truth_boundary"]
    assert truth["repository_X_identified_with_source_Y"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["complete_kokuno_composite_velocity"] is False
    assert truth["matched_global_pressure_available"] is False
    assert truth["restricted_forcing_composite_available"] is False
    assert truth["leading_only_ns_residual_assessed"] is False
    assert truth["leading_plus_oscillatory_ns_residual_assessed"] is False
    assert truth["after_correction_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False


def test_independent_audit_is_deterministic_for_frozen_seed() -> None:
    first = audit_leading_axis_profile(seed=9173441, sample_count=256)
    second = audit_leading_axis_profile(seed=9173441, sample_count=256)
    assert first == second
    assert first["receipt_sha256"] == second["receipt_sha256"]


def test_audit_source_does_not_call_agent1_axis_state_or_fraction_helpers() -> None:
    source = inspect.getsource(audit_module)
    assert ".axis_state(" not in source
    assert "from fractions import Fraction" not in source
    assert "bridge_post_j_replay" not in source


class _DerivativeMutation:
    def __init__(self, base: KokunoPA10LeadingAxisProfileContract) -> None:
        self.base = base

    def configuration(self):
        return self.base.configuration()

    def values(self, *args, **kwargs):
        return self.base.values(*args, **kwargs)

    def derivatives(self, *args, **kwargs):
        raw = self.base.derivatives(*args, **kwargs)
        out = {name: np.asarray(value).copy() for name, value in raw.items()}
        out["Z_star_eta"] *= 0.999
        return out


class _ValueMutation:
    def __init__(self, base: KokunoPA10LeadingAxisProfileContract) -> None:
        self.base = base

    def configuration(self):
        return self.base.configuration()

    def values(self, *args, **kwargs):
        raw = self.base.values(*args, **kwargs)
        out = {name: np.asarray(value).copy() for name, value in raw.items()}
        out["W_star"] *= 0.999
        return out

    def derivatives(self, *args, **kwargs):
        return self.base.derivatives(*args, **kwargs)


def test_mutated_analytic_derivative_fails_closed() -> None:
    mutated = _DerivativeMutation(KokunoPA10LeadingAxisProfileContract())
    report = audit_leading_axis_profile(mutated, seed=9173441, sample_count=256)
    assert report["independent_profile_contract_passed"] is False
    assert "fresh_analytic_derivative_replay" in report["failed_guards"]


def test_mutated_public_value_fails_closed() -> None:
    mutated = _ValueMutation(KokunoPA10LeadingAxisProfileContract())
    report = audit_leading_axis_profile(mutated, seed=9173441, sample_count=256)
    assert report["independent_profile_contract_passed"] is False
    assert "fresh_value_formula_replay" in report["failed_guards"]
