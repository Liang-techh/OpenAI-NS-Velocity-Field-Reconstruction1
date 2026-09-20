from __future__ import annotations

import inspect

import numpy as np

import openai_ns_reconstruction.kokuno_pa10_physical_center_independent_audit as audit_module
from openai_ns_reconstruction.kokuno_pa10_physical_center_independent_audit import (
    FD_STEPS,
    SIMPSON_LEVELS,
    audit_physical_center_profile,
)
from openai_ns_reconstruction.kokuno_pa10_physical_center_profile_contract import (
    KokunoPA10PhysicalCenterProfileContract,
)


class _MutatedF0Profile:
    def __init__(self, base: KokunoPA10PhysicalCenterProfileContract) -> None:
        self._base = base

    def __getattr__(self, name: str):
        return getattr(self._base, name)

    def values(self, X, eta):
        out = dict(self._base.values(X, eta))
        out["F_0"] = np.asarray(out["F_0"], dtype=float) * 0.999
        return out


def test_independent_physical_center_audit_is_fail_closed_and_detects_C_obstruction() -> None:
    report = audit_physical_center_profile(sample_count=256)
    assert report["audit_passed"] is True
    assert report["failed_guards"] == []
    assert report["mapping_formula_independently_consistent"] is True
    assert report["configured_C_real_axis_normalization_obstruction_detected"] is True
    assert report["normalization_margin_phi_minus_C"] > 0.0
    assert report["physical_center_profile_independently_admitted"] is False
    assert report["simpson_levels"] == list(SIMPSON_LEVELS)
    assert report["fd_steps"] == list(FD_STEPS)
    assert report["gates"]["final_momentum_max_L2"] == 1e-3
    assert report["gates"]["final_divergence_max_L2"] == 1e-5
    truth = report["truth_boundary"]
    assert truth["source_complex_C_normalization_certified"] is False
    assert truth["configured_C_real_axis_source_condition_admitted"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["cartesian_spacetime_velocity_materialized"] is False
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["pde_validated"] is False


def test_downward_public_F0_mutation_is_detected() -> None:
    mutated = _MutatedF0Profile(KokunoPA10PhysicalCenterProfileContract())
    report = audit_physical_center_profile(mutated, sample_count=64)
    assert report["guards"]["fresh_physical_mapping_replay"] is False
    assert report["mapping_formula_independently_consistent"] is False
    assert report["audit_passed"] is False


def test_audit_does_not_call_public_partial_fraction_or_self_identity_helpers() -> None:
    source = inspect.getsource(audit_module)
    assert ".zeta_primitive(" not in source
    assert ".phi_star(" not in source
    assert ".incompressibility_defect(" not in source
    assert "bridge_post_j_replay" not in source


def test_three_level_numerical_paths_are_preregistered() -> None:
    assert SIMPSON_LEVELS == (64, 128, 256)
    assert FD_STEPS == (0.004, 0.002, 0.001)
