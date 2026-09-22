from __future__ import annotations

import inspect
from pathlib import Path
import tempfile

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_relative_swirl_compensator import (
    KokunoPublicRelativeSwirlCompensator,
    RelativeSwirlCompensatorSolution,
)
from openai_ns_reconstruction import (
    kokuno_a4_relative_swirl_compensator_independent_audit as audit,
)


def _save_reload():
    original = KokunoPublicRelativeSwirlCompensator()
    with tempfile.TemporaryDirectory(prefix="k4_val_120_test_") as td:
        path = Path(td) / "candidate.json"
        original.save_configuration(path)
        loaded = KokunoPublicRelativeSwirlCompensator.load_configuration(path)
        return original, loaded


def test_exact_1154_save_reload_independent_audit_passes_frozen_scoped_gates():
    original, loaded = _save_reload()
    report = audit.audit_loaded_compensator(loaded, original)
    audit.enforce_preregistered_gates(report)
    assert report["audit_pass"] is True
    assert report["failures"] == []
    assert report["protocol"]["simpson_nodes"] == [2049, 4097, 8193]
    assert report["protocol"]["angular_targets"] == list(audit.ANGULAR_TARGETS)
    fine = report["resolution_reports"][-1]
    assert fine["row_metrics"]["relative_max"] <= audit.ROW_RELATIVE_MAX_GATE
    assert (
        fine["coefficient_metrics"]["relative_rms"]
        <= audit.COEFFICIENT_RELATIVE_RMS_GATE
    )
    assert (
        fine["coefficient_metrics"]["relative_max"]
        <= audit.COEFFICIENT_RELATIVE_MAX_GATE
    )


def test_materialize_receipt_executes_real_save_reload_and_keeps_pde_false():
    report = audit.materialize_receipt()
    audit.enforce_preregistered_gates(report)
    assert report["save_reload_coefficient_absolute_max"] <= 2.0e-12
    truth = report["truth_boundary"]
    assert truth["relative_swirl_algebra_scoped_assessed"] is True
    assert truth["current_lineage_angular_entry_I_materialized"] is False
    assert truth["current_cartesian_relative_swirl_composed"] is False
    assert truth["complete_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False


class _AngularRowCorruptionProxy:
    def __init__(self, base):
        self._base = base

    def __getattr__(self, name):
        return getattr(self._base, name)

    @property
    def angular_row_scaled(self):
        out = np.asarray(self._base.angular_row_scaled, float).copy()
        out[1] *= 1.02
        return out


def test_independent_reference_detects_production_angular_row_corruption():
    original, loaded = _save_reload()
    report = audit.audit_loaded_compensator(_AngularRowCorruptionProxy(loaded), original)
    assert report["audit_pass"] is False
    assert "fine_row_relative_max_le_2e-6" in report["failures"]


class _PressureSignCorruptionProxy:
    def __init__(self, base):
        self._base = base

    def __getattr__(self, name):
        return getattr(self._base, name)

    @property
    def pressure_linear_row_scaled(self):
        return -np.asarray(self._base.pressure_linear_row_scaled, float)


def test_independent_reference_detects_pressure_sign_corruption():
    original, loaded = _save_reload()
    report = audit.audit_loaded_compensator(_PressureSignCorruptionProxy(loaded), original)
    assert report["audit_pass"] is False
    assert "fine_row_relative_max_le_2e-6" in report["failures"]


class _BranchSwapProxy:
    def __init__(self, base):
        self._base = base

    def __getattr__(self, name):
        return getattr(self._base, name)

    def solve(self, angular_target):
        target = np.asarray(angular_target, float)
        sol = self._base.solve(target)
        mask = target != 0.0
        c1 = np.asarray(sol.c1, float).copy()
        c2 = np.asarray(sol.c2, float).copy()
        c1[mask] += 0.025
        c2[mask] -= 0.004
        return RelativeSwirlCompensatorSolution(
            c1=c1,
            c2=c2,
            dc1_dtarget=np.asarray(sol.dc1_dtarget, float),
            dc2_dtarget=np.asarray(sol.dc2_dtarget, float),
            angular_residual=np.asarray(sol.angular_residual, float),
            pressure_residual=np.asarray(sol.pressure_residual, float),
            discriminant=np.asarray(sol.discriminant, float),
        )


def test_wrong_branch_or_root_swap_is_not_accepted_as_same_evidence():
    original, loaded = _save_reload()
    report = audit.audit_loaded_compensator(_BranchSwapProxy(loaded), original)
    assert report["audit_pass"] is False
    assert (
        "fine_coefficient_relative_rms_le_5e-5" in report["failures"]
        or "fine_coefficient_relative_max_le_2e-4" in report["failures"]
    )


def test_geometry_and_provenance_drift_fail_closed():
    original = KokunoPublicRelativeSwirlCompensator()
    payload = original.configuration()
    payload["bump_width"] = 0.3001
    with pytest.raises(ValueError):
        KokunoPublicRelativeSwirlCompensator.from_configuration(payload)


@pytest.mark.parametrize(
    "forbidden",
    [
        "residual",
        "forcing",
        "pressure",
        "viscosity",
        "threshold",
        "tolerance",
        "step",
        "nodes",
        "target",
        "gain",
        "optimizer",
    ],
)
def test_public_scientific_entrypoint_has_no_tuning_knobs(forbidden):
    sig = inspect.signature(audit.audit_loaded_compensator)
    names = [name.lower() for name in sig.parameters]
    assert names == ["candidate", "pre_serialization_reference"]
    assert all(forbidden not in name for name in names)


def test_final_project_gates_and_truth_boundary_are_immutable():
    report = audit.materialize_receipt()
    assert report["gates"]["final_project_momentum_gate_unchanged"] == 1.0e-3
    assert report["gates"]["final_project_divergence_gate_unchanged"] == 1.0e-5
    truth = report["truth_boundary"]
    assert truth["leading_only_ns_residual_assessed"] is False
    assert truth["leading_plus_oscillatory_ns_residual_assessed"] is False
    assert truth["after_correction_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
